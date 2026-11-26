/*
 * dec.c - streaming decoder.
 *
 * Input is buffered until a whole block is available (its header gives
 * comp_len). Each compressed block is decoded in two phases:
 *   1. the four rANS streams (literals, ll/ml/of codes) are decoded into
 *      flat arrays with tight table-lookup loops;
 *   2. the sequence loop reads code + extra bits, copies the literal run
 *      and the match into the history buffer (16-byte wildcopy where safe).
 * The history buffer keeps the last 2^window_log bytes so matches can
 * reach into previous blocks.
 */
#include "lzr.h"
#include "rans.h"
#include "codes.h"
#include "crc32.h"

#define MAX_SEQ (LZR_BLOCK_SIZE / LZR_MIN_MATCH + 1)
#define MAX_COMP (2 * LZR_BLOCK_SIZE + 65536)   /* sanity cap on comp_len */

enum { ST_HEADER, ST_BLOCKS, ST_DONE };

struct lzr_dctx {
    int state;
    lzr_buf in;
    size_t ipos;
    uint8_t *hist;
    size_t hcap, wsize, op;
    uint32_t reps[LZR_NUM_REPS];
    uint32_t crc;
    uint64_t total_out;
    int64_t content_size;
    uint8_t *lits, *llc, *mlc, *ofc;
    rans_dslot *tbl;
    code_info ct[VAL_CODES];
};

lzr_dctx *lzr_dctx_new(void)
{
    lzr_dctx *d = (lzr_dctx *)calloc(1, sizeof(*d));
    if (!d) return NULL;
    d->reps[0] = 1; d->reps[1] = 4; d->reps[2] = 8;
    d->content_size = -1;
    d->lits = (uint8_t *)malloc(LZR_BLOCK_SIZE + LZR_SLACK);
    d->llc = (uint8_t *)malloc(MAX_SEQ);
    d->mlc = (uint8_t *)malloc(MAX_SEQ);
    d->ofc = (uint8_t *)malloc(MAX_SEQ);
    d->tbl = (rans_dslot *)malloc(RANS_M * sizeof(rans_dslot));
    if (!d->lits || !d->llc || !d->mlc || !d->ofc || !d->tbl) { lzr_dctx_free(d); return NULL; }
    code_table_init(d->ct);
    return d;
}

void lzr_dctx_free(lzr_dctx *d)
{
    if (!d) return;
    free(d->in.p); free(d->hist);
    free(d->lits); free(d->llc); free(d->mlc); free(d->ofc); free(d->tbl);
    free(d);
}

int lzr_frame_header(const uint8_t *src, size_t n, int *window_log, uint8_t *flags, int64_t *content_size)
{
    if (memcmp(src, LZR_MAGIC, n < 4 ? n : 4) != 0) return LZR_ERR_CORRUPT;
    if (n < 8) return 0;                          /* need more input */
    if (src[4] != LZR_VERSION) return LZR_ERR_UNSUPPORTED;
    if (flags) *flags = src[5];
    if (src[5] & ~LZR_KNOWN_FLAGS) return LZR_ERR_UNSUPPORTED;
    if (src[6] < LZR_WLOG_MIN || src[6] > LZR_WLOG_MAX || src[7]) return LZR_ERR_CORRUPT;
    if (window_log) *window_log = src[6];
    if (src[5] & LZR_FLAG_SIZE) {
        if (n < 16) return 0;
        uint64_t sz = rd64le(src + 8);
        if (sz > (uint64_t)INT64_MAX) return LZR_ERR_CORRUPT;
        if (content_size) *content_size = (int64_t)sz;
        return 16;
    }
    if (content_size) *content_size = -1;
    return 8;
}

/* Copy len bytes from op - dist to op; regions may overlap (dist < len).
 * May write up to 15 bytes past op + len (buffer has LZR_SLACK). */
static inline void copy_match(uint8_t *op, size_t dist, size_t len)
{
    const uint8_t *src = op - dist;
    if (dist >= 16) {
        for (size_t k = 0; k < len; k += 16) memcpy(op + k, src + k, 16);
    } else if (dist == 1) {
        memset(op, *src, len);
    } else {
        /* Short period: replicate one period, then double the span. */
        size_t k = 0;
        while (k < len && k < dist) { op[k] = src[k]; k++; }
        while (k < len) {
            size_t chunk = k < len - k ? k : len - k;
            memcpy(op + k, op, chunk);
            k += chunk;
        }
    }
}

static int decode_block(lzr_dctx *d, const uint8_t *in, size_t avail, uint8_t *dst, size_t raw)
{
    if (avail < 8) return LZR_ERR_CORRUPT;
    size_t ns = rd32le(in), nlit = rd32le(in + 4), k = 8, c;
    if (nlit > raw || ns > MAX_SEQ || ns * LZR_MIN_MATCH + nlit > raw) return LZR_ERR_CORRUPT;

    c = stream_decode(in + k, avail - k, d->lits, nlit, ALPHA_LIT, d->tbl);
    if (nlit && !c) return LZR_ERR_CORRUPT;
    k += c;
    uint8_t *codes[3] = {d->llc, d->mlc, d->ofc};
    for (int s = 0; s < 3; s++) {
        c = stream_decode(in + k, avail - k, codes[s], ns, ALPHA_CODE, d->tbl);
        if (ns && !c) return LZR_ERR_CORRUPT;
        k += c;
    }
    if (avail - k < 4) return LZR_ERR_CORRUPT;
    size_t xlen = rd32le(in + k);
    k += 4;
    if (xlen != avail - k) return LZR_ERR_CORRUPT;
    bitr br = {0, 0, in + k, in + k + xlen, 0};

    const code_info *ct = d->ct;
    const uint8_t *lp = d->lits, *lend = d->lits + nlit;
    uint8_t *o = dst, *const oend = dst + raw;
    uint32_t *reps = d->reps;

    for (size_t i = 0; i < ns; i++) {
        unsigned lc = d->llc[i], mc = d->mlc[i], oc = d->ofc[i];
        if (lc >= VAL_CODES || mc >= VAL_CODES || oc >= OF_CODES) return LZR_ERR_CORRUPT;
        size_t ll = ct[lc].base + br_get(&br, ct[lc].bits);
        size_t ml = ct[mc].base + br_get(&br, ct[mc].bits) + LZR_MIN_MATCH;
        uint32_t dist;
        if (oc >= LZR_NUM_REPS) {
            const code_info *e = &ct[oc - LZR_NUM_REPS];
            dist = e->base + br_get(&br, e->bits) + 1;
            reps[2] = reps[1]; reps[1] = reps[0]; reps[0] = dist;
        } else {
            dist = reps[oc];
            if (oc == 1) { reps[1] = reps[0]; reps[0] = dist; }
            else if (oc == 2) { reps[2] = reps[1]; reps[1] = reps[0]; reps[0] = dist; }
        }
        if (ll > (size_t)(lend - lp) || ll + ml > (size_t)(oend - o)) return LZR_ERR_CORRUPT;

        /* literal run (lits buffer and history both have slack) */
        if (ll <= 16) memcpy(o, lp, 16);
        else memcpy(o, lp, ll);
        o += ll;
        lp += ll;

        if (dist > (size_t)(o - d->hist) || dist > d->wsize) return LZR_ERR_CORRUPT;
        copy_match(o, dist, ml);
        o += ml;
    }
    size_t rest = (size_t)(lend - lp);
    if (rest != (size_t)(oend - o)) return LZR_ERR_CORRUPT;
    memcpy(o, lp, rest);
    if (br.err || br.p != br.end) return LZR_ERR_CORRUPT;
    return LZR_OK;
}

static int init_history(lzr_dctx *d, int wlog)
{
    d->wsize = (size_t)1 << wlog;
    d->hcap = 2 * d->wsize + LZR_BLOCK_SIZE;
    d->hist = (uint8_t *)big_alloc(d->hcap + LZR_SLACK, 0);
    return d->hist ? LZR_OK : LZR_ERR_NOMEM;
}

/* Process one unit at in.p + ipos. Returns 1 if consumed, 0 if more input
 * is needed, 2 at the verified end, or LZR_ERR_*. */
static int step(lzr_dctx *d, lzr_buf *out)
{
    const uint8_t *p = d->in.p + d->ipos;
    size_t avail = d->in.len - d->ipos;

    if (d->state == ST_HEADER) {
        int wlog;
        int r = lzr_frame_header(p, avail, &wlog, NULL, &d->content_size);
        if (r <= 0) return r;
        if (init_history(d, wlog)) return LZR_ERR_NOMEM;
        d->ipos += (size_t)r;
        d->state = ST_BLOCKS;
        return 1;
    }

    if (avail < LZR_BLOCK_HDR) return 0;
    uint8_t type = p[0];
    size_t raw = rd32le(p + 1), comp = rd32le(p + 5);
    if (type > 2 || raw > LZR_BLOCK_SIZE || comp > MAX_COMP) return LZR_ERR_CORRUPT;
    if (type == 0 && comp != raw) return LZR_ERR_CORRUPT;
    if (type == 2 && (raw != 0 || comp != 4)) return LZR_ERR_CORRUPT;
    if (type == 1 && raw == 0) return LZR_ERR_CORRUPT;
    if (avail - LZR_BLOCK_HDR < comp) return 0;
    const uint8_t *payload = p + LZR_BLOCK_HDR;

    if (type == 2) {
        if (d->content_size >= 0 && (uint64_t)d->content_size != d->total_out) return LZR_ERR_CORRUPT;
        d->ipos += LZR_BLOCK_HDR + comp;
        d->state = ST_DONE;
        return rd32le(payload) == d->crc ? 2 : LZR_ERR_CHECKSUM;
    }
    if (d->content_size >= 0 && d->total_out + raw > (uint64_t)d->content_size) return LZR_ERR_CORRUPT;

    /* Make room: keep exactly the last wsize bytes of history. */
    if (d->op + raw > d->hcap) {
        size_t delta = d->op - d->wsize;
        memmove(d->hist, d->hist + delta, d->wsize);
        d->op = d->wsize;
    }
    uint8_t *dst = d->hist + d->op;
    if (type == 0) {
        memcpy(dst, payload, raw);
    } else {
        int r = decode_block(d, payload, comp, dst, raw);
        if (r) return r;
    }
    d->crc = lzr_crc32(d->crc, dst, raw);
    d->op += raw;
    d->total_out += raw;
    d->ipos += LZR_BLOCK_HDR + comp;
    if (buf_append(out, dst, raw)) return LZR_ERR_NOMEM;
    return 1;
}

int lzr_dctx_write(lzr_dctx *d, const uint8_t *src, size_t n, lzr_buf *out, size_t *consumed)
{
    if (consumed) *consumed = 0;
    if (d->state == ST_DONE) return 1;
    if (n && buf_append(&d->in, src, n)) return LZR_ERR_NOMEM;

    int r;
    while ((r = step(d, out)) == 1) {}
    if (r < 0) return r;

    size_t left = d->in.len - d->ipos;
    if (r == 2) {                 /* bytes after the end block are not ours */
        if (consumed) *consumed = n - left;
        d->in.len = d->ipos = 0;
        return 1;
    }
    if (consumed) *consumed = n;
    if (d->ipos && (left == 0 || d->ipos >= (1u << 20))) {
        memmove(d->in.p, d->in.p + d->ipos, left);
        d->in.len = left;
        d->ipos = 0;
    }
    return LZR_OK;
}
