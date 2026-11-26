/*
 * dec.c - streaming decoder.
 *
 * Input is buffered until a whole block is available (its header gives
 * comp_len). Block type 5: the Huffman literals are decoded into an array
 * first; then one loop decodes each sequence from the tANS bit stream
 * (a few sequences ahead, prefetching the match source) and copies its
 * literal run and match with fixed-size chunks. Block type 4: one range
 * decoder loop. Output goes to a history buffer that keeps the last
 * 2^window_log bytes, or, in one-buffer mode, straight into the output.
 */
#include "lzr.h"
#include "codes.h"
#include "xxh64.h"
#include "arc.h"
#include "x86.h"
#include "huf.h"
#include "fse.h"

#define MAX_SEQ (LZR_BLOCK_SIZE / LZR_ML_BASE + 1)
#define MAX_COMP (2 * LZR_BLOCK_SIZE + 65536)   /* sanity cap on comp_len */

enum { ST_HEADER, ST_BLOCKS, ST_DONE };

typedef struct {
    uint32_t base;   /* value without its extra bits                 */
    uint16_t next;   /* tANS: next state = next + (nb bits)          */
    uint8_t nb, xb;  /* tANS bits, extra bits                        */
} seq_dent;

struct lzr_dctx {
    int state;
    lzr_buf in;
    size_t ipos;
    uint8_t *hist;     /* history: own buffer, or the output (direct mode) */
    lzr_buf *direct;   /* lzr_dctx_set_output(): decode into this buffer   */
    size_t hcap, wsize, op;
    uint32_t reps[LZR_NUM_REPS];
    xxh64_state xs;    /* checksum of the decoded content */
    uint64_t total_out;
    int64_t content_size;
    uint8_t *lits;     /* type 5 literals of the current block */
    code_info ct[VAL_CODES];
    tok_model tm;      /* block type 4 model, persists across blocks */
    huf_table huf;     /* block type 5 tables, kept for the repeat modes */
    fse_table fse[3];
    seq_dent sq[3][FSE_L];   /* per-block sequence decoding tables */
    lzr_x86 x86;
    uint8_t *x86tmp;   /* x86 filter work buffer, NULL = no filter */
};

lzr_dctx *lzr_dctx_new(void)
{
    lzr_dctx *d = (lzr_dctx *)calloc(1, sizeof(*d));
    if (d) xxh64_init(&d->xs);
    if (!d) return NULL;
    d->reps[0] = 1; d->reps[1] = 4; d->reps[2] = 8;
    d->content_size = -1;
    d->lits = (uint8_t *)malloc(LZR_BLOCK_SIZE + LZR_SLACK);
    if (!d->lits) { lzr_dctx_free(d); return NULL; }
    code_table_init(d->ct);
    tok_model_init(&d->tm);
    return d;
}

void lzr_dctx_free(lzr_dctx *d)
{
    if (!d) return;
    free(d->in.p);
    if (!d->direct) free(d->hist);
    free(d->x86tmp);
    free(d->lits);
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
 * Fixed-size chunks only (no library calls); may write up to 15 bytes
 * past op + len (buffers have LZR_SLACK). */
static inline void copy_match(uint8_t *op, size_t dist, size_t len)
{
    const uint8_t *src = op - dist;
    if (dist >= 16) {
        for (size_t k = 0; k < len; k += 16) memcpy(op + k, src + k, 16);
    } else if (dist >= 8) {
        for (size_t k = 0; k < len; k += 8) memcpy(op + k, src + k, 8);
    } else {
        /* one period byte by byte, then 8-byte chunks from a whole number
         * of periods back (8 <= period * m <= 14) */
        for (size_t k = 0; k < 8; k++) op[k] = src[k];
        size_t back = dist * ((8 + dist - 1) / dist);
        for (size_t k = 8; k < len; k += 8) memcpy(op + k, op + k - back, 8);
    }
}

/* Literal run: 16-byte chunks (source and destination have slack). */
static inline void copy_lits(uint8_t *o, const uint8_t *lp, size_t ll)
{
    memcpy(o, lp, 16);
    for (size_t k = 16; k < ll; k += 16) memcpy(o + k, lp + k, 16);
}

/* Block type 4: one range-coded token stream (arc.h), decoded straight
 * into the history. The coder state lives in locals (range, code, ip) so
 * it stays in registers; the macros mirror arce_bit() and the trees in
 * arc.h. */
#define RC_NORM() do { if (range < ARC_TOP) { range <<= 8; code = (code << 8) | *ip++; } } while (0)
#define RC_BIT(pr, b) do {                                                     \
        arc_p *pp_ = (pr);                                                     \
        uint32_t p_ = *pp_, bound_ = (range >> ARC_PBITS) * p_;                \
        if (code < bound_) { range = bound_; *pp_ = (arc_p)(p_ + ((ARC_ONE - p_) >> ARC_SHIFT)); b = 0; } \
        else { range -= bound_; code -= bound_; *pp_ = (arc_p)(p_ - (p_ >> ARC_SHIFT)); b = 1; } \
        RC_NORM();                                                             \
    } while (0)
#define RC_TREE(probs, nbits, out) do {                                        \
        unsigned m_ = 1, b_;                                                   \
        for (unsigned i_ = 0; i_ < (nbits); i_++) { RC_BIT(&(probs)[m_], b_); m_ = (m_ << 1) | b_; } \
        out = m_ - (1u << (nbits));                                            \
    } while (0)
#define RC_RTREE(probs, nbits, out) do {                                       \
        unsigned m_ = 1, b_, v_ = 0;                                           \
        for (unsigned i_ = 0; i_ < (nbits); i_++) { RC_BIT(&(probs)[m_], b_); m_ = (m_ << 1) | b_; v_ |= b_ << i_; } \
        out = v_;                                                              \
    } while (0)
#define RC_LEN(l, ps, out) do {                                                \
        unsigned b0_, v0_;                                                     \
        RC_BIT(&(l)->choice, b0_);                                             \
        if (!b0_) { RC_TREE((l)->low[ps], 3, v0_); out = 2 + v0_; break; }   \
        RC_BIT(&(l)->choice2, b0_);                                            \
        if (!b0_) { RC_TREE((l)->mid[ps], 3, v0_); out = 10 + v0_; break; }  \
        RC_TREE((l)->high, 8, v0_); out = 18 + v0_;                            \
    } while (0)

static int decode_tokens(lzr_dctx *d, const uint8_t *in, size_t avail, uint8_t *dst, size_t raw)
{
    arc_dec r;
    if (arcd_init(&r, in, avail)) return LZR_ERR_CORRUPT;
    uint32_t range = r.range, code = r.code;
    const uint8_t *ip = r.p, *const iend = r.end;
    tok_model *m = &d->tm;
    unsigned s = m->state;
    uint32_t rep0 = d->reps[0], rep1 = d->reps[1], rep2 = d->reps[2];
    const uint8_t *hist = d->hist;
    const size_t wsize = d->wsize;
    uint8_t *o = dst, *const oend = dst + raw;
    int err = 0;
    while (o < oend) {
        unsigned ps = (unsigned)(o - dst) & (TOK_POS - 1), bit;
        if (ip > iend) { err = 1; break; }
        RC_BIT(&m->is_match[s][ps], bit);
        if (!bit) {
            arc_p *probs = m->lit[(o > hist ? o[-1] : 0) >> 5];
            unsigned v = 1;
            if (s >= 7) {
                if (rep0 > (size_t)(o - hist)) { err = 1; break; }
                unsigned mbyte = o[-(ptrdiff_t)rep0];
                for (int i = 7; i >= 0; i--) {
                    unsigned mb = (mbyte >> i) & 1;
                    RC_BIT(&probs[256 + (mb << 8) + v], bit);
                    v = (v << 1) | bit;
                    if (bit != mb) break;
                }
            }
            while (v < 256) { RC_BIT(&probs[v], bit); v = (v << 1) | bit; }
            *o++ = (uint8_t)v;
            s = tok_after_lit(s);
            continue;
        }
        size_t len;
        uint32_t dist;
        RC_BIT(&m->is_rep[s], bit);
        if (!bit) {
            unsigned l, slot;
            RC_LEN(&m->mlen, ps, l);
            len = l;
            unsigned lctx = l - 2 < 3 ? l - 2 : 3;
            RC_TREE(m->slot[lctx], 6, slot);
            uint32_t d0 = slot;
            if (slot >= 4) {
                unsigned nb = (slot >> 1) - 1, v;
                d0 = (2u | (slot & 1)) << nb;
                if (slot < TOK_END_SPEC) {
                    RC_RTREE(m->spec[slot], nb, v);
                    d0 += v;
                } else {
                    uint32_t hi = 0;
                    for (unsigned k = nb - 4; k; k--) {
                        range >>= 1;
                        uint32_t t = code >= range;
                        code -= range & (0u - t);
                        hi = (hi << 1) | t;
                        RC_NORM();
                    }
                    RC_RTREE(m->align, 4, v);
                    d0 += (hi << 4) + v;
                }
            }
            if (d0 >= wsize) { err = 1; break; }
            dist = d0 + 1;
            rep2 = rep1; rep1 = rep0; rep0 = dist;
            s = tok_after_match(s);
        } else {
            RC_BIT(&m->is_rep0[s], bit);
            if (!bit) {
                dist = rep0;
                RC_BIT(&m->is_rep0_long[s][ps], bit);
                if (!bit) {   /* short rep */
                    if (dist > (size_t)(o - hist)) { err = 1; break; }
                    *o = o[-(ptrdiff_t)dist];
                    o++;
                    s = tok_after_short(s);
                    continue;
                }
            } else {
                RC_BIT(&m->is_rep1[s], bit);
                if (!bit) { dist = rep1; rep1 = rep0; }
                else { dist = rep2; rep2 = rep1; rep1 = rep0; }
                rep0 = dist;
            }
            unsigned l;
            RC_LEN(&m->rlen, ps, l);
            len = l;
            s = tok_after_rep(s);
        }
        if (len > (size_t)(oend - o) || dist > (size_t)(o - hist) || dist > wsize) { err = 1; break; }
        copy_match(o, dist, len);
        o += len;
    }
    m->state = (uint8_t)s;
    d->reps[0] = rep0; d->reps[1] = rep1; d->reps[2] = rep2;
    if (err || ip != iend) return LZR_ERR_CORRUPT;
    return LZR_OK;
}

/* Block type 5 (see emit_seqs5 in enc.c): Huffman literals into d->lits,
 * then one loop that decodes each sequence from the tANS bit stream and
 * executes it at once. */
static int decode_block5(lzr_dctx *d, const uint8_t *in, size_t avail, uint8_t *dst, size_t raw)
{
    if (avail < 8) return LZR_ERR_CORRUPT;
    size_t ns = rd32le(in), nlit = rd32le(in + 4), k = 8;
    if (nlit > raw || ns > MAX_SEQ || ns * LZR_ML_BASE + nlit > raw) return LZR_ERR_CORRUPT;
    if (nlit) {
        size_t c = huf_decode(in + k, avail - k, d->lits, nlit, &d->huf);
        if (!c) return LZR_ERR_CORRUPT;
        k += c;
    }
    const uint8_t *lp = d->lits, *const lend = d->lits + nlit;
    uint8_t *o = dst, *const oend = dst + raw;
    if (!ns) {
        if (k != avail || nlit != raw) return LZR_ERR_CORRUPT;
        memcpy(o, lp, nlit);
        return LZR_OK;
    }
    const unsigned alpha[3] = {VAL_CODES, VAL_CODES, OF_CODES};
    unsigned rsym[3] = {0, 0, 0};
    int rle[3];
    for (int t = 0; t < 3; t++) {
        if (k >= avail) return LZR_ERR_CORRUPT;
        unsigned mode = in[k++];
        fse_table *T = &d->fse[t];
        rle[t] = mode == 2;
        if (mode == 2) {
            if (k >= avail || in[k] >= alpha[t]) return LZR_ERR_CORRUPT;
            rsym[t] = in[k++];
        } else if (mode == 1) {
            if (!T->valid || T->alphabet != alpha[t]) return LZR_ERR_CORRUPT;
        } else if (mode == 0) {
            size_t c = fse_read_norm(in + k, avail - k, T->norm, alpha[t]);
            T->valid = 0;
            T->alphabet = alpha[t];
            if (!c || fse_build(T)) return LZR_ERR_CORRUPT;
            k += c;
        } else {
            return LZR_ERR_CORRUPT;
        }
    }
    /* Sequence tables: the tANS state step plus the code's base value and
     * extra-bit count in one entry. Offsets: value < 3 is a repeat index,
     * else distance + 3. A single-symbol stream is a one-state table. */
    seq_dent *tb[3] = {d->sq[0], d->sq[1], d->sq[2]};
    unsigned tlog[3];
    for (int t = 0; t < 3; t++) {
        const fse_table *T = &d->fse[t];
        unsigned n = rle[t] ? 1 : FSE_L;
        tlog[t] = rle[t] ? 0 : FSE_LOG;
        for (unsigned u = 0; u < n; u++) {
            unsigned sym = rle[t] ? rsym[t] : T->dt[u].sym;
            seq_dent *e = &tb[t][u];
            e->nb = rle[t] ? 0 : T->dt[u].nb;
            e->next = rle[t] ? 0 : T->dt[u].base;
            if (t < 2) { e->base = d->ct[sym].base + (t == 1 ? LZR_ML_BASE : 0); e->xb = d->ct[sym].bits; }
            else if (sym < LZR_NUM_REPS) { e->base = sym; e->xb = 0; }
            else { e->base = d->ct[sym - LZR_NUM_REPS].base + 1 + 3; e->xb = d->ct[sym - LZR_NUM_REPS].bits; }
        }
    }

    /* backward bit reader in locals (see fse.h) */
    const uint8_t *const bstart = in + k;
    const size_t bsize = avail - k;
    if (!bsize || !bstart[bsize - 1]) return LZR_ERR_CORRUPT;
    const uint8_t *bp = bstart + bsize - 8;
    uint64_t bc = fse_ld64(bp);
    unsigned used = (unsigned)__builtin_clzll(bc) + 1;
    int berr = 0;
#define BGET(nb) ((uint32_t)(((bc << (used & 63)) >> 1) >> (63 - (nb))))
#define BRELOAD()                                                              \
    do {                                                                       \
        if (bp >= bstart + 8) { bp -= used >> 3; used &= 7; }                  \
        else {                                                                 \
            size_t nby_ = used >> 3, room_ = bp > bstart ? (size_t)(bp - bstart) : 0; \
            if (nby_ > room_) nby_ = room_;                                    \
            bp -= nby_; used -= (unsigned)nby_ * 8;                            \
            if (used > 64) { berr = 1; used = 64; }                            \
        }                                                                      \
        bc = fse_ld64(bp);                                                     \
    } while (0)

    unsigned sll = BGET(tlog[0]); used += tlog[0];
    unsigned sml = BGET(tlog[1]); used += tlog[1];
    unsigned sof = BGET(tlog[2]); used += tlog[2];
    BRELOAD();
    const seq_dent *dll = tb[0], *dml = tb[1], *dof = tb[2];
    uint32_t rep0 = d->reps[0], rep1 = d->reps[1], rep2 = d->reps[2];
    const uint8_t *hist = d->hist;
    const size_t wsize = d->wsize;

    /* Decode a sequence (ll, ml, dist), advancing the bit stream. */
#define DECODE_SEQ(LL, ML, DIST)                                               \
    do {                                                                       \
        const seq_dent el_ = dll[sll], em_ = dml[sml], eo_ = dof[sof];         \
        uint32_t ov_ = eo_.base + BGET(eo_.xb);                                \
        used += eo_.xb;                                                        \
        if (ov_ >= LZR_NUM_REPS) { DIST = ov_ - 3; rep2 = rep1; rep1 = rep0; rep0 = DIST; } \
        else if (ov_ == 0) { DIST = rep0; }                                    \
        else if (ov_ == 1) { DIST = rep1; rep1 = rep0; rep0 = DIST; }          \
        else { DIST = rep2; rep2 = rep1; rep1 = rep0; rep0 = DIST; }           \
        BRELOAD();                                                             \
        ML = em_.base + BGET(em_.xb); used += em_.xb;                          \
        LL = el_.base + BGET(el_.xb); used += el_.xb;                          \
        BRELOAD();                                                             \
        sll = el_.next + BGET(el_.nb); used += el_.nb;                         \
        sml = em_.next + BGET(em_.nb); used += em_.nb;                         \
        sof = eo_.next + BGET(eo_.nb); used += eo_.nb;                         \
        BRELOAD();                                                             \
    } while (0)
    /* Execute it: literal run, then the match (checked here). */
#define EXEC_SEQ(LL, ML, DIST)                                                 \
    do {                                                                       \
        if ((LL) > (size_t)(lend - lp) || (LL) + (ML) > (size_t)(oend - o)) { berr = 1; goto done; } \
        copy_lits(o, lp, LL);                                                  \
        o += LL;                                                               \
        lp += LL;                                                              \
        if ((DIST) > (size_t)(o - hist) || (DIST) > wsize) { berr = 1; goto done; } \
        copy_match(o, DIST, ML);                                               \
        o += ML;                                                               \
    } while (0)

    /* Sequences are decoded SEQ_AHEAD ahead of execution, and the source of
     * each match is prefetched meanwhile: far matches miss the cache. */
    enum { SEQ_AHEAD = 8 };
    size_t rll[SEQ_AHEAD], rml[SEQ_AHEAD];
    uint32_t rdist[SEQ_AHEAD];
    size_t ahead = (size_t)(o - hist);   /* history position after the decoded sequences */
    size_t i = 0, lead = ns < SEQ_AHEAD ? ns : SEQ_AHEAD;
    for (; i < lead; i++) {
        size_t ll, ml;
        uint32_t dist;
        DECODE_SEQ(ll, ml, dist);
        ahead += ll;
        if (dist <= ahead) __builtin_prefetch(hist + ahead - dist);
        ahead += ml;
        rll[i] = ll; rml[i] = ml; rdist[i] = dist;
    }
    for (; i < ns; i++) {
        size_t ll, ml;
        uint32_t dist;
        DECODE_SEQ(ll, ml, dist);
        ahead += ll;
        if (dist <= ahead) __builtin_prefetch(hist + ahead - dist);
        ahead += ml;
        const size_t j = i % SEQ_AHEAD;
        EXEC_SEQ(rll[j], rml[j], rdist[j]);
        rll[j] = ll; rml[j] = ml; rdist[j] = dist;
    }
    for (size_t e = ns - lead; e < ns; e++) {
        const size_t j = e % SEQ_AHEAD;
        EXEC_SEQ(rll[j], rml[j], rdist[j]);
    }
done:
#undef DECODE_SEQ
#undef EXEC_SEQ
#undef BGET
#undef BRELOAD
    d->reps[0] = rep0; d->reps[1] = rep1; d->reps[2] = rep2;
    /* all bits consumed, exactly */
    if (berr || (size_t)(bstart + bsize - (bp + 8)) * 8 + used != bsize * 8) return LZR_ERR_CORRUPT;
    size_t tail = (size_t)(lend - lp);
    if (tail != (size_t)(oend - o)) return LZR_ERR_CORRUPT;
    memcpy(o, lp, tail);
    return LZR_OK;
}

static int init_history(lzr_dctx *d, int wlog)
{
    d->wsize = (size_t)1 << wlog;
    d->hcap = 2 * d->wsize + LZR_BLOCK_SIZE;
    d->hist = (uint8_t *)big_alloc(d->hcap + LZR_SLACK, 0);
    return d->hist ? LZR_OK : LZR_ERR_NOMEM;
}

typedef struct { lzr_dctx *d; lzr_buf *out; } out_ctx;

/* Deliver decoded (unfiltered) content; the CRC covers the original. */
static int emit_out(void *ctx, const uint8_t *p, size_t n)
{
    out_ctx *s = (out_ctx *)ctx;
    xxh64_update(&s->d->xs, p, n);
    return buf_append(s->out, p, n);
}

/* Process one unit at in.p + ipos. Returns 1 if consumed, 0 if more input
 * is needed, 2 at the verified end, or LZR_ERR_*. */
static int step(lzr_dctx *d, lzr_buf *out)
{
    const uint8_t *p = d->in.p + d->ipos;
    size_t avail = d->in.len - d->ipos;

    if (d->state == ST_HEADER) {
        int wlog;
        uint8_t flags;
        int r = lzr_frame_header(p, avail, &wlog, &flags, &d->content_size);
        if (r <= 0) return r;
        if (flags & LZR_FLAG_X86) d->direct = NULL;   /* history must hold filtered bytes */
        if (d->direct) d->wsize = (size_t)1 << wlog;
        else if (init_history(d, wlog)) return LZR_ERR_NOMEM;
        if ((flags & LZR_FLAG_X86) && !(d->x86tmp = (uint8_t *)malloc(X86_CHUNK + 4))) return LZR_ERR_NOMEM;
        d->ipos += (size_t)r;
        d->state = ST_BLOCKS;
        return 1;
    }

    if (avail < LZR_BLOCK_HDR) return 0;
    uint8_t type = p[0];
    size_t raw = rd32le(p + 1), comp = rd32le(p + 5);
    if (type > 5 || type == 1 || type == 3 || raw > LZR_BLOCK_SIZE || comp > MAX_COMP) return LZR_ERR_CORRUPT;
    if (type == 0 && comp != raw) return LZR_ERR_CORRUPT;
    if (type == 2 && (raw != 0 || comp != 4)) return LZR_ERR_CORRUPT;
    if (type >= 4 && raw == 0) return LZR_ERR_CORRUPT;
    if (avail - LZR_BLOCK_HDR < comp) return 0;
    const uint8_t *payload = p + LZR_BLOCK_HDR;

    if (type == 2) {
        if (d->content_size >= 0 && (uint64_t)d->content_size != d->total_out) return LZR_ERR_CORRUPT;
        if (d->x86tmp) {
            out_ctx s = {d, out};
            if (x86_finish(&d->x86, emit_out, &s)) return LZR_ERR_NOMEM;
        }
        d->ipos += LZR_BLOCK_HDR + comp;
        d->state = ST_DONE;
        return rd32le(payload) == (uint32_t)xxh64_digest(&d->xs) ? 2 : LZR_ERR_CHECKSUM;
    }
    if (d->content_size >= 0 && d->total_out + raw > (uint64_t)d->content_size) return LZR_ERR_CORRUPT;

    if (d->direct) {   /* the output itself is the history: no copy */
        if (out != d->direct) return LZR_ERR_PARAM;
        if (buf_reserve(out, raw + LZR_SLACK)) return LZR_ERR_NOMEM;
        d->hist = out->p;
        d->op = out->len;
    } else if (d->op + raw > d->hcap) {   /* make room: keep the last wsize bytes */
        size_t delta = d->op - d->wsize;
        memmove(d->hist, d->hist + delta, d->wsize);
        d->op = d->wsize;
    }
    uint8_t *dst = d->hist + d->op;
    if (type == 0) {
        memcpy(dst, payload, raw);
    } else {
        int r = type == 4 ? decode_tokens(d, payload, comp, dst, raw) : decode_block5(d, payload, comp, dst, raw);
        if (r) return r;
    }
    d->op += raw;
    d->total_out += raw;
    d->ipos += LZR_BLOCK_HDR + comp;
    if (d->direct) {
        xxh64_update(&d->xs, dst, raw);
        out->len += raw;
        return 1;
    }
    out_ctx s = {d, out};
    /* history keeps the filtered bytes (matches refer to them) */
    int r = d->x86tmp ? x86_run(&d->x86, dst, raw, 1, d->x86tmp, emit_out, &s) : emit_out(&s, dst, raw);
    return r ? LZR_ERR_NOMEM : 1;
}

int lzr_dctx_set_output(lzr_dctx *d, lzr_buf *out)
{
    if (d->state != ST_HEADER || d->in.len) return LZR_ERR_PARAM;
    d->direct = out;
    return LZR_OK;
}

int lzr_dctx_write(lzr_dctx *d, const uint8_t *src, size_t n, lzr_buf *out, size_t *consumed)
{
    if (consumed) *consumed = 0;
    if (d->state == ST_DONE) return 1;
    /* Keep ARC_SLACK readable bytes after the input (see arc.h). */
    if (n && (buf_reserve(&d->in, n + ARC_SLACK) || buf_append(&d->in, src, n))) return LZR_ERR_NOMEM;
    if (n) memset(d->in.p + d->in.len, 0, ARC_SLACK);

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
