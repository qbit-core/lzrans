/*
 * enc.c - streaming encoder.
 *
 * Structure:
 *   lzr_cctx          sliding input buffer (2*window + block), hash-chain
 *                     tables, repeat offsets, running CRC
 *   find_best()       repeat-offset probes + hash-chain walk, scored by gain
 *   parse_block()     lazy parse (depth 0..2) into sequences + literals
 *   emit_block()      sequences -> code streams + extra bits -> rANS
 *   slide()           drops history older than the window, rebases tables
 */
#include "lzr.h"
#include "rans.h"
#include "codes.h"
#include "crc32.h"

#include "enc_internal.h"

/* Levels 1-6: hash chains + lazy parsing (window, hash bits, lazy depth,
 * chain steps, nice length). Levels 7-9: binary tree + optimal parsing
 * (chain = tree depth, opt = number of parsing passes). */
static const level_params LEVELS[LZR_LEVEL_MAX + 1] = {
    {0, 0, 0, 0, 0, 0},
    {20, 16, 0, 4, 16, 0},    {21, 17, 1, 8, 24, 0},    {22, 18, 1, 16, 32, 0},
    {22, 19, 1, 24, 48, 0},   {23, 20, 2, 32, 64, 0},   {23, 20, 2, 64, 96, 0},
    {23, 20, 0, 16, 64, 1},   {24, 21, 0, 32, 128, 2},  {24, 22, 0, 48, 273, 2},
};


lzr_cctx *lzr_cctx_new(int level, int window_log, int64_t content_size)
{
    if (level < LZR_LEVEL_MIN || level > LZR_LEVEL_MAX) return NULL;
    lzr_cctx *c = (lzr_cctx *)calloc(1, sizeof(*c));
    if (!c) return NULL;
    c->lp = LEVELS[level];
    if (window_log) {
        if (window_log < LZR_WLOG_MIN || window_log > LZR_WLOG_MAX) { free(c); return NULL; }
        c->lp.wlog = (uint8_t)window_log;
    }
#ifdef LZR_TUNE   /* developer builds only: override level parameters */
    if (getenv("LZR_DEPTH")) c->lp.chain = (uint16_t)atoi(getenv("LZR_DEPTH"));
    if (getenv("LZR_NICE")) c->lp.nice = (uint16_t)atoi(getenv("LZR_NICE"));
    if (getenv("LZR_PASSES")) c->lp.opt = (uint8_t)atoi(getenv("LZR_PASSES"));
    if (getenv("LZR_WLOG")) c->lp.wlog = (uint8_t)atoi(getenv("LZR_WLOG"));
#endif
    /* A window larger than the known input only wastes memory. */
    while (content_size >= 0 && c->lp.wlog > LZR_WLOG_MIN && ((uint64_t)1 << (c->lp.wlog - 1)) >= (uint64_t)content_size)
        c->lp.wlog--;
    if (c->lp.hlog > c->lp.wlog) c->lp.hlog = c->lp.wlog;
    c->wsize = (size_t)1 << c->lp.wlog;
    c->wmask = c->wsize - 1;
    c->cap = 2 * c->wsize + LZR_BLOCK_SIZE;
    c->content_size = content_size;
    c->reps[0] = 1; c->reps[1] = 4; c->reps[2] = 8;

    c->buf = (uint8_t *)malloc(c->cap + LZR_SLACK);
    c->head = (uint32_t *)calloc((size_t)1 << c->lp.hlog, sizeof(uint32_t));
    if (!c->lp.opt) c->prev = (uint32_t *)calloc(c->wsize, sizeof(uint32_t));   /* chains: lazy levels only */
    c->seq_ll = (uint32_t *)malloc(MAX_SEQ * sizeof(uint32_t));
    c->seq_ml = (uint32_t *)malloc(MAX_SEQ * sizeof(uint32_t));
    c->seq_dist = (uint32_t *)malloc(MAX_SEQ * sizeof(uint32_t));
    c->lits = (uint8_t *)malloc(LZR_BLOCK_SIZE);
    c->llc = (uint8_t *)malloc(MAX_SEQ);
    c->mlc = (uint8_t *)malloc(MAX_SEQ);
    c->ofc = (uint8_t *)malloc(MAX_SEQ);
    c->extra = (uint8_t *)malloc(MAX_SEQ * 8 + 16);
    c->scratch = (uint8_t *)malloc(2 * LZR_BLOCK_SIZE + 16);
    c->blk = (uint8_t *)malloc(64 + stream_bound(LZR_BLOCK_SIZE, ALPHA_LIT) + 3 * stream_bound(MAX_SEQ, ALPHA_CODE) + MAX_SEQ * 8 + 16);
    if (c->lp.opt && opt_init(c)) { lzr_cctx_free(c); return NULL; }
    if (!c->buf || !c->head || (!c->prev && !c->lp.opt) || !c->seq_ll || !c->seq_ml || !c->seq_dist || !c->lits ||
        !c->llc || !c->mlc || !c->ofc || !c->extra || !c->scratch || !c->blk) {
        lzr_cctx_free(c);
        return NULL;
    }
    return c;
}

void lzr_cctx_free(lzr_cctx *c)
{
    if (!c) return;
    opt_free(c);
    free(c->buf); free(c->head); free(c->prev);
    free(c->seq_ll); free(c->seq_ml); free(c->seq_dist);
    free(c->lits); free(c->llc); free(c->mlc); free(c->ofc);
    free(c->extra); free(c->scratch); free(c->blk);
    free(c);
}

/* ---------------- match finding ---------------- */

/* Insert positions [c->ins, upto) into the hash chains. */
static inline void insert_upto(lzr_cctx *c, size_t upto)
{
    const unsigned hlog = c->lp.hlog;
    size_t limit = c->end >= LZR_HASH_LEN ? c->end - LZR_HASH_LEN + 1 : 0;  /* only hash real bytes */
    if (upto > limit) upto = limit;
    for (size_t i = c->ins; i < upto; i++) {
        uint32_t h = hash4(c->buf + i, hlog);
        c->prev[i & c->wmask] = c->head[h];
        c->head[h] = (uint32_t)(i + 1);
    }
    if (upto > c->ins) c->ins = upto;
}

typedef struct { uint32_t len, off; int gain; } match_t;  /* off: 0..2 rep idx, else dist + 2 */

/* Score in quarter-bits: each matched byte is worth ~4, the offset costs its log2. */
static inline int gain_of(uint32_t len, uint32_t off) { return (int)len * 4 - (int)highbit32(off + 1); }

static match_t find_best(lzr_cctx *c, size_t p, size_t maxlen, const uint32_t *reps)
{
    match_t best = {0, 0, 0};
    const uint8_t *cur = c->buf + p;
    const uint32_t first = rd32(cur);

    /* 1) repeat offsets: cheapest to code, try them first */
    for (uint32_t r = 0; r < LZR_NUM_REPS; r++) {
        uint32_t d = reps[r];
        if (d > p || d > c->wsize) continue;
        if (rd32(cur - d) != first) continue;
        uint32_t len = (uint32_t)match_len(cur - d, cur, maxlen);
        int g = gain_of(len, r);
        if (len >= LZR_MIN_MATCH && g > best.gain) { best.len = len; best.off = r; best.gain = g; }
    }

    if (best.len == maxlen || best.len >= c->lp.nice) return best;

    /* 2) hash chain: longest match within the window */
    uint32_t cand = c->head[hash4(cur, c->lp.hlog)];
    unsigned chain = c->lp.chain;
    size_t best_len = best.len >= LZR_MIN_MATCH ? best.len : LZR_MIN_MATCH - 1;
    while (cand && chain--) {
        size_t cp = cand - 1;
        size_t dist = p - cp;
        if (dist > c->wsize) break;
        const uint8_t *ref = c->buf + cp;
        if (ref[best_len] == cur[best_len] && rd32(ref) == first) {
            size_t len = match_len(ref, cur, maxlen);
            if (len > best_len) {
                int g = gain_of((uint32_t)len, (uint32_t)dist + 2);
                if (g > best.gain) {
                    best_len = len;
                    best.len = (uint32_t)len; best.off = (uint32_t)dist + 2; best.gain = g;
                    if (len == maxlen || len >= c->lp.nice) break;
                }
            }
        }
        uint32_t next = c->prev[cp & c->wmask];
        if (next == 0 || next - 1 >= cp) break;   /* recycled slot / end of chain */
        cand = next;
    }
    return best;
}

/* Lazy parse of buf[bs, be). Returns number of sequences; *nlit literals. */
static size_t parse_block(lzr_cctx *c, size_t bs, size_t be, size_t *nlit_out)
{
    size_t p = bs, lit_start = bs, ns = 0, nlit = 0;
    uint32_t *reps = c->reps;

    while (p + LZR_MIN_MATCH <= be) {
        insert_upto(c, p);
        match_t m = find_best(c, p, be - p, reps);
        if (!m.len) { p++; continue; }

        /* Lazy evaluation: would starting 1 byte later be clearly better? */
        for (unsigned depth = 1; depth <= c->lp.lazy && m.len < c->lp.nice && p + 1 + LZR_MIN_MATCH <= be; depth++) {
            insert_upto(c, p + 1);
            match_t m2 = find_best(c, p + 1, be - p - 1, reps);
            int threshold = m.gain + (depth == 1 ? 4 : 7);
            if (m2.len && m2.gain > threshold) { m = m2; p++; }
            else break;
        }

        /* Emit: literal run, then the match. */
        size_t ll = p - lit_start;
        memcpy(c->lits + nlit, c->buf + lit_start, ll);
        nlit += ll;
        c->seq_ll[ns] = (uint32_t)ll;
        c->seq_ml[ns] = m.len;
        c->seq_dist[ns] = commit_offset(reps, m.off);
        ns++;
        p += m.len;
        lit_start = p;
    }
    insert_upto(c, be);
    memcpy(c->lits + nlit, c->buf + lit_start, be - lit_start);
    nlit += be - lit_start;
    *nlit_out = nlit;
    return ns;
}

/* ---------------- block emission ---------------- */

static int write_header(lzr_cctx *c, lzr_buf *out)
{
    uint8_t h[16];
    memcpy(h, LZR_MAGIC, 4);
    h[4] = LZR_VERSION;
    h[5] = c->content_size >= 0 ? LZR_FLAG_SIZE : 0;
    h[6] = c->lp.wlog;
    h[7] = 0;
    size_t n = 8;
    if (c->content_size >= 0) { wr64le(h + 8, (uint64_t)c->content_size); n = 16; }
    c->header_done = 1;
    return buf_append(out, h, n);
}

static int emit_block(lzr_cctx *c, size_t bs, size_t be, lzr_buf *out)
{
    const size_t raw = be - bs;
    uint32_t saved_reps[LZR_NUM_REPS];
    memcpy(saved_reps, c->reps, sizeof(saved_reps));

    size_t nlit;
    size_t ns = c->lp.opt ? opt_parse_block(c, bs, be, &nlit) : parse_block(c, bs, be, &nlit);
    if (ns == (size_t)-1) return LZR_ERR_NOMEM;

    /* Sequences -> code symbols + raw extra bits (LL, ML, OF order). */
    bitw bw = {0, 0, c->extra};
    for (size_t i = 0; i < ns; i++) {
        unsigned nb;
        uint32_t ex;
        c->llc[i] = (uint8_t)val_code(c->seq_ll[i], &nb, &ex);
        bw_put(&bw, ex, nb);
        c->mlc[i] = (uint8_t)val_code(c->seq_ml[i] - LZR_MIN_MATCH, &nb, &ex);
        bw_put(&bw, ex, nb);
        uint32_t off = c->seq_dist[i];
        if (off < LZR_NUM_REPS) {
            c->ofc[i] = (uint8_t)off;
        } else {
            c->ofc[i] = (uint8_t)(LZR_NUM_REPS + val_code(off - 3, &nb, &ex));   /* dist - 1 */
            bw_put(&bw, ex, nb);
        }
    }
    bw_flush(&bw);
    size_t extra_len = (size_t)(bw.p - c->extra);

    uint8_t *b = c->blk;
    size_t k = LZR_BLOCK_HDR;
    wr32le(b + k, (uint32_t)ns); k += 4;
    wr32le(b + k, (uint32_t)nlit); k += 4;
    k += stream_encode(b + k, c->lits, nlit, ALPHA_LIT, c->scratch);
    k += stream_encode(b + k, c->llc, ns, ALPHA_CODE, c->scratch);
    k += stream_encode(b + k, c->mlc, ns, ALPHA_CODE, c->scratch);
    k += stream_encode(b + k, c->ofc, ns, ALPHA_CODE, c->scratch);
    wr32le(b + k, (uint32_t)extra_len); k += 4;
    memcpy(b + k, c->extra, extra_len); k += extra_len;

    size_t payload = k - LZR_BLOCK_HDR;
    if (payload < raw) {
        b[0] = 1;
        wr32le(b + 1, (uint32_t)raw);
        wr32le(b + 5, (uint32_t)payload);
        return buf_append(out, b, k);
    }
    /* Stored block; the decoder leaves its repeat offsets untouched. */
    memcpy(c->reps, saved_reps, sizeof(saved_reps));
    uint8_t h[LZR_BLOCK_HDR] = {0};
    wr32le(h + 1, (uint32_t)raw);
    wr32le(h + 5, (uint32_t)raw);
    if (buf_append(out, h, sizeof(h))) return LZR_ERR_NOMEM;
    return buf_append(out, c->buf + bs, raw);
}

/* Keep the last >= wsize bytes; shift by a multiple of wsize so the
 * prev[] ring indices (pos & wmask) stay valid. */
static void slide(lzr_cctx *c)
{
    if (c->pos <= c->wsize) return;
    size_t delta = ((c->pos - c->wsize) / c->wsize) * c->wsize;
    if (!delta) return;
    memmove(c->buf, c->buf + delta, c->end - delta);
    c->pos -= delta;
    c->end -= delta;
    c->ins = c->ins > delta ? c->ins - delta : 0;
    const uint32_t d = (uint32_t)delta;
    size_t hs = (size_t)1 << c->lp.hlog;
    for (size_t i = 0; i < hs; i++) c->head[i] = c->head[i] > d ? c->head[i] - d : 0;
    if (c->prev) for (size_t i = 0; i < c->wsize; i++) c->prev[i] = c->prev[i] > d ? c->prev[i] - d : 0;
    if (c->lp.opt) opt_slide(c, d);
}

int lzr_cctx_write(lzr_cctx *c, const uint8_t *src, size_t n, lzr_buf *out)
{
    if (c->finished) return LZR_ERR_PARAM;
    if (!c->header_done && write_header(c, out)) return LZR_ERR_NOMEM;
    c->crc = lzr_crc32(c->crc, src, n);
    c->total_in += n;
    while (n) {
        if (c->end == c->cap) slide(c);
        size_t take = c->cap - c->end;
        if (take > n) take = n;
        memcpy(c->buf + c->end, src, take);
        c->end += take;
        memset(c->buf + c->end, 0, 8);   /* deterministic hashing past the end */
        src += take;
        n -= take;
        /* Compress full blocks, but keep LZR_LOOKAHEAD bytes after the block
         * so hashing and tree comparisons at its last positions see real
         * data (and the output does not depend on input chunking). */
        while (c->end - c->pos >= LZR_BLOCK_SIZE + LZR_LOOKAHEAD) {
            int r = emit_block(c, c->pos, c->pos + LZR_BLOCK_SIZE, out);
            if (r) return r;
            c->pos += LZR_BLOCK_SIZE;
        }
    }
    return LZR_OK;
}

int lzr_cctx_end(lzr_cctx *c, lzr_buf *out)
{
    if (c->finished) return LZR_ERR_PARAM;
    if (!c->header_done && write_header(c, out)) return LZR_ERR_NOMEM;
    if (c->content_size >= 0 && (uint64_t)c->content_size != c->total_in) return LZR_ERR_PARAM;
    while (c->pos < c->end) {
        size_t be = c->end - c->pos > LZR_BLOCK_SIZE ? c->pos + LZR_BLOCK_SIZE : c->end;
        int r = emit_block(c, c->pos, be, out);
        if (r) return r;
        c->pos = be;
    }
    uint8_t t[LZR_BLOCK_HDR + 4] = {2};
    wr32le(t + 5, 4);
    wr32le(t + LZR_BLOCK_HDR, c->crc);
    c->finished = 1;
    return buf_append(out, t, sizeof(t));
}
