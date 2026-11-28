/*
 * enc.c - streaming encoder.
 *
 * Structure:
 *   lzr_cctx          sliding input buffer (2*window + block), hash-chain
 *                     tables, repeat offsets, running XXH64
 *   find_best()       repeat-offset probes + hash-chain walk, scored by gain
 *   parse_block()     lazy parse (depth 0..2) into sequences + literals;
 *                     level 0: parse_block_fast(), double-hash greedy
 *   emit_block()      sequences -> block type 5 (Huffman literals + tANS
 *                     sequences) or type 4 (range-coded tokens, level 11)
 *   slide()           drops history older than the window, rebases tables
 */
#include "core.h"
#include "codes.h"
#include "xxh64.h"
#include "x86.h"
#include "arc.h"

#include "enc_internal.h"
#if defined(__SSE2__) && !defined(LZR_NO_SIMD)
#include <emmintrin.h>
#endif
#if defined(LZR_TUNE) || defined(LZR_STATS)
#include <stdio.h>
#endif

/* Levels -3..-1: single-hash greedy finder (window, hash bits, accel;
 * see parse_block_fast1). Level 0: double-hash greedy finder (window,
 * hash bits; see parse_block_fast). Levels 1-6: lazy parsing (window, hash bits, lazy depth, search steps,
 * nice length) over rows of 2^row recent positions (hash bits = log2 of
 * all row entries; replaced hash chains in v0.6: 2-5x faster). Levels
 * 7-8: the same rows as the optimal parser's match source (chain = row
 * slots tried per position). Levels 9-11: binary tree + optimal parsing
 * (chain = tree depth, opt = number of
 * parsing passes; a second pass gained <0.5% ratio for ~40% more time, so
 * all use one). */
static const level_params LEVELS[LZR_LEVEL_MAX - LZR_LEVEL_MIN + 1] = {
    /* levels -3 .. -1: single-hash fast finder (window, hash bits, accel) */
    {19, 13, 0, 0, 0, 0, 0, 4, 0}, {20, 14, 0, 0, 0, 0, 0, 2, 0}, {20, 14, 0, 0, 0, 0, 0, 1, 0},
    {22, 16, 0, 0, 0, 0, 0, 0, 0},
    {20, 20, 0, 4, 16, 0, 4, 0, 0}, {21, 21, 1, 8, 24, 0, 4, 0, 0}, {22, 21, 1, 16, 32, 0, 4, 0, 0},
    {22, 22, 1, 32, 48, 0, 5, 0, 0}, {23, 23, 2, 48, 64, 0, 6, 0, 0}, {24, 23, 2, 64, 96, 0, 6, 0, 0},
    /* levels 7-8: row finder + light optimal parsing (lazy = inside-match
     * length, see row_gather; lite = OPT_LITE_* shortcuts) */
    {24, 22, 16, 16, 64, 1, 5, 0, OPT_LITE_NOCOMBO | OPT_LITE_NOSHORT | OPT_LITE_CUTLEN},
    {24, 23, 16, 32, 128, 1, 6, 0, OPT_LITE_NOSHORT | OPT_LITE_CUTLEN},
    /* levels 9-11: binary tree + optimal parsing */
    {25, 20, 0, 16, 128, 1, 0, 0, 0}, {26, 21, 0, 24, 273, 1, 0, 0, 0}, {26, 22, 0, 48, 273, 1, 0, 0, 0},
};


/* Match table entries: one per hash; level 0 (lp.row == 0) has a long-
 * and a short-hash table of that size. Row tables (rtag, rhead) are used
 * by levels 1-6 only (level 0 gets a small unused placeholder). */
static size_t head_entries(const level_params *lp)
{
    size_t n = (size_t)1 << lp->hlog;
    return !lp->opt && !lp->row && !lp->accel ? 2 * n : n;
}
static size_t rtag_size(const level_params *lp) { return lp->row ? (size_t)1 << lp->hlog : 64; }
static size_t rhead_size(const level_params *lp) { return lp->row ? (size_t)1 << (lp->hlog - lp->row) : 64; }

/* A window larger than the known input only wastes memory. */
static unsigned fit_wlog(const level_params *lp, int64_t content_size)
{
    /* the tree finder needs >= 2 blocks of slack, see opt_init() */
    unsigned min_wlog = lp->opt ? 20 : LZR_WLOG_MIN, w = lp->wlog;
    if (w < min_wlog) w = min_wlog;
    while (content_size >= 0 && w > min_wlog && ((uint64_t)1 << (w - 1)) >= (uint64_t)content_size) w--;
    return w;
}

lzr_cctx *lzr_cctx_new(int level, int window_log, int64_t content_size)
{
    if (level < LZR_LEVEL_MIN || level > LZR_LEVEL_MAX) return NULL;
    lzr_cctx *c = (lzr_cctx *)calloc(1, sizeof(*c));
    if (c) xxh64_init(&c->xs);
    if (!c) return NULL;
    c->lp = LEVELS[level - LZR_LEVEL_MIN];
    if (window_log) {
        if (window_log < LZR_WLOG_MIN || window_log > LZR_WLOG_MAX) { free(c); return NULL; }
        c->lp.wlog = (uint8_t)window_log;
    }
#ifdef LZR_TUNE   /* developer builds only: override level parameters */
    if (getenv("LZR_DEPTH")) c->lp.chain = (uint16_t)atoi(getenv("LZR_DEPTH"));
    if (getenv("LZR_NICE")) c->lp.nice = (uint16_t)atoi(getenv("LZR_NICE"));
    if (getenv("LZR_PASSES")) c->lp.opt = (uint8_t)atoi(getenv("LZR_PASSES"));
    if (getenv("LZR_WLOG")) c->lp.wlog = (uint8_t)atoi(getenv("LZR_WLOG"));
    if (getenv("LZR_HLOG")) c->lp.hlog = (uint8_t)atoi(getenv("LZR_HLOG"));
    if (getenv("LZR_ROW")) c->lp.row = (uint8_t)atoi(getenv("LZR_ROW"));
    if (getenv("LZR_LAZY")) c->lp.lazy = (uint8_t)atoi(getenv("LZR_LAZY"));
    if (getenv("LZR_ACCEL")) c->lp.accel = (uint8_t)atoi(getenv("LZR_ACCEL"));
    if (getenv("LZR_LITE")) c->lp.lite = (uint8_t)atoi(getenv("LZR_LITE"));
    if (getenv("LZR_THREADS")) c->threads = atoi(getenv("LZR_THREADS"));
#endif
    c->lp.wlog = (uint8_t)fit_wlog(&c->lp, content_size);
    if (c->lp.hlog > c->lp.wlog) c->lp.hlog = c->lp.wlog;
    c->wsize = (size_t)1 << c->lp.wlog;
    c->wmask = c->wsize - 1;
    c->cap = 2 * c->wsize + LZR_BLOCK_SIZE;
    c->content_size = content_size;
    c->reps[0] = 1; c->reps[1] = 4; c->reps[2] = 8;
    for (int k = 0; k < 8; k++) c->hcpos[k] = SIZE_MAX;

    c->buf = (uint8_t *)big_alloc(c->cap + LZR_SLACK, 0);
    c->head = (uint32_t *)big_alloc(head_entries(&c->lp) * sizeof(uint32_t), 1);
    if (!c->lp.opt || c->lp.row) {   /* row finder (lazy and mid levels) */
        c->rtag = (uint8_t *)big_alloc(rtag_size(&c->lp), 1);
        c->rhead = (uint8_t *)big_alloc(rhead_size(&c->lp), 1);
    }
    c->seq_ll = (uint32_t *)malloc(MAX_SEQ * sizeof(uint32_t));
    c->seq_ml = (uint32_t *)malloc(MAX_SEQ * sizeof(uint32_t));
    c->seq_dist = (uint32_t *)malloc(MAX_SEQ * sizeof(uint32_t));
    c->lits = (uint8_t *)malloc(LZR_BLOCK_SIZE + 32);
    /* block type 6: tree levels 9-10 only (levels 7-8 keep type 5 for decode speed) */
    if (c->lp.opt && !c->lp.row && level < LZR_ADAPTIVE_LEVEL) {
        c->lctx = (uint8_t *)malloc(LZR_BLOCK_SIZE + 32);
        c->ctmp = (uint8_t *)malloc(hufc_bound(LZR_BLOCK_SIZE));
        c->hufc = (hufc_table *)calloc(2, sizeof(hufc_table));
        if (c->hufc) c->hufc[0].work = (uint32_t *)malloc(HUFC_WORK * sizeof(uint32_t));
        c->ctx_try = 1;
        if (!c->lctx || !c->ctmp || !c->hufc || !c->hufc[0].work) { lzr_cctx_free(c); return NULL; }
    }
    c->llc = (uint8_t *)malloc(MAX_SEQ);
    c->mlc = (uint8_t *)malloc(MAX_SEQ);
    c->ofc = (uint8_t *)malloc(MAX_SEQ);
    c->xlm = (uint32_t *)malloc(MAX_SEQ * sizeof(uint32_t));
    c->xof = (uint32_t *)malloc(MAX_SEQ * sizeof(uint32_t));
    c->nlm = (uint8_t *)malloc(MAX_SEQ);
    c->nof = (uint8_t *)malloc(MAX_SEQ);
    {
        /* type 5: literals + tables + <= 11 bytes of bits per sequence;
         * type 4 stops at raw_len + one token (< 64 bytes) + pending bytes */
        size_t rb = huf_bound(LZR_BLOCK_SIZE) + 3 * (1 + 3 * FSE_MAXSYM) + 11 * (size_t)MAX_SEQ + 64;
        size_t ab = 2 * (size_t)LZR_BLOCK_SIZE + 4096;
        c->blk = (uint8_t *)malloc(64 + (rb > ab ? rb : ab));
    }
    c->adaptive = level >= LZR_ADAPTIVE_LEVEL;
#ifdef LZR_TUNE
    if (getenv("LZR_ADAPT")) c->adaptive = atoi(getenv("LZR_ADAPT"));
#endif
    c->tm = (tok_model *)malloc(2 * sizeof(tok_model));
    if (c->tm) tok_model_init(c->tm);
    c->huf = (huf_table *)calloc(2, sizeof(huf_table));
    c->fse = (fse_table *)calloc(6, sizeof(fse_table));
    if (c->lp.opt && opt_init(c)) { lzr_cctx_free(c); return NULL; }
    if (!c->buf || !c->head || ((!c->lp.opt || c->lp.row) && (!c->rtag || !c->rhead)) || !c->seq_ll || !c->seq_ml || !c->seq_dist || !c->lits ||
        !c->llc || !c->mlc || !c->ofc || !c->xlm || !c->xof || !c->nlm || !c->nof || !c->blk || !c->tm || !c->huf || !c->fse) {
        lzr_cctx_free(c);
        return NULL;
    }
    return c;
}

void lzr_cctx_set_threads(lzr_cctx *c, int threads)
{
    c->threads = threads < 1 ? 1 : (threads > LZR_MAX_THREADS ? LZR_MAX_THREADS : threads);
}

static void batch_free(struct lzr_batch *b);
static void seg_free(lzr_cctx *c);
static int seg_write(lzr_cctx *c, const uint8_t *src, size_t n, lzr_buf *out);
static int seg_end(lzr_cctx *c, lzr_buf *out);

void lzr_cctx_free(lzr_cctx *c)
{
    if (!c) return;
    if (c->sg) { seg_free(c); return; }
    if (c->bt) batch_free(c->bt);   /* lets an open batch's helpers finish first */
    if (c->wk) {
        for (int t = 0; t < c->threads; t++) {
            lzr_cctx *w = c->wk[t];
            if (!w) continue;
            w->buf = NULL;   /* shared */
            lzr_cctx_free(w);
        }
        for (int j = 0; j < c->jobs; j++) free(c->wout[j].p);
        free(c->wk);
        free(c->wout);
    }
    opt_free(c);
    free(c->buf); free(c->head); free(c->rtag); free(c->rhead);
    free(c->seq_ll); free(c->seq_ml); free(c->seq_dist);
    free(c->lits); free(c->llc); free(c->mlc); free(c->ofc); free(c->xlm); free(c->xof); free(c->nlm); free(c->nof);
    free(c->blk); free(c->tm); free(c->huf); free(c->fse); free(c->x86tmp);
    if (c->hufc) hufc_free(&c->hufc[0]);
    free(c->hufc); free(c->lctx); free(c->ctmp);
    free(c);
}

/* ---------------- match finding ---------------- */

/* ---- row match finder (levels 1-6) ----
 * The hash picks a row of R = 2^row slots and an 8-bit tag. A row keeps
 * the R most recent positions with that row hash, newest at rhead[row],
 * and their tags side by side, so one search reads one row of tags
 * (compared all at once), and the positions and data of the candidates
 * whose tag matches only. Unlike a chain walk, memory stays local. */
static inline uint32_t row_hash(const uint8_t *p, unsigned bits)
{
    uint64_t v;
    memcpy(&v, p, 8);
    return (uint32_t)(((v << (64 - 8 * LZR_HASH_LEN)) * 0x9FB21C651E98DF25ull) >> (64 - bits));
}

/* Row hash of buf[i], from the cache when present. */
static inline uint32_t row_hash_at(const lzr_cctx *c, size_t i, unsigned bits)
{
    return c->hcpos[i & 7] == i ? c->hc[i & 7] : row_hash(c->buf + i, bits);
}

static inline void hash_cache_clear(lzr_cctx *c)
{
    for (int k = 0; k < 8; k++) c->hcpos[k] = SIZE_MAX;
}

/* Bit k set where tags[k] == tag, k < R (R = 16, 32 or 64). */
static inline uint64_t row_match_mask(const uint8_t *tags, uint8_t tag, unsigned R)
{
#if defined(__SSE2__) && !defined(LZR_NO_SIMD)
    const __m128i t = _mm_set1_epi8((char)tag);
    uint64_t m = 0;
    for (unsigned k = 0; k < R; k += 16) {
        __m128i v = _mm_loadu_si128((const __m128i *)(tags + k));
        m |= (uint64_t)(uint16_t)_mm_movemask_epi8(_mm_cmpeq_epi8(v, t)) << k;
    }
    return m;
#else
    uint64_t m = 0;
    for (unsigned k = 0; k < R; k++) m |= (uint64_t)(tags[k] == tag) << k;
    return m;
#endif
}

#define LZR_AINLINE inline __attribute__((always_inline))
/* rlog is a constant in each specialized copy (parse_block_r4/5/6) */
static LZR_AINLINE void row_insert_upto_(lzr_cctx *c, size_t upto, const unsigned rlog)
{
    const unsigned bits = c->lp.hlog - rlog + 8, rmask = (1u << rlog) - 1;
    const size_t base = c->base, end = c->end;
    size_t limit = end >= LZR_HASH_LEN ? end - LZR_HASH_LEN + 1 : 0;
    if (upto > limit) upto = limit;
    /* locals (restrict): byte stores would otherwise force reloads */
    const uint8_t *restrict buf = c->buf;
    uint8_t *restrict rtag = c->rtag, *restrict rhead = c->rhead;
    uint32_t *restrict head = c->head, *restrict hc = c->hc;
    size_t *restrict hcpos = c->hcpos;
    for (size_t i = c->ins; i < upto; i++) {
        /* the hash of i is normally cached (computed 8 positions earlier);
         * a real branch, so the miss path costs nothing in the common case */
        uint32_t h;
        if (__builtin_expect(hcpos[i & 7] == i, 1)) h = hc[i & 7];
        else h = row_hash(buf + i, bits);
        const uint32_t row = h >> 8;
        const size_t ia = i + 8;   /* (takes the slot of i: read it first) */
        if (__builtin_expect(ia < limit, 1)) {   /* ia + LZR_HASH_LEN <= end */
            /* hash a position 8 ahead (real bytes only, so the result does
             * not depend on how input arrives), keep it, prefetch its row */
            uint32_t ha = row_hash(buf + ia, bits);
            hc[ia & 7] = ha; hcpos[ia & 7] = ia;
            __builtin_prefetch(rtag + ((size_t)(ha >> 8) << rlog));
            __builtin_prefetch(head + ((size_t)(ha >> 8) << rlog));
            __builtin_prefetch(rhead + (ha >> 8));
        }
        unsigned slot = (rhead[row] - 1) & rmask;
        rhead[row] = (uint8_t)slot;
        size_t e = ((size_t)row << rlog) + slot;
        rtag[e] = (uint8_t)h;
        head[e] = (uint32_t)(base + i + 1);
    }
    if (upto > c->ins) c->ins = upto;
}

static void row_insert_upto(lzr_cctx *c, size_t upto) { row_insert_upto_(c, upto, c->lp.row); }

/* ---- mid levels: row finder feeding the optimal parser ----
 * Candidates at p of strictly increasing length (first above out[0].len
 * when nonzero), at most `tries` row slots, then p is indexed. */
static LZR_AINLINE unsigned row_cands_(lzr_cctx *c, size_t p, size_t limit, mpair *out, unsigned max, unsigned tries,
                                       const unsigned rlog)
{
    const unsigned R = 1u << rlog;
    unsigned n = 0;
    if (out) {
        const uint8_t *cur = c->buf + p, *const buf = c->buf;
        const uint32_t h = row_hash_at(c, p, c->lp.hlog - rlog + 8), row = h >> 8;
        const size_t r0 = (size_t)row << rlog;
        const unsigned hd = c->rhead[row];
        uint64_t m = row_match_mask(c->rtag + r0, (uint8_t)h, R);
        if (hd) m = (m >> hd | m << (R - hd)) & (R == 64 ? ~0ull : (1ull << R) - 1);
        const uint32_t *const rpos = c->head + r0;
        const size_t base = c->base, wsize = c->wsize;
        const uint32_t first = rd32(cur);
        size_t best = out[0].len ? out[0].len : LZR_MIN_MATCH - 1;
        while (m && tries--) {
            unsigned k = (unsigned)__builtin_ctzll(m);
            m &= m - 1;
            size_t cp = (size_t)rpos[(hd + k) & (R - 1)] - 1 - base;
            size_t dist = p - cp;
            if (cp >= p || dist > wsize) break;
            const uint8_t *ref = buf + cp;
            if (ref[best] != cur[best] || rd32(ref) != first) continue;
            size_t len = match_len(ref, cur, limit);
            if (len > best) {
                best = len;
                out[n].len = (uint32_t)len; out[n].dist = (uint32_t)dist;
                if (++n == max || len == limit) break;
            }
        }
    }
    row_insert_upto_(c, p + 1, rlog);
    return n;
}

/* A match found at p - 1 continues at p one byte shorter, at no search
 * cost. Inside a match at least lp.lazy long (lazy is unused by the tree
 * levels), the row search makes only 2 tries: the carried match is
 * usually the best one there. */
static LZR_AINLINE int row_gather_(lzr_cctx *c, size_t bs, size_t be, size_t lim_end, mpair **pairs, size_t *pcap,
                                   uint32_t *mstart, unsigned max_cand, const unsigned rlog)
{
    const uint32_t inside = c->lp.lazy ? c->lp.lazy : UINT32_MAX;
    size_t np = 0;
    uint32_t cl = 0, cd = 0;
    if (c->ins < bs) row_insert_upto_(c, bs, rlog);
    for (size_t p = bs; p < be; p++) {
        mstart[p - bs] = (uint32_t)np;
        cl = cl > LZR_MIN_MATCH ? cl - 1 : 0;
        if (p + LZR_MIN_MATCH > lim_end) continue;
        if (np + max_cand + 1 > *pcap) {
            size_t ncap = *pcap * 2;
            mpair *o = (mpair *)realloc(*pairs, ncap * sizeof(mpair));
            if (!o) return 1;
            *pairs = o;
            *pcap = ncap;
        }
        mpair *out = *pairs + np;
        const uint32_t maxlen = (uint32_t)(be - p);
        size_t limit = lim_end - p;
        if (limit > c->lp.nice) limit = c->lp.nice;
        if (cl > maxlen) cl = maxlen;
        out[0].len = 0;
        unsigned n = row_cands_(c, p, limit, out, max_cand, cl >= inside ? 2 : c->lp.chain, rlog);
        uint32_t prev = 0, m = 0;
        for (unsigned k = 0; k < n; k++) {
            uint32_t len = out[k].len < maxlen ? out[k].len : maxlen;
            if (len <= prev) continue;    /* clamped duplicates */
            out[m].len = len; out[m].dist = out[k].dist; m++;
            prev = len;
        }
        if (cl > prev && cl >= LZR_MIN_MATCH) {   /* the carried match is the longest */
            out[m].len = cl; out[m].dist = cd; m++;
        }
        if (m) { cl = out[m - 1].len; cd = out[m - 1].dist; }
        np += m;
    }
    mstart[be - bs] = (uint32_t)np;
    row_insert_upto_(c, be, rlog);
    return 0;
}

int row_gather(lzr_cctx *c, size_t bs, size_t be, size_t lim_end, mpair **pairs, size_t *pcap, uint32_t *mstart,
               unsigned max_cand)
{
    switch (c->lp.row) {
    case 4: return row_gather_(c, bs, be, lim_end, pairs, pcap, mstart, max_cand, 4);
    case 5: return row_gather_(c, bs, be, lim_end, pairs, pcap, mstart, max_cand, 5);
    case 6: return row_gather_(c, bs, be, lim_end, pairs, pcap, mstart, max_cand, 6);
    default: return row_gather_(c, bs, be, lim_end, pairs, pcap, mstart, max_cand, c->lp.row);   /* LZR_TUNE */
    }
}

/* Insert positions [c->ins, upto). */
#define insert_upto(c, upto) row_insert_upto_(c, upto, rlog)

typedef struct { uint32_t len, off; int gain; } match_t;  /* off: 0..2 rep idx, else dist + 2 */

/* Score in quarter-bits: each matched byte is worth ~4, the offset costs its log2. */
static inline int gain_of(uint32_t len, uint32_t off) { return (int)len * 4 - (int)highbit32(off + 1); }

static LZR_AINLINE match_t find_best(lzr_cctx *c, size_t p, size_t maxlen, const uint32_t *reps, const unsigned rlog)
{
    match_t best = {0, 0, 0};
    const uint8_t *cur = c->buf + p;
    const uint32_t first = rd32(cur);

    /* 1) repeat offsets: cheapest to code, try them first. The three
     * 4-byte probes are made without branches (an unusable offset probes
     * cur itself and is masked out); only hits take the loop. */
    unsigned hits = 0;
    for (uint32_t r = 0; r < LZR_NUM_REPS; r++) {
        uint32_t d = reps[r];
        unsigned ok = d && d <= p && d <= c->wsize;   /* d 0: unknown (job start) */
        size_t dd = ok ? d : 0;
        hits |= (unsigned)(ok & (rd32(cur - dd) == first)) << r;
    }
    for (; hits; hits &= hits - 1) {
        uint32_t r = (uint32_t)__builtin_ctz(hits), d = reps[r];
        uint32_t len = (uint32_t)match_len(cur - d, cur, maxlen);
        int g = gain_of(len, r);
        if (len >= LZR_MIN_MATCH && g > best.gain) { best.len = len; best.off = r; best.gain = g; }
    }

    if (best.len == maxlen || best.len >= c->lp.nice) return best;

    /* 2) row finder: tag-matching slots, newest first */
    const unsigned R = 1u << rlog;
    const uint32_t h = row_hash_at(c, p, c->lp.hlog - rlog + 8), row = h >> 8;
    const size_t r0 = (size_t)row << rlog;
    const unsigned hd = c->rhead[row];
    uint64_t m = row_match_mask(c->rtag + r0, (uint8_t)h, R);
    /* rotate so that bit k is the k-th newest slot */
    if (hd) m = (m >> hd | m << (R - hd)) & (R == 64 ? ~0ull : (1ull << R) - 1);
    unsigned tries = c->lp.chain;
    size_t best_len = best.len >= LZR_MIN_MATCH ? best.len : LZR_MIN_MATCH - 1;
    const size_t base = c->base, wsize = c->wsize;
    const uint32_t *const rpos = c->head + r0;
    const uint8_t *const buf = c->buf;
    {   /* start loading the data of the first candidates: those loads
         * miss the cache and would otherwise run one after another */
        uint64_t q = m;
        for (int n = 0; n < 2; n++, q &= q - 1) {
            /* no branches: an empty mask or a dead slot only prefetches a
             * harmless address (prefetches never fault) */
            size_t cq = (size_t)rpos[(hd + (unsigned)__builtin_ctzll(q | 1ull << 63)) & (R - 1)] - 1 - base;
            __builtin_prefetch(buf + (cq < p ? cq : p));
        }
    }
    while (m && tries--) {
        unsigned k = (unsigned)__builtin_ctzll(m);
        m &= m - 1;
        /* slot -> buffer index; an empty (0) or slid-out slot wraps to a
         * huge index. Slots are newest first, so the first one out of
         * reach (also: dist > window) ends the search. Only earlier
         * positions are ever inserted, so cp < p for live slots. */
        size_t cp = (size_t)rpos[(hd + k) & (R - 1)] - 1 - base;
        size_t dist = p - cp;
        if (cp >= p || dist > wsize) break;
        const uint8_t *ref = buf + cp;
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
    }
    return best;
}

/* Lazy parse of buf[bs, be). Returns number of sequences; *nlit literals. */
static LZR_AINLINE size_t parse_block_(lzr_cctx *c, size_t bs, size_t be, size_t *nlit_out, const unsigned rlog)
{
    size_t p = bs, lit_start = bs, ns = 0, nlit = 0;
    uint32_t *reps = c->reps;

    while (p + LZR_MIN_MATCH <= be) {
        insert_upto(c, p);
        match_t m = find_best(c, p, be - p, reps, rlog);
        if (!m.len) { p++; continue; }

        /* Lazy evaluation: would starting 1 byte later be clearly better? */
        for (unsigned depth = 1; depth <= c->lp.lazy && m.len < c->lp.nice && p + 1 + LZR_MIN_MATCH <= be; depth++) {
            insert_upto(c, p + 1);
            match_t m2 = find_best(c, p + 1, be - p - 1, reps, rlog);
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

#undef insert_upto
static __attribute__((noinline)) size_t parse_block_r4(lzr_cctx *c, size_t bs, size_t be, size_t *n) { return parse_block_(c, bs, be, n, 4); }
static __attribute__((noinline)) size_t parse_block_r5(lzr_cctx *c, size_t bs, size_t be, size_t *n) { return parse_block_(c, bs, be, n, 5); }
static __attribute__((noinline)) size_t parse_block_r6(lzr_cctx *c, size_t bs, size_t be, size_t *n) { return parse_block_(c, bs, be, n, 6); }

/* ---- level 0: double-hash finder (zstd "dfast"-like) ----
 * Two direct-mapped tables of the newest position per hash: 8-byte
 * hashes (head[0, 2^hlog)) find long matches, 5-byte hashes (the next
 * 2^hlog entries) short ones. Repeat offsets never exceed the window
 * (they were all once explicit distances), so they need no window check. Only searched positions and a few per
 * match are indexed, and long literal runs are crossed with a growing
 * step, so the search does a fixed small amount of work per position. */
static inline uint64_t rd64(const uint8_t *p) { uint64_t v; memcpy(&v, p, 8); return v; }
static inline uint32_t hash8(uint64_t v, unsigned bits) { return (uint32_t)((v * 0xCF1BBCDCB7A56463ull) >> (64 - bits)); }
static inline uint32_t hash5(uint64_t v, unsigned bits) { return (uint32_t)(((v << 24) * 0x9FB21C651E98DF25ull) >> (64 - bits)); }
#define FAST_STEP_LOG 8   /* step + 1 every 256 literals */

static inline void fast_put(lzr_cctx *c, size_t i)
{
    const unsigned hl = c->lp.hlog;
    uint64_t v = rd64(c->buf + i);
    c->head[hash8(v, hl)] = (uint32_t)(c->base + i + 1);
    c->head[((size_t)1 << hl) + hash5(v, hl)] = (uint32_t)(c->base + i + 1);
}

/* Index [from, upto) (job history); positions need 8 real bytes. */
static void fast_fill(lzr_cctx *c, size_t from, size_t upto)
{
    size_t lim = c->end >= 8 ? c->end - 7 : 0;
    if (upto > lim) upto = lim;
    /* every third position (as zstd fills its tables from a dictionary):
     * the history only seeds matches, a third of the work */
    for (size_t i = from; i < upto; i += 3) fast_put(c, i);
}

static __attribute__((noinline)) size_t parse_block_fast(lzr_cctx *c, size_t bs, size_t be, size_t *nlit_out)
{
    uint32_t *restrict sll = c->seq_ll, *restrict sml = c->seq_ml, *restrict sdi = c->seq_dist;
    uint8_t *restrict lits = c->lits;
    const unsigned hl = c->lp.hlog;
    uint32_t *const tl = c->head, *const ts = c->head + ((size_t)1 << hl);
    const uint8_t *const buf = c->buf;
    const size_t base = c->base, wsize = c->wsize, lo = c->lo;
    uint32_t *reps = c->reps;
    size_t p = bs, anchor = bs, ns = 0, nlit = 0;
    /* searched positions read 8 bytes (+1 for the next position's long
     * hash); stay inside the block so the result does not depend on how
     * the input arrives */
    const size_t ilimit = be >= bs + 9 ? be - 9 : bs;

    while (p < ilimit) {
        const uint64_t v = rd64(buf + p);
        const uint32_t h8 = hash8(v, hl), h5 = hash5(v, hl);
        const size_t cl = (size_t)tl[h8] - 1 - base, cs = (size_t)ts[h5] - 1 - base;
        tl[h8] = ts[h5] = (uint32_t)(base + p + 1);
        size_t mstart, mlen, cand;
        uint32_t off;
        const uint32_t d0 = reps[0];
        if ((size_t)d0 - 1 < p + 1 && rd32(buf + p + 1 - d0) == rd32(buf + p + 1)) {
            /* repeat offset at p + 1 */
            mstart = p + 1;
            mlen = 4 + match_len(buf + mstart - d0 + 4, buf + mstart + 4, be - mstart - 4);
            cand = mstart - d0;
            off = 0;
        } else if (cl < p && p - cl <= wsize && rd64(buf + cl) == v) {
            mstart = p; cand = cl;
            mlen = 8 + match_len(buf + cl + 8, buf + p + 8, be - p - 8);
            off = (uint32_t)(p - cl) + 2;
        } else if (cs < p && p - cs <= wsize && rd32(buf + cs) == (uint32_t)v) {
            /* short match: a long match one byte later is usually better */
            const uint64_t v1 = rd64(buf + p + 1);
            const uint32_t h1 = hash8(v1, hl);
            const size_t c1 = (size_t)tl[h1] - 1 - base;
            tl[h1] = (uint32_t)(base + p + 2);
            if (c1 <= p && p + 1 - c1 <= wsize && rd64(buf + c1) == v1) {
                mstart = p + 1; cand = c1;
                mlen = 8 + match_len(buf + c1 + 8, buf + p + 9, be - p - 9);
            } else {
                mstart = p; cand = cs;
                mlen = 4 + match_len(buf + cs + 4, buf + p + 4, be - p - 4);
            }
            off = (uint32_t)(mstart - cand) + 2;
        } else {
            p += ((p - anchor) >> FAST_STEP_LOG) + 1;
            continue;
        }
        /* extend backwards over the pending literals */
        while (mstart > anchor && cand > lo && buf[mstart - 1] == buf[cand - 1]) {
            mstart--; cand--; mlen++;
        }
        size_t ll = mstart - anchor;
        for (size_t k = 0; k < ll; k += 16) memcpy(lits + nlit + k, buf + anchor + k, 16);   /* lits has slack */
        nlit += ll;
        sll[ns] = (uint32_t)ll;
        sml[ns] = (uint32_t)mlen;
        sdi[ns] = commit_offset(reps, off);
        ns++;
        p = mstart + mlen;
        anchor = p;
        /* index a few positions of the match */
        if (p < ilimit + 2) {
            fast_put(c, mstart + 2);
            fast_put(c, p - 2);
            ts[hash5(rd64(buf + p - 1), hl)] = (uint32_t)(base + p);
        }
        /* repeat offset 1 right after the match (literal length 0) */
        while (p < ilimit) {
            const uint32_t d1 = reps[1];
            if (!((size_t)d1 - 1 < p && rd32(buf + p - d1) == rd32(buf + p))) break;
            size_t l = 4 + match_len(buf + p - d1 + 4, buf + p + 4, be - p - 4);
            sll[ns] = 0;
            sml[ns] = (uint32_t)l;
            sdi[ns] = commit_offset(reps, 1);
            ns++;
            fast_put(c, p);
            p += l;
            anchor = p;
        }
    }
    memcpy(lits + nlit, buf + anchor, be - anchor);
    nlit += be - anchor;
    *nlit_out = nlit;
    return ns;
}

/* ---- levels -1 .. -3: single-hash finder (zstd "fast"-like) ----
 * One direct-mapped table of the newest position per 6-byte hash
 * (2^hlog entries, small enough to stay in L1/L2). A miss advances by
 * accel + (literals so far >> 8), so incompressible data is crossed
 * quickly; matches are indexed at two positions. */
static inline uint32_t hash6(uint64_t v, unsigned bits) { return (uint32_t)(((v << 16) * 0xCF1BBCDCB7A56463ull) >> (64 - bits)); }

static inline void fast1_put(lzr_cctx *c, size_t i)
{
    c->head[hash6(rd64(c->buf + i), c->lp.hlog)] = (uint32_t)(c->base + i + 1);
}

/* Index [from, upto) (job history), every third position. */
static void fast1_fill(lzr_cctx *c, size_t from, size_t upto)
{
    size_t lim = c->end >= 8 ? c->end - 7 : 0;
    if (upto > lim) upto = lim;
    for (size_t i = from; i < upto; i += 3) fast1_put(c, i);
}

static __attribute__((noinline)) size_t parse_block_fast1(lzr_cctx *c, size_t bs, size_t be, size_t *nlit_out)
{
    uint32_t *restrict sll = c->seq_ll, *restrict sml = c->seq_ml, *restrict sdi = c->seq_dist;
    uint8_t *restrict lits = c->lits;
    const unsigned hl = c->lp.hlog;
    const size_t accel = c->lp.accel;
    uint32_t *const t = c->head;
    const uint8_t *const buf = c->buf;
    const size_t base = c->base, wsize = c->wsize, lo = c->lo;
    uint32_t *reps = c->reps;
    size_t p = bs, anchor = bs, ns = 0, nlit = 0;
    const size_t ilimit = be >= bs + 9 ? be - 9 : bs;   /* as parse_block_fast */

    /* Two-stage pipeline: the table slot of the next position (and the
     * bytes of its candidate) are loaded before the current position is
     * resolved, so table and window loads overlap instead of queueing. */
    size_t ip = p;
    uint32_t h = hash6(rd64(buf + ip), hl), idx = t[h];
    while (ip < ilimit) {
        const size_t nip = ip + accel + ((ip - anchor) >> FAST_STEP_LOG);
        const size_t hp = nip < ilimit ? nip : ip;   /* stay in the block */
        const uint32_t nh = hash6(rd64(buf + hp), hl), nidx = t[nh];
        __builtin_prefetch(buf + (((size_t)nidx - 1 - base) < hp ? (size_t)nidx - 1 - base : hp));
        const uint32_t v = rd32(buf + ip);
        t[h] = (uint32_t)(base + ip + 1);
        const size_t cand0 = (size_t)idx - 1 - base;
        const int valid = cand0 < ip && ip - cand0 <= wsize;
        const int hit = valid & (rd32(buf + (valid ? cand0 : ip)) == v);
        size_t mstart, mlen, cand;
        uint32_t off;
        const uint32_t d0 = reps[0];
        if ((size_t)d0 - 1 < ip + 1 && rd32(buf + ip + 1 - d0) == rd32(buf + ip + 1)) {
            mstart = ip + 1;
            mlen = 4 + match_len(buf + mstart - d0 + 4, buf + mstart + 4, be - mstart - 4);
            cand = mstart - d0;
            off = 0;
        } else if (hit) {
            mstart = ip; cand = cand0;
            mlen = 4 + match_len(buf + cand + 4, buf + ip + 4, be - ip - 4);
            off = (uint32_t)(ip - cand) + 2;
        } else {
            ip = nip; h = nh; idx = nidx;
            continue;
        }
        while (mstart > anchor && cand > lo && buf[mstart - 1] == buf[cand - 1]) {
            mstart--; cand--; mlen++;
        }
        size_t ll = mstart - anchor;
        for (size_t k = 0; k < ll; k += 16) memcpy(lits + nlit + k, buf + anchor + k, 16);   /* lits has slack */
        nlit += ll;
        sll[ns] = (uint32_t)ll;
        sml[ns] = (uint32_t)mlen;
        sdi[ns] = commit_offset(reps, off);
        ns++;
        ip = mstart + mlen;
        anchor = ip;
        if (ip < ilimit + 2) {
            fast1_put(c, mstart + 2);
            fast1_put(c, ip - 2);
        }
        while (ip < ilimit) {
            const uint32_t d1 = reps[1];
            if (!((size_t)d1 - 1 < ip && rd32(buf + ip - d1) == rd32(buf + ip))) break;
            size_t l = 4 + match_len(buf + ip - d1 + 4, buf + ip + 4, be - ip - 4);
            sll[ns] = 0;
            sml[ns] = (uint32_t)l;
            sdi[ns] = commit_offset(reps, 1);
            ns++;
            fast1_put(c, ip);
            ip += l;
            anchor = ip;
        }
        if (ip < ilimit) { h = hash6(rd64(buf + ip), hl); idx = t[h]; }
    }
    memcpy(lits + nlit, buf + anchor, be - anchor);
    nlit += be - anchor;
    *nlit_out = nlit;
    return ns;
}

static size_t parse_block(lzr_cctx *c, size_t bs, size_t be, size_t *n)
{
    switch (c->lp.row) {
    case 0: return c->lp.accel ? parse_block_fast1(c, bs, be, n) : parse_block_fast(c, bs, be, n);
    case 4: return parse_block_r4(c, bs, be, n);
    case 5: return parse_block_r5(c, bs, be, n);
    case 6: return parse_block_r6(c, bs, be, n);
    default: return parse_block_(c, bs, be, n, c->lp.row);   /* LZR_TUNE row sizes */
    }
}

/* ---------------- block emission ---------------- */

static int write_header(lzr_cctx *c, lzr_buf *out)
{
    uint8_t h[16];
    memcpy(h, LZR_MAGIC, 4);
    h[4] = LZR_VERSION;
    h[5] = (uint8_t)((c->content_size >= 0 ? LZR_FLAG_SIZE : 0) | (c->x86tmp ? LZR_FLAG_X86 : 0));
    h[6] = c->lp.wlog;
    h[7] = 0;
    size_t n = 8;
    if (c->content_size >= 0) { wr64le(h + 8, (uint64_t)c->content_size); n = 16; }
    c->header_done = 1;
    return buf_append(out, h, n);
}

/* Block type 4 payload (arc.h): one range-coded token stream. The match
 * byte of a literal after a match is at distance rep0, so the repeat
 * offsets are replayed from their state at the start of the block.
 * Returns the payload size, or (size_t)-1 once it exceeds raw (the block
 * is then stored). */
static size_t encode_tokens(lzr_cctx *c, size_t bs, size_t be, size_t ns, const uint32_t *reps_in, uint8_t *out)
{
    tok_model *m = c->tm;
    arc_enc e;
    arce_init(&e, out);
    uint32_t reps[LZR_NUM_REPS];
    memcpy(reps, reps_in, sizeof(reps));
    const uint8_t *buf = c->buf;
    const size_t raw = be - bs;
    size_t p = bs;
    for (size_t i = 0; i <= ns; i++) {
        size_t n = i < ns ? c->seq_ll[i] : be - p;
        for (size_t t = 0; t < n; t++, p++)
            tok_enc_lit(&e, m, buf + p, p ? buf[p - 1] : 0, reps[0], (unsigned)(p - bs) & 3);
        if ((size_t)(e.p - out) + e.cache_size > raw) return (size_t)-1;
        if (i == ns) break;
        uint32_t off = c->seq_dist[i];
        commit_offset(reps, off);
        size_t ml = c->seq_ml[i];
        while (ml) {
            unsigned len = tok_piece(ml);
            tok_enc_match(&e, m, off, len, (unsigned)(p - bs) & 3);
            off = 0;
            p += len;
            ml -= len;
        }
    }
    size_t n = arce_finish(&e, out);
    return n < raw ? n : (size_t)-1;
}

/* Block type 5 payload into b + LZR_BLOCK_HDR; returns the block size.
 *   nseq u32 | nlit u32 | literal section (huf.h) |
 *   if nseq: 3 table modes (ll, ml, of: 0 new + counts, 1 previous,
 *   2 single symbol + that symbol) | bit stream (fse.h) to the block end */
/* Block type 6 only when its literal section is at least this many 1/1000
 * smaller than type 5's: decoding context literals is slower. */
#ifndef CTX_MIN_GAIN
#define CTX_MIN_GAIN 20
#endif
#ifndef CTX_MIN_LITS
#define CTX_MIN_LITS 1024
#endif

/* Type 5 (or 6, see *btype) block body for sequences [s0, s0 + ns) and
 * literals [l0, l0 + nlit). */
static size_t emit_seqs5(lzr_cctx *c, size_t s0, size_t ns, size_t l0, size_t nlit, uint8_t *b, uint8_t *btype)
{
    size_t k = LZR_BLOCK_HDR;
    wr32le(b + k, (uint32_t)ns); k += 4;
    wr32le(b + k, (uint32_t)nlit); k += 4;
    *btype = 5;
    if (c->hufc && c->ctx_try && nlit >= CTX_MIN_LITS) {
        /* both literal codings; type 6 leaves the type 5 table as it was */
        hufc_table hs = c->hufc[0];
        huf_table h5 = c->huf[0];
        size_t n6 = hufc_encode(c->ctmp, c->lits + l0, c->lctx + l0, nlit, &c->hufc[0]);
        size_t n5 = huf_encode(b + k, c->lits + l0, nlit, &c->huf[0]);
#ifdef LZR_STATS
        fprintf(stderr, "LZRCTX %zu %zu %zu %d %d\n", nlit, n5, n6, n6 ? c->ctmp[0] : -1, b[k]);
#endif
        if (n6 && n6 * 1000 < n5 * (1000 - CTX_MIN_GAIN)) {
            c->huf[0] = h5;
            memcpy(b + k, c->ctmp, n6);
            k += n6;
            *btype = 6;
        } else {
            hs.work = c->hufc[0].work;   /* (scratch, allocated once) */
            c->hufc[0] = hs;
            k += n5;
        }
    } else {
        k += huf_encode(b + k, c->lits + l0, nlit, &c->huf[0]);
    }
#ifdef LZR_STATS
    const size_t st_lit_end = k;
#endif
    if (!ns) {
#ifdef LZR_STATS
        fprintf(stderr, "LZRSTAT5 %zu %zu %zu 0 0\n", nlit, st_lit_end - LZR_BLOCK_HDR - 8, ns);
#endif
        return k;
    }

    uint8_t *codes[3] = {c->llc, c->mlc, c->ofc};
    const unsigned alpha[3] = {VAL_CODES, ML_CODES, OF_CODES};
    uint32_t cnt[3][FSE_MAXSYM];
    memset(cnt, 0, sizeof(cnt));
    /* codes, and the extra bits in the order they are written: ll and ml
     * together (<= 16 + 16 bits, block size 2^18), then the offset.
     * Arrays in restrict locals: byte stores would otherwise make the
     * compiler reload every c->... pointer after each one. */
    {
        const uint32_t *restrict sll = c->seq_ll + s0, *restrict sml = c->seq_ml + s0, *restrict sdi = c->seq_dist + s0;
        uint8_t *restrict llc = c->llc, *restrict mlc = c->mlc, *restrict ofc = c->ofc;
        uint8_t *restrict nlm = c->nlm, *restrict nof = c->nof;
        uint32_t *restrict xlm = c->xlm, *restrict xof = c->xof;
        for (size_t i = 0; i < ns; i++) {
            unsigned nl, nm, no;
            uint32_t el, em, eo;
            const unsigned lc = code_bf(sll[i], 4, &nl, &el);
            const unsigned mc = code_bf(sml[i] - LZR_ML_BASE, ML_DLOG, &nm, &em);
            const uint32_t off = sdi[i];
            const int rep = off < LZR_NUM_REPS;
            const unsigned oc0 = code_bf(rep ? 0 : off - LZR_NUM_REPS, 4, &no, &eo);
            const unsigned oc = rep ? off : LZR_NUM_REPS + oc0;
            llc[i] = (uint8_t)lc; mlc[i] = (uint8_t)mc; ofc[i] = (uint8_t)oc;
            xlm[i] = el | em << nl;
            nlm[i] = (uint8_t)(nl + nm);
            xof[i] = eo;
            nof[i] = (uint8_t)no;
            cnt[0][lc]++; cnt[1][mc]++; cnt[2][oc]++;
        }
    }
    int rle[3];
    for (int t = 0; t < 3; t++) {
        fse_table *T = &c->fse[t];
        unsigned distinct = 0, sym = 0;
        for (unsigned v = 0; v < alpha[t]; v++) if (cnt[t][v]) { distinct++; sym = v; }
        rle[t] = distinct == 1;
        if (rle[t]) { b[k++] = 2; b[k++] = (uint8_t)sym; continue; }
        uint16_t norm[FSE_MAXSYM];
        fse_normalize(cnt[t], alpha[t], (uint32_t)ns, norm);
        uint8_t hdr[3 * FSE_MAXSYM];
        size_t hl = fse_write_norm(hdr, norm, alpha[t]);
        uint64_t cn = fse_cost(cnt[t], norm, alpha[t]) + hl * 8 * 256;
        uint64_t cp = T->valid && T->alphabet == alpha[t] ? fse_cost(cnt[t], T->norm, alpha[t]) : UINT64_MAX;
        if (cp <= cn) { b[k++] = 1; continue; }
        b[k++] = 0;
        memcpy(b + k, hdr, hl); k += hl;
        memcpy(T->norm, norm, sizeof(norm));
        T->alphabet = alpha[t];
        fse_build(T);
    }
#ifdef LZR_STATS
    const size_t st_tab_end = k;
    uint64_t st_extra = 0;
    for (size_t i = 0; i < ns; i++) st_extra += (uint64_t)c->nlm[i] + c->nof[i];
#endif
    /* bit stream, backwards over the sequences (see dec.c for the order) */
    fbw w = {0, 0, b + k};
    uint32_t x[3] = {FSE_L, FSE_L, FSE_L};
    if (!rle[0] && !rle[1] && !rle[2]) {
        /* common case: the three state transitions in one write (<= 27 bits) */
        const fse_table *T0 = &c->fse[0], *T1 = &c->fse[1], *T2 = &c->fse[2];
        const uint8_t *restrict llc = c->llc, *restrict mlc = c->mlc, *restrict ofc = c->ofc;
        const uint8_t *restrict nlm = c->nlm, *restrict nof = c->nof;
        const uint32_t *restrict xlm = c->xlm, *restrict xof = c->xof;
        uint32_t x0 = x[0], x1 = x[1], x2 = x[2];
        fbw wl = w;   /* a local writer: its stores cannot alias the arrays */
        /* branch-free flushes: n <= 7 + 27, then <= 7 + 32 + 24 (offsets
         * < 2^27) bits between them */
        for (size_t i = ns; i-- > 0;) {
            unsigned n2, n1, n0;
            uint32_t b2 = fse_step(T2, &x2, ofc[i], &n2);
            uint32_t b1 = fse_step(T1, &x1, mlc[i], &n1);
            uint32_t b0 = fse_step(T0, &x0, llc[i], &n0);
            fbw_add(&wl, b2 | b1 << n2 | (uint64_t)b0 << (n2 + n1), n2 + n1 + n0);
            fbw_flush8(&wl);
            fbw_add(&wl, xlm[i], nlm[i]);
            fbw_add(&wl, xof[i], nof[i]);
            fbw_flush8(&wl);
        }
        w = wl;
        x[0] = x0; x[1] = x1; x[2] = x2;
    } else {
        for (size_t i = ns; i-- > 0;) {
            for (int t = 2; t >= 0; t--)
                if (!rle[t]) fse_put(&w, &c->fse[t], &x[t], codes[t][i]);
            fbw_put(&w, c->xlm[i], c->nlm[i]);
            fbw_put(&w, c->xof[i], c->nof[i]);
        }
    }
    for (int t = 2; t >= 0; t--)
        if (!rle[t]) fbw_put(&w, x[t] - FSE_L, FSE_LOG);
#ifdef LZR_STATS
    {
        size_t end = (size_t)(fbw_close(&w) - b);
        /* nlit, literal-section bytes, nseq, table bytes, bit-stream bytes, extra bits */
        fprintf(stderr, "LZRSTAT5 %zu %zu %zu %zu %zu %llu\n", nlit, st_lit_end - LZR_BLOCK_HDR - 8, ns,
                st_tab_end - st_lit_end, end - st_tab_end, (unsigned long long)st_extra);
        return end;
    }
#endif
    return (size_t)(fbw_close(&w) - b);
}

/* ---- cost-based block splitting (block type 5, optimal-parse levels) ----
 * A parsed block is also encoded as two halves (cut at the sequence
 * nearest the middle), recursively, and the halves are kept when their
 * own tables make them smaller in total. Nothing changes in the format:
 * each part is an ordinary type 5 block. Parts are never stored raw (a
 * stored block leaves the decoder's repeat offsets as they were, and the
 * sequences after it were coded with the advanced ones). */
#ifndef SPLIT_DEPTH
#define SPLIT_DEPTH 3
#endif
#ifndef SPLIT_MIN_SEQ
#define SPLIT_MIN_SEQ 128
#endif
typedef struct { huf_table huf; fse_table fse[3]; hufc_table hufc; } tab_state;
static void tab_save(const lzr_cctx *c, tab_state *t)
{
    t->huf = c->huf[0];
    memcpy(t->fse, c->fse, sizeof(t->fse));
    if (c->hufc) t->hufc = c->hufc[0];
}
static void tab_load(lzr_cctx *c, const tab_state *t)
{
    c->huf[0] = t->huf;
    memcpy(c->fse, t->fse, sizeof(t->fse));
    if (c->hufc) { uint32_t *w = c->hufc[0].work; c->hufc[0] = t->hufc; c->hufc[0].work = w; }
}

/* Encode sequences [s0, s1) and literals [l0, l1) (raw bytes [p0, p1)) as
 * type 5 block(s) into dst; returns the size, or (size_t)-1 when a part
 * would not be smaller than raw, or (size_t)-2 on allocation failure. */
static size_t split5(lzr_cctx *c, size_t s0, size_t s1, size_t l0, size_t l1, size_t p0, size_t p1,
                     unsigned depth, uint8_t *dst)
{
    tab_state t0, tw;
    tab_save(c, &t0);
    uint8_t bt;
    size_t k = emit_seqs5(c, s0, s1 - s0, l0, l1 - l0, dst, &bt);
    const size_t raw = p1 - p0;
    if (k - LZR_BLOCK_HDR >= raw) { tab_load(c, &t0); return (size_t)-1; }
    dst[0] = bt;
    wr32le(dst + 1, (uint32_t)raw);
    wr32le(dst + 5, (uint32_t)(k - LZR_BLOCK_HDR));
    if (!depth || s1 - s0 < 2 * SPLIT_MIN_SEQ) return k;

    /* cut at the sequence where half of the raw bytes are covered */
    size_t sm = s0, lm = l0, pm = p0;
    while (sm < s1 && pm - p0 < raw / 2) {
        pm += c->seq_ll[sm] + c->seq_ml[sm];
        lm += c->seq_ll[sm];
        sm++;
    }
    if (sm - s0 < SPLIT_MIN_SEQ || s1 - sm < SPLIT_MIN_SEQ) return k;
    tab_save(c, &tw);
    tab_load(c, &t0);
    uint8_t *tmp = (uint8_t *)malloc(2 * (k + 64) + 2 * (3 * (1 + 3 * FSE_MAXSYM) + 64) + 2 * hufc_bound(raw) + 11 * (s1 - s0));
    if (!tmp) return (size_t)-2;
    /* halves only try context literals when the whole part chose them */
    const int try6 = c->ctx_try;
    c->ctx_try = try6 && bt == 6;
    size_t a = split5(c, s0, sm, l0, lm, p0, pm, depth - 1, tmp);
    size_t b = a >= (size_t)-2 ? a : split5(c, sm, s1, lm, l1, pm, p1, depth - 1, tmp + a);
    c->ctx_try = try6;
    if (a == (size_t)-2 || b == (size_t)-2) { free(tmp); return (size_t)-2; }
    if (a != (size_t)-1 && b != (size_t)-1 && a + b < k) {
        memcpy(dst, tmp, a + b);
        free(tmp);
        return a + b;
    }
    free(tmp);
    tab_load(c, &tw);
    return k;
}

static size_t emit_split5(lzr_cctx *c, size_t s0, size_t s1, size_t l0, size_t l1, size_t p0, size_t p1,
                          unsigned depth, lzr_buf *out)
{
    size_t n = split5(c, s0, s1, l0, l1, p0, p1, depth, c->blk);
    if (n >= (size_t)-2) return n;
    return buf_append(out, c->blk, n) ? (size_t)-2 : n;
}

static int emit_block(lzr_cctx *c, size_t bs, size_t be, lzr_buf *out)
{
    const size_t raw = be - bs;
    uint32_t saved_reps[LZR_NUM_REPS];
    memcpy(saved_reps, c->reps, sizeof(saved_reps));

    size_t nlit;
    size_t ns = c->lp.opt ? opt_parse_block(c, bs, be, &nlit) : parse_block(c, bs, be, &nlit);
    if (ns == (size_t)-1) return LZR_ERR_NOMEM;

    uint8_t *b = c->blk;
    size_t k = LZR_BLOCK_HDR;
    uint8_t btype = 5;
    tok_model *tm = c->tm;
    if (!c->adaptive) {   /* save the repeatable tables for the stored fallback */
        c->huf[1] = c->huf[0];
        memcpy(c->fse + 3, c->fse, 3 * sizeof(fse_table));
    }
    if (c->adaptive) {
        btype = 4;
        memcpy(tm + 1, tm, sizeof(*tm));   /* for the stored fallback */
        size_t n = encode_tokens(c, bs, be, ns, saved_reps, b + k);
        k = n == (size_t)-1 ? (size_t)-1 : k + n;
    } else if (c->lp.opt) {
        /* optimal-parse levels: let the block be cut where separate tables pay */
        size_t n = emit_split5(c, 0, ns, 0, nlit, bs, be, SPLIT_DEPTH, out);
        if (n == (size_t)-2) return LZR_ERR_NOMEM;
        if (n != (size_t)-1) return LZR_OK;
        k = (size_t)-1;   /* stored below */
    } else {
        k = emit_seqs5(c, 0, ns, 0, nlit, b, &btype);
    }

    if (k != (size_t)-1 && k - LZR_BLOCK_HDR < raw) {
        b[0] = btype;
        wr32le(b + 1, (uint32_t)raw);
        wr32le(b + 5, (uint32_t)(k - LZR_BLOCK_HDR));
        return buf_append(out, b, k);
    }
    /* Stored block; the decoder leaves its repeat offsets and model untouched. */
    memcpy(c->reps, saved_reps, sizeof(saved_reps));
    if (btype == 4) memcpy(tm, tm + 1, sizeof(*tm));
    if (btype == 5) { c->huf[0] = c->huf[1]; memcpy(c->fse, c->fse + 3, 3 * sizeof(fse_table)); }
    uint8_t h[LZR_BLOCK_HDR] = {0};
    wr32le(h + 1, (uint32_t)raw);
    wr32le(h + 5, (uint32_t)raw);
    if (buf_append(out, h, sizeof(h))) return LZR_ERR_NOMEM;
    return buf_append(out, c->buf + bs, raw);
}

#ifndef LZR_REBASE_AT
#define LZR_REBASE_AT ((size_t)1 << 31)
#endif

/* Keep the last >= wsize bytes; shift by a multiple of wsize so the
 * tree indices (pos & wmask, levels 9-11) stay valid. */
static void slide(lzr_cctx *c)
{
    if (c->pos <= c->wsize) return;
    size_t delta = ((c->pos - c->wsize) / c->wsize) * c->wsize;
    if (!delta) return;
    if (c->opt) opt_sync(c);   /* a background gather reads buf */
    memmove(c->buf, c->buf + delta, c->end - delta);
    c->pos -= delta;
    c->end -= delta;
    c->ins = c->ins > delta ? c->ins - delta : 0;
    hash_cache_clear(c);   /* positions moved */
    size_t hs = head_entries(&c->lp);
    if (c->lp.opt && c->lp.row) opt_slide(c, (uint32_t)delta);   /* short-match tables */
    if (!c->lp.opt || c->lp.row) {
        /* Levels 1-6 store base + index, so sliding only moves the base
         * (row slots do not depend on positions).
         * The table is rewritten only before positions would pass 2^31. */
        c->base += delta;
        if (c->base + c->cap < LZR_REBASE_AT) return;
        const uint32_t d = (uint32_t)c->base;
        for (size_t i = 0; i < hs; i++) c->head[i] = c->head[i] > d ? c->head[i] - d : 0;
        c->base = 0;
        return;
    }
    const uint32_t d = (uint32_t)delta;
    for (size_t i = 0; i < hs; i++) c->head[i] = c->head[i] > d ? c->head[i] - d : 0;
    opt_slide(c, d);
}

/* ---------------- job mode (levels 1-6, threads > 1) ----------------
 * The input is cut into jobs of jsize bytes, independent of the thread
 * count. A worker compresses one job alone: it indexes the jover bytes
 * before the job (so matches still reach back), then parses and codes the
 * job's blocks with no state from earlier jobs. The format needs no change:
 *  - repeat offsets start unknown (0) and are used only once known; an
 *    explicit offset updates the decoder's repeat history the same way
 *    whatever it held, so after it the encoder knows the decoder's state;
 *  - the first block of a job builds new entropy tables (no repeat mode).
 * The main thread appends the jobs' blocks in order. Output is identical
 * for any thread count > 1; one thread keeps the sequential path, whose
 * ratio is slightly better (no job borders). */
#ifndef LZR_JOB_LOG
#define LZR_JOB_LOG 22        /* 4 MB jobs */
#endif
#if LZR_HAVE_THREADS
#include <pthread.h>
#endif

typedef struct { lzr_cctx *c; const uint8_t *buf; size_t span, js, je; int r; } job_arg;   /* buf[js, je), positions < span */

/* A batch: up to c->jobs jobs over buf[pos, cap). Jobs are queued as soon
 * as their input is buffered (batch_feed), and helper threads start on
 * them while the caller is still copying (and hashing) the rest; the
 * caller joins in when the batch is closed (batch_finish). Every job is
 * compressed the same way whichever thread takes it, so the output does
 * not depend on the thread count or on timing. */
struct lzr_batch {
    job_arg a[2 * LZR_MAX_THREADS];
    int nready, next, closed, nthr;
#if LZR_HAVE_THREADS
    pthread_mutex_t mu;
    pthread_cond_t cv;
    pthread_t tid[LZR_MAX_THREADS];
    struct { struct lzr_batch *b; lzr_cctx *w; lzr_buf *outs; } arg[LZR_MAX_THREADS];
#endif
};


static lzr_cctx *worker_new(const lzr_cctx *c)
{
    lzr_cctx *w = (lzr_cctx *)calloc(1, sizeof(*w));
    if (!w) return NULL;
    w->lp = c->lp;
    w->wsize = c->wsize; w->wmask = c->wmask; w->cap = c->cap;
    w->threads = 1;
    w->head = (uint32_t *)big_alloc(head_entries(&w->lp) * sizeof(uint32_t), 1);
    w->rtag = (uint8_t *)big_alloc(rtag_size(&w->lp), 1);
    w->rhead = (uint8_t *)big_alloc(rhead_size(&w->lp), 1);
    w->seq_ll = (uint32_t *)malloc(MAX_SEQ * sizeof(uint32_t));
    w->seq_ml = (uint32_t *)malloc(MAX_SEQ * sizeof(uint32_t));
    w->seq_dist = (uint32_t *)malloc(MAX_SEQ * sizeof(uint32_t));
    w->lits = (uint8_t *)malloc(LZR_BLOCK_SIZE + 32);
    w->llc = (uint8_t *)malloc(MAX_SEQ);
    w->mlc = (uint8_t *)malloc(MAX_SEQ);
    w->ofc = (uint8_t *)malloc(MAX_SEQ);
    w->xlm = (uint32_t *)malloc(MAX_SEQ * sizeof(uint32_t));
    w->xof = (uint32_t *)malloc(MAX_SEQ * sizeof(uint32_t));
    w->nlm = (uint8_t *)malloc(MAX_SEQ);
    w->nof = (uint8_t *)malloc(MAX_SEQ);
    w->blk = (uint8_t *)malloc(64 + huf_bound(LZR_BLOCK_SIZE) + 3 * (1 + 3 * FSE_MAXSYM) + 11 * (size_t)MAX_SEQ + 64);
    w->huf = (huf_table *)calloc(2, sizeof(huf_table));
    w->fse = (fse_table *)calloc(6, sizeof(fse_table));
    if (!w->head || !w->rtag || !w->rhead || !w->seq_ll || !w->seq_ml || !w->seq_dist || !w->lits || !w->llc ||
        !w->mlc || !w->ofc || !w->xlm || !w->xof || !w->nlm || !w->nof || !w->blk || !w->huf || !w->fse) {
        lzr_cctx_free(w);
        return NULL;
    }
    return w;
}

/* Before the first write: switch to job mode if it applies. The buffer
 * then holds jover bytes of history plus one batch (threads jobs). */
static int jobs_setup(lzr_cctx *c)
{
    if (c->lp.opt || c->threads < 2) return LZR_OK;
    /* history indexed before each job: the whole window up to 8 MB, else
     * half of it (measured: ratio within ~0.1% of one thread on our sets;
     * a shorter history lost >1% on data with far repeats) */
    c->jsize = (size_t)1 << LZR_JOB_LOG;
    c->jover = c->wsize <= ((size_t)8 << 20) ? c->wsize : c->wsize / 2;
#ifdef LZR_TUNE
    if (getenv("LZR_JOB")) c->jsize = (size_t)1 << atoi(getenv("LZR_JOB"));
    if (getenv("LZR_OVL")) c->jover = (size_t)1 << atoi(getenv("LZR_OVL"));
    if (c->jover > c->wsize) c->jover = c->wsize;
#endif
    /* a batch holds 2 jobs per thread: threads take jobs as they finish,
     * which evens out jobs of unequal difficulty */
    const int K = 2 * c->threads;
    size_t cap = c->jover + (size_t)K * c->jsize;
    uint8_t *nb = (uint8_t *)big_alloc(cap + LZR_SLACK, 0);
    lzr_cctx **wk = (lzr_cctx **)calloc((size_t)c->threads, sizeof(*wk));
    lzr_buf *wout = (lzr_buf *)calloc((size_t)K, sizeof(*wout));
    int ok = nb && wk && wout;
    for (int t = 0; ok && t < c->threads; t++) ok = (wk[t] = worker_new(c)) != NULL;
    if (!ok) {
        for (int t = 0; wk && t < c->threads; t++) lzr_cctx_free(wk[t]);
        free(wk); free(wout); free(nb);
        return LZR_ERR_NOMEM;
    }
    struct lzr_batch *bt = (struct lzr_batch *)calloc(1, sizeof(*bt));
    if (!bt) {
        for (int t = 0; t < c->threads; t++) lzr_cctx_free(wk[t]);
        free(wk); free(wout); free(nb);
        return LZR_ERR_NOMEM;
    }
#if LZR_HAVE_THREADS
    pthread_mutex_init(&bt->mu, NULL);
    pthread_cond_init(&bt->cv, NULL);
#endif
    c->bt = bt;
    c->wk = wk;
    c->wout = wout;
    c->jobs = K;
    free(c->buf);
    c->buf = nb;
    c->cap = cap;
    /* the main context's own match tables are not used in job mode */
    free(c->head); free(c->rtag); free(c->rhead);
    c->head = NULL; c->rtag = NULL; c->rhead = NULL;
    return LZR_OK;
}

/* Compress job j, buf[js, je), with worker w into wout[j]. */
static void job_compress(job_arg *a, lzr_cctx *w, lzr_buf *out)
{
    lzr_cctx *c = a->c;
    out->len = 0;
    w->buf = (uint8_t *)a->buf;   /* only read */
    /* forget everything from earlier jobs: positions below the new base
     * read as "slid out" (rows keep their slots, the search stops there) */
    w->base += a->span + 1;
    if (w->base + a->span >= LZR_REBASE_AT) {
        memset(w->head, 0, head_entries(&w->lp) * sizeof(uint32_t));
        w->base = 0;
    }
    w->end = a->je;
    w->ins = a->js > c->jover ? a->js - c->jover : 0;
    w->lo = w->ins;   /* the same bound wherever the job's bytes are */
    hash_cache_clear(w);
    if (w->lp.row) row_insert_upto(w, a->js);
    else if (w->lp.accel) fast1_fill(w, w->ins, a->js);
    else fast_fill(w, w->ins, a->js);
    w->reps[0] = w->reps[1] = w->reps[2] = 0;
    w->huf[0].valid = 0;
    w->fse[0].valid = w->fse[1].valid = w->fse[2].valid = 0;
    a->r = LZR_OK;
    for (size_t bs = a->js; bs < a->je && !a->r; bs += LZR_BLOCK_SIZE) {
        size_t be = a->je - bs > LZR_BLOCK_SIZE ? bs + LZR_BLOCK_SIZE : a->je;
        a->r = emit_block(w, bs, be, out);
    }
}

/* Take queued jobs in order until the batch is closed and drained. */
static void batch_work(struct lzr_batch *b, lzr_cctx *w, lzr_buf *outs)
{
#if LZR_HAVE_THREADS
    pthread_mutex_lock(&b->mu);
    for (;;) {
        while (b->next >= b->nready && !b->closed) pthread_cond_wait(&b->cv, &b->mu);
        if (b->next >= b->nready) break;
        int j = b->next++;
        pthread_mutex_unlock(&b->mu);
        job_compress(&b->a[j], w, &outs[j]);
        pthread_mutex_lock(&b->mu);
    }
    pthread_mutex_unlock(&b->mu);
#else
    while (b->next < b->nready) { int j = b->next++; job_compress(&b->a[j], w, &outs[j]); }
#endif
}

#if LZR_HAVE_THREADS
static void *batch_thread(void *p)
{
    struct { struct lzr_batch *b; lzr_cctx *w; lzr_buf *outs; } *t = p;
    batch_work(t->b, t->w, t->outs);
    return NULL;
}
#endif

static void batch_free(struct lzr_batch *b)
{
#if LZR_HAVE_THREADS
    if (b->nthr) {
        pthread_mutex_lock(&b->mu);
        b->closed = 1;
        pthread_cond_broadcast(&b->cv);
        pthread_mutex_unlock(&b->mu);
        for (int t = 1; t < b->nthr; t++) pthread_join(b->tid[t], NULL);
    }
    pthread_mutex_destroy(&b->mu);
    pthread_cond_destroy(&b->cv);
#endif
    free(b);
}

/* Queue job buf[js, je); the first job of a batch starts the helpers. */
static void batch_add_in(lzr_cctx *c, const uint8_t *buf, size_t span, size_t js, size_t je)
{
    struct lzr_batch *b = c->bt;
    job_arg *a = &b->a[b->nready];
    a->c = c; a->buf = buf; a->span = span; a->js = js; a->je = je; a->r = LZR_OK;
#if LZR_HAVE_THREADS
    pthread_mutex_lock(&b->mu);
    b->nready++;
    pthread_cond_signal(&b->cv);
    pthread_mutex_unlock(&b->mu);
    if (b->nthr == 0) {
        /* helpers wk[1..]: the caller is wk[0] and works once it closes */
        b->nthr = 1;
        for (int t = 1; t < c->threads; t++) {
            b->arg[t].b = b; b->arg[t].w = c->wk[t]; b->arg[t].outs = c->wout;
            if (pthread_create(&b->tid[t], NULL, batch_thread, &b->arg[t]) != 0) break;
            b->nthr = t + 1;
        }
    }
#else
    b->nready++;
#endif
}

static void batch_add(lzr_cctx *c, size_t js, size_t je) { batch_add_in(c, c->buf, c->cap, js, je); }

/* Queue every full job that is now buffered. */
static void batch_feed(lzr_cctx *c)
{
    struct lzr_batch *b = c->bt;
    while (b->nready < c->jobs) {
        size_t js = c->pos + (size_t)b->nready * c->jsize;
        if (c->end - js < c->jsize || c->end < js) break;
        batch_add(c, js, js + c->jsize);
    }
}

/* Queue the rest of buf[pos, end) (the last job may be short), run the
 * batch to the end, append the jobs in order and keep jover bytes of
 * history. */
static int batch_run(lzr_cctx *c, lzr_buf *out);
static int batch_finish(lzr_cctx *c, lzr_buf *out)
{
    struct lzr_batch *b = c->bt;
    batch_feed(c);
    size_t js = c->pos + (size_t)b->nready * c->jsize;
    if (b->nready < c->jobs && js < c->end) batch_add(c, js, c->end);
    int n = b->nready;
    int r = batch_run(c, out);
    if (r) return r;
    if (n) c->pos = b->a[n - 1].je;
    /* keep jover bytes of history */
    size_t keep_from = c->pos > c->jover ? c->pos - c->jover : 0;
    if (keep_from) {
        memmove(c->buf, c->buf + keep_from, c->end - keep_from);
        c->pos -= keep_from;
        c->end -= keep_from;
    }
    return LZR_OK;
}

/* Close the batch, work on it until it is done and append the jobs'
 * output in order. */
static int batch_run(lzr_cctx *c, lzr_buf *out)
{
    struct lzr_batch *b = c->bt;
    int n = b->nready;
#if LZR_HAVE_THREADS
    pthread_mutex_lock(&b->mu);
    b->closed = 1;
    pthread_cond_broadcast(&b->cv);
    pthread_mutex_unlock(&b->mu);
#endif
    batch_work(b, c->wk[0], c->wout);
#if LZR_HAVE_THREADS
    for (int t = 1; t < b->nthr; t++) pthread_join(b->tid[t], NULL);
#endif
    b->nready = b->next = b->closed = b->nthr = 0;
    int r = LZR_OK;
    for (int j = 0; j < n && !r; j++) {
        r = b->a[j].r;
        if (!r && buf_append(out, c->wout[j].p, c->wout[j].len)) r = LZR_ERR_NOMEM;
    }
    return r;
}

typedef struct { lzr_cctx *c; lzr_buf *out; int hash; } sink_ctx;   /* hash: src is the original content */

/* Append (already filtered) content and compress every full block. */
static int write_raw(void *ctx, const uint8_t *src, size_t n)
{
    lzr_cctx *c = ((sink_ctx *)ctx)->c;
    lzr_buf *out = ((sink_ctx *)ctx)->out;
    const int hash = ((sink_ctx *)ctx)->hash;
    if (c->jobs) {
        while (n) {
            /* copy at most one job at a time, so that each job is queued
             * (and compressed by the other threads) while the rest is copied */
            size_t take = c->cap - c->end;
            if (take > n) take = n;
            if (take > c->jsize) take = c->jsize;
            memcpy(c->buf + c->end, src, take);
            if (hash) xxh64_update(&c->xs, src, take);
            c->end += take;
            memset(c->buf + c->end, 0, 8);
            src += take;
            n -= take;
            batch_feed(c);
            if (c->end == c->cap) {   /* a full batch */
                int r = batch_finish(c, out);
                if (r) return r;
            }
        }
        return LZR_OK;
    }
    while (n) {
        if (c->end == c->cap) slide(c);
        size_t take = c->cap - c->end;
        if (take > n) take = n;
        memcpy(c->buf + c->end, src, take);
        if (hash) xxh64_update(&c->xs, src, take);
        c->end += take;
        memset(c->buf + c->end, 0, 8);   /* deterministic hashing past the end */
        src += take;
        n -= take;
        /* Compress full blocks, but keep LZR_LOOKAHEAD bytes after the block
         * so hashing and tree comparisons at its last positions see real
         * data (and the output does not depend on input chunking). With
         * match-finder threads, wait until the next block is buffered too,
         * so it can be gathered in the background (see opt_parse_block);
         * this changes only when blocks are written, not the output. */
        const size_t need = LZR_BLOCK_SIZE + LZR_LOOKAHEAD + (c->lp.opt && c->threads > 1 ? LZR_BLOCK_SIZE : 0);
        while (c->end - c->pos >= LZR_BLOCK_SIZE + LZR_LOOKAHEAD && (c->end - c->pos >= need || c->end == c->cap)) {
            int r = emit_block(c, c->pos, c->pos + LZR_BLOCK_SIZE, out);
            if (r) return r;
            c->pos += LZR_BLOCK_SIZE;
        }
    }
    return LZR_OK;
}

int lzr_cctx_set_filter(lzr_cctx *c, int filter)
{
    if (c->header_done || (filter != LZR_FILTER_NONE && filter != LZR_FILTER_X86)) return LZR_ERR_PARAM;
    free(c->x86tmp);
    c->x86tmp = NULL;
    if (filter == LZR_FILTER_X86) {
        c->x86tmp = (uint8_t *)malloc(X86_CHUNK + 4);
        if (!c->x86tmp) return LZR_ERR_NOMEM;
        memset(&c->x86, 0, sizeof(c->x86));
    }
    return LZR_OK;
}

/* One-shot input (the whole content in the first write, size known) in
 * job mode: the jobs read the caller's buffer instead of a copy, and the
 * caller hashes it while the helpers compress. Jobs are the same as when
 * copied (the same boundaries, history and bounds: see job_compress), so
 * the output is too. The last jobs, which read a little past their last
 * byte, are copied with their history into the job buffer and run in the
 * same final batch. Returns 1 if the input was taken (all of it). */
static int jobs_borrow(lzr_cctx *c, const uint8_t *src, size_t n, lzr_buf *out, int *err)
{
    size_t m = 0;
    while (m + c->jsize + 64 <= n) m += c->jsize;
    const size_t rest = n - m;   /* < 2 jobs + 64 */
    const int tail_jobs = (int)((rest + c->jsize - 1) / c->jsize);
    const size_t h = m < c->jover ? m : c->jover;
    if (!m || c->jobs < tail_jobs + 1 || h + rest > c->cap) return 0;
    int hashed = 0;
    size_t js = 0;
    for (;;) {
        /* borrowed jobs; the batch with the last of them has room for the
         * tail jobs too (c->jobs > tail_jobs, so a next batch always does) */
        const int R = (int)((m - js) / c->jsize), room = c->jobs;
        const int k = R + tail_jobs <= room ? R : (R - 1 < room ? R - 1 : room);
        for (int i = 0; i < k; i++) {
            batch_add_in(c, src, n, js, js + c->jsize);
            js += c->jsize;
        }
        if (js == m) {
            /* the tail: history and rest copied, as batch_finish keeps them */
            memcpy(c->buf, src + m - h, h + rest);
            c->pos = h;
            c->end = h + rest;
            memset(c->buf + c->end, 0, 8);
            for (size_t t = c->pos; t < c->end; t += c->jsize)
                batch_add(c, t, c->end - t > c->jsize ? t + c->jsize : c->end);
        }
        if (!hashed) { xxh64_update(&c->xs, src, n); hashed = 1; }
        int r = batch_run(c, out);
        if (r) { *err = r; return 0; }
        if (js == m) break;
    }
    c->pos = c->end;
    return 1;
}

int lzr_cctx_write(lzr_cctx *c, const uint8_t *src, size_t n, lzr_buf *out)
{
    if (c->finished) return LZR_ERR_PARAM;
    if (c->sg) return seg_write(c, src, n, out);
    if (!c->header_done && jobs_setup(c)) return LZR_ERR_NOMEM;
    if (!c->header_done && write_header(c, out)) return LZR_ERR_NOMEM;
    c->total_in += n;
    /* the checksum is of the original content: hashed here before the
     * filter, else by write_raw piece by piece as it is copied */
    sink_ctx s = {c, out, !c->x86tmp};
    if (c->jobs && !c->x86tmp && c->end == 0 && c->total_in == n && c->content_size == (int64_t)n &&
        n < ((size_t)1 << 30)) {
        int err = 0;
        if (jobs_borrow(c, src, n, out, &err) || err) return err;
    }
    if (c->x86tmp) {
        xxh64_update(&c->xs, src, n);
        return x86_run(&c->x86, src, n, 0, c->x86tmp, write_raw, &s);
    }
    return write_raw(&s, src, n);
}

int lzr_cctx_end(lzr_cctx *c, lzr_buf *out)
{
    if (c->finished) return LZR_ERR_PARAM;
    if (c->sg) return seg_end(c, out);
    if (!c->header_done && jobs_setup(c)) return LZR_ERR_NOMEM;
    if (!c->header_done && write_header(c, out)) return LZR_ERR_NOMEM;
    if (c->content_size >= 0 && (uint64_t)c->content_size != c->total_in) return LZR_ERR_PARAM;
    if (c->x86tmp) {
        sink_ctx s = {c, out, 0};
        int r = x86_finish(&c->x86, write_raw, &s);
        if (r) return r;
    }
    while (c->jobs && c->pos < c->end) {
        int r = batch_finish(c, out);
        if (r) return r;
    }
    while (c->pos < c->end) {
        size_t be = c->end - c->pos > LZR_BLOCK_SIZE ? c->pos + LZR_BLOCK_SIZE : c->end;
        int r = emit_block(c, c->pos, be, out);
        if (r) return r;
        c->pos = be;
    }
    uint8_t t[LZR_BLOCK_HDR + 4] = {2};
    wr32le(t + 5, 4);
    wr32le(t + LZR_BLOCK_HDR, (uint32_t)xxh64_digest(&c->xs));
    c->finished = 1;
    return buf_append(out, t, sizeof(t));
}

/* ---------------- segment mode: reset points + seek index ----------------
 * The content (after the x86 filter, if any) is cut at fixed offsets into
 * segments of sg->size bytes. Each is compressed by a fresh inner context
 * as a frame of its own; its blocks are copied into this frame and its
 * first block gets the reset bit (FORMAT.md 5.5). An inner frame starts in
 * the frame-start state, which is what a reset point restores, so the
 * blocks decode the same here. Segments do not depend on each other, so a
 * batch of them is compressed by parallel threads, and the output is the
 * same for any thread count. */
typedef struct { uint64_t foff, coff; } seg_ent;

struct lzr_seg {
    int level, wlog_req;
    uint64_t size;          /* segment size (reset interval)                 */
    int batch;              /* segments buffered before compressing (threads) */
    lzr_buf in;             /* pending content, < batch * size bytes         */
    seg_ent *ent;           /* seek index: one entry per segment             */
    size_t nent, cent;
    uint64_t foff, coff;    /* frame bytes written so far, content cut so far */
};

typedef struct {
    struct lzr_seg *sg;
    const uint8_t *src;     /* the batch's content                          */
    size_t len;             /* its length                                   */
    int k, mf;              /* segments, match-finder threads per segment   */
    lzr_buf *outs;          /* inner frames                                 */
    int *errs;
    int next;
#if LZR_HAVE_THREADS
    pthread_mutex_t mu;
#endif
} seg_job;

static int seg_one(seg_job *j, int i)
{
    const struct lzr_seg *sg = j->sg;
    size_t s0 = (size_t)i * sg->size, s1 = s0 + sg->size < j->len ? s0 + sg->size : j->len;
    lzr_cctx *w = lzr_cctx_new(sg->level, sg->wlog_req, (int64_t)(s1 - s0));
    if (!w) return LZR_ERR_NOMEM;
    lzr_cctx_set_threads(w, j->mf);
    lzr_buf *o = &j->outs[i];
    o->len = 0;
    int r = buf_reserve(o, (s1 - s0) / 2 + 64);
    if (!r) r = lzr_cctx_write(w, j->src + s0, s1 - s0, o);
    if (!r) r = lzr_cctx_end(w, o);
    lzr_cctx_free(w);
    return r;
}

static void *seg_thread(void *p)
{
    seg_job *j = (seg_job *)p;
    for (;;) {
        int i;
#if LZR_HAVE_THREADS
        pthread_mutex_lock(&j->mu);
        i = j->next++;
        pthread_mutex_unlock(&j->mu);
#else
        i = j->next++;
#endif
        if (i >= j->k) break;
        j->errs[i] = seg_one(j, i);
    }
    return NULL;
}

/* Compress src[0..len) (whole segments, the last one may be short) and
 * append the blocks to out. */
static int seg_batch(lzr_cctx *c, const uint8_t *src, size_t len, lzr_buf *out)
{
    struct lzr_seg *sg = c->sg;
    if (!len) return LZR_OK;
    seg_job j;
    memset(&j, 0, sizeof(j));
    j.sg = sg;
    j.src = src;
    j.len = len;
    j.k = (int)((len + sg->size - 1) / sg->size);
    int nt = sg->batch < j.k ? sg->batch : j.k;
    if (nt < 1) nt = 1;
    j.mf = c->threads / nt > 1 ? c->threads / nt : 1;   /* spare threads go to the match finders */
    j.outs = (lzr_buf *)calloc((size_t)j.k, sizeof(lzr_buf));
    j.errs = (int *)calloc((size_t)j.k, sizeof(int));
    if (!j.outs || !j.errs) { free(j.outs); free(j.errs); return LZR_ERR_NOMEM; }
#if LZR_HAVE_THREADS
    pthread_t tid[LZR_MAX_THREADS];
    int started = 0;
    pthread_mutex_init(&j.mu, NULL);
    for (int t = 1; t < nt; t++)
        if (pthread_create(&tid[started], NULL, seg_thread, &j) == 0) started++;
    seg_thread(&j);
    for (int t = 0; t < started; t++) pthread_join(tid[t], NULL);
    pthread_mutex_destroy(&j.mu);
#else
    seg_thread(&j);
#endif
    int r = LZR_OK;
    for (int i = 0; i < j.k && !r; i++) {
        r = j.errs[i];
        if (r) break;
        lzr_buf *o = &j.outs[i];
        int hl = lzr_frame_header(o->p, o->len, NULL, NULL, NULL);
        size_t body = o->len - (size_t)hl - (LZR_BLOCK_HDR + 4);   /* without header and end block */
        if (hl <= 0 || o->len < (size_t)hl + LZR_BLOCK_HDR + 4) { r = LZR_ERR_PARAM; break; }
        if (sg->nent == sg->cent) {
            size_t nc = sg->cent ? 2 * sg->cent : 16;
            seg_ent *ne = (seg_ent *)realloc(sg->ent, nc * sizeof(seg_ent));
            if (!ne) { r = LZR_ERR_NOMEM; break; }
            sg->ent = ne;
            sg->cent = nc;
        }
        sg->ent[sg->nent].foff = sg->foff;
        sg->ent[sg->nent].coff = sg->coff;
        if (sg->nent && body) o->p[hl] |= LZR_BT_RESET;   /* segment start = reset point */
        sg->nent++;
        r = buf_append(out, o->p + hl, body);
        sg->foff += body;
        size_t s0 = (size_t)i * sg->size;
        sg->coff += (s0 + sg->size < len ? sg->size : len - s0);
    }
    for (int i = 0; i < j.k; i++) free(j.outs[i].p);
    free(j.outs);
    free(j.errs);
    return r;
}

/* Levels 9-11 (32/64 MB windows, the slowest to decode): 48 MiB segments.
 * Measured ratio cost vs one segment: Silesia +0.28..0.38%, enwik8
 * +1.0..2.6% (long-range repeats); 32 MiB cost 2-3x more. Other levels
 * decode fast enough on one thread: no segments unless asked for. */
uint64_t lzr_default_segment(int level)
{
    return level >= 9 ? (uint64_t)48 << 20 : 0;
}

/* Memory of one inner context for a segment (big tables only). */
static uint64_t seg_ctx_mem(const level_params *lp)
{
    uint64_t w = (uint64_t)1 << lp->wlog, h = (uint64_t)1 << lp->hlog;
    uint64_t m = 2 * w + 4 * h;                       /* window buffer + heads */
    if (lp->opt && !lp->row) m += 8 * w;              /* binary tree           */
    else m += h;                                      /* row tags              */
    return m;
}
#ifndef LZR_SEG_MEM_BUDGET
#define LZR_SEG_MEM_BUDGET ((uint64_t)1 << 30)   /* inner contexts alive at once */
#endif

lzr_cctx *lzr_cctx_new_seg(int level, int window_log, int64_t content_size, uint64_t seg_size)
{
    if (!seg_size || (content_size >= 0 && (uint64_t)content_size <= seg_size)) return lzr_cctx_new(level, window_log, content_size);
    if (level < LZR_LEVEL_MIN || level > LZR_LEVEL_MAX || seg_size < LZR_SEG_MIN || seg_size > ((uint64_t)1 << 40)) return NULL;
    if (window_log && (window_log < LZR_WLOG_MIN || window_log > LZR_WLOG_MAX)) return NULL;
    lzr_cctx *c = (lzr_cctx *)calloc(1, sizeof(*c));
    struct lzr_seg *sg = (struct lzr_seg *)calloc(1, sizeof(*sg));
    if (!c || !sg) { free(c); free(sg); return NULL; }
    xxh64_init(&c->xs);
    c->sg = sg;
    c->threads = 1;
    c->content_size = content_size;
    c->lp = LEVELS[level - LZR_LEVEL_MIN];
    if (window_log) c->lp.wlog = (uint8_t)window_log;
    /* the window every inner frame fits in: the one for a whole segment */
    c->lp.wlog = (uint8_t)fit_wlog(&c->lp, (int64_t)seg_size);
    sg->level = level;
    sg->wlog_req = window_log;
    sg->size = seg_size;
    return c;
}

static void seg_free(lzr_cctx *c)
{
    free(c->sg->in.p);
    free(c->sg->ent);
    free(c->sg);
    free(c->x86tmp);
    free(c);
}

/* Filtered content in: compress whenever a whole batch is available. */
static int seg_raw(void *ctx, const uint8_t *src, size_t n)
{
    lzr_cctx *c = ((sink_ctx *)ctx)->c;
    lzr_buf *out = ((sink_ctx *)ctx)->out;
    struct lzr_seg *sg = c->sg;
    if (!sg->batch) {   /* segments compressed at once: threads, within the memory budget */
        uint64_t cap = LZR_SEG_MEM_BUDGET / seg_ctx_mem(&c->lp);
        sg->batch = c->threads < 1 ? 1 : c->threads;
        if ((uint64_t)sg->batch > cap) sg->batch = cap ? (int)cap : 1;
    }
    const size_t full = (size_t)sg->size * (size_t)sg->batch;
    while (n) {
        if (!sg->in.len && n >= full) {   /* straight from the caller's buffer */
            int r = seg_batch(c, src, full, out);
            if (r) return r;
            src += full;
            n -= full;
            continue;
        }
        size_t take = full - sg->in.len < n ? full - sg->in.len : n;
        if (buf_append(&sg->in, src, take)) return LZR_ERR_NOMEM;
        src += take;
        n -= take;
        if (sg->in.len == full) {
            int r = seg_batch(c, sg->in.p, full, out);
            sg->in.len = 0;
            if (r) return r;
        }
    }
    return LZR_OK;
}

static int seg_write(lzr_cctx *c, const uint8_t *src, size_t n, lzr_buf *out)
{
    if (!c->header_done) {
        size_t before = out->len;
        if (write_header(c, out)) return LZR_ERR_NOMEM;
        c->sg->foff = out->len - before;
    }
    if (!n) return LZR_OK;
    c->total_in += n;
    xxh64_update(&c->xs, src, n);
    sink_ctx s = {c, out, 0};
    return c->x86tmp ? x86_run(&c->x86, src, n, 0, c->x86tmp, seg_raw, &s) : seg_raw(&s, src, n);
}

/* Seek index payload (FORMAT.md 12.3): per segment u64 frame offset of its
 * first block, u64 content offset; then u64 content size, u32 count, "LZRX". */
static int seg_end(lzr_cctx *c, lzr_buf *out)
{
    struct lzr_seg *sg = c->sg;
    int r = seg_write(c, NULL, 0, out);
    if (r) return r;
    if (c->content_size >= 0 && (uint64_t)c->content_size != c->total_in) return LZR_ERR_PARAM;
    if (c->x86tmp) {
        sink_ctx s = {c, out, 0};
        if ((r = x86_finish(&c->x86, seg_raw, &s))) return r;
    }
    if ((r = seg_batch(c, sg->in.p, sg->in.len, out))) return r;
    sg->in.len = 0;
    if (sg->nent >= 2) {
        size_t pl = 16 * sg->nent + 16;
        uint8_t h[LZR_BLOCK_HDR] = {0x80};
        wr32le(h + 5, (uint32_t)pl);
        if (buf_reserve(out, LZR_BLOCK_HDR + pl)) return LZR_ERR_NOMEM;
        buf_append(out, h, LZR_BLOCK_HDR);
        uint8_t *q = out->p + out->len;
        for (size_t i = 0; i < sg->nent; i++, q += 16) {
            wr64le(q, sg->ent[i].foff);
            wr64le(q + 8, sg->ent[i].coff);
        }
        wr64le(q, sg->coff);
        wr32le(q + 8, (uint32_t)sg->nent);
        memcpy(q + 12, "LZRX", 4);
        out->len += pl;
    }
    uint8_t t[LZR_BLOCK_HDR + 4] = {2};
    wr32le(t + 5, 4);
    wr32le(t + LZR_BLOCK_HDR, (uint32_t)xxh64_digest(&c->xs));
    c->finished = 1;
    return buf_append(out, t, sizeof(t));
}
