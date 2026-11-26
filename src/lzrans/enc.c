/*
 * enc.c - streaming encoder.
 *
 * Structure:
 *   lzr_cctx          sliding input buffer (2*window + block), hash-chain
 *                     tables, repeat offsets, running XXH64
 *   find_best()       repeat-offset probes + hash-chain walk, scored by gain
 *   parse_block()     lazy parse (depth 0..2) into sequences + literals
 *   emit_block()      sequences -> block type 5 (Huffman literals + tANS
 *                     sequences) or type 4 (range-coded tokens, level 9)
 *   slide()           drops history older than the window, rebases tables
 */
#include "lzr.h"
#include "codes.h"
#include "xxh64.h"
#include "x86.h"
#include "arc.h"

#include "enc_internal.h"
#if defined(__SSE2__) && !defined(LZR_NO_SIMD)
#include <emmintrin.h>
#endif
#ifdef LZR_TUNE
#include <stdio.h>
#endif

/* Levels 1-6: lazy parsing (window, hash bits, lazy depth, search steps,
 * nice length) over rows of 2^row recent positions (hash bits = log2 of
 * all row entries; replaced hash chains in v0.6: 2-5x faster). Levels
 * 7-9: binary tree + optimal parsing (chain = tree depth, opt = number of
 * parsing passes; a second pass gained <0.5% ratio for ~40% more time, so
 * all use one). */
static const level_params LEVELS[LZR_LEVEL_MAX + 1] = {
    {0, 0, 0, 0, 0, 0, 0},
    {20, 20, 0, 4, 16, 0, 4}, {21, 21, 1, 8, 24, 0, 4}, {22, 21, 1, 16, 32, 0, 4},
    {22, 22, 1, 32, 48, 0, 5}, {23, 23, 2, 48, 64, 0, 6}, {24, 23, 2, 64, 96, 0, 6},
    {23, 20, 0, 16, 128, 1, 0}, {24, 21, 0, 24, 273, 1, 0}, {24, 22, 0, 48, 273, 1, 0},
};


lzr_cctx *lzr_cctx_new(int level, int window_log, int64_t content_size)
{
    if (level < LZR_LEVEL_MIN || level > LZR_LEVEL_MAX) return NULL;
    lzr_cctx *c = (lzr_cctx *)calloc(1, sizeof(*c));
    if (c) xxh64_init(&c->xs);
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
    if (getenv("LZR_HLOG")) c->lp.hlog = (uint8_t)atoi(getenv("LZR_HLOG"));
    if (getenv("LZR_ROW")) c->lp.row = (uint8_t)atoi(getenv("LZR_ROW"));
    if (getenv("LZR_LAZY")) c->lp.lazy = (uint8_t)atoi(getenv("LZR_LAZY"));
    if (getenv("LZR_THREADS")) c->threads = atoi(getenv("LZR_THREADS"));
#endif
    /* A window larger than the known input only wastes memory. */
    {
        /* the tree finder needs >= 2 blocks of slack, see opt_init() */
        unsigned min_wlog = c->lp.opt ? 20 : LZR_WLOG_MIN;
        if (c->lp.wlog < min_wlog) c->lp.wlog = (uint8_t)min_wlog;
        while (content_size >= 0 && c->lp.wlog > min_wlog && ((uint64_t)1 << (c->lp.wlog - 1)) >= (uint64_t)content_size)
            c->lp.wlog--;
    }
    if (c->lp.hlog > c->lp.wlog) c->lp.hlog = c->lp.wlog;
    c->wsize = (size_t)1 << c->lp.wlog;
    c->wmask = c->wsize - 1;
    c->cap = 2 * c->wsize + LZR_BLOCK_SIZE;
    c->content_size = content_size;
    c->reps[0] = 1; c->reps[1] = 4; c->reps[2] = 8;

    c->buf = (uint8_t *)big_alloc(c->cap + LZR_SLACK, 0);
    c->head = (uint32_t *)big_alloc(((size_t)1 << c->lp.hlog) * sizeof(uint32_t), 1);
    if (!c->lp.opt) {   /* row finder (lazy levels) */
        c->rtag = (uint8_t *)big_alloc((size_t)1 << c->lp.hlog, 1);
        c->rhead = (uint8_t *)big_alloc((size_t)1 << (c->lp.hlog - c->lp.row), 1);
    }
    c->seq_ll = (uint32_t *)malloc(MAX_SEQ * sizeof(uint32_t));
    c->seq_ml = (uint32_t *)malloc(MAX_SEQ * sizeof(uint32_t));
    c->seq_dist = (uint32_t *)malloc(MAX_SEQ * sizeof(uint32_t));
    c->lits = (uint8_t *)malloc(LZR_BLOCK_SIZE);
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
    if (!c->buf || !c->head || (!c->lp.opt && (!c->rtag || !c->rhead)) || !c->seq_ll || !c->seq_ml || !c->seq_dist || !c->lits ||
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

void lzr_cctx_free(lzr_cctx *c)
{
    if (!c) return;
    opt_free(c);
    free(c->buf); free(c->head); free(c->rtag); free(c->rhead);
    free(c->seq_ll); free(c->seq_ml); free(c->seq_dist);
    free(c->lits); free(c->llc); free(c->mlc); free(c->ofc); free(c->xlm); free(c->xof); free(c->nlm); free(c->nof);
    free(c->blk); free(c->tm); free(c->huf); free(c->fse); free(c->x86tmp);
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

static inline void row_insert_upto(lzr_cctx *c, size_t upto)
{
    const unsigned rlog = c->lp.row, bits = c->lp.hlog - rlog + 8, rmask = (1u << rlog) - 1;
    const size_t base = c->base;
    size_t limit = c->end >= LZR_HASH_LEN ? c->end - LZR_HASH_LEN + 1 : 0;
    if (upto > limit) upto = limit;
    for (size_t i = c->ins; i < upto; i++) {
        {   /* the row of a position a few bytes ahead: hides the miss */
            uint32_t ha = row_hash(c->buf + i + 8, bits) >> 8;
            __builtin_prefetch(c->rtag + ((size_t)ha << rlog));
            __builtin_prefetch(c->head + ((size_t)ha << rlog));
        }
        uint32_t h = row_hash(c->buf + i, bits), row = h >> 8;
        unsigned slot = (c->rhead[row] - 1) & rmask;
        c->rhead[row] = (uint8_t)slot;
        size_t e = ((size_t)row << rlog) + slot;
        c->rtag[e] = (uint8_t)h;
        c->head[e] = (uint32_t)(base + i + 1);
    }
    if (upto > c->ins) c->ins = upto;
}

/* Insert positions [c->ins, upto). */
static inline void insert_upto(lzr_cctx *c, size_t upto) { row_insert_upto(c, upto); }

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

    /* 2) row finder: tag-matching slots, newest first */
    const unsigned rlog = c->lp.row, R = 1u << rlog;
    const uint32_t h = row_hash(cur, c->lp.hlog - rlog + 8), row = h >> 8;
    const size_t r0 = (size_t)row << rlog;
    const unsigned hd = c->rhead[row];
    uint64_t m = row_match_mask(c->rtag + r0, (uint8_t)h, R);
    /* rotate so that bit k is the k-th newest slot */
    if (hd) m = (m >> hd | m << (R - hd)) & (R == 64 ? ~0ull : (1ull << R) - 1);
    unsigned tries = c->lp.chain;
    size_t best_len = best.len >= LZR_MIN_MATCH ? best.len : LZR_MIN_MATCH - 1;
    const size_t base = c->base;
    while (m && tries--) {
        unsigned k = (unsigned)__builtin_ctzll(m);
        m &= m - 1;
        uint32_t cand = c->head[r0 + ((hd + k) & (R - 1))];
        if (!cand || cand - 1 < base) break;          /* empty / slid out (older follow) */
        size_t cp = cand - 1 - base;
        if (cp >= p) continue;                          /* not yet valid (cannot happen) */
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
static size_t emit_seqs5(lzr_cctx *c, size_t ns, size_t nlit, uint8_t *b)
{
    size_t k = LZR_BLOCK_HDR;
    wr32le(b + k, (uint32_t)ns); k += 4;
    wr32le(b + k, (uint32_t)nlit); k += 4;
    k += huf_encode(b + k, c->lits, nlit, &c->huf[0]);
    if (!ns) return k;

    uint8_t *codes[3] = {c->llc, c->mlc, c->ofc};
    const unsigned alpha[3] = {VAL_CODES, ML_CODES, OF_CODES};
    uint32_t cnt[3][FSE_MAXSYM];
    memset(cnt, 0, sizeof(cnt));
    /* codes, and the extra bits in the order they are written: ll and ml
     * together (<= 16 + 16 bits, block size 2^18), then the offset */
    for (size_t i = 0; i < ns; i++) {
        unsigned nl, nm, no = 0;
        uint32_t el, em, eo = 0;
        c->llc[i] = (uint8_t)val_code(c->seq_ll[i], &nl, &el);
        c->mlc[i] = (uint8_t)ml_code(c->seq_ml[i] - LZR_ML_BASE, &nm, &em);
        uint32_t off = c->seq_dist[i];
        c->ofc[i] = (uint8_t)(off < LZR_NUM_REPS ? off : LZR_NUM_REPS + val_code(off - 3, &no, &eo));
        c->xlm[i] = el | em << nl;
        c->nlm[i] = (uint8_t)(nl + nm);
        c->xof[i] = eo;
        c->nof[i] = (uint8_t)no;
        cnt[0][c->llc[i]]++; cnt[1][c->mlc[i]]++; cnt[2][c->ofc[i]]++;
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
    /* bit stream, backwards over the sequences (see dec.c for the order) */
    fbw w = {0, 0, b + k};
    uint32_t x[3] = {FSE_L, FSE_L, FSE_L};
    if (!rle[0] && !rle[1] && !rle[2]) {
        /* common case: the three state transitions in one write (<= 27 bits) */
        const fse_table *T0 = &c->fse[0], *T1 = &c->fse[1], *T2 = &c->fse[2];
        uint32_t x0 = x[0], x1 = x[1], x2 = x[2];
        for (size_t i = ns; i-- > 0;) {
            unsigned n2, n1, n0;
            uint32_t b2 = fse_step(T2, &x2, c->ofc[i], &n2);
            uint32_t b1 = fse_step(T1, &x1, c->mlc[i], &n1);
            uint32_t b0 = fse_step(T0, &x0, c->llc[i], &n0);
            fbw_put(&w, b2 | b1 << n2 | (uint64_t)b0 << (n2 + n1), n2 + n1 + n0);
            fbw_put(&w, c->xlm[i], c->nlm[i]);
            fbw_put(&w, c->xof[i], c->nof[i]);
        }
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
    return (size_t)(fbw_close(&w) - b);
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
    } else {
        k = emit_seqs5(c, ns, nlit, b);
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
 * tree indices (pos & wmask, levels 7-9) stay valid. */
static void slide(lzr_cctx *c)
{
    if (c->pos <= c->wsize) return;
    size_t delta = ((c->pos - c->wsize) / c->wsize) * c->wsize;
    if (!delta) return;
    memmove(c->buf, c->buf + delta, c->end - delta);
    c->pos -= delta;
    c->end -= delta;
    c->ins = c->ins > delta ? c->ins - delta : 0;
    size_t hs = (size_t)1 << c->lp.hlog;
    if (!c->lp.opt) {
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

typedef struct { lzr_cctx *c; lzr_buf *out; } sink_ctx;

/* Append (already filtered) content and compress every full block. */
static int write_raw(void *ctx, const uint8_t *src, size_t n)
{
    lzr_cctx *c = ((sink_ctx *)ctx)->c;
    lzr_buf *out = ((sink_ctx *)ctx)->out;
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

int lzr_cctx_write(lzr_cctx *c, const uint8_t *src, size_t n, lzr_buf *out)
{
    if (c->finished) return LZR_ERR_PARAM;
    if (!c->header_done && write_header(c, out)) return LZR_ERR_NOMEM;
    xxh64_update(&c->xs, src, n);   /* of the original content */
    c->total_in += n;
    sink_ctx s = {c, out};
    if (c->x86tmp) return x86_run(&c->x86, src, n, 0, c->x86tmp, write_raw, &s);
    return write_raw(&s, src, n);
}

int lzr_cctx_end(lzr_cctx *c, lzr_buf *out)
{
    if (c->finished) return LZR_ERR_PARAM;
    if (!c->header_done && write_header(c, out)) return LZR_ERR_NOMEM;
    if (c->content_size >= 0 && (uint64_t)c->content_size != c->total_in) return LZR_ERR_PARAM;
    if (c->x86tmp) {
        sink_ctx s = {c, out};
        int r = x86_finish(&c->x86, write_raw, &s);
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
