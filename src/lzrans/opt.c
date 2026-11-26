/*
 * opt.c - binary-tree match finder + price-based optimal parser
 * (levels 7-9).
 *
 * Per 256 KB block:
 *   1. gather()        For every position, walk/insert the binary tree
 *                      (BT4: 4-byte hash heads, son[] holds left/right
 *                      children in a window-sized ring) and cache all
 *                      matches of increasing length. Inside a match of
 *                      length >= nice the positions are only inserted.
 *   2. parse_pass()    Forward dynamic programming over segments of up to
 *                      OPT_NUM positions. Each node keeps the cheapest
 *                      price to reach it, how it was reached (literal or
 *                      match length + offset) and the repeat offsets
 *                      (and token coder state) on that path. Block type 4
 *                      (levels 7-9 by default): prices from a shadow copy
 *                      of the adaptive token model, advanced over the
 *                      tokens chosen so far before each segment. Block
 *                      type 1: prices from the symbol statistics of the
 *                      previous pass/block (log2 of frequencies).
 *   3. With 2 passes (type 1), the first pass only gathers statistics for
 *      the second, which produces the real sequences.
 */
#include "enc_internal.h"
#include <math.h>
#if LZR_HAVE_THREADS
#include <pthread.h>
#endif

#ifndef OPT_NUM
#define OPT_NUM     4096
#endif
#define MAX_CAND    48
#define NODE_EXTRA  (600 + 16)    /* match + literal + rep0 reaches 2 * 273 + 1 past */
#define BITCOST     256                 /* price units per bit */
#define PRICE_INF   0x3FFFFFFF

typedef struct { uint32_t len, dist; } mpair;

typedef struct { uint32_t pos, len, dist; } mrec;

typedef struct {
    lzr_cctx *c;
    size_t from, bs, be, lim_end;
    unsigned part, nparts;
    mrec *out;
    size_t n, cap;
    int err;
} gjob;


typedef struct {
    int32_t price;
    uint32_t from, ml, off, litlen;
    uint32_t ml2;          /* != 0: match, 1 literal, then rep0 of this length */
    uint32_t reps[LZR_NUM_REPS];
    uint32_t state;        /* token model state (block type 4) on this path */
} onode;

typedef struct {
    uint32_t lit[256], llc[VAL_CODES], mlc[VAL_CODES], ofc[OF_CODES];
} prices_t;

#ifndef LZR_OPT_FULL
#define LZR_OPT_FULL 273
#endif
#ifndef SEG_MIN
#define SEG_MIN 512            /* shortest segment ended early (type 4) */
#endif
#ifndef LZR_REP_MIN
#define LZR_REP_MIN 2        /* repeat-offset matches may be this short */
#endif
#define ML_TABLE 512                     /* > max nice length */
#define LL_TABLE 4096

/* Match cache for one block (two of them, so the next block can be
 * gathered in the background while the current one is parsed). */
typedef struct {
    uint32_t *mstart;      /* [B+1] index of the first cached pair per position */
    mpair *pairs;
    size_t pcap;
    gjob jobs[LZR_MAX_THREADS];
} mcache;

typedef struct {
    mcache mc[2];
    int cur;               /* cache holding the block being parsed           */
#if LZR_HAVE_THREADS
    pthread_t pend_tid;    /* background gather of the next block            */
#endif
    int pending, pend_err;
    size_t pend_bs, pend_be, pend_lim;
    lzr_cctx *pend_c;
    onode *nodes;
    uint32_t (*path)[3];
    prices_t pr;
    int32_t ml_p[ML_TABLE];   /* match-length price incl. extra bits, per pass */
    int32_t ll_p[LL_TABLE];   /* literal-run price incl. extra bits           */
    /* block type 4: prices from a shadow copy of the token model, advanced
     * over the tokens chosen so far and re-tabulated after each segment */
    uint16_t bitc[ARC_ONE >> 4];   /* price of a 0 bit by probability / 16 */
    tok_model pm;                  /* shadow model                          */
    uint32_t dreps[LZR_NUM_REPS];  /* repeat offsets at dpos                 */
    size_t dpos, dseq;             /* shadow position / next sequence        */
    size_t dlit;                   /* start of that sequence's literal run   */
    int32_t litp[8][256];          /* plain literal, per context             */
    int32_t lenp[2][TOK_POS][TOK_MAXLEN + 1];   /* [0] new dist, [1] repeat  */
    int32_t slotp[4][TOK_DIST_SLOTS];
    int32_t distp[4][128];         /* slot + footer, dist0 < 128             */
    int32_t alignp[16];
    /* short new-distance matches (2-3 bytes, block type 4): last position
     * (+1) of each 2-byte value / 3-byte hash, and per block position the
     * nearest earlier one (0 = none) */
    uint32_t *last2, *last3, *sh2, *sh3;
} opt_state;
#define SH3_LOG 16
#ifndef SH2_MAXD
#define SH2_MAXD 256
#endif
#ifndef SH3_MAXD
#define SH3_MAXD 65536
#endif

int opt_init(lzr_cctx *c)
{
    opt_state *st = (opt_state *)calloc(1, sizeof(opt_state));
    c->opt = st;
    c->son = (uint32_t *)big_alloc(2 * c->wsize * sizeof(uint32_t), 1);
    /* son[] is a ring indexed by position & (wsize-1). Worker threads run
     * at different positions within a block (and one may gather the next
     * block), so a slot can be reused by position cm + wsize while another
     * thread still sees cm as in range. Treating nodes older than
     * wsize - 2 blocks as out of range makes that impossible, and makes the
     * result identical for every thread count. */
    c->bt_maxdist = c->wsize - 2 * (size_t)LZR_BLOCK_SIZE;
    if (!st || !c->son) return LZR_ERR_NOMEM;
    for (int k = 0; k < 2; k++) {
        mcache *m = &st->mc[k];
        m->pcap = (size_t)LZR_BLOCK_SIZE * 2;
        m->mstart = (uint32_t *)malloc((LZR_BLOCK_SIZE + 1) * sizeof(uint32_t));
        m->pairs = (mpair *)malloc(m->pcap * sizeof(mpair));
        if (!m->mstart || !m->pairs) return LZR_ERR_NOMEM;
    }
    st->nodes = (onode *)malloc((OPT_NUM + NODE_EXTRA) * sizeof(onode));
    st->path = malloc((OPT_NUM + NODE_EXTRA) * sizeof(*st->path));
    if (!st->nodes || !st->path) return LZR_ERR_NOMEM;
    st->last2 = (uint32_t *)calloc(65536, sizeof(uint32_t));
    st->last3 = (uint32_t *)calloc((size_t)1 << SH3_LOG, sizeof(uint32_t));
    st->sh2 = (uint32_t *)malloc(LZR_BLOCK_SIZE * sizeof(uint32_t));
    st->sh3 = (uint32_t *)malloc(LZR_BLOCK_SIZE * sizeof(uint32_t));
    if (!st->last2 || !st->last3 || !st->sh2 || !st->sh3) return LZR_ERR_NOMEM;
    for (unsigned i = 0; i < (ARC_ONE >> 4); i++)
        st->bitc[i] = (uint16_t)(-log2((i * 16 + 8) / (double)ARC_ONE) * BITCOST + 0.5);
    return LZR_OK;
}

static inline int32_t arc_bit_price(const opt_state *st, arc_p p, unsigned bit)
{
    return st->bitc[(bit ? ARC_ONE - p : p) >> 4];
}

/* out[sym] = price of sym in an nbits bit tree (reverse: low bit first). */
static void tree_prices(const opt_state *st, const arc_p *probs, unsigned nbits, int rev, int32_t *out)
{
    for (unsigned sym = 0; sym < (1u << nbits); sym++) {
        int32_t pr = 0;
        unsigned m = 1;
        for (unsigned i = 0; i < nbits; i++) {
            unsigned b = rev ? (sym >> i) & 1 : (sym >> (nbits - 1 - i)) & 1;
            pr += arc_bit_price(st, probs[m], b);
            m = (m << 1) | b;
        }
        out[sym] = pr;
    }
}

static void len_prices(const opt_state *st, const tok_len *l, int32_t (*out)[TOK_MAXLEN + 1])
{
    int32_t lo[8], mi[8], hi[256];
    int32_t c0 = arc_bit_price(st, l->choice, 0), c1 = arc_bit_price(st, l->choice, 1);
    int32_t d0 = arc_bit_price(st, l->choice2, 0), d1 = arc_bit_price(st, l->choice2, 1);
    tree_prices(st, l->high, 8, 0, hi);
    for (unsigned ps = 0; ps < TOK_POS; ps++) {
        tree_prices(st, l->low[ps], 3, 0, lo);
        tree_prices(st, l->mid[ps], 3, 0, mi);
        out[ps][0] = out[ps][1] = 0;
        for (unsigned L = 2; L <= TOK_MAXLEN; L++) {
            unsigned v = L - 2;
            out[ps][L] = v < 8 ? c0 + lo[v] : v < 16 ? c1 + d0 + mi[v - 8] : c1 + d1 + hi[v - 16];
        }
    }
}

/* Re-tabulate the prices from the shadow model. */
static void tok_prices(opt_state *st)
{
    const tok_model *m = &st->pm;
    for (unsigned c = 0; c < 8; c++) {
        int32_t node[256];
        const arc_p *pr = m->lit[c];
        node[1] = 0;
        for (unsigned k = 1; k < 128; k++) {
            node[2 * k] = node[k] + arc_bit_price(st, pr[k], 0);
            node[2 * k + 1] = node[k] + arc_bit_price(st, pr[k], 1);
        }
        for (unsigned b = 0; b < 256; b++) st->litp[c][b] = node[128 + (b >> 1)] + arc_bit_price(st, pr[128 + (b >> 1)], b & 1);
    }
    len_prices(st, &m->mlen, st->lenp[0]);
    len_prices(st, &m->rlen, st->lenp[1]);
    for (unsigned c = 0; c < 4; c++) tree_prices(st, m->slot[c], 6, 0, st->slotp[c]);
    tree_prices(st, m->align, 4, 1, st->alignp);
    int32_t spec[32];
    for (unsigned slot = 0; slot < TOK_END_SPEC; slot++) {
        unsigned nb = slot < 4 ? 0 : (slot >> 1) - 1;
        uint32_t base = slot < 4 ? slot : (2u | (slot & 1)) << nb;
        if (nb) tree_prices(st, m->spec[slot], nb, 1, spec);
        else spec[0] = 0;
        for (uint32_t r = 0; r < (1u << nb); r++)
            for (unsigned c = 0; c < 4; c++) st->distp[c][base + r] = st->slotp[c][slot] + spec[r];
    }
}

static inline int32_t dist_price(const opt_state *st, uint32_t dist0, unsigned lctx)
{
    if (dist0 < 128) return st->distp[lctx][dist0];
    unsigned slot = tok_slot(dist0), nb = (slot >> 1) - 1;
    return st->slotp[lctx][slot] + (int32_t)(nb - 4) * BITCOST + st->alignp[dist0 & 15];
}

static int32_t matched_price(const opt_state *st, const arc_p *pr, unsigned byte, unsigned mbyte)
{
    int32_t price = 0;
    unsigned m = 1, same = 1;
    for (int i = 7; i >= 0; i--) {
        unsigned b = (byte >> i) & 1, mb = (mbyte >> i) & 1;
        price += arc_bit_price(st, same ? pr[256 + (mb << 8) + m] : pr[m], b);
        m = (m << 1) | b;
        same &= (b == mb);
    }
    return price;
}

/* Advance the shadow model over the sequences emitted so far and the
 * literals up to buffer position upto (as encode_tokens() will code them). */
static void shadow_advance(lzr_cctx *c, opt_state *st, size_t bs, size_t ns, size_t upto)
{
    arc_enc e;
    arce_init(&e, NULL);
    e.dry = 1;
    tok_model *m = &st->pm;
    const uint8_t *buf = c->buf;
    size_t p = st->dpos;
    for (; st->dseq < ns; st->dseq++) {
        size_t i = st->dseq, mpos = st->dlit + c->seq_ll[i];   /* literals before dpos are done */
        for (; p < mpos; p++)
            tok_enc_lit(&e, m, buf + p, p ? buf[p - 1] : 0, st->dreps[0], (unsigned)(p - bs) & 3);
        uint32_t off = c->seq_dist[i];
        commit_offset(st->dreps, off);
        size_t ml = c->seq_ml[i];
        while (ml) {
            unsigned len = tok_piece(ml);
            tok_enc_match(&e, m, off, len, (unsigned)(p - bs) & 3);
            off = 0;
            p += len;
            ml -= len;
        }
        st->dlit = p;
    }
    for (; p < upto; p++) tok_enc_lit(&e, m, buf + p, p ? buf[p - 1] : 0, st->dreps[0], (unsigned)(p - bs) & 3);
    st->dpos = p;
}

void opt_free(lzr_cctx *c)
{
    opt_state *st = (opt_state *)c->opt;
    opt_sync(c);
    if (st) {
        for (int k = 0; k < 2; k++) {
            for (unsigned t = 0; t < LZR_MAX_THREADS; t++) free(st->mc[k].jobs[t].out);
            free(st->mc[k].mstart);
            free(st->mc[k].pairs);
        }
        free(st->nodes); free(st->path);
        free(st->last2); free(st->last3); free(st->sh2); free(st->sh3);
        free(st);
    }
    free(c->son);
    c->opt = NULL;
    c->son = NULL;
}

void opt_slide(lzr_cctx *c, uint32_t d)
{
    opt_state *st = (opt_state *)c->opt;
    opt_sync(c);                     /* the background gather reads buf/tree */
    if (st->pending) { st->pend_bs -= d; st->pend_be -= d; st->pend_lim -= d; }
    for (size_t i = 0; i < 2 * c->wsize; i++) c->son[i] = c->son[i] > d ? c->son[i] - d : 0;
    for (size_t i = 0; i < 65536; i++) st->last2[i] = st->last2[i] > d ? st->last2[i] - d : 0;
    for (size_t i = 0; i < ((size_t)1 << SH3_LOG); i++) st->last3[i] = st->last3[i] > d ? st->last3[i] - d : 0;
}

/* ---------------- binary tree ---------------- */

static inline uint32_t hash_bt(const uint8_t *p, unsigned hlog)
{
    return (rd32(p) * 2654435761u) >> (32 - hlog);
}

/* Insert position p into its tree; if `out` is given, also return matches of
 * strictly increasing length (up to MAX_CAND). `avail` = readable bytes. */
static unsigned bt_search(lzr_cctx *c, size_t p, size_t avail, uint32_t h, mpair *out)
{
    const uint8_t *buf = c->buf, *cur = buf + p;
    const size_t wmask = c->wmask;
    const size_t limit = avail < c->lp.nice ? avail : c->lp.nice;
    uint32_t cm1 = c->head[h];
    c->head[h] = (uint32_t)(p + 1);

    uint32_t *son = c->son;
    uint32_t *ptr0 = son + ((p & wmask) << 1) + 1;  /* right: greater suffixes */
    uint32_t *ptr1 = son + ((p & wmask) << 1);      /* left: smaller suffixes  */
    size_t len0 = 0, len1 = 0, best = LZR_MIN_MATCH - 1;
    unsigned n = 0, depth = c->lp.chain;

    for (;;) {
        if (!cm1 || depth-- == 0) { *ptr0 = *ptr1 = 0; break; }
        size_t cm = cm1 - 1, delta = p - cm;
        if (delta >= c->bt_maxdist) { *ptr0 = *ptr1 = 0; break; }
        uint32_t *pair = son + ((cm & wmask) << 1);
        const uint8_t *pb = buf + cm;
        size_t len = len0 < len1 ? len0 : len1;
#if 1
        /* The next node is one of the two children: start fetching both
         * (tree node and string) while this node is compared. */
        {
            uint32_t a = pair[0], b = pair[1];
            if (a) { __builtin_prefetch(son + (((a - 1) & wmask) << 1)); __builtin_prefetch(buf + a - 1 + len); }
            if (b) { __builtin_prefetch(son + (((b - 1) & wmask) << 1)); __builtin_prefetch(buf + b - 1 + len); }
        }
#endif
        if (pb[len] == cur[len]) {
            len++;
            len += match_len(pb + len, cur + len, limit - len);
            if (len > best) {
                best = len;
                if (out && n < MAX_CAND) { out[n].len = (uint32_t)len; out[n].dist = (uint32_t)delta; n++; }
                if (len == limit) { *ptr1 = pair[0]; *ptr0 = pair[1]; break; }
            }
        }
        if (pb[len] < cur[len]) { *ptr1 = cm1; ptr1 = pair + 1; cm1 = *ptr1; len1 = len; }
        else                    { *ptr0 = cm1; ptr0 = pair;     cm1 = *ptr0; len0 = len; }
    }
    return n;
}

/* ---------------- parallel match gathering ----------------
 *
 * Every position belongs to exactly one tree (the one rooted at its hash
 * bucket), and a tree only links positions of its own bucket. So buckets
 * can be split between threads: each worker walks all positions in order
 * but only handles those whose bucket it owns, touching disjoint head[] and
 * son[] entries. Each tree sees its positions in the same order whatever
 * the thread count, so the output is identical for 1 or N threads.
 */


#ifdef LZR_TUNE
#include <time.h>
double lzr_t_worker[LZR_MAX_THREADS];
static double now_w(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + t.tv_nsec * 1e-9; }
#endif
static void *gather_worker_(void *arg);
static void *gather_worker(void *arg)
{
#ifdef LZR_TUNE
    double t0 = now_w();
    void *r = gather_worker_(arg);
    lzr_t_worker[((gjob *)arg)->part] += now_w() - t0;
    return r;
#else
    return gather_worker_(arg);
#endif
}
static void *gather_worker_(void *arg)
{
    gjob *j = (gjob *)arg;
    lzr_cctx *c = j->c;
    const unsigned hlog = c->lp.hlog;
    mpair tmp[MAX_CAND];
    j->n = 0;
    for (size_t p = j->from; p < j->be; p++) {
        if (p + LZR_MIN_MATCH > j->lim_end) break;
        uint32_t h = hash_bt(c->buf + p, hlog);
        if (h % j->nparts != j->part) continue;
        if (p + 12 <= j->lim_end) __builtin_prefetch(c->head + hash_bt(c->buf + p + 8, hlog));
        unsigned n = bt_search(c, p, j->lim_end - p, h, p >= j->bs ? tmp : NULL);
        if (!n) continue;
        if (j->n + n > j->cap) {
            size_t ncap = j->cap ? j->cap * 2 : 65536;
            while (ncap < j->n + n) ncap *= 2;
            mrec *o = (mrec *)realloc(j->out, ncap * sizeof(mrec));
            if (!o) { j->err = 1; return NULL; }
            j->out = o;
            j->cap = ncap;
        }
        uint32_t maxlen = (uint32_t)(j->be - p), prev = 0;
        for (unsigned k = 0; k < n; k++) {
            uint32_t len = tmp[k].len < maxlen ? tmp[k].len : maxlen;
            /* Defensive re-check: never trust a tree-reported length. */
            len = (uint32_t)match_len(c->buf + p - tmp[k].dist, c->buf + p, len);
            if (len < LZR_MIN_MATCH || len <= prev) continue;   /* clamped duplicates */
            j->out[j->n].pos = (uint32_t)(p - j->bs);
            j->out[j->n].len = len;
            j->out[j->n].dist = tmp[k].dist;
            j->n++;
            prev = len;
        }
    }
    return NULL;
}

/* Fill the match cache for [bs, be). */
/* Bytes readable for tree comparisons: up to LZR_LOOKAHEAD past the block
 * (always present for non-final blocks). Keeping the compare limit at
 * `nice` everywhere except the very end of the input keeps the tree
 * ordering consistent, and output independent of chunking. */
static size_t gather_limit(const lzr_cctx *c, size_t be)
{
    return be + LZR_LOOKAHEAD < c->end ? be + LZR_LOOKAHEAD : c->end;
}

static int gather(lzr_cctx *c, mcache *st, size_t bs, size_t be, size_t lim_end)
{
    unsigned T = c->threads < 1 ? 1 : (c->threads > LZR_MAX_THREADS ? LZR_MAX_THREADS : c->threads);
    gjob *jobs = st->jobs;
    for (unsigned t = 0; t < T; t++) {
        jobs[t].c = c; jobs[t].from = c->ins; jobs[t].bs = bs; jobs[t].be = be;
        jobs[t].lim_end = lim_end; jobs[t].part = t; jobs[t].nparts = T; jobs[t].err = 0;
    }
#if LZR_HAVE_THREADS
    pthread_t tid[LZR_MAX_THREADS];
    unsigned started = 0;
    for (unsigned t = 1; t < T; t++) {
        if (pthread_create(&tid[t], NULL, gather_worker, &jobs[t]) != 0) break;
        started = t;
    }
    gather_worker(&jobs[0]);
    for (unsigned t = 1; t <= started; t++) pthread_join(tid[t], NULL);
    for (unsigned t = started + 1; t < T; t++) gather_worker(&jobs[t]);   /* thread start failed */
#else
    for (unsigned t = 0; t < T; t++) gather_worker(&jobs[t]);
#endif

    /* Merge the per-partition lists (each sorted by position). */
    size_t total = 0;
    for (unsigned t = 0; t < T; t++) { if (jobs[t].err) return LZR_ERR_NOMEM; total += jobs[t].n; }
    if (total > st->pcap) {
        mpair *np2 = (mpair *)realloc(st->pairs, total * sizeof(mpair));
        if (!np2) return LZR_ERR_NOMEM;
        st->pairs = np2;
        st->pcap = total;
    }
    size_t idx[LZR_MAX_THREADS] = {0}, np = 0;
    for (size_t p = 0; p < be - bs; p++) {
        st->mstart[p] = (uint32_t)np;
        for (unsigned t = 0; t < T; t++) {
            gjob *j = &jobs[t];
            while (idx[t] < j->n && j->out[idx[t]].pos == p) {
                st->pairs[np].len = j->out[idx[t]].len;
                st->pairs[np].dist = j->out[idx[t]].dist;
                np++;
                idx[t]++;
            }
        }
    }
    st->mstart[be - bs] = (uint32_t)np;
    c->ins = be;
    return LZR_OK;
}

/* ---------------- prices ---------------- */

static void prices_from_counts(uint32_t *price, const uint32_t *cnt, unsigned n)
{
    uint64_t total = 0;
    for (unsigned s = 0; s < n; s++) total += cnt[s];
    if (!total) { for (unsigned s = 0; s < n; s++) price[s] = 8 * BITCOST; return; }
    double lt = log2((double)total);
    for (unsigned s = 0; s < n; s++) {
        double bits = cnt[s] ? lt - log2((double)cnt[s]) : lt + 2.0;   /* unseen: expensive, finite */
        price[s] = (uint32_t)(bits * BITCOST + 0.5);
    }
}

static void default_prices(prices_t *pr)
{
    for (unsigned s = 0; s < 256; s++) pr->lit[s] = 6 * BITCOST;
    for (unsigned s = 0; s < VAL_CODES; s++) { pr->llc[s] = 3 * BITCOST; pr->mlc[s] = 4 * BITCOST; }
    pr->llc[0] = 2 * BITCOST;
    pr->ofc[0] = 2 * BITCOST; pr->ofc[1] = 3 * BITCOST; pr->ofc[2] = 3 * BITCOST;
    for (unsigned s = LZR_NUM_REPS; s < OF_CODES; s++) pr->ofc[s] = 6 * BITCOST;
}

/* Statistics of the sequences currently in c->seq_* / c->lits. */
static void prices_from_parse(lzr_cctx *c, prices_t *pr, size_t ns, size_t nlit)
{
    uint32_t lit[256] = {0}, llc[VAL_CODES] = {0}, mlc[VAL_CODES] = {0}, ofc[OF_CODES] = {0};
    unsigned nb;
    uint32_t ex;
    for (size_t i = 0; i < nlit; i++) lit[c->lits[i]]++;
    for (size_t i = 0; i < ns; i++) {
        llc[val_code(c->seq_ll[i], &nb, &ex)]++;
        if (c->seq_ml[i] >= LZR_ML_BASE) mlc[val_code(c->seq_ml[i] - LZR_ML_BASE, &nb, &ex)]++;   /* type 4: may be 1 */
        uint32_t off = c->seq_dist[i];
        ofc[off < LZR_NUM_REPS ? off : LZR_NUM_REPS + val_code(off - 3, &nb, &ex)]++;
    }
    prices_from_counts(pr->lit, lit, 256);
    prices_from_counts(pr->llc, llc, VAL_CODES);
    prices_from_counts(pr->mlc, mlc, VAL_CODES);
    prices_from_counts(pr->ofc, ofc, OF_CODES);
}

static inline int32_t ll_price(const prices_t *pr, uint32_t ll)
{
    unsigned nb; uint32_t ex;
    unsigned code = val_code(ll, &nb, &ex);
    return (int32_t)(pr->llc[code] + nb * BITCOST);
}

static inline int32_t of_price(const prices_t *pr, uint32_t off)
{
    if (off < LZR_NUM_REPS) return (int32_t)pr->ofc[off];
    unsigned nb; uint32_t ex;
    unsigned oc = val_code(off - 3, &nb, &ex);
    return (int32_t)(pr->ofc[LZR_NUM_REPS + oc] + nb * BITCOST);
}

static void build_tables(opt_state *st)
{
    for (uint32_t L = LZR_ML_BASE; L < ML_TABLE; L++) {
        unsigned nb; uint32_t ex;
        unsigned mc = val_code(L - LZR_ML_BASE, &nb, &ex);
        st->ml_p[L] = (int32_t)(st->pr.mlc[mc] + nb * BITCOST);
    }
    for (uint32_t L = 0; L < LL_TABLE; L++) st->ll_p[L] = ll_price(&st->pr, L);
}

static inline int32_t llp(const opt_state *st, uint32_t ll)
{
    return ll < LL_TABLE ? st->ll_p[ll] : ll_price(&st->pr, ll);
}

/* Nearest earlier 2- and 3-byte matches for every position of [bs, be). */
static void short_matches(lzr_cctx *c, opt_state *st, size_t bs, size_t be)
{
    const uint8_t *buf = c->buf;
    const size_t maxd = c->wsize;
    for (size_t p = bs; p < be; p++) {
        uint32_t d2 = 0, d3 = 0;
        if (p + 3 <= c->end) {
            unsigned v2 = buf[p] | (unsigned)buf[p + 1] << 8;
            uint32_t h3 = ((v2 | (uint32_t)buf[p + 2] << 16) * 2654435761u) >> (32 - SH3_LOG);
            uint32_t q = st->last2[v2];
            if (q && p - (q - 1) <= maxd) d2 = (uint32_t)(p - (q - 1));
            q = st->last3[h3];
            if (q && p - (q - 1) <= maxd && buf[q - 1] == buf[p] && buf[q] == buf[p + 1] && buf[q + 1] == buf[p + 2])
                d3 = (uint32_t)(p - (q - 1));
            st->last2[v2] = (uint32_t)(p + 1);
            st->last3[h3] = (uint32_t)(p + 1);
        }
        st->sh2[p - bs] = d2;
        st->sh3[p - bs] = d3;
    }
}

/* ---------------- dynamic programming ---------------- */

static inline void emit_seq(lzr_cctx *c, size_t *ns, size_t *nlit, size_t *lit_start, size_t mpos,
                            uint32_t ml, uint32_t off)
{
    size_t ll = mpos - *lit_start;
    memcpy(c->lits + *nlit, c->buf + *lit_start, ll);
    *nlit += ll;
    c->seq_ll[*ns] = (uint32_t)ll;
    c->seq_ml[*ns] = ml;
    c->seq_dist[*ns] = commit_offset(c->reps, off);
    (*ns)++;
    *lit_start = mpos + ml;
}

static size_t parse_pass(lzr_cctx *c, opt_state *st, size_t bs, size_t be, size_t *nlit_out)
{
    const mcache *mc = &st->mc[st->cur];
    const prices_t *pr = &st->pr;
    const uint8_t *buf = c->buf;
    const uint32_t nice = c->lp.nice;
    onode *N = st->nodes;
    size_t ns = 0, nlit = 0, p = bs, lit_start = bs;

    while (p + LZR_MIN_MATCH <= be) {
        size_t seg = be - LZR_MIN_MATCH + 1 - p;
        if (seg > OPT_NUM) seg = OPT_NUM;
        uint32_t ll0 = (uint32_t)(p - lit_start);
        const int adaptive = c->adaptive;
        if (adaptive) {
            shadow_advance(c, st, bs, ns, p);
            tok_prices(st);
            N[0].state = st->pm.state;
        }
        const tok_model *pm = &st->pm;
        N[0].price = adaptive ? 0 : llp(st, ll0);
        N[0].litlen = ll0;
        N[0].ml = 0;
        N[0].ml2 = 0;
        memcpy(N[0].reps, c->reps, sizeof(c->reps));
        size_t last = 0;
        int forced = 0;
        size_t f_i = 0;
        uint32_t f_len = 0, f_off = 0;

        for (size_t i = 0; i < seg; i++) {
            /* Every path passes through i: end the segment here, so the
             * next one is priced from a fresher model (block type 4). */
            if (adaptive && i >= SEG_MIN && last == i) break;
            onode *n = &N[i];
            /* Repeat offsets on the cheapest path to i are derived here,
             * once, rather than on every tentative update. */
            if (i) {
                const onode *f = &N[n->from];
                memcpy(n->reps, f->reps, sizeof(n->reps));
                if (n->ml) commit_offset(n->reps, n->off);
                if (n->ml2) n->state = tok_after_rep(tok_after_lit(!n->ml ? f->state : n->off < LZR_NUM_REPS ? tok_after_rep(f->state) : tok_after_match(f->state)));
                else n->state = !n->ml ? tok_after_lit(f->state)
                           : n->ml == 1 ? tok_after_short(f->state)
                           : n->off < LZR_NUM_REPS ? tok_after_rep(f->state) : tok_after_match(f->state);
            }
            const size_t cur = p + i;
            const uint8_t *cp = buf + cur;

            /* literal transition (charges the change in literal-run cost) */
            if (last < i + 1) { N[i + 1].price = PRICE_INF; last = i + 1; }
            const unsigned ps = (unsigned)(cur - bs) & (TOK_POS - 1), s = n->state;
            int32_t lp;
            if (adaptive) {
                unsigned ctx = (cur ? cp[-1] : 0) >> 5;
                lp = n->price + arc_bit_price(st, pm->is_match[s][ps], 0) +
                     (s >= 7 ? matched_price(st, pm->lit[ctx], *cp, cp[-(ptrdiff_t)n->reps[0]]) : st->litp[ctx][*cp]);
            } else {
                lp = n->price + (int32_t)pr->lit[*cp] + llp(st, n->litlen + 1) - llp(st, n->litlen);
            }
            if (lp < N[i + 1].price) {
                onode *t = &N[i + 1];
                t->price = lp; t->from = (uint32_t)i; t->ml = 0; t->off = 0; t->litlen = n->litlen + 1; t->ml2 = 0;
            }
            /* literal + rep0, for when the literal is not the cheapest way
             * to i + 1 */
            if (adaptive && cur + 3 <= be && n->reps[0] < cur + 1 &&
                cp[1] == cp[1 - (ptrdiff_t)n->reps[0]] && cp[2] == cp[2 - (ptrdiff_t)n->reps[0]]) {
                uint32_t d = n->reps[0], lim = (uint32_t)(be - cur - 1), len2;
                if (lim > TOK_MAXLEN) lim = TOK_MAXLEN;
                len2 = (uint32_t)match_len(cp + 1 - d, cp + 1, lim);
                unsigned s2 = tok_after_lit(s), ps2 = (ps + 1) & (TOK_POS - 1);
                int32_t price = lp + arc_bit_price(st, pm->is_match[s2][ps2], 1) + arc_bit_price(st, pm->is_rep[s2], 1) +
                                arc_bit_price(st, pm->is_rep0[s2], 0) + arc_bit_price(st, pm->is_rep0_long[s2][ps2], 1) +
                                st->lenp[1][ps2][len2];
                size_t t = i + 1 + len2;
                while (last < t) N[++last].price = PRICE_INF;
                if (price < N[t].price) {
                    onode *tn = &N[t];
                    tn->price = price; tn->from = (uint32_t)i; tn->ml = 0; tn->off = 0; tn->litlen = 0; tn->ml2 = len2;
                }
            }
            /* short rep: the byte at rep0, as a 1-byte match */
            if (adaptive && n->reps[0] <= cur && cp[-(ptrdiff_t)n->reps[0]] == cp[0]) {
                int32_t sp = n->price + arc_bit_price(st, pm->is_match[s][ps], 1) + arc_bit_price(st, pm->is_rep[s], 1) +
                             arc_bit_price(st, pm->is_rep0[s], 0) + arc_bit_price(st, pm->is_rep0_long[s][ps], 0);
                if (sp < N[i + 1].price) {
                    onode *t = &N[i + 1];
                    t->price = sp; t->from = (uint32_t)i; t->ml = 1; t->off = 0; t->litlen = 0; t->ml2 = 0;
                }
            }

            /* candidate matches: repeat offsets of this path, then the tree */
            uint32_t cl[LZR_NUM_REPS + 2 + MAX_CAND], co[LZR_NUM_REPS + 2 + MAX_CAND], cmin[LZR_NUM_REPS + 2 + MAX_CAND];
            unsigned nc = 0;
            const uint32_t maxlen = (uint32_t)(be - cur);
            uint32_t best_len = 0, best_off = 0;
            for (uint32_t r = 0; r < LZR_NUM_REPS; r++) {
                uint32_t d = n->reps[r];
                if (d > cur || d > c->wsize || cp[-(ptrdiff_t)d] != cp[0] || cp[1 - (ptrdiff_t)d] != cp[1]) continue;
                uint32_t len = (uint32_t)match_len(cp - d, cp, maxlen);
                if (len < LZR_REP_MIN) continue;
                cl[nc] = len; co[nc] = r; cmin[nc] = LZR_REP_MIN; nc++;
                if (len > best_len) { best_len = len; best_off = r; }
            }
            uint32_t prevlen = LZR_MIN_MATCH - 1;
            if (adaptive) {        /* 2- and 3-byte matches (tree finds >= 4) */
                uint32_t d2 = st->sh2[cur - bs], d3 = st->sh3[cur - bs];
                for (unsigned t = 0; t < 2; t++) {
                    uint32_t d = t ? d3 : d2;
                    if (!d || (t && d == d2) || cp[-(ptrdiff_t)d] != cp[0] || cp[1 - (ptrdiff_t)d] != cp[1]) continue;
                    uint32_t off = d + 2;
                    for (uint32_t r = 0; r < LZR_NUM_REPS; r++) if (n->reps[r] == d) { off = r; break; }
                    if (off < LZR_NUM_REPS) continue;   /* covered by the repeat probes */
                    /* short lengths only at short distances: a far 2-byte
                     * match rarely pays and evicts a repeat offset */
                    uint32_t len = (uint32_t)match_len(cp - d, cp, maxlen), lmin = d <= SH2_MAXD ? 2 : d <= SH3_MAXD ? 3 : 4;
                    if (len < lmin) continue;
                    cl[nc] = len; co[nc] = off; cmin[nc] = lmin; nc++;
                    if (len > best_len) { best_len = len; best_off = off; }
                }
            }
            for (uint32_t k = mc->mstart[cur - bs]; k < mc->mstart[cur - bs + 1]; k++) {
                uint32_t len = mc->pairs[k].len, d = mc->pairs[k].dist, off = d + 2;
                for (uint32_t r = 0; r < LZR_NUM_REPS; r++) if (n->reps[r] == d) { off = r; break; }
                cl[nc] = len; co[nc] = off; cmin[nc] = prevlen + 1; nc++;
                prevlen = len;
                if (len > best_len) { best_len = len; best_off = off; }
            }
            if (best_len >= nice) {       /* long match: take it, end segment */
                forced = 1; f_i = i; f_len = best_len; f_off = best_off;
                break;
            }
            int32_t base, rep_b[LZR_NUM_REPS], new_b = 0;
            if (adaptive) {
                base = n->price + arc_bit_price(st, pm->is_match[s][ps], 1);
                int32_t r1 = base + arc_bit_price(st, pm->is_rep[s], 1);
                new_b = base + arc_bit_price(st, pm->is_rep[s], 0);
                rep_b[0] = r1 + arc_bit_price(st, pm->is_rep0[s], 0) + arc_bit_price(st, pm->is_rep0_long[s][ps], 1);
                r1 += arc_bit_price(st, pm->is_rep0[s], 1);
                rep_b[1] = r1 + arc_bit_price(st, pm->is_rep1[s], 0);
                rep_b[2] = r1 + arc_bit_price(st, pm->is_rep1[s], 1);
            } else {
                base = n->price + st->ll_p[0];
            }
            for (unsigned k = 0; k < nc; k++) {
                int32_t bo;
                const int32_t *lenp = NULL;
                int32_t dpr[4] = {0, 0, 0, 0};
                if (!adaptive) {
                    bo = base + of_price(pr, co[k]);
                } else if (co[k] < LZR_NUM_REPS) {
                    bo = rep_b[co[k]];
                    lenp = st->lenp[1][ps];
                } else {
                    bo = new_b;
                    lenp = st->lenp[0][ps];
                    for (unsigned lc = 0; lc < 4; lc++) dpr[lc] = dist_price(st, co[k] - 3, lc);
                }
                /* Lengths past LZR_OPT_FULL are only tried at full length:
                 * cutting a long match short almost never pays, and this
                 * bounds the work per position. */
                uint32_t Lmax = cl[k], Lcut = cl[k] > LZR_OPT_FULL ? LZR_OPT_FULL : cl[k];
                int32_t full = 0;
                for (uint32_t L = cmin[k]; L <= Lmax; L = (L < Lcut || L == Lmax) ? L + 1 : Lmax) {
                    int32_t price;
                    if (!adaptive) price = bo + st->ml_p[L];
                    else if (L <= TOK_MAXLEN) price = bo + lenp[L] + dpr[L - 2 < 3 ? L - 2 : 3];
                    else price = bo + lenp[TOK_MAXLEN] + dpr[3] + st->lenp[1][0][L - TOK_MAXLEN < 2 ? 2 : (L - TOK_MAXLEN > TOK_MAXLEN ? TOK_MAXLEN : L - TOK_MAXLEN)];
                    size_t t = i + L;
                    while (last < t) N[++last].price = PRICE_INF;
                    if (price < N[t].price) {
                        onode *tn = &N[t];
                        tn->price = price; tn->from = (uint32_t)i; tn->ml = L; tn->off = co[k]; tn->litlen = 0; tn->ml2 = 0;
                    }
                    if (L == Lmax) full = price;
                }
                /* match + literal + rep0 (the literal breaks the match, the
                 * same distance continues after it) */
                if (adaptive && Lmax <= TOK_MAXLEN && Lmax + 3 <= maxlen) {
                    uint32_t d = co[k] < LZR_NUM_REPS ? n->reps[co[k]] : co[k] - 2;
                    const uint8_t *q = cp + Lmax + 1;
                    if (q[0] == q[-(ptrdiff_t)d] && q[1] == q[1 - (ptrdiff_t)d]) {
                        uint32_t lim = maxlen - Lmax - 1, len2;
                        if (lim > TOK_MAXLEN) lim = TOK_MAXLEN;
                        len2 = (uint32_t)match_len(q - d, q, lim);
                        unsigned s1 = co[k] < LZR_NUM_REPS ? tok_after_rep(s) : tok_after_match(s);
                        unsigned ps1 = (ps + Lmax) & (TOK_POS - 1), ps2 = (ps1 + 1) & (TOK_POS - 1);
                        unsigned s2 = tok_after_lit(s1);
                        int32_t price = full + arc_bit_price(st, pm->is_match[s1][ps1], 0) +
                                        matched_price(st, pm->lit[q[-2] >> 5], q[-1], q[-1 - (ptrdiff_t)d]) +
                                        arc_bit_price(st, pm->is_match[s2][ps2], 1) + arc_bit_price(st, pm->is_rep[s2], 1) +
                                        arc_bit_price(st, pm->is_rep0[s2], 0) + arc_bit_price(st, pm->is_rep0_long[s2][ps2], 1) +
                                        st->lenp[1][ps2][len2];
                        size_t t = i + Lmax + 1 + len2;
                        while (last < t) N[++last].price = PRICE_INF;
                        if (price < N[t].price) {
                            onode *tn = &N[t];
                            tn->price = price; tn->from = (uint32_t)i; tn->ml = Lmax; tn->off = co[k]; tn->litlen = 0;
                            tn->ml2 = len2;
                        }
                    }
                }
            }
        }

        size_t endi = forced ? f_i : last;
        while (!forced && endi > 0 && N[endi].price >= PRICE_INF) endi--;

        /* Backtrack the cheapest path, then emit it forwards. */
        size_t cnt = 0, j = endi;
        while (j > 0) {
            if (N[j].ml2) {
                st->path[cnt][0] = N[j].from + N[j].ml + 1; st->path[cnt][1] = N[j].ml2; st->path[cnt][2] = 0;
                cnt++;
            }
            if (N[j].ml || N[j].ml2) {
                if (N[j].ml) {
                    st->path[cnt][0] = N[j].from; st->path[cnt][1] = N[j].ml; st->path[cnt][2] = N[j].off;
                    cnt++;
                }
                j = N[j].from;
            } else {
                j--;
            }
        }
        while (cnt--) emit_seq(c, &ns, &nlit, &lit_start, p + st->path[cnt][0], st->path[cnt][1], st->path[cnt][2]);
        if (forced) {
            emit_seq(c, &ns, &nlit, &lit_start, p + f_i, f_len, f_off);
            p += f_i + f_len;
        } else {
            p += endi ? endi : 1;
        }
    }
    memcpy(c->lits + nlit, buf + lit_start, be - lit_start);
    nlit += be - lit_start;
    *nlit_out = nlit;
    return ns;
}

#ifdef LZR_TUNE
#include <time.h>
double lzr_t_gather, lzr_t_parse;
static double now_s(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + t.tv_nsec * 1e-9; }
#define TSTART double t0_ = now_s()
#define TACC(v) do { double t1_ = now_s(); v += t1_ - t0_; t0_ = t1_; } while (0)
#else
#define TSTART
#define TACC(v)
#endif

#if LZR_HAVE_THREADS
static void *gather_async(void *arg)
{
    opt_state *st = (opt_state *)arg;
    st->pend_err = gather(st->pend_c, &st->mc[st->cur ^ 1], st->pend_bs, st->pend_be, st->pend_lim);
    return NULL;
}
#endif

/* Wait for a background gather (keeps its result for the next block). */
void opt_sync(lzr_cctx *c)
{
#if LZR_HAVE_THREADS
    opt_state *st = (opt_state *)c->opt;
    if (st && st->pending == 1) {
        pthread_join(st->pend_tid, NULL);
        st->pending = 2;             /* finished, result not consumed yet */
    }
#else
    (void)c;
#endif
}

size_t opt_parse_block(lzr_cctx *c, size_t bs, size_t be, size_t *nlit_out)
{
    opt_state *st = (opt_state *)c->opt;
    TSTART;
    /* 1. matches for this block: from the background gather, or now */
    opt_sync(c);
    if (st->pending && st->pend_bs == bs && st->pend_be == be) {
        st->cur ^= 1;
        st->pending = 0;
        if (st->pend_err) return (size_t)-1;
    } else {
        if (st->pending) return (size_t)-1;      /* cannot happen: blocks are sequential */
        if (gather(c, &st->mc[st->cur], bs, be, gather_limit(c, be))) return (size_t)-1;
    }
    /* 2. with threads, gather the next block while this one is parsed.
     * Only when that block and its look-ahead are already buffered, so the
     * result is the same as gathering it later. */
#if LZR_HAVE_THREADS
    if (c->threads > 1 && c->end >= be + LZR_BLOCK_SIZE + LZR_LOOKAHEAD) {
        st->pend_c = c;
        st->pend_bs = be;
        st->pend_be = be + LZR_BLOCK_SIZE;
        st->pend_lim = gather_limit(c, st->pend_be);   /* read here: c->end changes */
        st->pend_err = 0;
        if (pthread_create(&st->pend_tid, NULL, gather_async, st) == 0) st->pending = 1;
    }
#endif
    if (c->adaptive) short_matches(c, st, bs, be);
    TACC(lzr_t_gather);
    if (!c->have_stats) default_prices(&st->pr);

    uint32_t saved[LZR_NUM_REPS];
    memcpy(saved, c->reps, sizeof(saved));
    size_t ns = 0;
    for (unsigned pass = 0; pass < c->lp.opt; pass++) {
        if (pass) {                       /* re-price with this block's statistics */
            prices_from_parse(c, &st->pr, ns, *nlit_out);
            memcpy(c->reps, saved, sizeof(saved));
        }
        build_tables(st);
        if (c->adaptive) {
            st->pm = *c->tm;
            memcpy(st->dreps, c->reps, sizeof(st->dreps));
            st->dpos = st->dlit = bs;
            st->dseq = 0;
        }
        ns = parse_pass(c, st, bs, be, nlit_out);
    }
    TACC(lzr_t_parse);
    prices_from_parse(c, &st->pr, ns, *nlit_out);   /* seed the next block */
    c->have_stats = 1;
    return ns;
}
