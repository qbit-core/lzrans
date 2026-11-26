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
 *                      match length + offset) and the repeat offsets on
 *                      that path. Prices come from the symbol statistics
 *                      of the previous pass/block (log2 of frequencies).
 *   3. With 2 passes, the first pass only gathers statistics for the
 *      second, which produces the real sequences.
 */
#include "enc_internal.h"
#include <math.h>
#if LZR_HAVE_THREADS
#include <pthread.h>
#endif

#define OPT_NUM     4096
#define MAX_CAND    48
#define NODE_EXTRA  (300 + 16)
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
    uint32_t reps[LZR_NUM_REPS];
} onode;

typedef struct {
    uint32_t lit[256], llc[VAL_CODES], mlc[VAL_CODES], ofc[OF_CODES];
} prices_t;

#ifndef LZR_OPT_FULL
#define LZR_OPT_FULL 32
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
} opt_state;

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
    return LZR_OK;
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
        free(st->nodes); free(st->path); free(st);
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
        mlc[val_code(c->seq_ml[i] - LZR_MIN_MATCH, &nb, &ex)]++;
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
    for (uint32_t L = LZR_MIN_MATCH; L < ML_TABLE; L++) {
        unsigned nb; uint32_t ex;
        unsigned mc = val_code(L - LZR_MIN_MATCH, &nb, &ex);
        st->ml_p[L] = (int32_t)(st->pr.mlc[mc] + nb * BITCOST);
    }
    for (uint32_t L = 0; L < LL_TABLE; L++) st->ll_p[L] = ll_price(&st->pr, L);
}

static inline int32_t llp(const opt_state *st, uint32_t ll)
{
    return ll < LL_TABLE ? st->ll_p[ll] : ll_price(&st->pr, ll);
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
        N[0].price = llp(st, ll0);
        N[0].litlen = ll0;
        N[0].ml = 0;
        memcpy(N[0].reps, c->reps, sizeof(c->reps));
        size_t last = 0;
        int forced = 0;
        size_t f_i = 0;
        uint32_t f_len = 0, f_off = 0;

        for (size_t i = 0; i < seg; i++) {
            onode *n = &N[i];
            /* Repeat offsets on the cheapest path to i are derived here,
             * once, rather than on every tentative update. */
            if (i) {
                const onode *f = &N[n->from];
                memcpy(n->reps, f->reps, sizeof(n->reps));
                if (n->ml) commit_offset(n->reps, n->off);
            }
            const size_t cur = p + i;
            const uint8_t *cp = buf + cur;

            /* literal transition (charges the change in literal-run cost) */
            if (last < i + 1) { N[i + 1].price = PRICE_INF; last = i + 1; }
            int32_t lp = n->price + (int32_t)pr->lit[*cp] + llp(st, n->litlen + 1) - llp(st, n->litlen);
            if (lp < N[i + 1].price) {
                onode *t = &N[i + 1];
                t->price = lp; t->from = (uint32_t)i; t->ml = 0; t->off = 0; t->litlen = n->litlen + 1;
            }

            /* candidate matches: repeat offsets of this path, then the tree */
            uint32_t cl[LZR_NUM_REPS + MAX_CAND], co[LZR_NUM_REPS + MAX_CAND], cmin[LZR_NUM_REPS + MAX_CAND];
            unsigned nc = 0;
            const uint32_t maxlen = (uint32_t)(be - cur);
            uint32_t best_len = 0, best_off = 0;
            for (uint32_t r = 0; r < LZR_NUM_REPS; r++) {
                uint32_t d = n->reps[r];
                if (d > cur || d > c->wsize || rd32(cp - d) != rd32(cp)) continue;
                uint32_t len = (uint32_t)match_len(cp - d, cp, maxlen);
                if (len < LZR_MIN_MATCH) continue;
                cl[nc] = len; co[nc] = r; cmin[nc] = LZR_MIN_MATCH; nc++;
                if (len > best_len) { best_len = len; best_off = r; }
            }
            uint32_t prevlen = LZR_MIN_MATCH - 1;
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
            int32_t base = n->price + st->ll_p[0];
            for (unsigned k = 0; k < nc; k++) {
                int32_t bo = base + of_price(pr, co[k]);
                /* Lengths past LZR_OPT_FULL are only tried at full length:
                 * cutting a long match short almost never pays, and this
                 * bounds the work per position. */
                uint32_t Lmax = cl[k], Lcut = cl[k] > LZR_OPT_FULL ? LZR_OPT_FULL : cl[k];
                for (uint32_t L = cmin[k]; L <= Lmax; L = (L < Lcut || L == Lmax) ? L + 1 : Lmax) {
                    int32_t price = bo + st->ml_p[L];
                    size_t t = i + L;
                    while (last < t) N[++last].price = PRICE_INF;
                    if (price < N[t].price) {
                        onode *tn = &N[t];
                        tn->price = price; tn->from = (uint32_t)i; tn->ml = L; tn->off = co[k]; tn->litlen = 0;
                    }
                }
            }
        }

        size_t endi = forced ? f_i : last;
        while (!forced && endi > 0 && N[endi].price >= PRICE_INF) endi--;

        /* Backtrack the cheapest path, then emit it forwards. */
        size_t cnt = 0, j = endi;
        while (j > 0) {
            if (N[j].ml) {
                st->path[cnt][0] = N[j].from; st->path[cnt][1] = N[j].ml; st->path[cnt][2] = N[j].off;
                cnt++;
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
        ns = parse_pass(c, st, bs, be, nlit_out);
    }
    TACC(lzr_t_parse);
    prices_from_parse(c, &st->pr, ns, *nlit_out);   /* seed the next block */
    c->have_stats = 1;
    return ns;
}
