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

#define OPT_NUM     4096
#define MAX_CAND    48
#define NODE_EXTRA  (300 + 16)
#define BITCOST     256                 /* price units per bit */
#define PRICE_INF   0x3FFFFFFF

typedef struct { uint32_t len, dist; } mpair;

typedef struct {
    int32_t price;
    uint32_t from, ml, off, litlen;
    uint32_t reps[LZR_NUM_REPS];
} onode;

typedef struct {
    uint32_t lit[256], llc[VAL_CODES], mlc[VAL_CODES], ofc[OF_CODES];
} prices_t;

#define ML_TABLE 512                     /* > max nice length */
#define LL_TABLE 4096

typedef struct {
    uint32_t *mstart;      /* [B+1] index of the first cached pair per position */
    mpair *pairs;
    size_t pcap;
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
    c->son = (uint32_t *)calloc(2 * c->wsize, sizeof(uint32_t));
    if (!st || !c->son) return LZR_ERR_NOMEM;
    st->pcap = (size_t)LZR_BLOCK_SIZE * 4;
    st->mstart = (uint32_t *)malloc((LZR_BLOCK_SIZE + 1) * sizeof(uint32_t));
    st->pairs = (mpair *)malloc(st->pcap * sizeof(mpair));
    st->nodes = (onode *)malloc((OPT_NUM + NODE_EXTRA) * sizeof(onode));
    st->path = malloc((OPT_NUM + NODE_EXTRA) * sizeof(*st->path));
    if (!st->mstart || !st->pairs || !st->nodes || !st->path) return LZR_ERR_NOMEM;
    return LZR_OK;
}

void opt_free(lzr_cctx *c)
{
    opt_state *st = (opt_state *)c->opt;
    if (st) { free(st->mstart); free(st->pairs); free(st->nodes); free(st->path); free(st); }
    free(c->son);
    c->opt = NULL;
    c->son = NULL;
}

void opt_slide(lzr_cctx *c, uint32_t d)
{
    for (size_t i = 0; i < 2 * c->wsize; i++) c->son[i] = c->son[i] > d ? c->son[i] - d : 0;
}

/* ---------------- binary tree ---------------- */

static inline uint32_t hash_bt(const uint8_t *p, unsigned hlog)
{
    return (rd32(p) * 2654435761u) >> (32 - hlog);
}

/* Insert position p into its tree; if `out` is given, also return matches of
 * strictly increasing length (up to MAX_CAND). `avail` = readable bytes. */
static unsigned bt_search(lzr_cctx *c, size_t p, size_t avail, mpair *out)
{
    const uint8_t *buf = c->buf, *cur = buf + p;
    const size_t wmask = c->wmask;
    const size_t limit = avail < c->lp.nice ? avail : c->lp.nice;
    uint32_t h = hash_bt(cur, c->lp.hlog);
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
        if (delta >= c->wsize) { *ptr0 = *ptr1 = 0; break; }
        uint32_t *pair = son + ((cm & wmask) << 1);
        const uint8_t *pb = buf + cm;
        size_t len = len0 < len1 ? len0 : len1;
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

/* Fill the match cache for [bs, be). */
static int gather(lzr_cctx *c, opt_state *st, size_t bs, size_t be)
{
    /* Bytes readable for tree comparisons: up to LZR_LOOKAHEAD past the
     * block (always present for non-final blocks). Keeping the compare
     * limit at `nice` everywhere except the very end of the input keeps
     * the tree ordering consistent, and output independent of chunking. */
    const size_t lim_end = be + LZR_LOOKAHEAD < c->end ? be + LZR_LOOKAHEAD : c->end;
    mpair tmp[MAX_CAND];

    for (size_t q = c->ins; q < bs; q++)                    /* catch up */
        if (q + LZR_MIN_MATCH <= lim_end) bt_search(c, q, lim_end - q, NULL);

    size_t np = 0, p = bs;
    while (p < be) {
        st->mstart[p - bs] = (uint32_t)np;
        if (p + LZR_MIN_MATCH > lim_end) { p++; continue; }
        unsigned n = bt_search(c, p, lim_end - p, tmp);
        if (np + n > st->pcap) {
            size_t ncap = st->pcap * 2;
            mpair *np2 = (mpair *)realloc(st->pairs, ncap * sizeof(mpair));
            if (!np2) return LZR_ERR_NOMEM;
            st->pairs = np2;
            st->pcap = ncap;
        }
        uint32_t maxlen = (uint32_t)(be - p), prev = 0, longest = 0;
        for (unsigned k = 0; k < n; k++) {
            uint32_t len = tmp[k].len < maxlen ? tmp[k].len : maxlen;
            /* Defensive re-check: never trust a tree-reported length. */
            len = (uint32_t)match_len(c->buf + p - tmp[k].dist, c->buf + p, len);
            if (len < LZR_MIN_MATCH || len <= prev) continue;   /* clamped duplicates */
            st->pairs[np].len = len;
            st->pairs[np].dist = tmp[k].dist;
            np++;
            prev = len;
        }
        if (n) longest = tmp[n - 1].len;
        if (longest >= c->lp.nice) {
            /* Long match: the parser will take it, so only insert. */
            size_t stop = p + (prev ? prev : 1);
            for (size_t q = p + 1; q < stop; q++) {
                st->mstart[q - bs] = (uint32_t)np;
                if (q + LZR_MIN_MATCH <= lim_end) bt_search(c, q, lim_end - q, NULL);
            }
            p = stop;
        } else {
            p++;
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
            const onode *n = &N[i];
            const size_t cur = p + i;
            const uint8_t *cp = buf + cur;

            /* literal transition (charges the change in literal-run cost) */
            if (last < i + 1) { N[i + 1].price = PRICE_INF; last = i + 1; }
            int32_t lp = n->price + (int32_t)pr->lit[*cp] + llp(st, n->litlen + 1) - llp(st, n->litlen);
            if (lp < N[i + 1].price) {
                onode *t = &N[i + 1];
                t->price = lp; t->from = (uint32_t)i; t->ml = 0; t->off = 0; t->litlen = n->litlen + 1;
                memcpy(t->reps, n->reps, sizeof(t->reps));
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
            for (uint32_t k = st->mstart[cur - bs]; k < st->mstart[cur - bs + 1]; k++) {
                uint32_t len = st->pairs[k].len, d = st->pairs[k].dist, off = d + 2;
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
                for (uint32_t L = cmin[k]; L <= cl[k]; L++) {
                    int32_t price = bo + st->ml_p[L];
                    size_t t = i + L;
                    while (last < t) N[++last].price = PRICE_INF;
                    if (price < N[t].price) {
                        onode *tn = &N[t];
                        tn->price = price; tn->from = (uint32_t)i; tn->ml = L; tn->off = co[k]; tn->litlen = 0;
                        memcpy(tn->reps, n->reps, sizeof(tn->reps));
                        commit_offset(tn->reps, co[k]);
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

size_t opt_parse_block(lzr_cctx *c, size_t bs, size_t be, size_t *nlit_out)
{
    opt_state *st = (opt_state *)c->opt;
    if (gather(c, st, bs, be)) return (size_t)-1;
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
    prices_from_parse(c, &st->pr, ns, *nlit_out);   /* seed the next block */
    c->have_stats = 1;
    return ns;
}
