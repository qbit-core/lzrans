/*
 * huf.c - length-limited canonical Huffman coding of literals, 4 streams
 * (format in huf.h).
 */
#include "huf.h"
#include <stdlib.h>
#include <string.h>

/* ---- code construction ------------------------------------------- */

/* Code lengths (<= HUF_MAXLEN) for the counts; absent symbols get 0. At
 * least one symbol must be present. */
static int cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}

void huf_lengths(const uint32_t *cnt, uint8_t *len)
{
    uint16_t sym[256];
    unsigned n = 0;
    memset(len, 0, 256);
    for (unsigned s = 0; s < 256; s++) if (cnt[s]) sym[n++] = (uint16_t)s;
    if (n == 1) { len[sym[0]] = 1; return; }
    /* sort by count ascending, ties by symbol: keys count << 8 | symbol
     * are distinct, so any sort gives this order */
    {
        uint64_t key[256];
        for (unsigned i = 0; i < n; i++) key[i] = (uint64_t)cnt[sym[i]] << 8 | sym[i];
        qsort(key, n, sizeof(key[0]), cmp_u64);
        for (unsigned i = 0; i < n; i++) sym[i] = (uint16_t)(key[i] & 255);
    }
    /* two-queue Huffman: leaves 0..n-1 (sorted), internal nodes n.. */
    uint64_t w[512];
    uint16_t parent[512];
    for (unsigned i = 0; i < n; i++) w[i] = cnt[sym[i]];
    unsigned leaf = 0, node = n, next = n;
    for (unsigned k = 0; k < n - 1; k++) {
        unsigned pick[2];
        for (int t = 0; t < 2; t++) {
            if (leaf < n && (node >= next || w[leaf] <= w[node])) pick[t] = leaf++;
            else pick[t] = node++;
        }
        w[next] = w[pick[0]] + w[pick[1]];
        parent[pick[0]] = parent[pick[1]] = (uint16_t)next;
        next++;
    }
    /* depths: root is next - 1 */
    uint8_t depth[512];
    depth[next - 1] = 0;
    for (unsigned i = next - 1; i-- > 0;) depth[i] = (uint8_t)(depth[parent[i]] + 1 > 255 ? 255 : depth[parent[i]] + 1);
    /* limit: clamp, then lengthen the rarest codes until Kraft holds */
    const uint32_t full = 1u << HUF_MAXLEN;
    uint32_t kraft = 0;
    for (unsigned i = 0; i < n; i++) {
        unsigned d = depth[i] > HUF_MAXLEN ? HUF_MAXLEN : depth[i];
        len[sym[i]] = (uint8_t)d;
        kraft += full >> d;
    }
    while (kraft > full) {
        for (unsigned i = 0; i < n && kraft > full; i++) {
            unsigned s = sym[i];
            if (len[s] < HUF_MAXLEN) { kraft -= full >> (len[s] + 1); len[s]++; }
        }
    }
    /* use leftover code space: shorten the most frequent codes */
    for (unsigned i = n; i-- > 0;) {
        unsigned s = sym[i];
        while (len[s] > 1 && kraft + (full >> len[s]) <= full) { kraft += full >> len[s]; len[s]--; }
    }
    /* format v7 requires a complete code (Kraft sum exactly 2^11). The
     * shortfall is a multiple of the unit of the longest code, so
     * shortening a longest code always fits (rarely needed). */
    while (kraft < full) {
        unsigned best = sym[0];
        for (unsigned i = 1; i < n; i++) if (len[sym[i]] > len[best]) best = sym[i];
        kraft += full >> len[best];
        len[best]--;
    }
}

/* Canonical codes, bit-reversed for LSB-first output. */
static void huf_codes(huf_table *t)
{
    unsigned cnt[HUF_MAXLEN + 2] = {0};
    for (unsigned s = 0; s < 256; s++) cnt[t->len[s]]++;
    cnt[0] = 0;
    unsigned start[HUF_MAXLEN + 2], code = 0;
    for (unsigned l = 1; l <= HUF_MAXLEN; l++) { code = (code + cnt[l - 1]) << 1; start[l] = code; }
    for (unsigned s = 0; s < 256; s++) {
        unsigned l = t->len[s];
        if (!l) continue;
        unsigned c = start[l]++, r = 0;
        for (unsigned b = 0; b < l; b++) r |= ((c >> b) & 1) << (l - 1 - b);
        t->code[s] = (uint16_t)r;
    }
}

/* Decoder table from lengths. Returns 0, or -1 if the lengths overflow
 * the code space (an incomplete code is accepted: unused slots decode as
 * a full-length symbol 0 and the stream-end check catches them). */
static int huf_build_dt(huf_table *t, int strict)
{
    uint32_t kraft = 0;
    for (unsigned s = 0; s < 256; s++) {
        if (t->len[s] > HUF_MAXLEN) return -1;
        if (t->len[s]) kraft += HUF_TSIZE >> t->len[s];
    }
    if (kraft > HUF_TSIZE || !kraft || (strict && kraft != HUF_TSIZE)) return -1;
    huf_codes(t);
    for (unsigned i = 0; i < HUF_TSIZE; i++) t->dt[i] = (uint16_t)(HUF_MAXLEN << 8);
    for (unsigned s = 0; s < 256; s++) {
        unsigned l = t->len[s];
        if (!l) continue;
        for (unsigned j = t->code[s]; j < HUF_TSIZE; j += 1u << l) t->dt[j] = (uint16_t)(s | l << 8);
    }
    t->valid = 1;
    return 0;
}

/* Mode 2 header of a literal section at in (mode byte first): build its
 * table into t. Returns the header bytes (mode, last, nibbles), 0 if
 * malformed. */
size_t huf_read_table(const uint8_t *in, size_t avail, huf_table *t, int v7)
{
    if (avail < 2 || in[0] != 2) return 0;
    unsigned last = in[1];
    size_t hb = (last + 2) / 2;
    if (avail - 2 < hb) return 0;
    memset(t->len, 0, sizeof(t->len));
    for (unsigned s = 0; s <= last; s++) t->len[s] = (uint8_t)((in[2 + s / 2] >> (4 * (s & 1))) & 15);
    t->valid = 0;
    if (huf_build_dt(t, v7)) return 0;
    return 2 + hb;
}

/* ---- encoder ------------------------------------------------------ */

typedef struct { uint64_t acc; unsigned n; uint8_t *p; } hbw;

static inline void hbw_put(hbw *w, uint32_t v, unsigned len)
{
    w->acc |= (uint64_t)v << w->n;
    w->n += len;
    if (w->n >= 32) { wr32le(w->p, (uint32_t)w->acc); w->p += 4; w->acc >>= 32; w->n -= 32; }
}

static inline void hbw_flush(hbw *w)
{
    while (w->n > 0) { *w->p++ = (uint8_t)w->acc; w->acc >>= 8; w->n = w->n > 8 ? w->n - 8 : 0; }
}

size_t huf_encode(uint8_t *out, const uint8_t *lits, size_t n, huf_table *prev)
{
    if (!n) return 0;
    /* 4 histograms: consecutive equal bytes do not wait on each other */
    uint32_t cnt[256], c4[3][256];
    memset(cnt, 0, sizeof(cnt));
    memset(c4, 0, sizeof(c4));
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        cnt[lits[i]]++; c4[0][lits[i + 1]]++; c4[1][lits[i + 2]]++; c4[2][lits[i + 3]]++;
    }
    for (; i < n; i++) cnt[lits[i]]++;
    for (unsigned s = 0; s < 256; s++) cnt[s] += c4[0][s] + c4[1][s] + c4[2][s];
    unsigned distinct = 0, last = 0;
    for (unsigned s = 0; s < 256; s++) if (cnt[s]) { distinct++; last = s; }
    if (distinct == 1) { out[0] = 1; out[1] = lits[0]; return 2; }

    huf_table nt;
    huf_lengths(cnt, nt.len);
    uint64_t bits_new = 0, bits_prev = 0;
    int prev_ok = prev->valid;
    for (unsigned s = 0; s < 256; s++) {
        if (!cnt[s]) continue;
        bits_new += (uint64_t)cnt[s] * nt.len[s];
        if (!prev->len[s]) prev_ok = 0;
        bits_prev += (uint64_t)cnt[s] * prev->len[s];
    }
    size_t hdr_new = 1 + (last + 2) / 2;
    uint64_t cost_new = bits_new + 8 * hdr_new, cost_prev = prev_ok ? bits_prev : UINT64_MAX;
    if ((cost_new < cost_prev ? cost_new : cost_prev) / 8 + 13 + 8 >= n) {   /* store raw */
        out[0] = 0;
        memcpy(out + 1, lits, n);
        return n + 1;
    }
    size_t k = 0;
    const huf_table *t;
    if (cost_prev <= cost_new) {
        out[k++] = 3;
        t = prev;
    } else {
        out[k++] = 2;
        out[k++] = (uint8_t)last;
        for (unsigned s = 0; s <= last; s += 2)
            out[k++] = (uint8_t)(nt.len[s] | (s + 1 <= last ? nt.len[s + 1] << 4 : 0));
        huf_codes(&nt);
        nt.valid = 1;
        /* the decoder builds the same table (the encoder needs no dt) */
        memcpy(prev->len, nt.len, sizeof(nt.len));
        memcpy(prev->code, nt.code, sizeof(nt.code));
        prev->valid = 1;
        t = prev;
    }
    uint32_t ce[256];   /* code << 4 | length */
    for (unsigned s = 0; s < 256; s++) ce[s] = (uint32_t)t->code[s] << 4 | t->len[s];
    size_t q = (n + 3) / 4, szpos = k;
    k += 16;
    for (unsigned st = 0; st < 4; st++) {
        size_t a = st * q < n ? st * q : n, b = (st + 1) * q < n ? (st + 1) * q : n;
        /* branch-free: 8-byte writes, 4 codes (<= 44 bits) per flush; the
         * bytes are those of hbw_put (bytes past the end are zero, then
         * overwritten by the next stream; out has huf_bound slack) */
        uint64_t acc = 0;
        unsigned nb = 0;
        uint8_t *p = out + k;
        size_t i = a;
        for (; i + 4 <= b; i += 4) {
            uint32_t e0 = ce[lits[i]], e1 = ce[lits[i + 1]], e2 = ce[lits[i + 2]], e3 = ce[lits[i + 3]];
            acc |= (uint64_t)(e0 >> 4) << nb; nb += e0 & 15;
            acc |= (uint64_t)(e1 >> 4) << nb; nb += e1 & 15;
            acc |= (uint64_t)(e2 >> 4) << nb; nb += e2 & 15;
            acc |= (uint64_t)(e3 >> 4) << nb; nb += e3 & 15;
            wr64le(p, acc);
            unsigned by = nb >> 3;
            p += by;
            acc = (acc >> (by * 4)) >> (by * 4);
            nb &= 7;
        }
        hbw w = {acc, nb, p};
        for (; i < b; i++) hbw_put(&w, t->code[lits[i]], t->len[lits[i]]);
        hbw_flush(&w);
        size_t sz = (size_t)(w.p - (out + k));
        wr32le(out + szpos + 4 * st, (uint32_t)sz);
        k += sz;
    }
    return k;
}

/* ---- decoder ------------------------------------------------------ */

static inline uint64_t ld64(const uint8_t *p)
{
    uint64_t v;
    memcpy(&v, p, 8);
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    v = __builtin_bswap64(v);
#endif
    return v;
}

static inline __attribute__((always_inline)) size_t huf_decode_impl(const uint8_t *in, size_t avail, uint8_t *out, size_t n, huf_table *t, int v7)
{
    if (!n) return 0;
    if (avail < 2) return 0;
    size_t k = 1;
    if (in[0] == 0) {
        if (avail - 1 < n) return 0;
        memcpy(out, in + 1, n);
        return n + 1;
    }
    if (in[0] == 1) { memset(out, in[1], n); return 2; }
    if (in[0] == 2) {
        unsigned last = in[1];
        size_t hb = (last + 2) / 2;
        if (avail - 2 < hb) return 0;
        memset(t->len, 0, sizeof(t->len));
        for (unsigned s = 0; s <= last; s++) t->len[s] = (uint8_t)((in[2 + s / 2] >> (4 * (s & 1))) & 15);
        t->valid = 0;
        if (huf_build_dt(t, v7)) return 0;
        k = 2 + hb;
    } else if (in[0] != 3 || !t->valid) {
        return 0;
    }
    /* v7: 4 stream sizes; v5/v6: 3, and stream 3 ends where its last
     * literal ends */
    const size_t nsz = v7 ? 4 : 3;
    if (avail - k < 4 * nsz) return 0;
    size_t sz[4], total = 0;
    for (size_t s = 0; s < nsz; s++) { sz[s] = rd32le(in + k + 4 * s); total += sz[s]; }
    k += 4 * nsz;
    if (total > avail - k) return 0;
    if (!v7) sz[3] = avail - k - (sz[0] + sz[1] + sz[2]);

    const uint8_t *p[4], *end[4];
    unsigned pos[4] = {0, 0, 0, 0};
    const uint8_t *q0 = in + k;
    for (int s = 0; s < 4; s++) { p[s] = q0; q0 += s < 3 ? sz[s] : 0; }
    /* stream 3's size is not stored: it is found from where decoding stops */
    for (int s = 0; s < 4; s++) end[s] = p[s] + sz[s];

    const size_t q = (n + 3) / 4;
    size_t lo[4], hi[4];
    for (int s = 0; s < 4; s++) { lo[s] = (size_t)s * q < n ? (size_t)s * q : n; hi[s] = (size_t)(s + 1) * q < n ? (size_t)(s + 1) * q : n; }
    const uint16_t *dt = t->dt;
    const uint64_t mask = HUF_TSIZE - 1;

    uint8_t *o[4];
    for (int s = 0; s < 4; s++) o[s] = out + lo[s];
#define HSTEP(S)                                                               \
    do {                                                                       \
        uint16_t e_ = dt[(ld64(p[S]) >> pos[S]) & mask];                       \
        *o[S]++ = (uint8_t)e_;                                                 \
        pos[S] += e_ >> 8;                                                     \
    } while (0)
#define HREFILL(S) do { p[S] += pos[S] >> 3; pos[S] &= 7; } while (0)

    /* main loop: 4 symbols per stream per round (<= 44 + 7 bits, so one
     * 64-bit load per stream and round), streams interleaved for
     * instruction-level parallelism. Everything lives in scalars: a byte
     * store may alias anything, so values in memory would be reloaded. */
    size_t rounds = (hi[3] - lo[3]) / 4, i;
    {
        const uint8_t *p0 = p[0], *p1 = p[1], *p2 = p[2], *p3 = p[3];
        unsigned s0 = 0, s1 = 0, s2 = 0, s3 = 0;
        uint8_t *o0 = o[0], *o1 = o[1], *o2 = o[2], *o3 = o[3];
#define HSYM(C, S, O, K)                                                       \
    do {                                                                       \
        uint16_t e_ = dt[(C >> S) & mask];                                     \
        O[K] = (uint8_t)e_;                                                    \
        S += e_ >> 8;                                                          \
    } while (0)
        for (i = 0; i < rounds; i++) {
            const uint64_t c0 = ld64(p0), c1 = ld64(p1), c2 = ld64(p2), c3 = ld64(p3);
            HSYM(c0, s0, o0, 0); HSYM(c1, s1, o1, 0); HSYM(c2, s2, o2, 0); HSYM(c3, s3, o3, 0);
            HSYM(c0, s0, o0, 1); HSYM(c1, s1, o1, 1); HSYM(c2, s2, o2, 1); HSYM(c3, s3, o3, 1);
            HSYM(c0, s0, o0, 2); HSYM(c1, s1, o1, 2); HSYM(c2, s2, o2, 2); HSYM(c3, s3, o3, 2);
            HSYM(c0, s0, o0, 3); HSYM(c1, s1, o1, 3); HSYM(c2, s2, o2, 3); HSYM(c3, s3, o3, 3);
            o0 += 4; o1 += 4; o2 += 4; o3 += 4;
            p0 += s0 >> 3; s0 &= 7; p1 += s1 >> 3; s1 &= 7;
            p2 += s2 >> 3; s2 &= 7; p3 += s3 >> 3; s3 &= 7;
            if (p0 > end[0] || p1 > end[1] || p2 > end[2] || p3 > end[3]) return 0;
        }
#undef HSYM
        p[0] = p0; p[1] = p1; p[2] = p2; p[3] = p3;
        pos[0] = s0; pos[1] = s1; pos[2] = s2; pos[3] = s3;
        o[0] = o0; o[1] = o1; o[2] = o2; o[3] = o3;
    }
    for (int s = 0; s < 4; s++) {
        while (o[s] < out + hi[s]) {
            HSTEP(s);
            HREFILL(s);
            if (p[s] > end[s]) return 0;
        }
    }
#undef HSTEP
#undef HREFILL
    /* each stream must end within its last byte */
    for (int s = 0; s < 4; s++) {
        const uint8_t *e = p[s] + (pos[s] + 7) / 8;
        if ((s < 3 || v7) ? e != end[s] : e > end[s]) return 0;
    }
    return (size_t)(p[3] + (pos[3] + 7) / 8 - in);
}

static __attribute__((noinline)) size_t huf_decode_generic(const uint8_t *in, size_t avail, uint8_t *out, size_t n, huf_table *t, int v7)
{
    return huf_decode_impl(in, avail, out, n, t, v7);
}

#ifdef LZR_BMI2_DISPATCH
static __attribute__((noinline)) LZR_TARGET_BMI2 size_t huf_decode_bmi2(const uint8_t *in, size_t avail, uint8_t *out, size_t n, huf_table *t, int v7)
{
    return huf_decode_impl(in, avail, out, n, t, v7);
}
#endif

size_t huf_decode(const uint8_t *in, size_t avail, uint8_t *out, size_t n, huf_table *t, int v7)
{
#ifdef LZR_BMI2_DISPATCH
    if (LZR_HAS_BMI2()) return huf_decode_bmi2(in, avail, out, n, t, v7);
#endif
    return huf_decode_generic(in, avail, out, n, t, v7);
}

/* ---- context literals (block type 6) ------------------------------ */

#include <math.h>

#ifndef HUFC_SMOOTH
#define HUFC_SMOOTH 1
#endif

void hufc_free(hufc_table *st)
{
    free(st->work);
    st->work = NULL;
}

/* Bytes of a cluster table in the section: last symbol + nibbles. */
static unsigned hufc_tab_bytes(const uint8_t *len)
{
    unsigned last = 0;
    for (unsigned s = 0; s < 256; s++) if (len[s]) last = s;
    return 1 + (last + 2) / 2;
}

size_t hufc_encode(uint8_t *out, const uint8_t *lits, const uint8_t *prev, size_t n, hufc_table *st)
{
    if (!n) return 0;
    /* work: counts by (previous byte, literal), then the nonzero ones of
     * each previous byte as (symbol, count) lists */
    if (!st->work && !(st->work = (uint32_t *)malloc(HUFC_WORK * sizeof(uint32_t)))) return 0;
    uint32_t *H = st->work, *ls = st->work + 65536, *lc = st->work + 2 * 65536;
    memset(H, 0, 65536 * sizeof(uint32_t));
    uint32_t pc[256] = {0}, tot[256] = {0};
    for (size_t i = 0; i < n; i++) H[(unsigned)prev[i] << 8 | lits[i]]++;
    uint32_t off[257];
    unsigned np = 0;
    off[0] = 0;
    for (unsigned p = 0; p < 256; p++) {
        uint32_t m = off[p];
        for (unsigned s = 0; s < 256; s++) {
            uint32_t v = H[p << 8 | s];
            if (!v) continue;
            ls[m] = s; lc[m] = v; m++;
            pc[p] += v; tot[s] += v;
        }
        off[p + 1] = m;
        if (pc[p]) np++;
    }
    /* a single literal value: type 5 codes it as a run (and a one-symbol
     * Huffman code would be incomplete, which format v7 forbids) */
    unsigned ndistinct = 0;
    for (unsigned s = 0; s < 256; s++) ndistinct += tot[s] != 0;
    if (ndistinct < 2) return 0;
    /* one table (what type 5 would use, without its 4-stream sizes) */
    uint8_t len1[256];
    huf_lengths(tot, len1);
    uint64_t cost1 = 8ull * hufc_tab_bytes(len1);
    for (unsigned s = 0; s < 256; s++) cost1 += (uint64_t)tot[s] * len1[s];
    if (np < 2) return 0;

    /* candidates: K = 2, 4, 8 clusters by a few rounds of assigning each
     * previous byte to the cluster whose statistics code it cheapest,
     * starting from the top bits of the previous byte */
    uint8_t best_map[256], best_len[HUFC_MAXK][256], best_k = 0;
    uint64_t best = UINT64_MAX;
    for (unsigned cand = 0; cand < 4; cand++) {
        const unsigned kl = cand < 3 ? cand + 1 : 3, K = 1u << kl;
        uint8_t a[256];
        for (unsigned p = 0; p < 256; p++) a[p] = (uint8_t)(p >> (8 - kl));
        if (cand == 3)   /* text: start from character classes */
            for (unsigned p = 0; p < 256; p++)
                a[p] = p == ' ' ? 0 : p >= 'a' && p <= 'z' ? 1 : p >= 'A' && p <= 'Z' ? 2 : p >= '0' && p <= '9' ? 3
                     : p == '\n' || p == '\r' ? 4 : p >= 128 ? 5 : p < 32 ? 6 : 7;
        static const unsigned ROUNDS = 4;
        uint32_t C[HUFC_MAXK][256], N[HUFC_MAXK];
        for (unsigned r = 0; r <= ROUNDS; r++) {
            memset(C, 0, sizeof(C));
            memset(N, 0, sizeof(N));
            for (unsigned p = 0; p < 256; p++)
                for (uint32_t m = off[p]; m < off[p + 1]; m++) C[a[p]][ls[m]] += lc[m];
            for (unsigned k = 0; k < K; k++) for (unsigned s = 0; s < 256; s++) N[k] += C[k][s];
            if (r == ROUNDS) break;
            float L[HUFC_MAXK][256];
            for (unsigned k = 0; k < K; k++) {
                float den = (float)N[k] + 0.25f * 256;
                for (unsigned s = 0; s < 256; s++) L[k][s] = -log2f(((float)C[k][s] + 0.25f) / den);
            }
            int moved = 0;
            for (unsigned p = 0; p < 256; p++) {
                if (!pc[p]) continue;
                float bc = 0;
                unsigned bk = a[p];
                for (uint32_t m = off[p]; m < off[p + 1]; m++) bc += (float)lc[m] * L[bk][ls[m]];
                for (unsigned k = 0; k < K; k++) {
                    if (k == a[p]) continue;
                    float c = 0;
                    for (uint32_t m = off[p]; m < off[p + 1] && c < bc; m++) c += (float)lc[m] * L[k][ls[m]];
                    if (c < bc) { bc = c; bk = k; }
                }
                if (bk != a[p]) { a[p] = (uint8_t)bk; moved = 1; }
            }
            if (!moved) r = ROUNDS - 1;
        }
        /* drop empty clusters, exact cost */
        uint8_t ren[HUFC_MAXK], kk = 0, lens[HUFC_MAXK][256];
        uint64_t cost = 8ull * (2 + 128);
        for (unsigned k = 0; k < K; k++) {
            if (!N[k]) { ren[k] = 0; continue; }
            ren[k] = kk;
            /* every symbol of the block gets a code in every cluster, so
             * the tables can be repeated by the next blocks */
            uint32_t Cs[256];
            for (unsigned s = 0; s < 256; s++) Cs[s] = C[k][s] ? C[k][s] : (tot[s] ? HUFC_SMOOTH : 0);
            huf_lengths(Cs, lens[kk]);
            cost += 8ull * hufc_tab_bytes(lens[kk]);
            for (unsigned s = 0; s < 256; s++) cost += (uint64_t)C[k][s] * lens[kk][s];
            kk++;
        }
        if (cost < best) {
            best = cost;
            best_k = kk;
            for (unsigned p = 0; p < 256; p++) best_map[p] = ren[a[p]];
            memcpy(best_len, lens, sizeof(lens));
        }
    }
    /* the tables the decoder holds */
    int rep = 0;
    if (st->valid) {
        uint64_t cost = 8;
        for (unsigned p = 0; p < 256 && cost != UINT64_MAX; p++) {
            const uint8_t *l = st->len[st->cmap[p]];
            for (uint32_t m = off[p]; m < off[p + 1]; m++) {
                if (!l[ls[m]]) { cost = UINT64_MAX; break; }
                cost += (uint64_t)lc[m] * l[ls[m]];
            }
        }
        if (cost <= best) { best = cost; rep = 1; }
    }
    if (best >= cost1) return 0;

    size_t k = 0;
    if (rep) {
        out[k++] = 1;
    } else {
        out[k++] = 0;
        out[k++] = best_k;
        for (unsigned p = 0; p < 256; p += 2) out[k++] = (uint8_t)(best_map[p] | best_map[p + 1] << 4);
        for (unsigned c = 0; c < best_k; c++) {
            unsigned last = 0;
            for (unsigned s = 0; s < 256; s++) if (best_len[c][s]) last = s;
            out[k++] = (uint8_t)last;
            for (unsigned s = 0; s <= last; s += 2)
                out[k++] = (uint8_t)(best_len[c][s] | (s + 1 <= last ? best_len[c][s + 1] << 4 : 0));
            huf_table t;
            memcpy(t.len, best_len[c], 256);
            huf_codes(&t);
            memcpy(st->len[c], t.len, 256);
            memcpy(st->code[c], t.code, sizeof(t.code));
        }
        memcpy(st->cmap, best_map, 256);
        st->k = best_k;
        st->valid = 1;
    }
    size_t szpos = k;
    k += 4;
    hbw w = {0, 0, out + k};
    for (size_t i = 0; i < n; i++) {
        unsigned c = st->cmap[prev[i]];
        hbw_put(&w, st->code[c][lits[i]], st->len[c][lits[i]]);
    }
    hbw_flush(&w);
    size_t sz = (size_t)(w.p - (out + k));
    wr32le(out + szpos, (uint32_t)sz);
    return k + sz;
}

size_t hufc_read(const uint8_t *in, size_t avail, hufc_dtable *t, int strict)
{
    if (avail < 1) return 0;
    if (in[0] == 1) return t->valid ? 1 : 0;
    if (in[0] != 0 || avail < 2 + 128) return 0;
    unsigned K = in[1];
    if (!K || K > HUFC_MAXK) return 0;
    t->valid = 0;
    for (unsigned p = 0; p < 256; p++) {
        unsigned c = (in[2 + p / 2] >> (4 * (p & 1))) & 15;
        if (c >= K) return 0;
        t->cmap[p] = (uint8_t)c;
    }
    size_t k = 2 + 128;
    huf_table h;
    for (unsigned c = 0; c < K; c++) {
        if (k >= avail) return 0;
        unsigned last = in[k];
        size_t hb = (last + 2) / 2;
        if (avail - k - 1 < hb) return 0;
        memset(h.len, 0, sizeof(h.len));
        for (unsigned s = 0; s <= last; s++) h.len[s] = (uint8_t)((in[k + 1 + s / 2] >> (4 * (s & 1))) & 15);
        if (huf_build_dt(&h, strict)) return 0;
        k += 1 + hb;
        uint16_t *dt = t->dt + ((size_t)c << HUF_MAXLEN);
        for (unsigned j = 0; j < HUF_TSIZE; j++) {
            unsigned e = h.dt[j];
            dt[j] = (uint16_t)(e | (unsigned)t->cmap[e & 255] << 12);
        }
    }
    t->valid = 1;
    return k;
}
