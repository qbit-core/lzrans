/*
 * huf.c - length-limited canonical Huffman coding of literals, 4 streams
 * (format in huf.h).
 */
#include "huf.h"

/* ---- code construction ------------------------------------------- */

/* Code lengths (<= HUF_MAXLEN) for the counts; absent symbols get 0. At
 * least one symbol must be present. */
void huf_lengths(const uint32_t *cnt, uint8_t *len)
{
    uint16_t sym[256];
    unsigned n = 0;
    memset(len, 0, 256);
    for (unsigned s = 0; s < 256; s++) if (cnt[s]) sym[n++] = (uint16_t)s;
    if (n == 1) { len[sym[0]] = 1; return; }
    /* sort by count ascending (insertion sort: n <= 256, once per block) */
    for (unsigned i = 1; i < n; i++) {
        uint16_t v = sym[i];
        unsigned j = i;
        while (j && cnt[sym[j - 1]] > cnt[v]) { sym[j] = sym[j - 1]; j--; }
        sym[j] = v;
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
static int huf_build_dt(huf_table *t)
{
    uint32_t kraft = 0;
    for (unsigned s = 0; s < 256; s++) {
        if (t->len[s] > HUF_MAXLEN) return -1;
        if (t->len[s]) kraft += HUF_TSIZE >> t->len[s];
    }
    if (kraft > HUF_TSIZE || !kraft) return -1;
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
    uint32_t cnt[256] = {0};
    for (size_t i = 0; i < n; i++) cnt[lits[i]]++;
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
        *prev = nt;       /* the decoder builds the same table */
        huf_build_dt(prev);
        t = prev;
    }
    size_t q = (n + 3) / 4, szpos = k;
    k += 12;
    for (unsigned st = 0; st < 4; st++) {
        size_t a = st * q < n ? st * q : n, b = (st + 1) * q < n ? (st + 1) * q : n;
        hbw w = {0, 0, out + k};
        for (size_t i = a; i < b; i++) hbw_put(&w, t->code[lits[i]], t->len[lits[i]]);
        hbw_flush(&w);
        size_t sz = (size_t)(w.p - (out + k));
        if (st < 3) wr32le(out + szpos + 4 * st, (uint32_t)sz);
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

static inline __attribute__((always_inline)) size_t huf_decode_impl(const uint8_t *in, size_t avail, uint8_t *out, size_t n, huf_table *t)
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
        if (huf_build_dt(t)) return 0;
        k = 2 + hb;
    } else if (in[0] != 3 || !t->valid) {
        return 0;
    }
    if (avail - k < 12) return 0;
    size_t sz[4], total = 0;
    for (int s = 0; s < 3; s++) { sz[s] = rd32le(in + k + 4 * s); total += sz[s]; }
    k += 12;
    if (total > avail - k) return 0;
    sz[3] = avail - k - total;   /* provisional: the last stream ends at the section end */

    const uint8_t *p[4], *end[4];
    unsigned pos[4] = {0, 0, 0, 0};
    const uint8_t *q0 = in + k;
    for (int s = 0; s < 4; s++) { p[s] = q0; q0 += s < 3 ? sz[s] : 0; }
    /* stream 3's size is not stored: it is found from where decoding stops */
    for (int s = 0; s < 3; s++) end[s] = p[s] + sz[s];
    end[3] = in + avail;

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
        if (s < 3 ? e != end[s] : e > end[s]) return 0;
    }
    return (size_t)(p[3] + (pos[3] + 7) / 8 - in);
}

static __attribute__((noinline)) size_t huf_decode_generic(const uint8_t *in, size_t avail, uint8_t *out, size_t n, huf_table *t)
{
    return huf_decode_impl(in, avail, out, n, t);
}

#ifdef LZR_BMI2_DISPATCH
static __attribute__((noinline)) LZR_TARGET_BMI2 size_t huf_decode_bmi2(const uint8_t *in, size_t avail, uint8_t *out, size_t n, huf_table *t)
{
    return huf_decode_impl(in, avail, out, n, t);
}
#endif

size_t huf_decode(const uint8_t *in, size_t avail, uint8_t *out, size_t n, huf_table *t)
{
#ifdef LZR_BMI2_DISPATCH
    if (LZR_HAS_BMI2()) return huf_decode_bmi2(in, avail, out, n, t);
#endif
    return huf_decode_generic(in, avail, out, n, t);
}
