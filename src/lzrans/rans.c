/*
 * rans.c - frequency normalisation, table (de)serialisation, and the
 * 2-way interleaved rANS stream encoder/decoder.
 */
#include "rans.h"
#include <math.h>

/* ---- frequency normalisation ------------------------------------- */

/* Scale counts to sum exactly RANS_M, keeping every present symbol >= 1.
 * Rounding errors are fixed greedily by the adjustment that costs the
 * fewest bits (cnt * log2 ratio), which is close to optimal. Requires at
 * least two distinct symbols, so every freq ends up < RANS_M. */
static void normalize(const uint32_t *cnt, unsigned alphabet, uint64_t total, uint32_t *freq)
{
    int64_t sum = 0;
    for (unsigned s = 0; s < alphabet; s++) {
        if (!cnt[s]) { freq[s] = 0; continue; }
        uint32_t f = (uint32_t)(((uint64_t)cnt[s] * RANS_M + total / 2) / total);
        freq[s] = f ? f : 1;
        sum += freq[s];
    }
    while (sum != RANS_M) {
        int best = -1;
        double best_cost = 0;
        if (sum > RANS_M) {          /* remove one unit where it hurts least */
            for (unsigned s = 0; s < alphabet; s++) {
                if (freq[s] <= 1) continue;
                double c = cnt[s] * log2((double)freq[s] / (freq[s] - 1));
                if (best < 0 || c < best_cost) { best = (int)s; best_cost = c; }
            }
            freq[best]--;
            sum--;
        } else {                     /* add one unit where it helps most */
            for (unsigned s = 0; s < alphabet; s++) {
                if (!cnt[s]) continue;
                double c = cnt[s] * log2((double)(freq[s] + 1) / freq[s]);
                if (best < 0 || c > best_cost) { best = (int)s; best_cost = c; }
            }
            freq[best]++;
            sum++;
        }
    }
}

/* ---- table serialisation (zero-run + varint) --------------------- */

static size_t put_varint(uint8_t *o, uint32_t v)
{
    size_t k = 0;
    while (v >= 0x80) { o[k++] = (uint8_t)(v | 0x80); v >>= 7; }
    o[k++] = (uint8_t)v;
    return k;
}

static size_t get_varint(const uint8_t *in, size_t avail, uint32_t *v)
{
    uint32_t r = 0;
    for (size_t k = 0; k < avail && k < 3; k++) {
        r |= (uint32_t)(in[k] & 0x7F) << (7 * k);
        if (!(in[k] & 0x80)) { *v = r; return k + 1; }
    }
    return 0;
}

/* f > 0: frequency of the next symbol. 0, r: r+1 absent symbols. */
static size_t write_freqs(uint8_t *out, const uint32_t *freq, unsigned alphabet)
{
    size_t k = 0;
    unsigned s = 0;
    while (s < alphabet) {
        if (freq[s]) { k += put_varint(out + k, freq[s]); s++; continue; }
        unsigned r = 0;
        while (s + r < alphabet && !freq[s + r]) r++;
        out[k++] = 0;
        k += put_varint(out + k, r - 1);
        s += r;
    }
    return k;
}

static size_t read_freqs(const uint8_t *in, size_t avail, unsigned alphabet, rans_dslot *tbl)
{
    size_t k = 0;
    unsigned s = 0;
    uint32_t cum = 0;
    while (s < alphabet) {
        uint32_t v, r;
        size_t c = get_varint(in + k, avail - k, &v);
        if (!c) return 0;
        k += c;
        if (!v) {
            c = get_varint(in + k, avail - k, &r);
            if (!c || (uint64_t)s + r + 1 > alphabet) return 0;
            k += c;
            s += r + 1;
            continue;
        }
        if (v >= RANS_M || cum + v > RANS_M) return 0;
        for (uint32_t i = 0; i < v; i++) tbl[cum + i] = s | (v << 8) | (i << 20);
        cum += v;
        s++;
    }
    return cum == RANS_M ? k : 0;
}

/* ---- stream encode ------------------------------------------------ */

size_t stream_encode(uint8_t *out, const uint8_t *syms, size_t n, unsigned alphabet, uint8_t *scratch)
{
    if (!n) return 0;
    uint32_t cnt[256] = {0};
    for (size_t i = 0; i < n; i++) cnt[syms[i]]++;
    unsigned distinct = 0;
    for (unsigned s = 0; s < alphabet; s++) distinct += cnt[s] != 0;
    if (distinct == 1) { out[0] = 1; out[1] = syms[0]; return 2; }

    uint32_t freq[256], cum[257];
    normalize(cnt, alphabet, n, freq);
    cum[0] = 0;
    for (unsigned s = 0; s < alphabet; s++) cum[s + 1] = cum[s] + freq[s];

    /* Encode backwards with two states: symbol i uses state i & 1. */
    uint8_t *end = scratch + 2 * n + 16, *ptr = end;
    uint32_t x[2] = {RANS_L, RANS_L};
    for (size_t i = n; i-- > 0;) {
        uint32_t *xs = &x[i & 1];
        unsigned s = syms[i];
        uint32_t f = freq[s], v = *xs;
        uint32_t x_max = ((RANS_L >> RANS_SCALE_BITS) << 8) * f;
        while (v >= x_max) { *--ptr = (uint8_t)v; v >>= 8; }
        *xs = ((v / f) << RANS_SCALE_BITS) + (v % f) + cum[s];
    }
    ptr -= 4; wr32le(ptr, x[1]);
    ptr -= 4; wr32le(ptr, x[0]);
    size_t plen = (size_t)(end - ptr);

    if (plen + 8 >= n && alphabet == 256) {   /* not worth it: store raw */
        out[0] = 0;
        memcpy(out + 1, syms, n);
        return n + 1;
    }
    size_t k = 0;
    out[k++] = 2;
    k += write_freqs(out + k, freq, alphabet);
    wr32le(out + k, (uint32_t)plen);
    k += 4;
    memcpy(out + k, ptr, plen);
    return k + plen;
}

/* ---- stream decode ------------------------------------------------ */

size_t stream_decode(const uint8_t *in, size_t avail, uint8_t *out, size_t n, unsigned alphabet,
                     rans_dslot *tbl)
{
    if (!n) return 0;
    if (avail < 2) return 0;
    if (in[0] == 0) {
        if (avail - 1 < n) return 0;
        memcpy(out, in + 1, n);
        for (size_t i = 0; i < n; i++) if (out[i] >= alphabet) return 0;
        return n + 1;
    }
    if (in[0] == 1) {
        if (in[1] >= alphabet) return 0;
        memset(out, in[1], n);
        return 2;
    }
    if (in[0] != 2) return 0;

    size_t k = 1 + read_freqs(in + 1, avail - 1, alphabet, tbl);
    if (k == 1 || avail - k < 4) return 0;
    uint32_t plen = rd32le(in + k);
    k += 4;
    if (plen < 8 || plen > avail - k) return 0;
    const uint8_t *p = in + k, *end = p + plen;
    uint32_t x0 = rd32le(p), x1 = rd32le(p + 4);
    p += 8;
    if (x0 < RANS_L || x1 < RANS_L) return 0;

    /* Hot loop: one table lookup per symbol, two independent states.
     * Renormalisation reads are bounds-checked; an overrun feeds zeros and
     * is caught by the final-state check below. */
#define RANS_STEP(X, I)                                                     \
    do {                                                                    \
        rans_dslot e = tbl[(X) & (RANS_M - 1)];                             \
        out[I] = (uint8_t)e;                                                \
        X = ((e >> 8) & 0xFFF) * ((X) >> RANS_SCALE_BITS) + (e >> 20);      \
        while (X < RANS_L) X = (X << 8) | (p < end ? *p++ : (bad = 1, 0));  \
    } while (0)

    int bad = 0;
    size_t i = 0;
    for (; i + 1 < n; i += 2) {
        RANS_STEP(x0, i);
        RANS_STEP(x1, i + 1);
    }
    if (i < n) RANS_STEP(x0, i);
#undef RANS_STEP
    if (bad || p != end || x0 != RANS_L || x1 != RANS_L) return 0;
    return k + plen;
}
