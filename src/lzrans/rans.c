/*
 * rans.c - frequency modelling, table (de)serialisation, rANS encoder.
 * The decoder hot loop lives in frame.c so the three streams can be
 * consumed in token order without function-call overhead.
 */
#include "rans.h"
#include <stdlib.h>

void rans_build_freqs(const uint16_t *syms, size_t n, unsigned alphabet, uint32_t *freq)
{
    uint64_t cnt[SYM_LITLEN];
    memset(cnt, 0, sizeof(cnt));
    memset(freq, 0, alphabet * sizeof(uint32_t));
    for (size_t i = 0; i < n; i++) cnt[syms[i]]++;
    if (n == 0) return;

    uint32_t sum = 0;
    unsigned maxi = 0;
    for (unsigned s = 0; s < alphabet; s++) {
        if (!cnt[s]) continue;
        uint32_t f = (uint32_t)((cnt[s] * RANS_M) / n);
        if (f == 0) f = 1;
        freq[s] = f;
        sum += f;
        if (cnt[s] > cnt[maxi]) maxi = s;
    }
    if (sum < RANS_M) {
        freq[maxi] += RANS_M - sum;          /* give slack to the most likely symbol */
    } else {
        /* Over budget because of the f>=1 floor: shave from the largest. */
        while (sum > RANS_M) {
            unsigned big = 0;
            for (unsigned s = 1; s < alphabet; s++) if (freq[s] > freq[big]) big = s;
            uint32_t take = sum - RANS_M;
            if (take > freq[big] - 1) take = freq[big] - 1;
            freq[big] -= take;
            sum -= take;
        }
    }
}

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
    for (size_t k = 0; k < avail && k < 3; k++) {   /* values <= 4096 need <= 2 bytes */
        r |= (uint32_t)(in[k] & 0x7F) << (7 * k);
        if (!(in[k] & 0x80)) { *v = r; return k + 1; }
    }
    return 0;
}

/* Format: sequence of varints. f>0 => that frequency for the next symbol.
 * 0 followed by varint(r) => r+1 consecutive zero-frequency symbols. */
size_t rans_write_freqs(uint8_t *out, const uint32_t *freq, unsigned alphabet)
{
    size_t k = 0;
    unsigned s = 0;
    while (s < alphabet) {
        if (freq[s]) { k += put_varint(out + k, freq[s]); s++; continue; }
        unsigned r = 0;
        while (s + r < alphabet && freq[s + r] == 0) r++;
        out[k++] = 0;
        k += put_varint(out + k, r - 1);
        s += r;
    }
    return k;
}

size_t rans_read_freqs(const uint8_t *in, size_t avail, unsigned alphabet, rans_slot *slots)
{
    size_t k = 0;
    unsigned s = 0;
    uint32_t cum = 0;
    while (s < alphabet) {
        uint32_t v;
        size_t c = get_varint(in + k, avail - k, &v);
        if (!c) return 0;
        k += c;
        if (v == 0) {
            uint32_t r;
            c = get_varint(in + k, avail - k, &r);
            if (!c) return 0;
            k += c;
            if ((uint64_t)s + r + 1 > alphabet) return 0;
            s += r + 1;
        } else {
            if (v > RANS_M || cum + v > RANS_M) return 0;
            for (uint32_t i = 0; i < v; i++) {
                rans_slot *e = &slots[cum + i];
                e->sym = (uint16_t)s;
                e->freq = (uint16_t)v;
                e->bias = (uint16_t)i;   /* slot - cum */
                e->pad = 0;
            }
            cum += v;
            s++;
        }
    }
    return cum == RANS_M ? k : 0;
}

size_t rans_encode(const uint16_t *syms, size_t n, const uint32_t *freq, unsigned alphabet,
                   uint8_t *scratch, size_t scratch_cap, uint8_t *out)
{
    uint32_t cum[SYM_LITLEN + 1];
    cum[0] = 0;
    for (unsigned s = 0; s < alphabet; s++) cum[s + 1] = cum[s] + freq[s];

    /* rANS is LIFO: encode backwards, writing bytes from the end of scratch. */
    uint8_t *ptr = scratch + scratch_cap;
    uint32_t x = RANS_L;
    for (size_t i = n; i-- > 0;) {
        unsigned s = syms[i];
        uint32_t f = freq[s];
        uint32_t x_max = ((RANS_L >> RANS_SCALE_BITS) << 8) * f;
        while (x >= x_max) { *--ptr = (uint8_t)x; x >>= 8; }
        x = ((x / f) << RANS_SCALE_BITS) + (x % f) + cum[s];
    }
    ptr -= 4;
    wr32le(ptr, x);
    size_t len = (size_t)(scratch + scratch_cap - ptr);
    memcpy(out, ptr, len);
    return len;
}
