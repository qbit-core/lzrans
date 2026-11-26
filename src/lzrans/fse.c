/*
 * fse.c - tANS tables: normalisation, serialisation, construction.
 */
#include "fse.h"
#include <math.h>

void fse_normalize(const uint32_t *cnt, unsigned alphabet, uint32_t total, uint16_t *norm)
{
    int64_t sum = 0;
    for (unsigned s = 0; s < alphabet; s++) {
        if (!cnt[s]) { norm[s] = 0; continue; }
        uint32_t f = (uint32_t)(((uint64_t)cnt[s] * FSE_L + total / 2) / total);
        norm[s] = (uint16_t)(f ? f : 1);
        sum += norm[s];
    }
    /* fix rounding where it costs the fewest bits */
    while (sum != FSE_L) {
        int best = -1;
        double best_cost = 0;
        for (unsigned s = 0; s < alphabet; s++) {
            if (!cnt[s]) continue;
            double c;
            if (sum > FSE_L) {
                if (norm[s] <= 1) continue;
                c = cnt[s] * log2((double)norm[s] / (norm[s] - 1));
                if (best < 0 || c < best_cost) { best = (int)s; best_cost = c; }
            } else {
                c = cnt[s] * log2((double)(norm[s] + 1) / norm[s]);
                if (best < 0 || c > best_cost) { best = (int)s; best_cost = c; }
            }
        }
        if (sum > FSE_L) { norm[best]--; sum--; }
        else { norm[best]++; sum++; }
    }
}

int fse_build(fse_table *t)
{
    uint8_t spread[FSE_L];
    unsigned pos = 0, total = 0;
    const unsigned step = (FSE_L >> 1) + (FSE_L >> 3) + 3, mask = FSE_L - 1;
    if (t->alphabet > FSE_MAXSYM) return -1;
    for (unsigned s = 0; s < t->alphabet; s++) {
        total += t->norm[s];
        if (total > FSE_L) return -1;
        for (unsigned i = 0; i < t->norm[s]; i++) { spread[pos] = (uint8_t)s; pos = (pos + step) & mask; }
    }
    if (total != FSE_L) return -1;
    uint16_t next[FSE_MAXSYM];
    unsigned cum = 0;
    for (unsigned s = 0; s < t->alphabet; s++) {
        next[s] = t->norm[s];
        t->first[s] = (uint16_t)cum;
        t->hbf[s] = t->norm[s] ? (uint8_t)(31 - __builtin_clz(t->norm[s])) : 0;
        cum += t->norm[s];
    }
    for (unsigned u = 0; u < FSE_L; u++) {
        unsigned s = spread[u], x = next[s]++;
        unsigned nb = FSE_LOG - (31 - __builtin_clz(x));
        t->dt[u].sym = (uint8_t)s;
        t->dt[u].nb = (uint8_t)nb;
        t->dt[u].base = (uint16_t)((x << nb) - FSE_L);
        t->st[t->first[s] + x - t->norm[s]] = (uint16_t)(FSE_L + u);
    }
    t->valid = 1;
    return 0;
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
    for (size_t k = 0; k < avail && k < 3; k++) {
        r |= (uint32_t)(in[k] & 0x7F) << (7 * k);
        if (!(in[k] & 0x80)) { *v = r; return k + 1; }
    }
    return 0;
}

/* f > 0: count of the next symbol; 0, r: r+1 absent symbols. */
size_t fse_write_norm(uint8_t *out, const uint16_t *norm, unsigned alphabet)
{
    size_t k = 0;
    unsigned s = 0;
    while (s < alphabet) {
        if (norm[s]) { k += put_varint(out + k, norm[s]); s++; continue; }
        unsigned r = 0;
        while (s + r < alphabet && !norm[s + r]) r++;
        out[k++] = 0;
        k += put_varint(out + k, r - 1);
        s += r;
    }
    return k;
}

size_t fse_read_norm(const uint8_t *in, size_t avail, uint16_t *norm, unsigned alphabet)
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
            for (uint32_t i = 0; i <= r; i++) norm[s + i] = 0;
            s += r + 1;
            continue;
        }
        if (v > FSE_L || cum + v > FSE_L) return 0;
        norm[s++] = (uint16_t)v;
        cum += v;
    }
    return cum == FSE_L ? k : 0;
}

uint64_t fse_cost(const uint32_t *cnt, const uint16_t *norm, unsigned alphabet)
{
    double bits = 0;
    for (unsigned s = 0; s < alphabet; s++) {
        if (!cnt[s]) continue;
        if (!norm[s]) return UINT64_MAX;
        bits += cnt[s] * (FSE_LOG - log2((double)norm[s]));
    }
    return (uint64_t)(bits * 256);
}
