/*
 * lz77.c - 4-byte hash-chain LZ77 over a 64 KB sliding window.
 *
 * Structure:
 *   hash4()       multiplicative hash of the next 4 bytes -> HASH_BITS (16)
 *   insert()      push a position onto its hash chain
 *   find_match()  walk the chain, compare 32 bits at a time, keep longest
 *   lz_parse()    greedy parse with optional one-step lazy evaluation
 */
#include "lz77.h"

#define WMASK (LZ_WINDOW - 1)

static inline uint32_t hash4(const uint8_t *p)
{
    return (rd32(p) * 2654435761u) >> (32 - LZ_HASH_BITS);
}

void lz_init(lz_state *st, int level)
{
    for (size_t i = 0; i < LZ_HASH_SIZE; i++) st->head[i] = -1;
    /* prev[] entries are always written before being followed, and chain
     * walks are bounded by the window check, so no init is needed;
     * clear anyway for deterministic behaviour under sanitizers. */
    for (size_t i = 0; i < LZ_WINDOW; i++) st->prev[i] = -1;
    if (level < 1) level = 1;
    if (level > 9) level = 9;
    static const int chains[10] = {0, 4, 8, 16, 32, 64, 128, 256, 1024, 4096};
    st->max_chain = chains[level];
    st->lazy = level >= 4;
}

static inline void insert(lz_state *st, const uint8_t *buf, int64_t pos)
{
    uint32_t h = hash4(buf + pos);
    st->prev[pos & WMASK] = st->head[h];
    st->head[h] = pos;
}

/* Longest match for `pos` limited to `maxlen`. Returns length (0 if < MIN). */
static inline size_t find_match(const lz_state *st, const uint8_t *buf, int64_t pos,
                                size_t maxlen, size_t *out_dist)
{
    if (maxlen < LZ_MIN_MATCH) return 0;
    const uint8_t *cur = buf + pos;
    const uint32_t first = rd32(cur);
    size_t best = LZ_MIN_MATCH - 1;
    int64_t cand = st->head[hash4(cur)];
    int chain = st->max_chain;

    while (cand >= 0 && pos - cand <= (int64_t)LZ_WINDOW && chain-- > 0) {
        const uint8_t *ref = buf + cand;
        /* Cheap reject: byte just past current best must match, then the
         * first 32-bit word (also filters hash collisions). */
        if (ref[best] == cur[best] && rd32(ref) == first) {
            size_t len = 4;
            while (len + 4 <= maxlen && rd32(ref + len) == rd32(cur + len)) len += 4;
            while (len < maxlen && ref[len] == cur[len]) len++;
            if (len > best) {
                best = len;
                *out_dist = (size_t)(pos - cand);
                if (len == maxlen) break;
            }
        }
        int64_t next = st->prev[cand & WMASK];
        if (next >= cand) break;   /* slot was recycled: chain ends here */
        cand = next;
    }
    return best >= LZ_MIN_MATCH ? best : 0;
}

size_t lz_parse(lz_state *st, const uint8_t *buf, size_t n, size_t start, size_t end,
                uint16_t *sym, uint16_t *dist, size_t *nmatch)
{
    size_t ns = 0, nm = 0;
    size_t pos = start;
    /* Hashing needs 4 readable bytes; matches must stop at `end`. */
    const size_t hash_end = n >= 4 ? n - 3 : 0;   /* pos < hash_end => pos+4 <= n */

    while (pos < end) {
        size_t maxlen = end - pos;
        if (maxlen > LZ_MAX_MATCH) maxlen = LZ_MAX_MATCH;
        size_t d = 0, len = 0;
        if (pos < hash_end && maxlen >= LZ_MIN_MATCH) {
            len = find_match(st, buf, (int64_t)pos, maxlen, &d);
            insert(st, buf, (int64_t)pos);

            /* Lazy: if the next position yields a strictly longer match,
             * emit a literal now and take that one instead. */
            while (st->lazy && len && len < LZ_MAX_MATCH && pos + 1 < hash_end && pos + 1 < end) {
                size_t maxlen2 = end - (pos + 1);
                if (maxlen2 > LZ_MAX_MATCH) maxlen2 = LZ_MAX_MATCH;
                size_t d2 = 0;
                size_t len2 = find_match(st, buf, (int64_t)pos + 1, maxlen2, &d2);
                if (len2 <= len) break;
                sym[ns++] = buf[pos];
                pos++;
                insert(st, buf, (int64_t)pos);
                len = len2; d = d2;
            }
        }

        if (len) {
            sym[ns++] = (uint16_t)(SYM_BYTE + (len - LZ_MIN_MATCH));
            dist[nm++] = (uint16_t)(d - 1);
            /* Index the positions covered by the match (pos already done). */
            size_t stop = pos + len;
            for (size_t p = pos + 1; p < stop && p < hash_end; p++) insert(st, buf, (int64_t)p);
            pos = stop;
        } else {
            sym[ns++] = buf[pos];
            pos++;
        }
    }
    *nmatch = nm;
    return ns;
}
