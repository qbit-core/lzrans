#ifndef LZRANS_LZ77_H
#define LZRANS_LZ77_H
#include "common.h"

/*
 * Hash-chain LZ77 match finder state.
 *  head[h]            most recent absolute position whose 4-byte hash is h
 *  prev[pos & mask]   previous position with the same hash (chain link)
 * Positions are absolute offsets into the whole input, so the 64 KB window
 * naturally spans block boundaries.
 */
typedef struct {
    int64_t head[LZ_HASH_SIZE];
    int64_t prev[LZ_WINDOW];
    int max_chain;   /* chain depth limit (speed/ratio knob)              */
    int lazy;        /* 1 = one-step lazy matching                        */
} lz_state;

void lz_init(lz_state *st, int level);

/*
 * Parse buf[start, end) into tokens. Matches never cross `end`, but may
 * reference data before `start` (up to LZ_WINDOW back). `n` is the total
 * buffer length (hash reads near `end` may look past it, never past n).
 *
 * Outputs:
 *   sym[i]   0..255 literal, 256..511 match of length (sym-256+LZ_MIN_MATCH)
 *   dist[j]  (distance-1) for the j-th match, 0..65535
 * Returns number of symbols; *nmatch receives number of matches.
 * Caller sizes sym[] >= end-start and dist[] >= (end-start)/LZ_MIN_MATCH + 1.
 */
size_t lz_parse(lz_state *st, const uint8_t *buf, size_t n, size_t start, size_t end,
                uint16_t *sym, uint16_t *dist, size_t *nmatch);

#endif
