/*
 * codes.h - value <-> (code, extra bits) mapping.
 *
 * Values 0..15 are their own code. Larger values use two codes per power
 * of two: code = 16 + 2*(nb-4) + (second-highest bit), and the remaining
 * nb-1 low bits are sent as raw extra bits (nb = index of top bit).
 * Values up to 2^27-1 map to codes 0..61.
 *
 * Sequence fields:
 *   LL  literal run length          -> ll code
 *   ML  match length - LZR_ML_BASE  -> ml code
 *   OF  0,1,2 = repeat offset #0..2; otherwise 3 + code(distance - 1)
 */
#ifndef LZRANS_CODES_H
#define LZRANS_CODES_H
#include "common.h"

#define VAL_CODES 62
#define OF_CODES  (3 + VAL_CODES)

static inline unsigned val_code(uint32_t v, unsigned *nbits, uint32_t *extra)
{
    if (v < 16) { *nbits = 0; *extra = 0; return v; }
    unsigned nb = highbit32(v);
    *nbits = nb - 1;
    *extra = v & ((1u << (nb - 1)) - 1);
    return 16 + (nb - 4) * 2 + ((v >> (nb - 1)) & 1);
}

typedef struct { uint32_t base; uint8_t bits; } code_info;

/* Fill table for codes 0..VAL_CODES-1. */
static inline void code_table_init(code_info *t)
{
    for (unsigned c = 0; c < VAL_CODES; c++) {
        if (c < 16) { t[c].base = c; t[c].bits = 0; continue; }
        unsigned k = c - 16, nb = 4 + k / 2, hb = k & 1;
        t[c].base = (1u << nb) | (hb << (nb - 1));
        t[c].bits = (uint8_t)(nb - 1);
    }
}

/* ---- LSB-first raw bit writer / reader ---- */
#endif
