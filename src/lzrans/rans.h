/*
 * rans.h - static-model byte-wise rANS, 32-bit state, 12-bit scale,
 * two interleaved states per stream (independent dependency chains give
 * the CPU instruction-level parallelism when decoding).
 *
 * A "stream" is n symbols (n known to the decoder) coded in one of:
 *   mode 0  raw bytes
 *   mode 1  run: every symbol is the same (1 byte)
 *   mode 2  rANS: freq table | u32 payload length | payload
 */
#ifndef LZRANS_RANS_H
#define LZRANS_RANS_H
#include "common.h"

/* Decode slot: sym | freq << 8 | bias << 20  (freq, bias < 4096). */
typedef uint32_t rans_dslot;

/* Upper bound on an encoded stream (header + payload) of n symbols. */
static inline size_t stream_bound(size_t n, unsigned alphabet)
{
    return 1 + (size_t)alphabet * 3 + 8 + 2 * n + 16;
}

/* Encode syms[0..n) as a stream at `out`. `scratch` must hold
 * 2*n + 16 bytes. Returns bytes written (0 if n == 0). */
size_t stream_encode(uint8_t *out, const uint8_t *syms, size_t n, unsigned alphabet, uint8_t *scratch);

/* Decode a stream of exactly n symbols into out[]. `tbl` is a work area of
 * RANS_M slots. Returns bytes consumed, or 0 on malformed input (n > 0). */
size_t stream_decode(const uint8_t *in, size_t avail, uint8_t *out, size_t n, unsigned alphabet,
                     rans_dslot *tbl);

#endif
