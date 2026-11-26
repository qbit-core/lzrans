#ifndef LZRANS_RANS_H
#define LZRANS_RANS_H
#include "common.h"

/*
 * Byte-wise rANS with a 32-bit state and 12-bit frequency scale (M = 4096).
 * Each of the three token streams is coded independently with its own
 * static frequency table (stored in the block header).
 */

/* Decoder slot entry: x' = freq * (x >> 12) + bias, where bias = slot - cum. */
typedef struct {
    uint16_t sym;
    uint16_t freq;
    uint16_t bias;
    uint16_t pad;
} rans_slot;

/* Count + normalise symbol frequencies so that they sum to RANS_M.
 * Every symbol that occurs gets freq >= 1. */
void rans_build_freqs(const uint16_t *syms, size_t n, unsigned alphabet, uint32_t *freq);

/* Serialise a frequency table (zero-run + varint). Returns bytes written. */
size_t rans_write_freqs(uint8_t *out, const uint32_t *freq, unsigned alphabet);

/* Parse a frequency table and build the O(1) decode slot table.
 * Returns bytes consumed, or 0 on malformed input. */
size_t rans_read_freqs(const uint8_t *in, size_t avail, unsigned alphabet, rans_slot *slots);

/* Upper bound on encoded size for n symbols. */
static inline size_t rans_bound(size_t n) { return 2 * n + 8; }

/* Encode syms[0..n) and write the payload to `out`. Returns size. */
size_t rans_encode(const uint16_t *syms, size_t n, const uint32_t *freq, unsigned alphabet,
                   uint8_t *scratch, size_t scratch_cap, uint8_t *out);

/* Max serialised frequency table size. */
#define RANS_FREQ_MAX(alphabet) ((size_t)(alphabet) * 3 + 8)

#endif
