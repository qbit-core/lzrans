/*
 * huf.h - literal coding of block type 5: canonical Huffman codes of at
 * most HUF_MAXLEN bits, 4 independent streams decoded interleaved by table
 * lookup.
 *
 * Literal section (n literals, n known to the decoder):
 *   mode u8   0 raw bytes | 1 run (1 byte) | 2 new table | 3 previous table
 *   mode 2    table: last symbol u8, then one nibble (length 0..11) per
 *             symbol 0..last, low nibble first
 *   mode 2/3  sizes of streams 1-3 (u32 each), then the 4 streams; stream
 *             k holds literals [k*q, min((k+1)*q, n)), q = ceil(n / 4)
 * Streams are written LSB-first; codes are stored bit-reversed so the
 * decoder can index its table with the low bits of its bit buffer.
 */
#ifndef LZRANS_HUF_H
#define LZRANS_HUF_H
#include "common.h"

#define HUF_MAXLEN 11
#define HUF_TSIZE (1u << HUF_MAXLEN)

typedef struct {
    uint8_t len[256];          /* code lengths, 0 = absent                 */
    uint16_t code[256];        /* bit-reversed canonical codes (encoder)    */
    uint16_t dt[HUF_TSIZE];    /* decoder table: sym | len << 8             */
    int valid;                 /* a table exists (for mode 3)               */
} huf_table;

/* Upper bound of an encoded literal section. */
static inline size_t huf_bound(size_t n) { return 1 + 129 + 12 + n + n / 4 + 64; }

/* Encode lits[0..n) at out. prev: the table the decoder holds (mode 3),
 * updated to the table used. Returns bytes written (0 if n == 0). */
size_t huf_encode(uint8_t *out, const uint8_t *lits, size_t n, huf_table *prev);

/* Decode exactly n literals. Returns bytes consumed, 0 on malformed input.
 * Reads up to 8 bytes past the section (the caller keeps slack). */
size_t huf_decode(const uint8_t *in, size_t avail, uint8_t *out, size_t n, huf_table *t);

/* Price (in 1/256 bits) of each byte value under the code that huf_encode
 * would build for these counts. */
void huf_lengths(const uint32_t *cnt, uint8_t *len);

#endif
