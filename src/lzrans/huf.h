/*
 * huf.h - literal coding of block type 5: canonical Huffman codes of at
 * most HUF_MAXLEN bits, 4 independent streams decoded interleaved by table
 * lookup.
 *
 * Literal section (n literals, n known to the decoder):
 *   mode u8   0 raw bytes | 1 run (1 byte) | 2 new table | 3 previous table
 *   mode 2    table: last symbol u8, then one nibble (length 0..11) per
 *             symbol 0..last, low nibble first
 *   mode 2/3  sizes of streams 0-3 (u32 each; frame versions 5-6: 0-2 only,
 *             stream 3 ending where its last literal ends), then the 4 streams; stream
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
static inline size_t huf_bound(size_t n) { return 1 + 129 + 16 + n + n / 4 + 64; }

/* Encode lits[0..n) at out. prev: the table the decoder holds (mode 3),
 * updated to the table used. Returns bytes written (0 if n == 0). */
size_t huf_encode(uint8_t *out, const uint8_t *lits, size_t n, huf_table *prev);

/* Decode exactly n literals. Returns bytes consumed, 0 on malformed input.
 * Reads up to 8 bytes past the section (the caller keeps slack). */
size_t huf_decode(const uint8_t *in, size_t avail, uint8_t *out, size_t n, huf_table *t, int v7);

/* Code lengths (<= HUF_MAXLEN, 0 = absent) that huf_encode would build
 * for these counts; at least one count must be nonzero. */
void huf_lengths(const uint32_t *cnt, uint8_t *len);

/* ---- context literals (block type 6) ----
 * Each literal is coded with one of K <= HUFC_MAXK Huffman tables, picked
 * by the byte before it in the output (0 at the start of the content)
 * through a 256-entry map. One forward LSB-first stream, decoded while
 * the sequences are executed.
 *
 * Section:  mode u8   0 new tables | 1 previous tables
 *   mode 0  K u8, map: 256 nibbles (cluster of each previous byte, low
 *           nibble first), then per cluster: last symbol u8 + one nibble
 *           (length) per symbol 0..last as in mode 2 above
 *   then    stream size u32, stream
 */
#define HUFC_MAXK 8
#define HUFC_WORK (3 * 65536)

typedef struct {
    uint8_t k;                       /* clusters                              */
    uint8_t cmap[256];               /* previous byte -> cluster              */
    uint8_t len[HUFC_MAXK][256];
    uint16_t code[HUFC_MAXK][256];
    int valid;
    uint32_t *work;                  /* encoder scratch, HUFC_WORK entries    */
} hufc_table;

typedef struct {
    /* entry: symbol | length << 8 | cluster of the symbol (as the next
     * previous byte) << 12, indexed by cluster << HUF_MAXLEN | next bits */
    uint16_t dt[HUFC_MAXK << HUF_MAXLEN];
    uint8_t cmap[256];
    int valid;
} hufc_dtable;

/* Bound of a context literal section of n literals. */
static inline size_t hufc_bound(size_t n) { return 2 + 128 + HUFC_MAXK * 129 + 4 + n + n / 4 + 64; }

/* Encode lits[0..n) (prev[i]: the byte before lits[i] in the content) at
 * out, n > 0. st: the tables the decoder holds, updated. Returns bytes
 * written, or 0 when a single table would do as well (use type 5). */
size_t hufc_encode(uint8_t *out, const uint8_t *lits, const uint8_t *prev, size_t n, hufc_table *st);
void hufc_free(hufc_table *st);

/* Read the section header (mode and tables) into t. Returns bytes read,
 * 0 on malformed input. */
size_t hufc_read(const uint8_t *in, size_t avail, hufc_dtable *t, int strict);

#endif
