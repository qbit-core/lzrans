/*
 * lzr.h - public C API of the lzrans core (frame format v3, streaming).
 *
 * Frame layout (little-endian):
 *   "LZRN" | version u8 (=3) | flags u8 | window_log u8 | reserved u8
 *   [content_size u64]                       if flags & LZR_FLAG_SIZE
 *   block*                                   until an end block
 *
 *   block header: type u8 | raw_len u32 | comp_len u32 | payload[comp_len]
 *     type 0  stored:      payload = raw_len bytes
 *     type 1  compressed:  payload = nseq u32 | nlit u32 |
 *                          literals | ll codes | ml codes | of codes
 *                          (rANS streams, see rans.h) | extra_len u32 | extra bits
 *     type 2  end:         raw_len = 0, payload = CRC32 of all content (u32)
 *     type 3  compressed, adaptive literals: as type 1, but the literals are
 *                          nplain u32 | 5 x (len u32 | range-coded stream):
 *                          matched literals, then 4 lanes of plain ones (arc.h)
 *   Match lengths are >= LZR_ML_BASE (2); the ml code holds length - 2.
 *
 * Matches may reach back up to 2^window_log bytes, across block borders.
 * The CRC sits at the end so a stream can be produced without knowing its
 * size in advance.
 */
#ifndef LZRANS_LZR_H
#define LZRANS_LZR_H
#include "common.h"

#define LZR_MAGIC          "LZRN"
#define LZR_VERSION        3
#define LZR_FLAG_ENCRYPTED 0x01u   /* reserved: future AEAD envelope       */
#define LZR_FLAG_SIZE      0x02u   /* header carries the content size       */
#define LZR_FLAG_X86       0x04u   /* content passed the x86 filter (x86.h) */
#define LZR_KNOWN_FLAGS    (LZR_FLAG_SIZE | LZR_FLAG_X86)
#define LZR_FILTER_NONE    0
#define LZR_FILTER_X86     1
#define LZR_BLOCK_HDR      9
#define LZR_LEVEL_MIN      1
#define LZR_LEVEL_MAX      9

/* ---- streaming compressor ---- */
typedef struct lzr_cctx lzr_cctx;

/* window_log 0 = level default. content_size < 0 = unknown. */
lzr_cctx *lzr_cctx_new(int level, int window_log, int64_t content_size);
void lzr_cctx_free(lzr_cctx *c);
/* Worker threads for the level 7-9 match finder (1 = none). The output is
 * identical for any thread count. Call before the first write. */
void lzr_cctx_set_threads(lzr_cctx *c, int threads);
/* Pre-filter: LZR_FILTER_NONE or LZR_FILTER_X86 (better ratio on x86
 * machine code). Call before the first write. */
int lzr_cctx_set_filter(lzr_cctx *c, int filter);
/* Feed input; complete blocks are appended to out. */
int lzr_cctx_write(lzr_cctx *c, const uint8_t *src, size_t n, lzr_buf *out);
/* Compress buffered input and write the end block. The context is then done. */
int lzr_cctx_end(lzr_cctx *c, lzr_buf *out);

/* ---- streaming decompressor ---- */
typedef struct lzr_dctx lzr_dctx;

lzr_dctx *lzr_dctx_new(void);
void lzr_dctx_free(lzr_dctx *d);
/* Feed compressed bytes; decoded bytes are appended to out. Returns LZR_OK
 * (more input welcome), 1 when the end block was verified, or LZR_ERR_*.
 * After the end, *consumed (if given) tells how many bytes of src belonged
 * to the frame. */
int lzr_dctx_write(lzr_dctx *d, const uint8_t *src, size_t n, lzr_buf *out, size_t *consumed);

/* Parse the frame header only. Returns header length, or LZR_ERR_*. */
int lzr_frame_header(const uint8_t *src, size_t n, int *window_log, uint8_t *flags, int64_t *content_size);

#endif
