/*
 * lzr.h - public C API of the lzrans core (frame format v5, streaming).
 *
 * Frame layout (little-endian):
 *   "LZRN" | version u8 (=5) | flags u8 | window_log u8 | reserved u8
 *   [content_size u64]                       if flags & LZR_FLAG_SIZE
 *   block*                                   until an end block
 *
 *   block header: type u8 | raw_len u32 | comp_len u32 | payload[comp_len]
 *     type 0  stored:      payload = raw_len bytes
 *     type 2  end:         raw_len = 0, payload = low 32 bits of the XXH64
 *                          (seed 0) of all content (u32)
 *     type 4  token coding (level 11): one adaptive range-coded stream of
 *                          literal / match tokens (arc.h); its model, state
 *                          and the repeat offsets carry over between blocks
 *     type 5  fast (levels 1-8): nseq u32 | nlit u32 | literal section
 *                          (huf.h: 4 Huffman streams) | if nseq: ll/ml/of
 *                          table modes + one backward tANS bit stream of
 *                          states and extra bits (fse.h, dec.c). Tables and
 *                          repeat offsets carry over between blocks.
 *   Match lengths are >= LZR_ML_BASE (2); the ml code holds length - 2.
 *   Types 1 and 3 (older formats) are rejected.
 *
 * Matches may reach back up to 2^window_log bytes, across block borders.
 * The checksum sits at the end so a stream can be produced without knowing
 * its size in advance.
 */
#ifndef LZRANS_LZR_H
#define LZRANS_LZR_H
#include "common.h"

#define LZR_MAGIC          "LZRN"
#define LZR_VERSION        5
#define LZR_FLAG_ENCRYPTED 0x01u   /* reserved: future AEAD envelope       */
#define LZR_FLAG_SIZE      0x02u   /* header carries the content size       */
#define LZR_FLAG_X86       0x04u   /* content passed the x86 filter (x86.h) */
#define LZR_KNOWN_FLAGS    (LZR_FLAG_SIZE | LZR_FLAG_X86)
#define LZR_FILTER_NONE    0
#define LZR_FILTER_X86     1
#define LZR_BLOCK_HDR      9
#define LZR_LEVEL_MIN     -3
#define LZR_LEVEL_MAX     11

/* ---- streaming compressor ---- */
typedef struct lzr_cctx lzr_cctx;

/* window_log 0 = level default. content_size < 0 = unknown. */
lzr_cctx *lzr_cctx_new(int level, int window_log, int64_t content_size);
void lzr_cctx_free(lzr_cctx *c);
/* Worker threads (1 = none). Levels 9-11: match finder threads, output
 * identical for any thread count. Levels 1-6: parallel 4 MB jobs, output
 * identical for any thread count > 1 (and ~0.1% larger than with one).
 * Call before the first write. */
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
/* One-buffer mode: the caller passes this same `out` to every write and
 * leaves it alone until the end. Output is then decoded in place (the
 * output is the history): no copy, no history allocation. Call before the
 * first write. Ignored for x86-filtered frames. */
int lzr_dctx_set_output(lzr_dctx *d, lzr_buf *out);

/* Parse the frame header only. Returns header length, or LZR_ERR_*. */
int lzr_frame_header(const uint8_t *src, size_t n, int *window_log, uint8_t *flags, int64_t *content_size);

#endif
