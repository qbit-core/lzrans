/*
 * lzr.h - public C API of the lzrans core (frame format v7 = LZR 1.0,
 * streaming; full specification: docs/FORMAT.md; v5 and v6 frames still
 * decode). Summary of v7 below; see FORMAT.md for v5/v6 differences.
 *
 * Frame layout (little-endian):
 *   "LZRN" | version u8 (=7) | flags u8 | window_log u8 | reserved u8
 *   [content_size u64]                       if flags & LZR_FLAG_SIZE
 *   block*                                   until an end block
 *
 *   block header: type u8 | raw_len u32 | comp_len u32 | payload[comp_len]
 *     type byte: bit 7 skippable (raw_len 0, payload ignored), else bits
 *     0-5 the type and bit 6 a reset point; types 7-63 unsupported
 *     type 0  stored:      payload = raw_len bytes
 *     type 2  end:         raw_len = 0, payload = low 32 bits of the XXH64
 *                          (seed 0) of all content (u32)
 *     type 4  token coding (level 11): one adaptive range-coded stream of
 *                          literal / match tokens (arc.h); its model, state
 *                          and the repeat offsets carry over between blocks
 *     type 5  fast (levels -3..10): nseq u32 | nlit u32 | literal section
 *                          (huf.h: 4 Huffman streams) | if nseq: ll/ml/of
 *                          table modes + one backward tANS bit stream of
 *                          states and extra bits (fse.h, dec.c). Tables and
 *                          repeat offsets carry over between blocks.
 *     type 6  context literals (levels 9-10, frame version >= 6): as type
 *                          5, but the literal section is a context literal
 *                          section (huf.h): literals are coded with one of
 *                          up to 8 Huffman tables chosen by the byte before
 *                          them, and decoded while the sequences run. Its
 *                          tables carry over between type 6 blocks.
 *   window_log is 18..31; decoders refuse windows above their memory limit.
 *   Match lengths are >= LZR_ML_BASE (2); the ml code holds length - 2.
 *   Types 1 and 3 (older formats) are rejected.
 *
 *   type byte 0x80 = seek index (skippable): per segment u64 frame offset,
 *                          u64 content offset; u64 content size, u32 count,
 *                          "LZRX" (FORMAT.md 12.3); written by segment mode
 *
 * Matches may reach back up to 2^window_log bytes, across block borders.
 * The checksum sits at the end so a stream can be produced without knowing
 * its size in advance.
 */
#ifndef LZRANS_LZR_H
#define LZRANS_LZR_H
#include "common.h"

#define LZR_MAGIC          "LZRN"
#define LZR_VERSION        7   /* format 1.0 */
#define LZR_VERSION_MIN    5   /* oldest frame version still decoded */
#define LZR_FLAG_ENCRYPTED 0x01u   /* reserved: future AEAD envelope       */
#define LZR_FLAG_SIZE      0x02u   /* header carries the content size       */
#define LZR_FLAG_X86       0x04u   /* content passed the x86 filter (x86.h) */
#define LZR_KNOWN_FLAGS    (LZR_FLAG_SIZE | LZR_FLAG_X86)
#define LZR_FILTER_NONE    0
#define LZR_FILTER_X86     1
#define LZR_BLOCK_HDR      9
/* Block type byte, frame version >= 7: bit 7 = skippable block (raw_len 0,
 * payload ignored), else bits 0-5 = type, bit 6 = reset point (decoder
 * state back to frame-start values; nothing before the block is used). */
#define LZR_BT_SKIP        0x80u
#define LZR_BT_RESET       0x40u
#define LZR_BT_TYPE        0x3Fu
#define LZR_LEVEL_MIN     -3
#define LZR_LEVEL_MAX     11

/* ---- streaming compressor ---- */
typedef struct lzr_cctx lzr_cctx;

/* window_log 0 = level default. content_size < 0 = unknown. */
lzr_cctx *lzr_cctx_new(int level, int window_log, int64_t content_size);
void lzr_cctx_free(lzr_cctx *c);
/* Segment mode: like lzr_cctx_new, but the content is cut into segments of
 * seg_size bytes (the last one shorter), each compressed independently
 * and starting at a reset point, and a seek index (skippable block 0x80)
 * is written before the end block. Segments decode in parallel
 * (lzr_decompress_mt). With threads, up to `threads` segments are
 * compressed at once; the output does not depend on the thread count.
 * seg_size 0, or a known content_size <= seg_size, gives exactly the
 * output of lzr_cctx_new. Output is produced a batch of segments at a time. */
lzr_cctx *lzr_cctx_new_seg(int level, int window_log, int64_t content_size, uint64_t seg_size);
#define LZR_SEG_MIN ((uint64_t)1 << 20)   /* smallest segment size accepted */
/* The segment size a level uses by default (0 = no segments). */
uint64_t lzr_default_segment(int level);
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

/* One-shot decode of a whole frame src[0..n) into out, with up to
 * `threads` threads decoding segments (reset points) in parallel. Frames
 * without reset points, legacy versions and x86-filtered frames are
 * decoded by one thread. Same result and error kinds as lzr_dctx_write in
 * one-buffer mode; *consumed = frame length. Returns 1 at the verified end,
 * 0 if the frame is truncated, or LZR_ERR_* (LZR_ERR_MEMLIMIT if the frame's
 * window is above 2^max_window_log). */
int lzr_decompress_mt(const uint8_t *src, size_t n, lzr_buf *out, int threads, int max_window_log,
                      size_t *consumed);

/* Decoder memory limit: frames declaring a window above 2^max_window_log
 * bytes fail with LZR_ERR_MEMLIMIT (default LZR_WLOG_DEC_DEFAULT = 27;
 * 18..31). Call before the first write. */
int lzr_dctx_set_max_window_log(lzr_dctx *d, int max_window_log);
/* The check itself, for decoders that parse the header on their own. */
static inline int lzr_window_ok(int window_log, int max_window_log)
{
    return window_log <= max_window_log ? LZR_OK : LZR_ERR_MEMLIMIT;
}

/* Parse the frame header only. Returns header length, or LZR_ERR_*. */
int lzr_frame_header(const uint8_t *src, size_t n, int *window_log, uint8_t *flags, int64_t *content_size);

#endif
