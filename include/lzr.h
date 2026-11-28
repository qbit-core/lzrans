/*
 * lzr.h - public API of liblzr, the lzrans compression library.
 *
 * Writes LZR 1.0 frames (frame format version 7, docs/FORMAT.md) and reads
 * versions 5-7. The output of every function here is byte-identical to the
 * Python package (lzrans.compress / lzrans.Compressor) for the same
 * parameters.
 *
 *   One-shot:   lzr_compress / lzr_compress_cs, lzr_decompress / lzr_decompress_ds
 *   Streaming:  lzr_cstream_* + lzr_compress_stream,
 *               lzr_dstream_* + lzr_decompress_stream (bounded memory)
 *
 * Every function returns LZR_OK (0) or a negative LZR_ERR_* code unless
 * noted; lzr_error_string() turns a code into a message. Contexts are not
 * thread-safe, but separate contexts may be used from different threads at
 * the same time. Link with -llzr -lpthread -lm (static) or -llzr (shared).
 */
#ifndef LZR_PUBLIC_H
#define LZR_PUBLIC_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN32) && defined(LZR_DLL_EXPORT)
#define LZR_API __declspec(dllexport)
#elif defined(__GNUC__) && __GNUC__ >= 4
#define LZR_API __attribute__((visibility("default")))
#else
#define LZR_API
#endif

/* ---- versions ---- */
#define LZR_LIB_VERSION_MAJOR 0
#define LZR_LIB_VERSION_MINOR 15
#define LZR_LIB_VERSION_PATCH 0
#define LZR_LIB_VERSION_NUMBER (LZR_LIB_VERSION_MAJOR * 10000 + LZR_LIB_VERSION_MINOR * 100 + LZR_LIB_VERSION_PATCH)
#define LZR_FORMAT_VERSION 7          /* frames written: LZR 1.0 */

LZR_API unsigned lzr_version_number(void);     /* LZR_LIB_VERSION_NUMBER of the linked library */
LZR_API const char *lzr_version_string(void);  /* e.g. "0.15.0" */

/* ---- error codes ---- */
enum {
    LZR_OK = 0,
    LZR_ERR_NOMEM = -1,
    LZR_ERR_CORRUPT = -2,     /* malformed stream (or trailing data after a frame) */
    LZR_ERR_CHECKSUM = -3,    /* content checksum mismatch after decompression     */
    LZR_ERR_UNSUPPORTED = -4, /* unknown frame version, flags or block type        */
    LZR_ERR_PARAM = -5,       /* bad argument or API misuse                        */
    LZR_ERR_MEMLIMIT = -6,    /* frame window above the decoder's max window log   */
    LZR_ERR_DST_SIZE = -7,    /* destination buffer too small (one-shot functions) */
    LZR_ERR_TRUNCATED = -8,   /* input ended before the end of the frame           */
};
LZR_API const char *lzr_error_string(int code);

/* ---- limits and defaults ---- */
#define LZR_LEVEL_MIN     -3          /* fastest */
#define LZR_LEVEL_MAX     11          /* best ratio (xz -9 class) */
#define LZR_LEVEL_DEFAULT 6
#define LZR_WINDOW_LOG_MIN 18
#define LZR_WINDOW_LOG_ENC_MAX 26     /* largest window the encoder uses (64 MiB) */
#define LZR_WINDOW_LOG_MAX 31         /* largest window a frame may declare */
#define LZR_MAX_WINDOW_LOG_DEFAULT 27 /* decoder memory limit by default (128 MiB) */
#define LZR_SEGMENT_MIN ((uint64_t)1 << 20)
#define LZR_THREADS_MAX 16
#define LZR_FILTER_NONE 0
#define LZR_FILTER_X86  1             /* x86/x86-64 CALL filter (executables) */

/* ---- compression ---- */
typedef struct lzr_cstream lzr_cstream;

typedef enum {
    LZR_C_LEVEL = 1,        /* LZR_LEVEL_MIN..LZR_LEVEL_MAX, default 6 */
    LZR_C_WINDOW_LOG,       /* 0 = level default, else 18..26 */
    LZR_C_THREADS,          /* worker threads, default 1; 0 = one per online CPU (max 16).
                               Levels 7-11 and segment mode: output identical for any count;
                               levels -3..6: identical for any count > 1 (~0.1% larger than 1). */
    LZR_C_SEGMENT,          /* reset interval in bytes for parallel decompression:
                               -1 = level default (48 MiB at levels 9-11, else off),
                               0 = off, else >= LZR_SEGMENT_MIN */
    LZR_C_FILTER,           /* LZR_FILTER_NONE or LZR_FILTER_X86 */
    LZR_C_CONTENT_SIZE,     /* pledged content size, -1 = unknown (default). When known it is
                               stored in the frame header and the input must match it exactly. */
} lzr_cparam;

/* A new context with default parameters. NULL when out of memory. */
LZR_API lzr_cstream *lzr_cstream_create(void);
LZR_API void lzr_cstream_free(lzr_cstream *cs);
/* Set a parameter. Allowed before the first lzr_compress_stream call of a
 * frame (or after lzr_cstream_reset). LZR_C_CONTENT_SIZE is cleared to -1
 * by reset; the others persist. */
LZR_API int lzr_cstream_set(lzr_cstream *cs, lzr_cparam param, int64_t value);
/* Abandon the current frame (if any) and prepare for a new one. */
LZR_API int lzr_cstream_reset(lzr_cstream *cs);

typedef struct { const void *src; size_t size; size_t pos; } lzr_inbuf;
typedef struct { void *dst; size_t size; size_t pos; } lzr_outbuf;

typedef enum { LZR_CONTINUE = 0, LZR_FINISH = 1 } lzr_end_op;

/* Consume input from in->src[in->pos..size) and write compressed bytes to
 * out->dst[out->pos..size), advancing both positions. With LZR_FINISH the
 * frame is completed once all input is consumed.
 * Returns:
 *   1   output is pending: call again with more room in out
 *   0   all input consumed and all output written; with LZR_FINISH the
 *       frame is complete (call lzr_cstream_reset to start another)
 *   <0  LZR_ERR_*; the context must then be reset
 * The encoder works in blocks, so output appears in bursts. Passing the
 * whole input in one call with LZR_FINISH and a pledged content size lets
 * levels -3..6 with threads read it in place (faster); the bytes produced
 * are the same either way. */
LZR_API int lzr_compress_stream(lzr_cstream *cs, lzr_outbuf *out, lzr_inbuf *in, lzr_end_op op);

/* Largest compressed size of n input bytes (any level and parameters). */
LZR_API size_t lzr_compress_bound(size_t n);
/* One-shot: compress src[0..n) into dst (dst_cap >= lzr_compress_bound(n)
 * always suffices) at `level`, one thread, level-default segments.
 * *dst_len receives the frame length. LZR_ERR_DST_SIZE if it does not fit. */
LZR_API int lzr_compress(void *dst, size_t dst_cap, size_t *dst_len, const void *src, size_t n, int level);
/* One-shot with the parameters of cs (the content size is set to n). Resets cs. */
LZR_API int lzr_compress_cs(lzr_cstream *cs, void *dst, size_t dst_cap, size_t *dst_len, const void *src, size_t n);

/* ---- decompression ---- */
typedef struct lzr_dstream lzr_dstream;

typedef enum {
    LZR_D_MAX_WINDOW_LOG = 1,   /* 18..31, default 27: frames with a larger window fail
                                   with LZR_ERR_MEMLIMIT (decoder memory ~ 2^value) */
    LZR_D_THREADS,              /* one-shot lzr_decompress_ds only: default 1, 0 = one per
                                   online CPU (max 16). Segments and blocks of levels -3..10
                                   decode in parallel; the output is the same. */
} lzr_dparam;

LZR_API lzr_dstream *lzr_dstream_create(void);
LZR_API void lzr_dstream_free(lzr_dstream *ds);
LZR_API int lzr_dstream_set(lzr_dstream *ds, lzr_dparam param, int64_t value);
/* Prepare to decode a new frame (e.g. the next of concatenated frames). */
LZR_API int lzr_dstream_reset(lzr_dstream *ds);

#define LZR_FRAME_END 1
/* Consume compressed input and write decoded bytes, advancing both
 * positions. Returns:
 *   LZR_FRAME_END (1)  the frame ended, its checksum matched and all of its
 *                      output was written; in->pos is just past the frame.
 *                      Bytes after it are left unread (call
 *                      lzr_dstream_reset to decode a following frame).
 *   0                  progress: call again with more input, or more room
 *                      in out when out->pos == out->size
 *   <0                 LZR_ERR_*; the context must then be reset
 * Memory stays bounded by the frame's window (about 2x 2^window_log). */
LZR_API int lzr_decompress_stream(lzr_dstream *ds, lzr_outbuf *out, lzr_inbuf *in);

/* One-shot: decode exactly one frame that fills src[0..n) into dst with
 * default parameters. LZR_ERR_TRUNCATED if the frame is cut short,
 * LZR_ERR_CORRUPT if bytes follow it, LZR_ERR_DST_SIZE if dst is too small. */
LZR_API int lzr_decompress(void *dst, size_t dst_cap, size_t *dst_len, const void *src, size_t n);
/* One-shot with the parameters of ds (threads, max window log): decodes the
 * frame at the start of src; trailing bytes are allowed and *src_used
 * (if not NULL) receives the frame length. Fastest when the frame header
 * carries the content size and dst_cap >= size + LZR_DECODE_MARGIN (the
 * decoder then writes in place). */
#define LZR_DECODE_MARGIN 64
LZR_API int lzr_decompress_ds(lzr_dstream *ds, void *dst, size_t dst_cap, size_t *dst_len,
                              const void *src, size_t n, size_t *src_used);

/* ---- frame inspection ---- */
typedef struct {
    int64_t content_size;   /* -1 when the header does not carry it */
    int window_log;
    int version;            /* frame version (7 = LZR 1.0) */
    unsigned flags;         /* LZR_FRAME_FLAG_* */
    size_t header_size;
} lzr_frame_info;
#define LZR_FRAME_FLAG_SIZE 0x02u
#define LZR_FRAME_FLAG_X86  0x04u
#define LZR_FRAME_MAGIC_SIZE 4   /* "LZRN" */
/* Parse a frame header at src. Returns LZR_OK, LZR_ERR_TRUNCATED when more
 * bytes are needed, or another LZR_ERR_* (LZR_ERR_CORRUPT: not a frame). */
LZR_API int lzr_get_frame_info(const void *src, size_t n, lzr_frame_info *info);

#ifdef __cplusplus
}
#endif
#endif
