/*
 * common.h - shared constants, little-endian helpers, growable buffers,
 * and error codes for the lzrans C core.
 */
#ifndef LZRANS_COMMON_H
#define LZRANS_COMMON_H

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef _MSC_VER
#include <intrin.h>
static inline unsigned highbit32(uint32_t v) { unsigned long r; _BitScanReverse(&r, v); return (unsigned)r; }
#else
/* Index of the highest set bit; v must be non-zero. */
static inline unsigned highbit32(uint32_t v) { return 31u - (unsigned)__builtin_clz(v); }
#endif

/* ---- frame / block parameters ---- */
#define LZR_BLOCK_LOG   18
#define LZR_BLOCK_SIZE  (1u << LZR_BLOCK_LOG)  /* 256 KB of input per block      */
#define LZR_WLOG_MIN    18                     /* window >= block size           */
#define LZR_WLOG_MAX    26                     /* 64 MB                          */
#define LZR_MIN_MATCH   4                      /* new-offset matches */
#define LZR_ML_BASE     2                      /* smallest codable match (repeat offsets) */
#define LZR_NUM_REPS    3
#define LZR_SLACK       32                     /* over-allocation for wildcopy   */

/* ---- rANS parameters ---- */
#define RANS_SCALE_BITS 12
#define RANS_M          (1u << RANS_SCALE_BITS)   /* 4096 */
#define RANS_L          (1u << 23)                /* lower state bound */

/* Alphabets: literal bytes, and the three sequence code streams. */
#define ALPHA_LIT   256
#define ALPHA_CODE  72

/* ---- error codes ---- */
enum {
    LZR_OK = 0,
    LZR_ERR_NOMEM = -1,
    LZR_ERR_CORRUPT = -2,     /* malformed / truncated stream              */
    LZR_ERR_CHECKSUM = -3,    /* CRC32 mismatch after decompression        */
    LZR_ERR_UNSUPPORTED = -4, /* unknown version or flags (e.g. encrypted) */
    LZR_ERR_PARAM = -5,       /* bad argument / API misuse                 */
};

static inline uint32_t rd32(const uint8_t *p) { uint32_t v; memcpy(&v, p, 4); return v; }
static inline void wr32le(uint8_t *p, uint32_t v) {
    p[0] = (uint8_t)v; p[1] = (uint8_t)(v >> 8); p[2] = (uint8_t)(v >> 16); p[3] = (uint8_t)(v >> 24);
}
static inline uint32_t rd32le(const uint8_t *p) {
    return (uint32_t)p[0] | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
}
static inline void wr64le(uint8_t *p, uint64_t v) { wr32le(p, (uint32_t)v); wr32le(p + 4, (uint32_t)(v >> 32)); }
static inline uint64_t rd64le(const uint8_t *p) { return (uint64_t)rd32le(p) | ((uint64_t)rd32le(p + 4) << 32); }

/* Large tables (match-finder trees, windows) are accessed randomly; on
 * Linux, back them with transparent huge pages to cut TLB misses. */
#if defined(__linux__)
#include <sys/mman.h>
#endif
static inline void *big_alloc(size_t n, int zero)
{
#if defined(__linux__) && defined(MADV_HUGEPAGE)
    if (n >= (1u << 21)) {
        void *p = NULL;
        if (posix_memalign(&p, 1u << 21, n) != 0) return NULL;
        madvise(p, n, MADV_HUGEPAGE);
        if (zero) memset(p, 0, n);
        return p;
    }
#endif
    return zero ? calloc(1, n) : malloc(n);
}

/* Growable output buffer shared by encoder and decoder. */
typedef struct {
    uint8_t *p;
    size_t len, cap;
} lzr_buf;

static inline int buf_reserve(lzr_buf *b, size_t extra)
{
    if (b->cap - b->len >= extra) return LZR_OK;
    size_t need = b->len + extra, cap = b->cap ? b->cap : 4096;
    while (cap < need) cap *= 2;
    uint8_t *np = (uint8_t *)realloc(b->p, cap);
    if (!np) return LZR_ERR_NOMEM;
    b->p = np;
    b->cap = cap;
    return LZR_OK;
}

static inline int buf_append(lzr_buf *b, const void *src, size_t n)
{
    if (buf_reserve(b, n)) return LZR_ERR_NOMEM;
    memcpy(b->p + b->len, src, n);
    b->len += n;
    return LZR_OK;
}

#endif
