/*
 * common.h - shared constants, little-endian helpers and error codes.
 *
 * All multi-byte integers in the container format are little-endian.
 * The helpers use memcpy so they are alignment-safe and compile to a
 * single load/store on x86-64 / AArch64.
 */
#ifndef LZRANS_COMMON_H
#define LZRANS_COMMON_H

#include <stddef.h>
#include <stdint.h>
#include <string.h>

/* ---- LZ77 parameters (from the project spec) ---- */
#define LZ_HASH_BITS   16
#define LZ_HASH_SIZE   (1u << LZ_HASH_BITS)
#define LZ_WINDOW      65536u          /* 64 KB sliding window            */
#define LZ_MIN_MATCH   4               /* 4-byte hash => min match 4      */
#define LZ_MAX_MATCH   (LZ_MIN_MATCH + 255) /* 259: fits symbol 256..511  */

/* ---- rANS parameters ---- */
#define RANS_SCALE_BITS 12
#define RANS_M          (1u << RANS_SCALE_BITS)   /* 4096 */
#define RANS_L          (1u << 23)                /* lower state bound */

/* Stream 1 alphabet: 0..255 literal bytes, 256..511 match lengths. */
#define SYM_LITLEN      512
#define SYM_BYTE        256

/* ---- error codes returned by the C core ---- */
enum {
    LZR_OK = 0,
    LZR_ERR_NOMEM = -1,
    LZR_ERR_CORRUPT = -2,     /* malformed / truncated stream            */
    LZR_ERR_CHECKSUM = -3,    /* CRC32 mismatch after decompression      */
    LZR_ERR_UNSUPPORTED = -4, /* unknown version or flags (e.g. encrypted) */
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

#endif
