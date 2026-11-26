/*
 * xxh64.h - streaming XXH64 (the xxHash 64-bit algorithm, seed 0), used
 * as the frame checksum: the end block stores its low 32 bits. Several
 * times faster than CRC-32, so it stays out of the way of fast decoding.
 */
#ifndef LZRANS_XXH64_H
#define LZRANS_XXH64_H
#include "common.h"

#define XXH_P1 0x9E3779B185EBCA87ull
#define XXH_P2 0xC2B2AE3D27D4EB4Full
#define XXH_P3 0x165667B19E3779F9ull
#define XXH_P4 0x85EBCA77C2B2AE63ull
#define XXH_P5 0x27D4EB2F165667C5ull

typedef struct {
    uint64_t v[4], total;
    uint8_t buf[32];
    unsigned nbuf;
} xxh64_state;

static inline uint64_t xxh_rotl(uint64_t x, int r) { return (x << r) | (x >> (64 - r)); }
static inline uint64_t xxh_rd64(const uint8_t *p) { return rd64le(p); }

static inline uint64_t xxh_round(uint64_t acc, uint64_t in)
{
    acc += in * XXH_P2;
    acc = xxh_rotl(acc, 31);
    return acc * XXH_P1;
}

static inline uint64_t xxh_merge(uint64_t acc, uint64_t v)
{
    acc ^= xxh_round(0, v);
    return acc * XXH_P1 + XXH_P4;
}

static inline void xxh64_init(xxh64_state *s)
{
    s->v[0] = XXH_P1 + XXH_P2;
    s->v[1] = XXH_P2;
    s->v[2] = 0;
    s->v[3] = 0 - XXH_P1;
    s->total = 0;
    s->nbuf = 0;
}

static inline void xxh64_update(xxh64_state *s, const uint8_t *p, size_t n)
{
    s->total += n;
    if (s->nbuf) {
        size_t take = 32 - s->nbuf < n ? 32 - s->nbuf : n;
        memcpy(s->buf + s->nbuf, p, take);
        s->nbuf += (unsigned)take;
        p += take;
        n -= take;
        if (s->nbuf < 32) return;
        for (int i = 0; i < 4; i++) s->v[i] = xxh_round(s->v[i], xxh_rd64(s->buf + 8 * i));
        s->nbuf = 0;
    }
    uint64_t v0 = s->v[0], v1 = s->v[1], v2 = s->v[2], v3 = s->v[3];
    while (n >= 32) {
        v0 = xxh_round(v0, xxh_rd64(p));
        v1 = xxh_round(v1, xxh_rd64(p + 8));
        v2 = xxh_round(v2, xxh_rd64(p + 16));
        v3 = xxh_round(v3, xxh_rd64(p + 24));
        p += 32;
        n -= 32;
    }
    s->v[0] = v0; s->v[1] = v1; s->v[2] = v2; s->v[3] = v3;
    memcpy(s->buf, p, n);
    s->nbuf = (unsigned)n;
}

static inline uint64_t xxh64_digest(const xxh64_state *s)
{
    uint64_t h;
    if (s->total >= 32) {
        h = xxh_rotl(s->v[0], 1) + xxh_rotl(s->v[1], 7) + xxh_rotl(s->v[2], 12) + xxh_rotl(s->v[3], 18);
        for (int i = 0; i < 4; i++) h = xxh_merge(h, s->v[i]);
    } else {
        h = s->v[2] + XXH_P5;
    }
    h += s->total;
    const uint8_t *p = s->buf;
    size_t n = s->nbuf;
    while (n >= 8) {
        h ^= xxh_round(0, xxh_rd64(p));
        h = xxh_rotl(h, 27) * XXH_P1 + XXH_P4;
        p += 8;
        n -= 8;
    }
    if (n >= 4) {
        h ^= (uint64_t)rd32le(p) * XXH_P1;
        h = xxh_rotl(h, 23) * XXH_P2 + XXH_P3;
        p += 4;
        n -= 4;
    }
    while (n--) {
        h ^= (*p++) * XXH_P5;
        h = xxh_rotl(h, 11) * XXH_P1;
    }
    h ^= h >> 33;
    h *= XXH_P2;
    h ^= h >> 29;
    h *= XXH_P3;
    h ^= h >> 32;
    return h;
}

#endif
