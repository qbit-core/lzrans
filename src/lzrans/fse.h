/*
 * fse.h - table-driven ANS (tANS / "FSE") for the sequence codes of block
 * type 5, and the backward bit stream they share with the extra bits.
 *
 * A table of L = 2^FSE_LOG states. Decoding: e = T[state]; symbol = e.sym;
 * state = e.base + (next e.nb bits). Encoding runs backwards over the
 * symbols; a symbol of normalised count f maps the state X in [L, 2L) to
 * enc[s][(X >> nb) - f] after writing the low nb bits of X, where nb is the
 * one shift that puts X >> nb in [f, 2f).
 *
 * Bit stream: the encoder writes LSB-first, forwards, and ends with a
 * single 1 bit; the decoder starts from that marker and reads backwards,
 * so it sees the last-written bits first.
 */
#ifndef LZRANS_FSE_H
#define LZRANS_FSE_H
#include "common.h"

#define FSE_LOG 9
#define FSE_L (1u << FSE_LOG)
#define FSE_MAXSYM 72

typedef struct { uint8_t sym, nb; uint16_t base; } fse_dent;

typedef struct {
    uint16_t norm[FSE_MAXSYM];   /* normalised counts (sum FSE_L)       */
    unsigned alphabet;
    int valid;
    fse_dent dt[FSE_L];          /* decoder                             */
    uint16_t st[FSE_L];          /* encoder: next states, by symbol run */
    uint16_t first[FSE_MAXSYM];  /* encoder: start of each symbol's run */
    uint8_t hbf[FSE_MAXSYM];     /* encoder: floor(log2(norm))          */
    uint32_t dnb[FSE_MAXSYM];    /* encoder: nb = (X + dnb) >> 16       */
    int32_t dfs[FSE_MAXSYM];     /* encoder: first - norm               */
} fse_table;

/* Normalise counts (at least one present) to sum FSE_L. */
void fse_normalize(const uint32_t *cnt, unsigned alphabet, uint32_t total, uint16_t *norm);
/* Build both directions from t->norm / t->alphabet. Returns 0 or -1. */
int fse_build(fse_table *t);
size_t fse_write_norm(uint8_t *out, const uint16_t *norm, unsigned alphabet);
/* Returns bytes read, 0 on malformed input. */
size_t fse_read_norm(const uint8_t *in, size_t avail, uint16_t *norm, unsigned alphabet);
/* Cost in 1/256 bits of coding cnt[] with norm[] (UINT64_MAX if a used
 * symbol has no state). */
uint64_t fse_cost(const uint32_t *cnt, const uint16_t *norm, unsigned alphabet);

/* ---- forward bit writer (LSB-first) ---- */
typedef struct { uint64_t acc; unsigned n; uint8_t *p; } fbw;

static inline void fbw_put(fbw *w, uint64_t v, unsigned len)
{
    w->acc |= v << w->n;
    w->n += len;
    if (w->n >= 32) { wr32le(w->p, (uint32_t)w->acc); w->p += 4; w->acc >>= 32; w->n -= 32; }
}

/* Writes the end marker and flushes; returns the end pointer. */
static inline uint8_t *fbw_close(fbw *w)
{
    fbw_put(w, 1, 1);
    while (w->n > 0) { *w->p++ = (uint8_t)w->acc; w->acc >>= 8; w->n = w->n > 8 ? w->n - 8 : 0; }
    return w->p;
}

/* Encode symbol s: state *x in [L, 2L). */
static inline void fse_put(fbw *w, const fse_table *t, uint32_t *x, unsigned s)
{
    uint32_t X = *x, f = t->norm[s];
    unsigned nb = (unsigned)(31 - __builtin_clz(X)) - t->hbf[s];
    if ((X >> nb) < f) nb--;   /* X >> nb is in [2^hbf, 2^(hbf+1)): fix to [f, 2f) */
    fbw_put(w, X & ((1u << nb) - 1), nb);
    *x = t->st[t->first[s] + (X >> nb) - f];
}

/* fse_put without the write: returns the bits (count in *nb). */
static inline uint32_t fse_step(const fse_table *t, uint32_t *x, unsigned s, unsigned *nbp)
{
    /* nb as in fse_put, from one add and shift (dnb: see fse_build) */
    uint32_t X = *x;
    unsigned nb = (X + t->dnb[s]) >> 16;
    *nbp = nb;
    *x = t->st[(int32_t)(X >> nb) + t->dfs[s]];
    return X & ((1u << nb) - 1);
}

/* ---- backward bit reader ---- */
typedef struct {
    uint64_t c;
    unsigned used;          /* bits of c already consumed, from the top */
    const uint8_t *p, *start;
    int err;
} bbr;

static inline uint64_t fse_ld64(const uint8_t *p)
{
    uint64_t v;
    memcpy(&v, p, 8);
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    v = __builtin_bswap64(v);
#endif
    return v;
}

/* The 8 bytes before `start` must be readable (their bits are never
 * used by a valid stream). Returns -1 if the stream has no end marker. */
static inline int bbr_init(bbr *r, const uint8_t *start, size_t size)
{
    r->start = start; r->err = 0;
    if (!size || !start[size - 1]) return -1;
    r->p = start + size - 8;
    r->c = fse_ld64(r->p);
    r->used = (unsigned)__builtin_clzll(r->c) + 1;
    return 0;
}

static inline uint32_t bbr_get(bbr *r, unsigned nb)
{
    /* (used & 63): only a corrupt stream reads past 64 bits; no UB then */
    uint32_t v = (uint32_t)(((r->c << (r->used & 63)) >> 1) >> (63 - nb));
    r->used += nb;
    return v;
}

static inline void bbr_reload(bbr *r)
{
    size_t nbytes = r->used >> 3;
    size_t room = r->p > r->start ? (size_t)(r->p - r->start) : 0;
    if (nbytes > room) nbytes = room;
    r->p -= nbytes;
    r->used -= (unsigned)nbytes * 8;
    if (r->used > 64) { r->err = 1; r->used = 64; }   /* overrun: corrupt */
    r->c = fse_ld64(r->p);
}

/* All bits consumed, exactly? */
static inline int bbr_done(const bbr *r, size_t size)
{
    return !r->err && (size_t)(r->start + size - (r->p + 8)) * 8 + r->used == size * 8;
}


#endif
