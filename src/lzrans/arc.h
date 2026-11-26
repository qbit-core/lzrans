/*
 * arc.h - adaptive binary range coder for the literals of compressed
 * block type 3.
 *
 * Coder: LZMA-style range coder, 32-bit range, 12-bit probabilities of a
 * zero bit, updated with shift ARC_SHIFT after every decision. The encoder
 * flushes 5 bytes; the decoder reads 5 bytes up front and must end with
 * its input pointer exactly at the end of the stream.
 *
 * Model: an 8-bit binary tree per context, context = previous byte >> 5.
 * The first literal after a match is coded with the byte at distance rep0
 * (the "match byte") as extra context while their bits agree, as in LZMA:
 * it is known to differ from that byte. The model persists across the
 * blocks of a frame; stored and type 1 blocks leave it untouched.
 *
 * Why only literals: adaptive coding of the ll/ml/of code streams gained
 * another 1-2.5% but costs ~12 serial binary decisions per sequence, which
 * made decoding no faster than xz on text. Those stay static rANS.
 */
#ifndef LZRANS_ARC_H
#define LZRANS_ARC_H
#include "common.h"

#define ARC_PBITS 12
#define ARC_ONE (1u << ARC_PBITS)
#ifndef ARC_SHIFT
#define ARC_SHIFT 5
#endif
#define ARC_TOP (1u << 24)

typedef uint16_t arc_p;

typedef struct {
    arc_p lit[8][768];   /* [1..255] tree, [256..767] matched mode */
} arc_model;

static inline void arc_model_init(arc_model *m)
{
    arc_p *p = &m->lit[0][0];
    for (size_t i = 0; i < sizeof(m->lit) / sizeof(arc_p); i++) p[i] = ARC_ONE / 2;
}

/* ---------------- encoder ---------------- */
typedef struct {
    uint64_t low;
    uint32_t range;
    uint8_t cache;
    uint64_t cache_size;
    uint8_t *p;
} arc_enc;

static inline void arce_init(arc_enc *e, uint8_t *out)
{
    e->low = 0; e->range = 0xFFFFFFFFu; e->cache = 0; e->cache_size = 1; e->p = out;
}

static inline void arce_shift_low(arc_enc *e)
{
    if ((uint32_t)e->low < 0xFF000000u || (e->low >> 32)) {
        uint8_t carry = (uint8_t)(e->low >> 32);
        uint8_t t = e->cache;
        do { *e->p++ = (uint8_t)(t + carry); t = 0xFF; } while (--e->cache_size);
        e->cache = (uint8_t)(e->low >> 24);
    }
    e->cache_size++;
    e->low = (e->low & 0x00FFFFFFu) << 8;
}

static inline void arce_bit(arc_enc *e, arc_p *pr, unsigned bit)
{
    uint32_t bound = (e->range >> ARC_PBITS) * *pr;
    if (!bit) {
        e->range = bound;
        *pr = (arc_p)(*pr + ((ARC_ONE - *pr) >> ARC_SHIFT));
    } else {
        e->low += bound;
        e->range -= bound;
        *pr = (arc_p)(*pr - (*pr >> ARC_SHIFT));
    }
    while (e->range < ARC_TOP) { e->range <<= 8; arce_shift_low(e); }
}

/* Returns the number of bytes written since arce_init. */
static inline size_t arce_finish(arc_enc *e, uint8_t *start)
{
    for (int i = 0; i < 5; i++) arce_shift_low(e);
    return (size_t)(e->p - start);
}

/* MSB-first bit tree of nbits over probs[1 .. 2^nbits - 1]. */
static inline void arce_tree(arc_enc *e, arc_p *probs, unsigned nbits, unsigned sym)
{
    unsigned m = 1;
    for (int i = (int)nbits - 1; i >= 0; i--) {
        unsigned b = (sym >> i) & 1;
        arce_bit(e, &probs[m], b);
        m = (m << 1) | b;
    }
}

static inline void arce_literal(arc_enc *e, arc_p *probs, unsigned byte)
{
    arce_tree(e, probs, 8, byte);
}

/* First literal after a match: the match byte's bits select the
 * probabilities until the first differing bit. */
static inline void arce_literal_matched(arc_enc *e, arc_p *probs, unsigned byte, unsigned mbyte)
{
    unsigned m = 1, same = 1;
    for (int i = 7; i >= 0; i--) {
        unsigned b = (byte >> i) & 1, mb = (mbyte >> i) & 1;
        arce_bit(e, same ? &probs[256 + (mb << 8) + m] : &probs[m], b);
        m = (m << 1) | b;
        same &= (b == mb);
    }
}

/* ---------------- decoder ---------------- */
typedef struct {
    uint32_t range, code;
    const uint8_t *p, *end;   /* a valid stream ends with p == end */
    int err;
} arc_dec;

static inline uint8_t arcd_byte(arc_dec *d)
{
    if (d->p < d->end) return *d->p++;
    d->err = 1;
    return 0;
}

static inline void arcd_init(arc_dec *d, const uint8_t *in, size_t n)
{
    d->p = in; d->end = in + n; d->err = 0;
    d->range = 0xFFFFFFFFu; d->code = 0;
    if (arcd_byte(d) != 0) d->err = 1;
    for (int i = 0; i < 4; i++) d->code = (d->code << 8) | arcd_byte(d);
}

static inline unsigned arcd_bit(arc_dec *d, arc_p *pr)
{
    /* Branchy on purpose: a branchless variant measured ~15% slower (the
     * bit is predicted often enough, and it needs twice the instructions). */
    uint32_t p = *pr;
    uint32_t bound = (d->range >> ARC_PBITS) * p;
    unsigned bit;
    if (d->code < bound) {
        d->range = bound;
        *pr = (arc_p)(p + ((ARC_ONE - p) >> ARC_SHIFT));
        bit = 0;
    } else {
        d->code -= bound;
        d->range -= bound;
        *pr = (arc_p)(p - (p >> ARC_SHIFT));
        bit = 1;
    }
    if (d->range < ARC_TOP) { d->range <<= 8; d->code = (d->code << 8) | arcd_byte(d); }
    return bit;
}

static inline unsigned arcd_tree(arc_dec *d, arc_p *probs, unsigned nbits)
{
    unsigned m = 1;
    for (unsigned i = 0; i < nbits; i++) m = (m << 1) | arcd_bit(d, &probs[m]);
    return m - (1u << nbits);
}

static inline unsigned arcd_literal_matched(arc_dec *d, arc_p *probs, unsigned mbyte)
{
    unsigned m = 1;
    for (int i = 7; i >= 0; i--) {
        unsigned mb = (mbyte >> i) & 1;
        unsigned b = arcd_bit(d, &probs[256 + (mb << 8) + m]);
        m = (m << 1) | b;
        if (b != mb) {
            while (m < 256) m = (m << 1) | arcd_bit(d, &probs[m]);
            break;
        }
    }
    return m & 0xFF;
}

#endif
