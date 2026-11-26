/*
 * arc.h - adaptive binary range coder for the literals of compressed
 * block type 3.
 *
 * Coder: LZMA-style range coder, 32-bit range, 12-bit probabilities of a
 * zero bit, updated with shift ARC_SHIFT after every decision. The encoder
 * flushes 5 bytes; the decoder reads 5 bytes up front and must end with
 * its input pointer exactly at the end of the stream.
 *
 * The literals of a block are split by role:
 *   matched  the first literal after a match. Coded with the byte at
 *            distance rep0 (the "match byte") as extra context while their
 *            bits agree, as in LZMA (it is known to differ from that byte),
 *            plus context = previous byte >> 5. One stream, decoded inside
 *            the sequence loop because both contexts need the output.
 *   plain    all other literals, in order. Dealt round-robin to ARC_LANES
 *            streams with their own order-0 bit-tree models, so the decoder
 *            decodes them up front with the lanes' dependency chains
 *            interleaved. A single adaptive coder is latency bound at
 *            ~15-20 cycles per bit, slower than xz on literal-heavy data.
 * Splitting costs ~0.5% of literal size against one coder with full
 * contexts. The models persist across the blocks of a frame; stored and
 * type 1 blocks leave them untouched.
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
#define ARC_LANES 4
#define ARC_PLAIN_SLACK 16   /* readable bytes required after each plain lane */

typedef uint16_t arc_p;

typedef struct {
    arc_p matched[8][768];         /* [1..255] tree, [256..767] while matching */
    arc_p plain[ARC_LANES][256];
} arc_model;

static inline void arc_model_init(arc_model *m)
{
    for (int c = 0; c < 8; c++)
        for (int i = 0; i < 768; i++) m->matched[c][i] = ARC_ONE / 2;
    for (int l = 0; l < ARC_LANES; l++)
        for (int i = 0; i < 256; i++) m->plain[l][i] = ARC_ONE / 2;
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

/* Probabilities stay in [31, ARC_ONE - 31], so a decision costs at most
 * ~7.05 bits and one byte of renormalisation always suffices. */
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

static inline void arce_literal(arc_enc *e, arc_p *probs, unsigned byte)
{
    unsigned m = 1;
    for (int i = 7; i >= 0; i--) {
        unsigned b = (byte >> i) & 1;
        arce_bit(e, &probs[m], b);
        m = (m << 1) | b;
    }
}

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

/* Branchy on purpose: branchless variants measured slower here (the bit is
 * predicted often enough, and they need about twice the instructions). */
static inline unsigned arcd_bit(arc_dec *d, arc_p *pr)
{
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

/* Decode n plain literals; literal j uses lane j % ARC_LANES. The lanes
 * are independent, so groups of ARC_LANES are decoded bit-interleaved and
 * the CPU overlaps their chains. Coder state is kept in locals: through a
 * pointer, every output byte store (char may alias anything) would force
 * it back to memory. Returns nonzero if a lane ran past its end. */
static inline int arcd_plain(arc_dec *ln, arc_model *m, uint8_t *out, size_t n)
{
    uint32_t r0 = ln[0].range, r1 = ln[1].range, r2 = ln[2].range, r3 = ln[3].range;
    uint32_t c0 = ln[0].code, c1 = ln[1].code, c2 = ln[2].code, c3 = ln[3].code;
    const uint8_t *q0 = ln[0].p, *q1 = ln[1].p, *q2 = ln[2].p, *q3 = ln[3].p;
    const uint8_t *e0 = ln[0].end, *e1 = ln[1].end, *e2 = ln[2].end, *e3 = ln[3].end;
    arc_p *p0 = m->plain[0], *p1 = m->plain[1], *p2 = m->plain[2], *p3 = m->plain[3];
    int err = 0;
/* Branchless: plain literals are close to random, so a branch on the bit
 * mispredicts about half the time, and mispredictions also defeat the
 * interleaving. Renormalisation reads *Q unconditionally: the caller
 * guarantees ARC_PLAIN_SLACK readable bytes after every lane, and the
 * overrun check after each group of four literals stops a corrupt stream
 * before it gets further. */
#define ARCD_BIT(R, C, Q, E, PR, SYM)                                              \
    do {                                                                           \
        arc_p *pr_ = (PR);                                                         \
        uint32_t pv_ = *pr_, bound_ = ((R) >> ARC_PBITS) * pv_;                    \
        uint32_t b_ = (C) >= bound_, m_ = 0u - b_;                                 \
        (C) -= bound_ & m_;                                                        \
        (R) = bound_ + (((R) - 2 * bound_) & m_);                                  \
        uint32_t p0_ = pv_ + ((ARC_ONE - pv_) >> ARC_SHIFT), p1_ = pv_ - (pv_ >> ARC_SHIFT); \
        *pr_ = (arc_p)(p0_ ^ ((p0_ ^ p1_) & m_));                                  \
        (SYM) = ((SYM) << 1) | b_;                                                 \
        uint32_t n_ = (R) < ARC_TOP, sh_ = n_ << 3;                                \
        (R) <<= sh_;                                                               \
        (C) = ((C) << sh_) | (*(Q) & (0u - n_));                                   \
        (Q) += n_;                                                                 \
    } while (0)
    size_t j = 0;
    for (; j + 4 <= n; j += 4) {
        unsigned s0 = 1, s1 = 1, s2 = 1, s3 = 1;
        for (int b = 0; b < 8; b++) {
            ARCD_BIT(r0, c0, q0, e0, &p0[s0], s0);
            ARCD_BIT(r1, c1, q1, e1, &p1[s1], s1);
            ARCD_BIT(r2, c2, q2, e2, &p2[s2], s2);
            ARCD_BIT(r3, c3, q3, e3, &p3[s3], s3);
        }
        out[j] = (uint8_t)s0; out[j + 1] = (uint8_t)s1; out[j + 2] = (uint8_t)s2; out[j + 3] = (uint8_t)s3;
        if ((q0 > e0) | (q1 > e1) | (q2 > e2) | (q3 > e3)) { err = 1; break; }
    }
    /* tail: lanes 0, 1, 2 in order */
    if (!err && j < n) { unsigned s = 1; for (int b = 0; b < 8; b++) ARCD_BIT(r0, c0, q0, e0, &p0[s], s); out[j++] = (uint8_t)s; }
    if (j < n) { unsigned s = 1; for (int b = 0; b < 8; b++) ARCD_BIT(r1, c1, q1, e1, &p1[s], s); out[j++] = (uint8_t)s; }
    if (j < n) { unsigned s = 1; for (int b = 0; b < 8; b++) ARCD_BIT(r2, c2, q2, e2, &p2[s], s); out[j++] = (uint8_t)s; }
#undef ARCD_BIT
    ln[0].range = r0; ln[1].range = r1; ln[2].range = r2; ln[3].range = r3;
    ln[0].code = c0; ln[1].code = c1; ln[2].code = c2; ln[3].code = c3;
    ln[0].p = q0; ln[1].p = q1; ln[2].p = q2; ln[3].p = q3;
    return err | (q0 > e0) | (q1 > e1) | (q2 > e2) | (q3 > e3);
}

#endif
