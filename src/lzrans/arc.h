/*
 * arc.h - adaptive binary range coder and the token model of compressed
 * block type 4 (levels 6-9).
 *
 * Coder: LZMA-style range coder, 32-bit range, 12-bit probabilities of a
 * zero bit, updated with shift ARC_SHIFT after every decision. The encoder
 * flushes 5 bytes; the decoder reads 5 bytes up front and must end with
 * its input pointer exactly at the end of the block payload.
 *
 * Token model (the coding layer of LZMA, on our sequences):
 *   every position is a literal or a match, chosen by is_match[state][pos&3].
 *   state (12 values) remembers the last few token kinds: literals are
 *   coded against the "match byte" at rep0 right after a match, and the
 *   rep / new-distance choice depends on what came before.
 *   literal       8-bit tree, context = previous byte >> 5
 *   match         is_rep[state]: 0 = new distance, 1 = repeat offset
 *     rep         is_rep0[state], then is_rep1[state]: which of 3 offsets;
 *                 for rep0, is_rep0_long[state][pos&3] = 0 means a 1-byte
 *                 match ("short rep") with no length
 *     new dist    slot (6-bit tree, context = min(len - 2, 3)); slots >= 4
 *                 carry nb = slot/2 - 1 footer bits: reverse trees for
 *                 slots < 14, else direct bits + 4 adaptive "align" bits
 *   length        2..273: choice bits select 3-bit low / 3-bit mid trees
 *                 (context pos&3) or an 8-bit high tree; separate coders
 *                 for new-distance and repeat matches
 * Longer matches are emitted as a match followed by rep0 matches.
 * The model, the state and the repeat offsets persist across the blocks of
 * a frame; stored and type 1 blocks leave the model untouched.
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

#define TOK_STATES 12
#define TOK_PB 2                     /* position bits for contexts */
#define TOK_POS (1 << TOK_PB)
#define TOK_MAXLEN 273               /* longest single match token */
#define TOK_DIST_SLOTS 64
#define TOK_END_SPEC 14              /* slots below this use reverse trees */

typedef uint16_t arc_p;

typedef struct {
    arc_p choice, choice2;
    arc_p low[TOK_POS][8], mid[TOK_POS][8], high[256];
} tok_len;

typedef struct {
    arc_p is_match[TOK_STATES][TOK_POS];
    arc_p is_rep[TOK_STATES], is_rep0[TOK_STATES], is_rep1[TOK_STATES];
    arc_p is_rep0_long[TOK_STATES][TOK_POS];   /* 0 = "short rep": 1 byte at rep0 */
    arc_p lit[8][0x300];             /* [1..255] tree, [256..767] matched mode */
    arc_p slot[4][TOK_DIST_SLOTS];
    arc_p spec[TOK_END_SPEC][32];    /* reverse footer trees, per slot */
    arc_p align[16];
    tok_len mlen, rlen;
    uint8_t state;
} tok_model;

static inline void tok_model_init(tok_model *m)
{
    arc_p *p = (arc_p *)m;
    for (size_t i = 0; i < offsetof(tok_model, state) / sizeof(arc_p); i++) p[i] = ARC_ONE / 2;
    m->state = 0;
}

static inline unsigned tok_after_lit(unsigned s) { return s < 4 ? 0 : (s < 10 ? s - 3 : s - 6); }
static inline unsigned tok_after_match(unsigned s) { return s < 7 ? 7 : 10; }
static inline unsigned tok_after_rep(unsigned s) { return s < 7 ? 8 : 11; }
static inline unsigned tok_after_short(unsigned s) { return s < 7 ? 9 : 11; }

/* dist0 = distance - 1 -> slot (2 per power of two, 0..3 direct). */
static inline unsigned tok_slot(uint32_t dist0)
{
    if (dist0 < 4) return dist0;
    unsigned nb = highbit32(dist0);
    return 2 * nb + ((dist0 >> (nb - 1)) & 1);
}

/* ---------------- encoder ---------------- */
typedef struct {
    uint64_t low;
    uint32_t range;
    uint8_t cache;
    uint64_t cache_size;
    uint8_t *p;
    int dry;          /* only update the model (parser price tracking) */
} arc_enc;

static inline void arce_init(arc_enc *e, uint8_t *out)
{
    e->low = 0; e->range = 0xFFFFFFFFu; e->cache = 0; e->cache_size = 1; e->p = out; e->dry = 0;
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

/* Probabilities stay in [31, ARC_ONE - 31]: a decision costs at most
 * ~7.05 bits, and one byte of renormalisation always suffices. */
static inline void arce_bit(arc_enc *e, arc_p *pr, unsigned bit)
{
    if (e->dry) {
        if (!bit) *pr = (arc_p)(*pr + ((ARC_ONE - *pr) >> ARC_SHIFT));
        else *pr = (arc_p)(*pr - (*pr >> ARC_SHIFT));
        return;
    }
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

static inline void arce_direct(arc_enc *e, uint32_t v, unsigned nb)
{
    if (e->dry) return;
    while (nb--) {
        e->range >>= 1;
        if ((v >> nb) & 1) e->low += e->range;
        while (e->range < ARC_TOP) { e->range <<= 8; arce_shift_low(e); }
    }
}

static inline size_t arce_finish(arc_enc *e, uint8_t *start)
{
    for (int i = 0; i < 5; i++) arce_shift_low(e);
    return (size_t)(e->p - start);
}

static inline void arce_tree(arc_enc *e, arc_p *probs, unsigned nbits, unsigned sym)
{
    unsigned m = 1;
    for (int i = (int)nbits - 1; i >= 0; i--) {
        unsigned b = (sym >> i) & 1;
        arce_bit(e, &probs[m], b);
        m = (m << 1) | b;
    }
}

static inline void arce_rtree(arc_enc *e, arc_p *probs, unsigned nbits, unsigned sym)
{
    unsigned m = 1;
    for (unsigned i = 0; i < nbits; i++) {
        unsigned b = (sym >> i) & 1;
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

static inline void arce_len(arc_enc *e, tok_len *l, unsigned len, unsigned ps)
{
    len -= 2;
    if (len < 8) { arce_bit(e, &l->choice, 0); arce_tree(e, l->low[ps], 3, len); return; }
    arce_bit(e, &l->choice, 1);
    if (len < 16) { arce_bit(e, &l->choice2, 0); arce_tree(e, l->mid[ps], 3, len - 8); return; }
    arce_bit(e, &l->choice2, 1);
    arce_tree(e, l->high, 8, len - 16);
}

static inline void arce_dist(arc_enc *e, tok_model *m, uint32_t dist0, unsigned len)
{
    unsigned slot = tok_slot(dist0), lctx = len - 2 < 3 ? len - 2 : 3;
    arce_tree(e, m->slot[lctx], 6, slot);
    if (slot < 4) return;
    unsigned nb = (slot >> 1) - 1;
    uint32_t base = (2u | (slot & 1)) << nb, rest = dist0 - base;
    if (slot < TOK_END_SPEC) { arce_rtree(e, m->spec[slot], nb, rest); return; }
    arce_direct(e, rest >> 4, nb - 4);
    arce_rtree(e, m->align, 4, rest & 15);
}

/* One literal at cp (prev = the byte before it, 0 at the stream start). */
static inline void tok_enc_lit(arc_enc *e, tok_model *m, const uint8_t *cp, unsigned prev, uint32_t rep0, unsigned ps)
{
    unsigned s = m->state;
    arce_bit(e, &m->is_match[s][ps], 0);
    arc_p *probs = m->lit[prev >> 5];
    if (s >= 7) arce_literal_matched(e, probs, cp[0], cp[-(ptrdiff_t)rep0]);
    else arce_tree(e, probs, 8, cp[0]);
    m->state = (uint8_t)tok_after_lit(s);
}

/* One match token, len 2..TOK_MAXLEN (or 1 with off 0); off: 0..2 = repeat index, else
 * distance + 2 (the convention of commit_offset). */
static inline void tok_enc_match(arc_enc *e, tok_model *m, uint32_t off, unsigned len, unsigned ps)
{
    unsigned s = m->state;
    arce_bit(e, &m->is_match[s][ps], 1);
    if (off >= 3) {
        arce_bit(e, &m->is_rep[s], 0);
        arce_len(e, &m->mlen, len, ps);
        arce_dist(e, m, off - 3, len);
        m->state = (uint8_t)tok_after_match(s);
    } else {
        arce_bit(e, &m->is_rep[s], 1);
        arce_bit(e, &m->is_rep0[s], off != 0);
        if (off) {
            arce_bit(e, &m->is_rep1[s], off - 1);
        } else {
            arce_bit(e, &m->is_rep0_long[s][ps], len > 1);
            if (len == 1) { m->state = (uint8_t)tok_after_short(s); return; }
        }
        arce_len(e, &m->rlen, len, ps);
        m->state = (uint8_t)tok_after_rep(s);
    }
}

/* Length of the first token of a match of length ml (the rest follows as
 * rep0 tokens); every token is >= 2. */
static inline unsigned tok_piece(size_t ml)
{
    if (ml <= TOK_MAXLEN) return (unsigned)ml;
    return ml - TOK_MAXLEN < 2 ? TOK_MAXLEN - 2 : TOK_MAXLEN;
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

static inline uint32_t arcd_direct(arc_dec *d, unsigned nb)
{
    uint32_t v = 0;
    while (nb--) {
        d->range >>= 1;
        uint32_t t = d->code >= d->range;
        d->code -= d->range & (0u - t);
        v = (v << 1) | t;
        if (d->range < ARC_TOP) { d->range <<= 8; d->code = (d->code << 8) | arcd_byte(d); }
    }
    return v;
}

static inline unsigned arcd_tree(arc_dec *d, arc_p *probs, unsigned nbits)
{
    unsigned m = 1;
    for (unsigned i = 0; i < nbits; i++) m = (m << 1) | arcd_bit(d, &probs[m]);
    return m - (1u << nbits);
}

static inline unsigned arcd_rtree(arc_dec *d, arc_p *probs, unsigned nbits)
{
    unsigned m = 1, v = 0;
    for (unsigned i = 0; i < nbits; i++) {
        unsigned b = arcd_bit(d, &probs[m]);
        m = (m << 1) | b;
        v |= b << i;
    }
    return v;
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

static inline unsigned arcd_len(arc_dec *d, tok_len *l, unsigned ps)
{
    if (!arcd_bit(d, &l->choice)) return 2 + arcd_tree(d, l->low[ps], 3);
    if (!arcd_bit(d, &l->choice2)) return 10 + arcd_tree(d, l->mid[ps], 3);
    return 18 + arcd_tree(d, l->high, 8);
}

static inline uint32_t arcd_dist(arc_dec *d, tok_model *m, unsigned len)
{
    unsigned lctx = len - 2 < 3 ? len - 2 : 3;
    unsigned slot = arcd_tree(d, m->slot[lctx], 6);
    if (slot < 4) return slot;
    unsigned nb = (slot >> 1) - 1;
    uint32_t base = (2u | (slot & 1)) << nb;
    if (slot < TOK_END_SPEC) return base + arcd_rtree(d, m->spec[slot], nb);
    uint32_t hi = arcd_direct(d, nb - 4);
    return base + (hi << 4) + arcd_rtree(d, m->align, 4);
}

#endif
