/*
 * dec.c - streaming decoder.
 *
 * Input is buffered until a whole block is available (its header gives
 * comp_len). Block type 5: the Huffman literals are decoded into an array
 * first; then one loop decodes each sequence from the tANS bit stream
 * and copies its literal run and match with fixed-size chunks. Block type 4: one range
 * decoder loop. Output goes to a history buffer that keeps the last
 * 2^window_log bytes, or, in one-buffer mode, straight into the output.
 */
#include "lzr.h"
#include "codes.h"
#include "xxh64.h"
#include "arc.h"
#include "x86.h"
#include "huf.h"
#include "fse.h"

#define MAX_SEQ (LZR_BLOCK_SIZE / LZR_ML_BASE + 1)
#define MAX_COMP (2 * LZR_BLOCK_SIZE + 65536)   /* sanity cap on comp_len */

enum { ST_HEADER, ST_BLOCKS, ST_DONE };

typedef struct {
    uint32_t base;   /* value without its extra bits                 */
    uint16_t next;   /* tANS: next state = next + (nb bits)          */
    uint8_t nb, xb;  /* tANS bits, extra bits                        */
} seq_dent;

struct lzr_dctx {
    int state;
    lzr_buf in;
    size_t ipos;
    uint8_t *hist;     /* history: own buffer, or the output (direct mode) */
    lzr_buf *direct;   /* lzr_dctx_set_output(): decode into this buffer   */
    size_t hcap, wsize, op;
    uint32_t reps[LZR_NUM_REPS];
    xxh64_state xs;    /* checksum of the decoded content */
    uint64_t total_out;
    int64_t content_size;
    uint8_t *lits;     /* type 5 literals of the current block */
    code_info ct[VAL_CODES], mt[ML_CODES];
    tok_model tm;      /* block type 4 model, persists across blocks */
    huf_table huf;     /* block type 5 tables, kept for the repeat modes */
    fse_table fse[3];
    seq_dent sq[3][FSE_L];   /* per-block sequence decoding tables */
    lzr_x86 x86;
    uint8_t *x86tmp;   /* x86 filter work buffer, NULL = no filter */
};

lzr_dctx *lzr_dctx_new(void)
{
    lzr_dctx *d = (lzr_dctx *)calloc(1, sizeof(*d));
    if (d) xxh64_init(&d->xs);
    if (!d) return NULL;
    d->reps[0] = 1; d->reps[1] = 4; d->reps[2] = 8;
    d->content_size = -1;
    d->lits = (uint8_t *)malloc(LZR_BLOCK_SIZE + LZR_SLACK);
    if (!d->lits) { lzr_dctx_free(d); return NULL; }
    code_table_init(d->ct);
    ml_table_init(d->mt);
    tok_model_init(&d->tm);
    return d;
}

void lzr_dctx_free(lzr_dctx *d)
{
    if (!d) return;
    free(d->in.p);
    if (!d->direct) free(d->hist);
    free(d->x86tmp);
    free(d->lits);
    free(d);
}

int lzr_frame_header(const uint8_t *src, size_t n, int *window_log, uint8_t *flags, int64_t *content_size)
{
    if (memcmp(src, LZR_MAGIC, n < 4 ? n : 4) != 0) return LZR_ERR_CORRUPT;
    if (n < 8) return 0;                          /* need more input */
    if (src[4] != LZR_VERSION) return LZR_ERR_UNSUPPORTED;
    if (flags) *flags = src[5];
    if (src[5] & ~LZR_KNOWN_FLAGS) return LZR_ERR_UNSUPPORTED;
    if (src[6] < LZR_WLOG_MIN || src[6] > LZR_WLOG_MAX || src[7]) return LZR_ERR_CORRUPT;
    if (window_log) *window_log = src[6];
    if (src[5] & LZR_FLAG_SIZE) {
        if (n < 16) return 0;
        uint64_t sz = rd64le(src + 8);
        if (sz > (uint64_t)INT64_MAX) return LZR_ERR_CORRUPT;
        if (content_size) *content_size = (int64_t)sz;
        return 16;
    }
    if (content_size) *content_size = -1;
    return 8;
}

/* Copy len bytes from op - dist to op; regions may overlap (dist < len).
 * Fixed-size chunks only (no library calls); may write up to 15 bytes
 * past op + len (buffers have LZR_SLACK). */
static inline void copy_match(uint8_t *op, size_t dist, size_t len)
{
    const uint8_t *src = op - dist;
    if (dist >= 16) {
        for (size_t k = 0; k < len; k += 16) memcpy(op + k, src + k, 16);
    } else if (dist >= 8) {
        for (size_t k = 0; k < len; k += 8) memcpy(op + k, src + k, 8);
    } else {
        /* one period byte by byte, then 8-byte chunks from a whole number
         * of periods back (8 <= period * m <= 14) */
        for (size_t k = 0; k < 8; k++) op[k] = src[k];
        size_t back = dist * ((8 + dist - 1) / dist);
        for (size_t k = 8; k < len; k += 8) memcpy(op + k, op + k - back, 8);
    }
}

/* Literal run: 16-byte chunks (source and destination have slack). */
static inline void copy_lits(uint8_t *o, const uint8_t *lp, size_t ll)
{
    memcpy(o, lp, 16);
    for (size_t k = 16; k < ll; k += 16) memcpy(o + k, lp + k, 16);
}

/* Block type 4: one range-coded token stream (arc.h), decoded straight
 * into the history. The coder state lives in locals (range, code, ip) so
 * it stays in registers; the macros mirror arce_bit() and the trees in
 * arc.h. */
#define RC_NORM() do { if (range < ARC_TOP) { range <<= 8; code = (code << 8) | *ip++; } } while (0)
#define RC_BIT(pr, b) do {                                                     \
        arc_p *pp_ = (pr);                                                     \
        uint32_t p_ = *pp_, bound_ = (range >> ARC_PBITS) * p_;                \
        if (code < bound_) { range = bound_; *pp_ = (arc_p)(p_ + ((ARC_ONE - p_) >> ARC_SHIFT)); b = 0; } \
        else { range -= bound_; code -= bound_; *pp_ = (arc_p)(p_ - (p_ >> ARC_SHIFT)); b = 1; } \
        RC_NORM();                                                             \
    } while (0)
/* Branchless decisions for tree bits (literals, lengths, distances),
 * whose outcomes are hard to predict. Same arithmetic as RC_BIT.
 * RC_SEL sets m = all ones if the bit is 0 (code < bound), else 0, and
 * range/code to bound/code or range - bound/code - bound. RC_PUPD then
 * updates p: bit 0 p += (ONE - p) >> SHIFT, bit 1 p -= p >> SHIFT; both
 * are p - ((p - k) >> SHIFT) with an arithmetic shift, k = ARC_K0 for
 * bit 0 and 0 for bit 1 (checked for every p in 0..ONE). */
#define ARC_K0 (ARC_ONE - (1u << ARC_SHIFT) + 1)
#if defined(__x86_64__) && defined(__GNUC__) && !defined(LZR_NO_ASM)
/* The compiler does not reliably emit conditional moves here. */
#define RC_SEL(m, bound) do {                                                  \
        uint32_t rg_ = (bound), r1_ = range - (bound), c1_ = code - (bound);   \
        __asm__("cmpl %[bd], %[cd]\n\t"                                        \
                "cmovael %[xr], %[rg]\n\t"                                     \
                "cmovael %[xc], %[cd]\n\t"                                     \
                "sbbl %[mm], %[mm]"                                            \
                : [rg] "+r"(rg_), [cd] "+r"(code), [mm] "=r"(m)               \
                : [bd] "r"(bound), [xr] "r"(r1_), [xc] "r"(c1_)               \
                : "cc");                                                       \
        range = rg_;                                                           \
    } while (0)
#else
#define RC_SEL(m, bound) do {                                                  \
        m = 0u - (uint32_t)(code < (bound));                                   \
        range = (bound) + (~m & (range - 2 * (bound)));                        \
        code -= ~m & (bound);                                                  \
    } while (0)
#endif
#define RC_PUPD(p, m) ((arc_p)((p) - (uint32_t)((int32_t)((p) - ((m) & ARC_K0)) >> ARC_SHIFT)))
/* Bit trees. The next node is probs[2m] or probs[2m + 1]; both are loaded
 * while the current bit is decoded and the right one is picked after, so
 * the load is off the bit-to-bit dependency chain. rev: the bits are the
 * value LSB first (reverse tree); else MSB first. */
#define RC_TREE_(probs, nbits, out, rev) do {                                  \
        arc_p *pb_ = (probs);                                                  \
        unsigned m_ = 1, v_ = 0;                                               \
        uint32_t pc_ = pb_[1];                                                 \
        for (unsigned i_ = 0; i_ < (nbits); i_++) {                            \
            uint32_t bd_ = (range >> ARC_PBITS) * pc_, mm_, pl_ = 0, ph_ = 0;  \
            if (i_ + 1 < (nbits)) { pl_ = pb_[2 * m_]; ph_ = pb_[2 * m_ + 1]; } \
            RC_SEL(mm_, bd_);                                                  \
            pb_[m_] = RC_PUPD(pc_, mm_);                                       \
            unsigned b_ = mm_ + 1;                                             \
            if (rev) v_ |= b_ << i_;                                           \
            m_ = 2 * m_ + b_;                                                  \
            pc_ = ph_ ^ ((pl_ ^ ph_) & mm_);                                   \
            RC_NORM();                                                         \
        }                                                                      \
        out = (rev) ? v_ : m_ - (1u << (nbits));                               \
    } while (0)
#define RC_TREE(probs, nbits, out) RC_TREE_(probs, nbits, out, 0)
#define RC_RTREE(probs, nbits, out) RC_TREE_(probs, nbits, out, 1)
#define RC_LEN(l, ps, out) do {                                                \
        unsigned b0_, v0_;                                                     \
        RC_BIT(&(l)->choice, b0_);                                             \
        if (!b0_) { RC_TREE((l)->low[ps], 3, v0_); out = 2 + v0_; break; }   \
        RC_BIT(&(l)->choice2, b0_);                                            \
        if (!b0_) { RC_TREE((l)->mid[ps], 3, v0_); out = 10 + v0_; break; }  \
        RC_TREE((l)->high, 8, v0_); out = 18 + v0_;                            \
    } while (0)

static __attribute__((noinline)) int decode_tokens(lzr_dctx *d, const uint8_t *in, size_t avail, uint8_t *dst, size_t raw)
{
    arc_dec r;
    if (arcd_init(&r, in, avail)) return LZR_ERR_CORRUPT;
    uint32_t range = r.range, code = r.code;
    const uint8_t *ip = r.p, *const iend = r.end;
    tok_model *m = &d->tm;
    unsigned s = m->state;
    uint32_t rep0 = d->reps[0], rep1 = d->reps[1], rep2 = d->reps[2];
    const uint8_t *hist = d->hist;
    const size_t wsize = d->wsize;
    uint8_t *o = dst, *const oend = dst + raw;
    int err = 0;
    while (o < oend) {
        unsigned ps = (unsigned)(o - dst) & (TOK_POS - 1), bit;
        if (ip > iend) { err = 1; break; }
        RC_BIT(&m->is_match[s][ps], bit);
        if (!bit) {
            arc_p *probs = m->lit[(o > hist ? o[-1] : 0) >> 5];
            unsigned v = 1;
            if (s >= 7) {
                /* matched literal: probs[0x100 + (match bit << 8) + v]
                 * while the decoded bits equal the match byte's, then
                 * probs[v]; offs drops to 0 at the first difference */
                if (rep0 > (size_t)(o - hist)) { err = 1; break; }
                unsigned mbyte = (unsigned)o[-(ptrdiff_t)rep0] << 1, offs = 0x100;
                unsigned idx = offs + (mbyte & offs) + 1;
                uint32_t pc = probs[idx];
                for (int i = 0; i < 8; i++) {
                    /* the next index for either outcome, loaded early */
                    unsigned mb = mbyte & offs, mn = mbyte << 1;
                    unsigned o0 = offs & ~mb, o1 = offs & mb;
                    uint32_t pl = 0, ph = 0;
                    if (i < 7) {   /* index < 0x300 */
                        pl = probs[o0 + (mn & o0) + 2 * v];
                        ph = probs[o1 + (mn & o1) + 2 * v + 1];
                    }
                    uint32_t bd = (range >> ARC_PBITS) * pc, mm;
                    RC_SEL(mm, bd);
                    probs[idx] = RC_PUPD(pc, mm);
                    bit = mm + 1;
                    v = 2 * v + bit;
                    offs = o1 ^ ((o0 ^ o1) & mm);
                    idx = offs + (mn & offs) + v;
                    pc = ph ^ ((pl ^ ph) & mm);
                    mbyte = mn;
                    RC_NORM();
                }
            } else {
                RC_TREE(probs, 8, v);
            }
            v &= 0xFF;
            *o++ = (uint8_t)v;
            s = tok_after_lit(s);
            continue;
        }
        size_t len;
        uint32_t dist;
        RC_BIT(&m->is_rep[s], bit);
        if (!bit) {
            unsigned l, slot;
            RC_LEN(&m->mlen, ps, l);
            len = l;
            unsigned lctx = l - 2 < 3 ? l - 2 : 3;
            RC_TREE(m->slot[lctx], 6, slot);
            uint32_t d0 = slot;
            if (slot >= 4) {
                unsigned nb = (slot >> 1) - 1, v;
                d0 = (2u | (slot & 1)) << nb;
                if (slot < TOK_END_SPEC) {
                    RC_RTREE(m->spec[slot], nb, v);
                    d0 += v;
                } else {
                    uint32_t hi = 0;
                    for (unsigned k = nb - 4; k; k--) {
                        range >>= 1;
                        uint32_t t = code >= range;
                        code -= range & (0u - t);
                        hi = (hi << 1) | t;
                        RC_NORM();
                    }
                    RC_RTREE(m->align, 4, v);
                    d0 += (hi << 4) + v;
                }
            }
            if (d0 >= wsize) { err = 1; break; }
            dist = d0 + 1;
            rep2 = rep1; rep1 = rep0; rep0 = dist;
            s = tok_after_match(s);
        } else {
            RC_BIT(&m->is_rep0[s], bit);
            if (!bit) {
                dist = rep0;
                RC_BIT(&m->is_rep0_long[s][ps], bit);
                if (!bit) {   /* short rep */
                    if (dist > (size_t)(o - hist)) { err = 1; break; }
                    *o = o[-(ptrdiff_t)dist];
                    o++;
                    s = tok_after_short(s);
                    continue;
                }
            } else {
                RC_BIT(&m->is_rep1[s], bit);
                if (!bit) { dist = rep1; rep1 = rep0; }
                else { dist = rep2; rep2 = rep1; rep1 = rep0; }
                rep0 = dist;
            }
            unsigned l;
            RC_LEN(&m->rlen, ps, l);
            len = l;
            s = tok_after_rep(s);
        }
        if (len > (size_t)(oend - o) || dist > (size_t)(o - hist) || dist > wsize) { err = 1; break; }
        copy_match(o, dist, len);
        o += len;
    }
    m->state = (uint8_t)s;
    d->reps[0] = rep0; d->reps[1] = rep1; d->reps[2] = rep2;
    if (err || ip != iend) return LZR_ERR_CORRUPT;
    return LZR_OK;
}

/* Block type 5 (see emit_seqs5 in enc.c): Huffman literals into d->lits,
 * then one loop that decodes each sequence from the tANS bit stream and
 * executes it at once. */
static inline __attribute__((always_inline)) int decode_block5_impl(lzr_dctx *d, const uint8_t *in, size_t avail, uint8_t *dst, size_t raw)
{
    if (avail < 8) return LZR_ERR_CORRUPT;
    size_t ns = rd32le(in), nlit = rd32le(in + 4), k = 8;
    if (nlit > raw || ns > MAX_SEQ || ns * LZR_ML_BASE + nlit > raw) return LZR_ERR_CORRUPT;
    if (nlit) {
        size_t c = huf_decode(in + k, avail - k, d->lits, nlit, &d->huf);
        if (!c) return LZR_ERR_CORRUPT;
        k += c;
    }
    const uint8_t *lp = d->lits, *const lend = d->lits + nlit;
    uint8_t *o = dst, *const oend = dst + raw;
    if (!ns) {
        if (k != avail || nlit != raw) return LZR_ERR_CORRUPT;
        memcpy(o, lp, nlit);
        return LZR_OK;
    }
    const unsigned alpha[3] = {VAL_CODES, ML_CODES, OF_CODES};
    unsigned rsym[3] = {0, 0, 0};
    int rle[3];
    for (int t = 0; t < 3; t++) {
        if (k >= avail) return LZR_ERR_CORRUPT;
        unsigned mode = in[k++];
        fse_table *T = &d->fse[t];
        rle[t] = mode == 2;
        if (mode == 2) {
            if (k >= avail || in[k] >= alpha[t]) return LZR_ERR_CORRUPT;
            rsym[t] = in[k++];
        } else if (mode == 1) {
            if (!T->valid || T->alphabet != alpha[t]) return LZR_ERR_CORRUPT;
        } else if (mode == 0) {
            size_t c = fse_read_norm(in + k, avail - k, T->norm, alpha[t]);
            T->valid = 0;
            T->alphabet = alpha[t];
            if (!c || fse_build(T)) return LZR_ERR_CORRUPT;
            k += c;
        } else {
            return LZR_ERR_CORRUPT;
        }
    }
    /* Sequence tables: the tANS state step plus the code's base value and
     * extra-bit count in one entry. Offsets: value < 3 is a repeat index,
     * else distance + 3. A single-symbol stream is a one-state table. */
    seq_dent *tb[3] = {d->sq[0], d->sq[1], d->sq[2]};   /* contiguous */
    unsigned tlog[3];
    for (int t = 0; t < 3; t++) {
        const fse_table *T = &d->fse[t];
        unsigned n = rle[t] ? 1 : FSE_L;
        tlog[t] = rle[t] ? 0 : FSE_LOG;
        for (unsigned u = 0; u < n; u++) {
            unsigned sym = rle[t] ? rsym[t] : T->dt[u].sym;
            seq_dent *e = &tb[t][u];
            e->nb = rle[t] ? 0 : T->dt[u].nb;
            /* states index the three tables as one array */
            e->next = (uint16_t)((rle[t] ? 0 : T->dt[u].base) + t * FSE_L);
            if (t == 0) { e->base = d->ct[sym].base; e->xb = d->ct[sym].bits; }
            else if (t == 1) { e->base = d->mt[sym].base + LZR_ML_BASE; e->xb = d->mt[sym].bits; }
            else if (sym < LZR_NUM_REPS) { e->base = sym; e->xb = 0; }
            else { e->base = d->ct[sym - LZR_NUM_REPS].base + 1 + 3; e->xb = d->ct[sym - LZR_NUM_REPS].bits; }
        }
    }

    /* backward bit reader in locals (see fse.h) */
    const uint8_t *const bstart = in + k;
    const size_t bsize = avail - k;
    if (!bsize || !bstart[bsize - 1]) return LZR_ERR_CORRUPT;
    const uint8_t *bp = bstart + bsize - 8;
    uint64_t bc = fse_ld64(bp);
    unsigned used = (unsigned)__builtin_clzll(bc) + 1;
    int berr = 0;
#define BGET(nb) ((uint32_t)(((bc << (used & 63)) >> 1) >> (63 - (nb))))
#define BGET64(nb) (((bc << (used & 63)) >> 1) >> (63 - (nb)))
#define BRELOAD()                                                              \
    do {                                                                       \
        if (bp >= bstart + 8) { bp -= used >> 3; used &= 7; }                  \
        else {                                                                 \
            size_t nby_ = used >> 3, room_ = bp > bstart ? (size_t)(bp - bstart) : 0; \
            if (nby_ > room_) nby_ = room_;                                    \
            bp -= nby_; used -= (unsigned)nby_ * 8;                            \
            if (used > 64) { berr = 1; used = 64; }                            \
        }                                                                      \
        bc = fse_ld64(bp);                                                     \
    } while (0)

    unsigned sll = BGET(tlog[0]); used += tlog[0];
    unsigned sml = BGET(tlog[1]) + FSE_L; used += tlog[1];
    unsigned sof = BGET(tlog[2]) + 2 * FSE_L; used += tlog[2];
    BRELOAD();
    const seq_dent *const sq = tb[0];
    uint32_t rep0 = d->reps[0], rep1 = d->reps[1], rep2 = d->reps[2];
    const uint8_t *hist = d->hist;
    const size_t wsize = d->wsize;

    /* Decode a sequence (ll, ml, dist), advancing the bit stream. */
#define DECODE_SEQ(LL, ML, DIST)                                               \
    do {                                                                       \
        const seq_dent el_ = sq[sll], em_ = sq[sml], eo_ = sq[sof];                 \
        /* extra bits (ll, ml, of written in that order), then the three      \
         * state updates (of, ml, ll): one read and one reload when they fit  \
         * in 56 bits (nearly always), else two of each */                    \
        const unsigned xl_ = el_.xb, xm_ = em_.xb, xo_ = eo_.xb;               \
        const unsigned no_ = eo_.nb, nm_ = em_.nb, nl_ = el_.nb;               \
        const unsigned xs_ = xl_ + xm_ + xo_, ss_ = no_ + nm_ + nl_;           \
        uint64_t v_;                                                           \
        uint32_t w_;                                                           \
        if (LZR_LIKELY(xs_ + ss_ <= 56)) {                                     \
            uint64_t t_ = BGET64(xs_ + ss_); used += xs_ + ss_;                \
            v_ = t_ >> ss_;                                                    \
            w_ = (uint32_t)(t_ & ((1u << ss_) - 1));                           \
            BRELOAD();                                                         \
        } else {                                                               \
            v_ = BGET64(xs_); used += xs_;                                     \
            BRELOAD();                                                         \
            w_ = BGET(ss_); used += ss_;                                       \
            BRELOAD();                                                         \
        }                                                                      \
        LL = el_.base + (uint32_t)(v_ & ((1u << xl_) - 1));                    \
        ML = em_.base + (uint32_t)((v_ >> xl_) & ((1u << xm_) - 1));           \
        uint32_t ov_ = eo_.base + (uint32_t)(v_ >> (xl_ + xm_));               \
        /* repeat index 0-2 or distance + 3; written for conditional moves */ \
        uint32_t r_ = ov_ == 0 ? rep0 : ov_ == 1 ? rep1 : rep2;                \
        DIST = ov_ >= LZR_NUM_REPS ? ov_ - 3 : r_;                             \
        rep2 = ov_ >= 2 ? rep1 : rep2;                                         \
        rep1 = ov_ >= 1 ? rep0 : rep1;                                         \
        rep0 = DIST;                                                           \
        sof = eo_.next + (w_ & ((1u << no_) - 1));                             \
        sml = em_.next + ((w_ >> no_) & ((1u << nm_) - 1));                    \
        sll = el_.next + (w_ >> (no_ + nm_));                                  \
    } while (0)
    /* Execute it: literal run, then the match (checked here). */
#define EXEC_SEQ(LL, ML, DIST)                                                 \
    do {                                                                       \
        if ((LL) > (size_t)(lend - lp) || (LL) + (ML) > (size_t)(oend - o)) { berr = 1; goto done; } \
        copy_lits(o, lp, LL);                                                  \
        o += LL;                                                               \
        lp += LL;                                                              \
        if ((DIST) > (size_t)(o - hist) || (DIST) > wsize) { berr = 1; goto done; } \
        copy_match(o, DIST, ML);                                               \
        o += ML;                                                               \
    } while (0)

    /* (Decoding a few sequences ahead with prefetching measured no faster
     * here, even past the cache size, so the loop stays simple.) */
    for (size_t i = 0; i < ns; i++) {
        size_t ll, ml;
        uint32_t dist;
        DECODE_SEQ(ll, ml, dist);
        EXEC_SEQ(ll, ml, dist);
    }
done:
#undef DECODE_SEQ
#undef EXEC_SEQ
#undef BGET
#undef BRELOAD
    d->reps[0] = rep0; d->reps[1] = rep1; d->reps[2] = rep2;
    /* all bits consumed, exactly */
    if (berr || (size_t)(bstart + bsize - (bp + 8)) * 8 + used != bsize * 8) return LZR_ERR_CORRUPT;
    size_t tail = (size_t)(lend - lp);
    if (tail != (size_t)(oend - o)) return LZR_ERR_CORRUPT;
    memcpy(o, lp, tail);
    return LZR_OK;
}

static __attribute__((noinline)) int decode_block5_generic(lzr_dctx *d, const uint8_t *in, size_t avail, uint8_t *dst, size_t raw)
{
    return decode_block5_impl(d, in, avail, dst, raw);
}

#ifdef LZR_BMI2_DISPATCH
static __attribute__((noinline)) LZR_TARGET_BMI2 int decode_block5_bmi2(lzr_dctx *d, const uint8_t *in, size_t avail, uint8_t *dst, size_t raw)
{
    return decode_block5_impl(d, in, avail, dst, raw);
}
#endif

static int decode_block5(lzr_dctx *d, const uint8_t *in, size_t avail, uint8_t *dst, size_t raw)
{
#ifdef LZR_BMI2_DISPATCH
    if (LZR_HAS_BMI2()) return decode_block5_bmi2(d, in, avail, dst, raw);
#endif
    return decode_block5_generic(d, in, avail, dst, raw);
}

static int init_history(lzr_dctx *d, int wlog)
{
    d->wsize = (size_t)1 << wlog;
    d->hcap = 2 * d->wsize + LZR_BLOCK_SIZE;
    d->hist = (uint8_t *)big_alloc(d->hcap + LZR_SLACK, 0);
    return d->hist ? LZR_OK : LZR_ERR_NOMEM;
}

typedef struct { lzr_dctx *d; lzr_buf *out; } out_ctx;

/* Deliver decoded (unfiltered) content; the CRC covers the original. */
static int emit_out(void *ctx, const uint8_t *p, size_t n)
{
    out_ctx *s = (out_ctx *)ctx;
    xxh64_update(&s->d->xs, p, n);
    return buf_append(s->out, p, n);
}

/* Process one unit at in.p + ipos. Returns 1 if consumed, 0 if more input
 * is needed, 2 at the verified end, or LZR_ERR_*. */
static int step(lzr_dctx *d, lzr_buf *out)
{
    const uint8_t *p = d->in.p + d->ipos;
    size_t avail = d->in.len - d->ipos;

    if (d->state == ST_HEADER) {
        int wlog;
        uint8_t flags;
        int r = lzr_frame_header(p, avail, &wlog, &flags, &d->content_size);
        if (r <= 0) return r;
        if (flags & LZR_FLAG_X86) d->direct = NULL;   /* history must hold filtered bytes */
        if (d->direct) d->wsize = (size_t)1 << wlog;
        else if (init_history(d, wlog)) return LZR_ERR_NOMEM;
        if ((flags & LZR_FLAG_X86) && !(d->x86tmp = (uint8_t *)malloc(X86_CHUNK + 4))) return LZR_ERR_NOMEM;
        d->ipos += (size_t)r;
        d->state = ST_BLOCKS;
        return 1;
    }

    if (avail < LZR_BLOCK_HDR) return 0;
    uint8_t type = p[0];
    size_t raw = rd32le(p + 1), comp = rd32le(p + 5);
    if (type > 5 || type == 1 || type == 3 || raw > LZR_BLOCK_SIZE || comp > MAX_COMP) return LZR_ERR_CORRUPT;
    if (type == 0 && comp != raw) return LZR_ERR_CORRUPT;
    if (type == 2 && (raw != 0 || comp != 4)) return LZR_ERR_CORRUPT;
    if (type >= 4 && raw == 0) return LZR_ERR_CORRUPT;
    if (avail - LZR_BLOCK_HDR < comp) return 0;
    const uint8_t *payload = p + LZR_BLOCK_HDR;

    if (type == 2) {
        if (d->content_size >= 0 && (uint64_t)d->content_size != d->total_out) return LZR_ERR_CORRUPT;
        if (d->x86tmp) {
            out_ctx s = {d, out};
            if (x86_finish(&d->x86, emit_out, &s)) return LZR_ERR_NOMEM;
        }
        d->ipos += LZR_BLOCK_HDR + comp;
        d->state = ST_DONE;
        return rd32le(payload) == (uint32_t)xxh64_digest(&d->xs) ? 2 : LZR_ERR_CHECKSUM;
    }
    if (d->content_size >= 0 && d->total_out + raw > (uint64_t)d->content_size) return LZR_ERR_CORRUPT;

    if (d->direct) {   /* the output itself is the history: no copy */
        if (out != d->direct) return LZR_ERR_PARAM;
        if (buf_reserve(out, raw + LZR_SLACK)) return LZR_ERR_NOMEM;
        d->hist = out->p;
        d->op = out->len;
    } else if (d->op + raw > d->hcap) {   /* make room: keep the last wsize bytes */
        size_t delta = d->op - d->wsize;
        memmove(d->hist, d->hist + delta, d->wsize);
        d->op = d->wsize;
    }
    uint8_t *dst = d->hist + d->op;
    if (type == 0) {
        memcpy(dst, payload, raw);
    } else {
        int r = type == 4 ? decode_tokens(d, payload, comp, dst, raw) : decode_block5(d, payload, comp, dst, raw);
        if (r) return r;
    }
    d->op += raw;
    d->total_out += raw;
    d->ipos += LZR_BLOCK_HDR + comp;
    if (d->direct) {
        xxh64_update(&d->xs, dst, raw);
        out->len += raw;
        return 1;
    }
    out_ctx s = {d, out};
    /* history keeps the filtered bytes (matches refer to them) */
    int r = d->x86tmp ? x86_run(&d->x86, dst, raw, 1, d->x86tmp, emit_out, &s) : emit_out(&s, dst, raw);
    return r ? LZR_ERR_NOMEM : 1;
}

int lzr_dctx_set_output(lzr_dctx *d, lzr_buf *out)
{
    if (d->state != ST_HEADER || d->in.len) return LZR_ERR_PARAM;
    d->direct = out;
    return LZR_OK;
}

int lzr_dctx_write(lzr_dctx *d, const uint8_t *src, size_t n, lzr_buf *out, size_t *consumed)
{
    if (consumed) *consumed = 0;
    if (d->state == ST_DONE) return 1;
    /* Keep ARC_SLACK readable bytes after the input (see arc.h). */
    if (n && (buf_reserve(&d->in, n + ARC_SLACK) || buf_append(&d->in, src, n))) return LZR_ERR_NOMEM;
    if (n) memset(d->in.p + d->in.len, 0, ARC_SLACK);

    int r;
    while ((r = step(d, out)) == 1) {}
    if (r < 0) return r;

    size_t left = d->in.len - d->ipos;
    if (r == 2) {                 /* bytes after the end block are not ours */
        if (consumed) *consumed = n - left;
        d->in.len = d->ipos = 0;
        return 1;
    }
    if (consumed) *consumed = n;
    if (d->ipos && (left == 0 || d->ipos >= (1u << 20))) {
        memmove(d->in.p, d->in.p + d->ipos, left);
        d->in.len = left;
        d->ipos = 0;
    }
    return LZR_OK;
}
