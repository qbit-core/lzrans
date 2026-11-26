/*
 * x86.h - x86 branch filter (BCJ) for frame flag LZR_FLAG_X86.
 *
 * In x86 code, the rel32 operand of CALL (E8) is relative to the next
 * instruction, so calls to one function look different from each
 * place they are made. Converting it to an absolute target (rel + offset of
 * the next instruction) makes repeated calls repeat as bytes, which the LZ
 * stage can match.
 *
 * Rule (both directions scan the same opcode positions):
 *   scan the stream; at an E8 byte at offset i with 4 more bytes
 *   available, take v = the 4 bytes (LE) and skip them. Only "small"
 *   displacements, v in [-2^25, 2^25) (top 7 bits equal), are converted,
 *   to 26 bits: a = (v + i + 5) mod 2^26, written sign-extended from bit
 *   25. Every converted value passes the same test and no unconverted one
 *   does, so the decoder tests what it reads and inverts with
 *   v = (a - i - 5) mod 2^26, sign-extended.
 *   Opcode bytes are never changed, so the decoder (scanning the filtered
 *   bytes) stops at the same positions. The mapping is a bijection, so a
 *   false positive in data costs a little ratio but never correctness.
 *   Offsets are positions in the content, mod 2^32.
 *
 * Streaming: an opcode within 4 bytes of the input seen so far is held
 * back (at most 4 bytes) until more input arrives or the stream ends; at
 * the end, a truncated operand is left unchanged on both sides.
 */
#ifndef LZRANS_X86_H
#define LZRANS_X86_H
#include "common.h"

typedef struct {
    uint32_t pos;     /* content offset of hold[0] */
    uint8_t hold[4];
    unsigned nhold;
} lzr_x86;

/* Tuned on a 16 MB x86-64 shared object: 26 bits beat 20-30, and leaving
 * JMP (E9) alone beat converting it (+0.6%). */
#define X86_BITS 26   /* converted displacements are X86_BITS-bit signed */
#define X86_MASK ((1u << X86_BITS) - 1)

static inline int x86_small(uint32_t v)
{
    uint32_t top = v >> (X86_BITS - 1);
    return top == 0 || top == (0xFFFFFFFFu >> (X86_BITS - 1));
}

/* Filter buf[0..n) in place (buf[0] is at content offset pos). Returns how
 * many bytes are final; the rest (< 5, starting at an opcode) must be
 * filtered again together with following input. */
static inline size_t x86_filter(uint8_t *buf, size_t n, uint32_t pos, int decode)
{
    size_t i = 0;
    while (i < n) {
        if (buf[i] != 0xE8) { i++; continue; }
        if (n - i < 5) return i;
        uint32_t v = rd32le(buf + i + 1);
        if (x86_small(v)) {
            uint32_t next = pos + (uint32_t)i + 5;
            uint32_t a = (decode ? v - next : v + next) & X86_MASK;
            if (a >> (X86_BITS - 1)) a |= ~X86_MASK;
            wr32le(buf + i + 1, a);
        }
        /* Skip the operand even when unconverted: the decoder's test at i
         * must see bytes no later conversion has touched. */
        i += 5;
    }
    return n;
}

/* Run src[0..n) through the filter. Final bytes are passed to emit(ctx,
 * p, len) (which may be called several times); tmp holds X86_CHUNK + 4. */
#define X86_CHUNK 65536
typedef int (*x86_sink)(void *ctx, const uint8_t *p, size_t n);

static inline int x86_run(lzr_x86 *f, const uint8_t *src, size_t n, int decode, uint8_t *tmp, x86_sink emit,
                          void *ctx)
{
    while (n) {
        size_t take = n < X86_CHUNK ? n : X86_CHUNK;
        memcpy(tmp, f->hold, f->nhold);
        memcpy(tmp + f->nhold, src, take);
        size_t len = f->nhold + take;
        size_t done = x86_filter(tmp, len, f->pos, decode);
        f->nhold = (unsigned)(len - done);
        memcpy(f->hold, tmp + done, f->nhold);
        if (done) {
            int r = emit(ctx, tmp, done);
            if (r) return r;
        }
        f->pos += (uint32_t)done;
        src += take;
        n -= take;
    }
    return 0;
}

/* End of stream: the held bytes (a truncated operand) pass unchanged. */
static inline int x86_finish(lzr_x86 *f, x86_sink emit, void *ctx)
{
    int r = f->nhold ? emit(ctx, f->hold, f->nhold) : 0;
    f->pos += f->nhold;
    f->nhold = 0;
    return r;
}

#endif
