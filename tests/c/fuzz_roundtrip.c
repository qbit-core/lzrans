/*
 * fuzz_roundtrip.c - standalone ASan/UBSan stress test for the C core.
 * Build and run:  make -C tests/c run
 *   1. round-trips generated inputs at every level (with and without the x86
 *      filter), feeding the encoder and
 *      decoder in random chunk sizes (streaming path);
 *   2. flips / truncates bytes of valid frames; the decoder must return an
 *      error (or the exact original) and never touch memory out of bounds.
 */
#include <stdio.h>
#include "../../src/lzrans/core.h"

static uint64_t rng = 88172645463325252ull;
static uint32_t rnd(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return (uint32_t)rng; }

static void gen(uint8_t *b, size_t n, int kind)
{
    for (size_t i = 0; i < n; i++) {
        switch (kind) {
        case 0: b[i] = (uint8_t)rnd(); break;                         /* random       */
        case 1: b[i] = 'a'; break;                                    /* constant     */
        case 2: b[i] = (uint8_t)("abcabd"[i % 6]); break;             /* short period */
        case 3: b[i] = (uint8_t)(rnd() % 4); break;                   /* low entropy  */
        case 5: b[i] = (rnd() % 6) ? (uint8_t)(rnd() % 3 ? 0 : rnd()) : 0xE8; break;   /* call-dense */
        case 6: b[i] = (uint8_t)(i ? 32 + (b[i - 1] * 7 + rnd() % 5) % 90 : 32); break;   /* context literals */
        default: /* text-like with copies from near and far back */
            b[i] = (i > 300000 && rnd() % 8) ? b[i - 1 - rnd() % (i > 2000000 ? 2000000 : 300000)]
                                              : (uint8_t)(32 + rnd() % 64);
        }
    }
}

static int compress_chunked(const uint8_t *src, size_t n, int level, lzr_buf *out)
{
    lzr_cctx *c = lzr_cctx_new(level, 0, (rnd() & 1) ? (int64_t)n : -1);
    if (!c) return -1;
    lzr_cctx_set_threads(c, 1 + (int)(rnd() % 4));
    if (rnd() & 1) lzr_cctx_set_filter(c, LZR_FILTER_X86);
    size_t i = 0;
    int r = 0;
    while (!r && i < n) {
        size_t k = rnd() % 400000;
        if (k > n - i) k = n - i;
        r = lzr_cctx_write(c, src + i, k, out);
        i += k;
    }
    if (!r) r = lzr_cctx_end(c, out);
    lzr_cctx_free(c);
    return r;
}

/* Returns 1 = verified end, <0 error, 0 = truncated. Every other call
 * uses one-buffer mode (decoding straight into out). */
static int decompress_chunked(const uint8_t *src, size_t n, lzr_buf *out)
{
    static unsigned calls;
    lzr_dctx *d = lzr_dctx_new();
    if (calls++ & 1) lzr_dctx_set_output(d, out);
    size_t i = 0;
    int r = 0;
    while (r == 0 && i < n) {
        size_t k = 1 + rnd() % 100000;
        if (k > n - i) k = n - i;
        r = lzr_dctx_write(d, src + i, k, out, NULL);
        i += k;
    }
    lzr_dctx_free(d);
    return r;
}

int main(int argc, char **argv)
{
    int iters = argc > 1 ? atoi(argv[1]) : 200;
    size_t maxn = 5 * LZR_BLOCK_SIZE * 2 + 123;
    uint8_t *src = malloc(maxn);
    long fails = 0, rejected = 0;
    for (int it = 0; it < iters; it++) {
        size_t n = (it % 10 == 0) ? (rnd() % maxn) : (rnd() % 5000);
        int kind = rnd() % 7, level = (int)(rnd() % 15) - 3;   /* LZR_LEVEL_MIN..MAX */
        gen(src, n, kind);
        lzr_buf comp = {0}, dec = {0};
        if (compress_chunked(src, n, level, &comp) || decompress_chunked(comp.p, comp.len, &dec) != 1 ||
            dec.len != n || (n && memcmp(src, dec.p, n))) {
            printf("ROUNDTRIP FAIL it=%d n=%zu kind=%d level=%d\n", it, n, kind, level);
            fails++;
        }
        for (int m = 0; m < 20 && !fails; m++) {
            uint8_t *mut = malloc(comp.len);
            memcpy(mut, comp.p, comp.len);
            size_t mlen = comp.len;
            if (m % 4 == 0) mlen = rnd() % mlen;
            else for (int f = 0; f < 1 + m % 3; f++) mut[rnd() % mlen] ^= (uint8_t)(1 + rnd() % 255);
            lzr_buf o = {0};
            int r = decompress_chunked(mut, mlen, &o);
            if (r != 1) rejected++;
            else if (o.len != n || (n && memcmp(o.p, src, n))) { printf("SILENT CORRUPTION it=%d\n", it); fails++; }
            free(o.p);
            free(mut);
        }
        free(comp.p);
        free(dec.p);
    }
    /* segment mode (reset points + seek index) and the parallel decoder:
     * same content as the serial decoder; damaged frames rejected */
    size_t segn = 3 * (1u << 20) + 77;
    uint8_t *big = malloc(segn);
    for (int it = 0; it < iters / 20 + 2 && !fails; it++) {
        size_t n = (1u << 20) + rnd() % (segn - (1u << 20));
        int kind = rnd() % 7, level = (int)(rnd() % 15) - 3;
        gen(big, n, kind);
        lzr_cctx *c = lzr_cctx_new_seg(level, 0, (rnd() & 1) ? (int64_t)n : -1, (1u << 20) + rnd() % 300000);
        lzr_buf comp = {0};
        lzr_cctx_set_threads(c, 1 + (int)(rnd() % 4));
        if (rnd() & 1) lzr_cctx_set_filter(c, LZR_FILTER_X86);
        int r = lzr_cctx_write(c, big, n, &comp);
        if (!r) r = lzr_cctx_end(c, &comp);
        lzr_cctx_free(c);
        for (int m = -1; m < 30 && !fails; m++) {
            uint8_t *mut = malloc(comp.len);
            memcpy(mut, comp.p, comp.len);
            size_t mlen = comp.len;
            if (m >= 0 && m % 5 == 0) mlen = rnd() % mlen;
            else if (m >= 0) for (int f = 0; f < 1 + m % 3; f++) mut[rnd() % mlen] ^= (uint8_t)(1 + rnd() % 255);
            lzr_buf o = {0}, o2 = {0};
            size_t used = 0;
            int r1 = lzr_decompress_mt(mut, mlen, &o, 2 + (int)(rnd() % 3), LZR_WLOG_DEC_DEFAULT, &used);
            int r2 = decompress_chunked(mut, mlen, &o2);
            if (m < 0 && (r || r1 != 1 || o.len != n || memcmp(o.p, big, n))) {
                printf("SEGMENT ROUNDTRIP FAIL it=%d n=%zu kind=%d level=%d r=%d r1=%d\n", it, n, kind, level, r, r1);
                fails++;
            }
            if (r1 == 1 && (o.len != n || memcmp(o.p, big, n))) { printf("SEGMENT SILENT CORRUPTION it=%d\n", it); fails++; }
            if ((r1 == 1) != (r2 == 1 && used == mlen) && !(r1 == 1 && used < mlen)) {
                printf("SEGMENT MISMATCH it=%d m=%d r1=%d r2=%d\n", it, m, r1, r2);
                fails++;
            }
            if (r1 != 1) rejected++;
            free(o.p);
            free(o2.p);
            free(mut);
        }
        free(comp.p);
    }
    free(big);
    free(src);
    printf("iters=%d fails=%ld corrupt_rejected=%ld\n", iters, fails, rejected);
    return fails != 0;
}
