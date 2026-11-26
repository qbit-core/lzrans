/*
 * fuzz_roundtrip.c - standalone ASan/UBSan stress test for the C core.
 * Build: make -C tests/c   (or see Makefile). Not part of the Python build.
 *   1. round-trips many generated inputs at every level
 *   2. flips/truncates bytes of valid frames; decoder must return an
 *      error (or the exact original) and never touch memory out of bounds.
 */
#include <stdio.h>
#include <stdlib.h>
#include "../../src/lzrans/frame.h"

static uint64_t rng = 88172645463325252ull;
static uint32_t rnd(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return (uint32_t)rng; }

static void gen(uint8_t *b, size_t n, int kind)
{
    for (size_t i = 0; i < n; i++) {
        switch (kind) {
        case 0: b[i] = (uint8_t)rnd(); break;                       /* random      */
        case 1: b[i] = 'a'; break;                                  /* constant    */
        case 2: b[i] = (uint8_t)("abcabd"[i % 6]); break;           /* short period */
        case 3: b[i] = (uint8_t)(rnd() % 4); break;                 /* low entropy */
        default: /* text-like with copies from far back */
            b[i] = (i > 70000 && rnd() % 8) ? b[i - 1 - rnd() % 65536] : (uint8_t)(32 + rnd() % 64);
        }
    }
}

int main(int argc, char **argv)
{
    int iters = argc > 1 ? atoi(argv[1]) : 200;
    size_t maxn = 3 * LZR_BLOCK_SIZE + 123;
    uint8_t *src = malloc(maxn), *dst = malloc(maxn), *cmp = malloc(lzr_compress_bound(maxn));
    long fails = 0, rejected = 0;
    for (int it = 0; it < iters; it++) {
        size_t n = (it % 10 == 0) ? (rnd() % maxn) : (rnd() % 5000);
        int kind = rnd() % 5, level = 1 + rnd() % 9;
        gen(src, n, kind);
        int64_t c = lzr_compress(src, n, cmp, lzr_compress_bound(n), level);
        if (c < 0 || lzr_decompress(cmp, (size_t)c, dst, n) != 0 || memcmp(src, dst, n)) {
            printf("ROUNDTRIP FAIL it=%d n=%zu kind=%d level=%d\n", it, n, kind, level);
            fails++;
            continue;
        }
        /* Mutation: corrupt a few bytes or truncate, decode into a fresh
         * exactly-sized heap buffer so ASan catches any overflow. */
        for (int m = 0; m < 20; m++) {
            uint8_t *mut = malloc((size_t)c);
            memcpy(mut, cmp, (size_t)c);
            size_t mlen = (size_t)c;
            if (m % 4 == 0 && mlen > 1) mlen = rnd() % mlen;
            else for (int f = 0; f < 1 + m % 3; f++) mut[rnd() % mlen] ^= (uint8_t)(1 + rnd() % 255);
            uint64_t orig;
            if (lzr_frame_info(mut, mlen, &orig, NULL, NULL) == 0 && orig < 64u * 1024 * 1024) {
                uint8_t *o = malloc(orig ? orig : 1);
                int r = lzr_decompress(mut, mlen, o, orig);
                if (r != 0) rejected++;
                else if (orig != n || memcmp(o, src, n)) { printf("SILENT CORRUPTION it=%d\n", it); fails++; }
                free(o);
            } else rejected++;
            free(mut);
        }
    }
    free(src); free(dst); free(cmp);
    printf("iters=%d fails=%ld corrupt_rejected=%ld\n", iters, fails, rejected);
    return fails != 0;
}
