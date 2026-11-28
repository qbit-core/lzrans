/*
 * ldm.c - long-range dedup prototype (measurement tool, not part of the
 * library). Finds long matches over the whole input with a content-defined
 * rolling hash (zstd LDM style) and splits the input into
 *   residual: the bytes not covered by a long match (to be compressed by a
 *             normal lzrans level), and
 *   refs:     one record per long match (literal run, distance, length).
 * `ldm d` rebuilds the input from both, to prove the split is lossless.
 *
 *   ldm e IN RESIDUAL REFS [minlen=256] [min_dist=0] [anchor_log=6]
 *   ldm d RESIDUAL REFS OUT
 *
 * Hash: Gear (h = 2h + G[b]) whose value depends on the last 64 bytes. A
 * position is an anchor when the top `anchor_log` bits of h are zero
 * (content defined, so inserted / deleted bytes shift anchors with the
 * data). Anchors index a table of 2^tbl_log buckets x 4 ways (position +
 * 32-bit check); a lookup verifies bytes, extends forward and backward
 * (back into the pending literals) and keeps the longest candidate.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define WAYS 4
typedef struct { uint64_t pos; uint32_t chk; uint32_t pad; } slot_t;

static uint64_t G[256];
static void gear_init(void) {
    uint64_t s = 0x9E3779B97F4A7C15ull;
    for (int i = 0; i < 256; i++) {           /* splitmix64 */
        uint64_t z = (s += 0x9E3779B97F4A7C15ull);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        G[i] = z ^ (z >> 31);
    }
}

static uint8_t *readall(const char *p, size_t *n) {
    FILE *f = fopen(p, "rb");
    if (!f) { perror(p); exit(1); }
    fseek(f, 0, SEEK_END); *n = (size_t)ftell(f); fseek(f, 0, SEEK_SET);
    uint8_t *b = malloc(*n + 1);
    if (*n && fread(b, 1, *n, f) != *n) { perror("read"); exit(1); }
    fclose(f);
    return b;
}

static void put_varint(FILE *f, uint64_t v) {
    while (v >= 0x80) { fputc((int)(v & 0x7f) | 0x80, f); v >>= 7; }
    fputc((int)v, f);
}
static int get_varint(const uint8_t **p, const uint8_t *e, uint64_t *v) {
    uint64_t r = 0; int sh = 0;
    while (*p < e) {
        uint8_t b = *(*p)++;
        r |= (uint64_t)(b & 0x7f) << sh;
        if (!(b & 0x80)) { *v = r; return 1; }
        sh += 7;
    }
    return 0;
}

#define W 64   /* gear window */

static int encode(int argc, char **argv) {
    size_t n; uint8_t *in = readall(argv[2], &n);
    size_t minlen = argc > 5 ? strtoull(argv[5], 0, 10) : 256;
    uint64_t min_dist = argc > 6 ? strtoull(argv[6], 0, 10) : 0;
    int alog = argc > 7 ? atoi(argv[7]) : 6;
    int tlog = 22;
    while (tlog < 26 && ((size_t)1 << (tlog + alog)) < n) tlog++;
    size_t nb = (size_t)1 << tlog;
    slot_t *tab = calloc(nb * WAYS, sizeof(slot_t));
    FILE *fr = fopen(argv[3], "wb"), *fm = fopen(argv[4], "wb");
    if (!tab || !fr || !fm) { perror("alloc/open"); return 1; }
    clock_t t0 = clock();
    const uint64_t amask = alog ? (~0ull << (64 - alog)) : 0;
    uint64_t h = 0;
    size_t lit = 0;        /* start of pending literals */
    size_t nmatch = 0, covered = 0;
    size_t i = 0;
    while (i < n) {
        h = (h << 1) + G[in[i]];
        i++;
        if (i < W || (h & amask)) continue;
        /* anchor: window [i-W, i) */
        size_t s = i - W;
        uint32_t chk = (uint32_t)(h >> 8);
        slot_t *b = tab + ((h >> 24) & (nb - 1)) * WAYS;
        size_t best_len = 0, best_start = 0, best_src = 0;
        if (s >= lit) {
            for (int w = 0; w < WAYS; w++) {
                if (b[w].chk != chk || !b[w].pos) continue;
                size_t c = (size_t)b[w].pos - 1;      /* window start of candidate */
                if (s - c < min_dist || c >= s) continue;
                if (memcmp(in + c, in + s, W)) continue;
                size_t f = W;
                while (s + f < n && in[c + f] == in[s + f]) f++;
                size_t bk = 0;
                while (bk < c && s - bk > lit && in[c - bk - 1] == in[s - bk - 1]) bk++;
                if (f + bk > best_len) { best_len = f + bk; best_start = s - bk; best_src = c - bk; }
            }
        }
        /* insert (newest first: shift ways) */
        memmove(b + 1, b, (WAYS - 1) * sizeof(slot_t));
        b[0].pos = s + 1; b[0].chk = chk;
        if (best_len >= minlen) {
            put_varint(fm, best_start - lit);
            put_varint(fm, best_start - best_src);
            put_varint(fm, best_len - minlen);
            fwrite(in + lit, 1, best_start - lit, fr);
            nmatch++; covered += best_len;
            size_t e = best_start + best_len;
            /* index anchors inside the match too (keeps newest copies) but
             * do not search: roll the hash through the match. */
            while (i < e) {
                h = (h << 1) + G[in[i]];
                i++;
                if (!(h & amask) && i >= W) {
                    slot_t *bb = tab + ((h >> 24) & (nb - 1)) * WAYS;
                    memmove(bb + 1, bb, (WAYS - 1) * sizeof(slot_t));
                    bb[0].pos = i - W + 1; bb[0].chk = (uint32_t)(h >> 8);
                }
            }
            lit = e;
        }
    }
    fwrite(in + lit, 1, n - lit, fr);
    double dt = (double)(clock() - t0) / CLOCKS_PER_SEC;
    fprintf(stderr, "ldm: %zu bytes, %zu matches, covered %.2f%%, residual %zu, %.2fs (%.0f MB/s), table 2^%dx%d\n",
            n, nmatch, 100.0 * covered / (n ? n : 1), n - covered, dt, n / 1e6 / (dt > 0 ? dt : 1e-9), tlog, WAYS);
    fclose(fr); fclose(fm);
    /* minlen is needed by the decoder: write it as a trailer file name? keep it simple: env */
    return 0;
}

static int decode(int argc, char **argv) {
    (void)argc;
    size_t nr, nm; uint8_t *r = readall(argv[2], &nr), *m = readall(argv[3], &nm);
    size_t minlen = getenv("LDM_MINLEN") ? strtoull(getenv("LDM_MINLEN"), 0, 10) : 256;
    size_t cap = nr * 4 + (1 << 20), o = 0, rp = 0;
    uint8_t *out = malloc(cap);
    const uint8_t *p = m, *e = m + nm;
    uint64_t ll, d, ml;
    while (p < e) {
        if (!get_varint(&p, e, &ll) || !get_varint(&p, e, &d) || !get_varint(&p, e, &ml)) return 2;
        ml += minlen;
        while (o + ll + ml > cap) { cap *= 2; out = realloc(out, cap); }
        if (rp + ll > nr || d > o + ll) return 3;
        memcpy(out + o, r + rp, ll); o += ll; rp += ll;
        for (uint64_t k = 0; k < ml; k++) out[o + k] = out[o + k - d];
        o += ml;
    }
    while (o + (nr - rp) > cap) { cap *= 2; out = realloc(out, cap); }
    memcpy(out + o, r + rp, nr - rp); o += nr - rp;
    FILE *f = fopen(argv[4], "wb"); fwrite(out, 1, o, f); fclose(f);
    return 0;
}

int main(int argc, char **argv) {
    gear_init();
    if (argc >= 5 && argv[1][0] == 'e') return encode(argc, argv);
    if (argc >= 5 && argv[1][0] == 'd') return decode(argc, argv);
    fprintf(stderr, "usage: ldm e IN RES REFS [minlen] [min_dist] [anchor_log] | ldm d RES REFS OUT\n");
    return 1;
}
