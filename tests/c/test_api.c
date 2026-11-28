/*
 * test_api.c - tests of the public liblzr API (include/lzr.h).
 *   make check                      (against build/c/liblzr.so)
 *   make -C tests/c api             (same, built from source under ASan/UBSan)
 *
 * Checks: one-shot round trips at every level; streaming with tiny and
 * odd-sized buffers gives the same bytes as one-shot; compress bound;
 * concatenated frames; parallel one-shot decode; error codes for truncated,
 * corrupt, trailing data, small buffers, memory limit and API misuse.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../include/lzr.h"

static int failures = 0;
#define CHECK(cond)                                                                    \
    do {                                                                               \
        if (!(cond)) {                                                                 \
            fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond);   \
            failures++;                                                                \
        }                                                                              \
    } while (0)
#define CHECK_RC(expr, want)                                                           \
    do {                                                                               \
        int rc_ = (expr);                                                              \
        if (rc_ != (want)) {                                                           \
            fprintf(stderr, "%s:%d: %s = %d (%s), want %d\n", __FILE__, __LINE__, #expr, \
                    rc_, lzr_error_string(rc_), (want));                               \
            failures++;                                                                \
        }                                                                              \
    } while (0)

static uint64_t rng = 0x9E3779B97F4A7C15ull;
static uint32_t rnd(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return (uint32_t)rng; }

/* text-like data with near and far repeats */
static uint8_t *gen(size_t n, int random_only)
{
    uint8_t *b = malloc(n ? n : 1);
    for (size_t i = 0; i < n; i++) {
        if (random_only) b[i] = (uint8_t)rnd();
        else if (i > 1000 && rnd() % 8) b[i] = b[i - 1 - rnd() % (i > 200000 ? 200000 : 1000)];
        else b[i] = (uint8_t)("etaoin shrdlu\n"[rnd() % 14]);
    }
    return b;
}

static uint8_t *oneshot(const uint8_t *src, size_t n, int level, int threads, int64_t segment, size_t *clen)
{
    size_t cap = lzr_compress_bound(n);
    uint8_t *dst = malloc(cap);
    lzr_cstream *cs = lzr_cstream_create();
    CHECK_RC(lzr_cstream_set(cs, LZR_C_LEVEL, level), LZR_OK);
    CHECK_RC(lzr_cstream_set(cs, LZR_C_THREADS, threads), LZR_OK);
    CHECK_RC(lzr_cstream_set(cs, LZR_C_SEGMENT, segment), LZR_OK);
    CHECK_RC(lzr_compress_cs(cs, dst, cap, clen, src, n), LZR_OK);
    lzr_cstream_free(cs);
    return dst;
}

/* Streaming compress with random input pieces and output room. */
static uint8_t *streamed(const uint8_t *src, size_t n, int level, int threads, int64_t segment, int pledge,
                         size_t maxio, size_t *clen)
{
    size_t cap = lzr_compress_bound(n), got = 0;
    uint8_t *dst = malloc(cap);
    lzr_cstream *cs = lzr_cstream_create();
    CHECK_RC(lzr_cstream_set(cs, LZR_C_LEVEL, level), LZR_OK);
    CHECK_RC(lzr_cstream_set(cs, LZR_C_THREADS, threads), LZR_OK);
    CHECK_RC(lzr_cstream_set(cs, LZR_C_SEGMENT, segment), LZR_OK);
    if (pledge) CHECK_RC(lzr_cstream_set(cs, LZR_C_CONTENT_SIZE, (int64_t)n), LZR_OK);
    size_t ip = 0;
    for (;;) {
        size_t ik = 1 + rnd() % maxio, ok = 1 + rnd() % maxio;
        if (ik > n - ip) ik = n - ip;
        if (ok > cap - got) ok = cap - got;
        lzr_inbuf in = {src + ip, ik, 0};
        lzr_outbuf out = {dst + got, ok, 0};
        int last = ip + ik == n;
        int r = lzr_compress_stream(cs, &out, &in, last ? LZR_FINISH : LZR_CONTINUE);
        ip += in.pos;
        got += out.pos;
        if (r < 0) { CHECK_RC(r, 0); break; }
        if (last && in.pos == in.size && r == 0) break;
    }
    lzr_cstream_free(cs);
    *clen = got;
    return dst;
}

/* Streaming decompress; returns 1 if it ended with LZR_FRAME_END for every
 * frame and the whole input was used. */
static int stream_decode(const uint8_t *src, size_t n, uint8_t *dst, size_t cap, size_t maxio, size_t *dlen)
{
    lzr_dstream *ds = lzr_dstream_create();
    size_t ip = 0, got = 0;
    int ok = 1, ended = 0;
    while (ip < n || !ended) {
        size_t ik = 1 + rnd() % maxio, oroom = 1 + rnd() % maxio;
        if (ik > n - ip) ik = n - ip;
        if (oroom > cap - got) oroom = cap - got;
        lzr_inbuf in = {src + ip, ik, 0};
        lzr_outbuf out = {dst + got, oroom, 0};
        int r = lzr_decompress_stream(ds, &out, &in);
        ip += in.pos;
        got += out.pos;
        if (r < 0) { ok = 0; break; }
        if (r == LZR_FRAME_END) {
            ended = 1;
            if (ip < n) { lzr_dstream_reset(ds); ended = 0; }
            continue;
        }
        if (ip == n && out.pos < out.size) { ok = 0; break; }   /* truncated */
        if (ip == n && !in.pos && !out.pos) { ok = 0; break; }  /* no progress: dst full */
    }
    lzr_dstream_free(ds);
    *dlen = got;
    return ok;
}

static void test_basics(void)
{
    CHECK(lzr_version_number() == LZR_LIB_VERSION_NUMBER);
    CHECK(strlen(lzr_version_string()) >= 5);
    for (int c = 0; c >= -8; c--) CHECK(strcmp(lzr_error_string(c), "unknown error") != 0);
    lzr_cstream *cs = lzr_cstream_create();
    CHECK_RC(lzr_cstream_set(cs, LZR_C_LEVEL, 12), LZR_ERR_PARAM);
    CHECK_RC(lzr_cstream_set(cs, LZR_C_LEVEL, -4), LZR_ERR_PARAM);
    CHECK_RC(lzr_cstream_set(cs, LZR_C_WINDOW_LOG, 17), LZR_ERR_PARAM);
    CHECK_RC(lzr_cstream_set(cs, LZR_C_WINDOW_LOG, 27), LZR_ERR_PARAM);
    CHECK_RC(lzr_cstream_set(cs, LZR_C_SEGMENT, 1000), LZR_ERR_PARAM);
    CHECK_RC(lzr_cstream_set(cs, LZR_C_THREADS, 17), LZR_ERR_PARAM);
    CHECK_RC(lzr_cstream_set(cs, LZR_C_FILTER, 2), LZR_ERR_PARAM);
    CHECK_RC(lzr_cstream_set(cs, (lzr_cparam)99, 0), LZR_ERR_PARAM);
    /* parameters are frozen once a frame started */
    uint8_t buf[256];
    lzr_inbuf in = {"abc", 3, 0};
    lzr_outbuf out = {buf, sizeof buf, 0};
    CHECK_RC(lzr_compress_stream(cs, &out, &in, LZR_CONTINUE), 0);
    CHECK_RC(lzr_cstream_set(cs, LZR_C_LEVEL, 1), LZR_ERR_PARAM);
    CHECK_RC(lzr_compress_stream(cs, &out, &in, LZR_FINISH), 0);
    CHECK_RC(lzr_cstream_reset(cs), LZR_OK);
    /* pledged size must match */
    CHECK_RC(lzr_cstream_set(cs, LZR_C_CONTENT_SIZE, 2), LZR_OK);
    in = (lzr_inbuf){"abc", 3, 0};
    out = (lzr_outbuf){buf, sizeof buf, 0};
    CHECK(lzr_compress_stream(cs, &out, &in, LZR_FINISH) == LZR_ERR_PARAM);
    CHECK(lzr_compress_stream(cs, &out, &in, LZR_FINISH) == LZR_ERR_PARAM);   /* stays failed */
    CHECK_RC(lzr_cstream_reset(cs), LZR_OK);
    CHECK_RC(lzr_cstream_set(cs, LZR_C_CONTENT_SIZE, 4), LZR_OK);
    in = (lzr_inbuf){"abc", 3, 0};
    out = (lzr_outbuf){buf, sizeof buf, 0};
    CHECK(lzr_compress_stream(cs, &out, &in, LZR_FINISH) == LZR_ERR_PARAM);
    lzr_cstream_free(cs);

    lzr_dstream *ds = lzr_dstream_create();
    CHECK_RC(lzr_dstream_set(ds, LZR_D_MAX_WINDOW_LOG, 17), LZR_ERR_PARAM);
    CHECK_RC(lzr_dstream_set(ds, LZR_D_MAX_WINDOW_LOG, 32), LZR_ERR_PARAM);
    CHECK_RC(lzr_dstream_set(ds, LZR_D_THREADS, 4), LZR_OK);
    lzr_dstream_free(ds);

    /* empty input */
    size_t clen = 0, dlen = 1;
    CHECK_RC(lzr_compress(buf, sizeof buf, &clen, NULL, 0, 6), LZR_OK);
    CHECK_RC(lzr_decompress(buf + 128, 64, &dlen, buf, clen), LZR_OK);
    CHECK(dlen == 0);
    lzr_frame_info fi;
    CHECK_RC(lzr_get_frame_info(buf, clen, &fi), LZR_OK);
    CHECK(fi.content_size == 0 && fi.version == LZR_FORMAT_VERSION && (fi.flags & LZR_FRAME_FLAG_SIZE));
    CHECK_RC(lzr_get_frame_info(buf, 3, &fi), LZR_ERR_TRUNCATED);
    CHECK_RC(lzr_get_frame_info("XYZW1234", 8, &fi), LZR_ERR_CORRUPT);
}

static void test_levels(void)
{
    const size_t n = 3 << 20;
    uint8_t *src = gen(n, 0), *dec = malloc(n + LZR_DECODE_MARGIN);
    for (int level = LZR_LEVEL_MIN; level <= LZR_LEVEL_MAX; level++) {
        size_t clen = 0, slen = 0, dlen = 0;
        uint8_t *a = oneshot(src, n, level, 1, -1, &clen);
        CHECK_RC(lzr_decompress(dec, n + LZR_DECODE_MARGIN, &dlen, a, clen), LZR_OK);
        CHECK(dlen == n && memcmp(dec, src, n) == 0);
        /* tight destination: decoded through a scratch buffer */
        memset(dec, 0, n);
        CHECK_RC(lzr_decompress(dec, n, &dlen, a, clen), LZR_OK);
        CHECK(dlen == n && memcmp(dec, src, n) == 0);
        CHECK_RC(lzr_decompress(dec, n - 1, &dlen, a, clen), LZR_ERR_DST_SIZE);
        /* the simple entry point is the same as compress_cs with defaults */
        uint8_t *b = malloc(lzr_compress_bound(n));
        CHECK_RC(lzr_compress(b, lzr_compress_bound(n), &slen, src, n, level), LZR_OK);
        CHECK(slen == clen && memcmp(a, b, clen) == 0);
        CHECK_RC(lzr_compress(b, clen - 1, &slen, src, n, level), LZR_ERR_DST_SIZE);
        free(b);
        /* streaming, pledged size: same bytes */
        uint8_t *s = streamed(src, n, level, 1, -1, 1, level >= 9 ? 300000 : 70000, &slen);
        CHECK(slen == clen && memcmp(s, a, clen) == 0);
        free(s);
        /* streaming decode with small pieces */
        CHECK(stream_decode(a, clen, dec, n, 5000, &dlen));
        CHECK(dlen == n && memcmp(dec, src, n) == 0);
        free(a);
        if (failures) fprintf(stderr, "  (level %d)\n", level);
    }
    free(src);
    free(dec);
}

static void test_stream_tiny_buffers(void)
{
    const size_t n = 600000;
    uint8_t *src = gen(n, 0), *dec = malloc(n);
    size_t clen, slen, dlen;
    uint8_t *a = oneshot(src, n, 6, 1, -1, &clen);
    uint8_t *s = streamed(src, n, 6, 1, -1, 1, 7, &slen);   /* 1..7 byte pieces */
    CHECK(slen == clen && memcmp(a, s, clen) == 0);
    free(s);
    /* unknown size: valid frame, no size in header */
    s = streamed(src, n, 6, 1, -1, 0, 4096, &slen);
    lzr_frame_info fi;
    CHECK_RC(lzr_get_frame_info(s, slen, &fi), LZR_OK);
    CHECK(fi.content_size == -1);
    CHECK(stream_decode(s, slen, dec, n, 3, &dlen));
    CHECK(dlen == n && memcmp(dec, src, n) == 0);
    free(s);
    free(a);
    free(src);
    free(dec);
}

static void test_threads_segments(void)
{
    const size_t n = 9 << 20;
    uint8_t *src = gen(n, 0), *dec = malloc(n + LZR_DECODE_MARGIN);
    struct { int level, threads; int64_t seg; } cfg[] = {{1, 3, -1}, {6, 4, -1}, {0, 2, -1}, {9, 4, 2 << 20}, {10, 1, 1 << 20}, {3, 2, 3 << 20}};
    for (size_t k = 0; k < sizeof cfg / sizeof cfg[0]; k++) {
        size_t clen, slen, dlen = 0, used = 0;
        uint8_t *a = oneshot(src, n, cfg[k].level, cfg[k].threads, cfg[k].seg, &clen);
        uint8_t *s = streamed(src, n, cfg[k].level, cfg[k].threads, cfg[k].seg, 1, 1 << 20, &slen);
        CHECK(slen == clen && memcmp(a, s, clen) == 0);
        free(s);
        lzr_dstream *ds = lzr_dstream_create();
        CHECK_RC(lzr_dstream_set(ds, LZR_D_THREADS, 4), LZR_OK);
        CHECK_RC(lzr_decompress_ds(ds, dec, n + LZR_DECODE_MARGIN, &dlen, a, clen, &used), LZR_OK);
        CHECK(dlen == n && used == clen && memcmp(dec, src, n) == 0);
        lzr_dstream_free(ds);
        free(a);
        if (failures) fprintf(stderr, "  (config %zu)\n", k);
    }
    free(src);
    free(dec);
}

static void test_bound(void)
{
    const size_t sizes[] = {0, 1, 100, 262144, 262145, 3000000, 5 << 20};
    for (size_t k = 0; k < sizeof sizes / sizeof sizes[0]; k++) {
        size_t n = sizes[k], clen;
        uint8_t *src = gen(n, 1);
        int64_t segs[] = {0, 1 << 20, (1 << 20) + 1};
        for (int j = 0; j < 3; j++) {
            for (int level = -3; level <= 11; level += 7) {
                uint8_t *a = oneshot(src, n, level, 2, segs[j], &clen);
                CHECK(clen <= lzr_compress_bound(n));
                free(a);
            }
        }
        free(src);
    }
}

static void test_errors_and_frames(void)
{
    const size_t n = 400000;
    uint8_t *src = gen(n, 0), *dec = malloc(2 * n + 64);
    size_t alen, blen, dlen, used;
    uint8_t *a = oneshot(src, n, 11, 1, -1, &alen);
    uint8_t *b = oneshot(src, n / 2, 1, 1, -1, &blen);

    /* truncated */
    CHECK_RC(lzr_decompress(dec, n + 64, &dlen, a, alen - 1), LZR_ERR_TRUNCATED);
    CHECK_RC(lzr_decompress(dec, n + 64, &dlen, a, 10), LZR_ERR_TRUNCATED);
    CHECK(!stream_decode(a, alen - 1, dec, n, 1000, &dlen));
    /* corrupt */
    a[alen / 2] ^= 0x10;
    int r = lzr_decompress(dec, n + 64, &dlen, a, alen);
    CHECK(r == LZR_ERR_CORRUPT || r == LZR_ERR_CHECKSUM);
    a[alen / 2] ^= 0x10;
    a[alen - 1] ^= 1;   /* checksum byte */
    CHECK_RC(lzr_decompress(dec, n + 64, &dlen, a, alen), LZR_ERR_CHECKSUM);
    a[alen - 1] ^= 1;
    /* not a frame */
    CHECK_RC(lzr_decompress(dec, n, &dlen, src, 100), LZR_ERR_CORRUPT);

    /* concatenated frames: one-shot reads one at a time, streaming all */
    uint8_t *cat = malloc(alen + blen + 5);
    memcpy(cat, a, alen);
    memcpy(cat + alen, b, blen);
    CHECK_RC(lzr_decompress(dec, 2 * n, &dlen, cat, alen + blen), LZR_ERR_CORRUPT);   /* trailing */
    lzr_dstream *ds = lzr_dstream_create();
    CHECK_RC(lzr_decompress_ds(ds, dec, 2 * n, &dlen, cat, alen + blen, &used), LZR_OK);
    CHECK(used == alen && dlen == n);
    CHECK_RC(lzr_decompress_ds(ds, dec, 2 * n, &dlen, cat + used, blen, &used), LZR_OK);
    CHECK(used == blen && dlen == n / 2 && memcmp(dec, src, n / 2) == 0);
    CHECK(stream_decode(cat, alen + blen, dec, 2 * n, 9000, &dlen));
    CHECK(dlen == n + n / 2 && memcmp(dec, src, n) == 0 && memcmp(dec + n, src, n / 2) == 0);
    /* streaming stops at the frame end and leaves the rest unread */
    lzr_dstream_reset(ds);
    memcpy(cat + alen, "junk!", 5);
    lzr_inbuf in = {cat, alen + 5, 0};
    lzr_outbuf out = {dec, 2 * n, 0};
    int rr;
    do rr = lzr_decompress_stream(ds, &out, &in); while (rr == 0 && in.pos < in.size);
    CHECK_RC(rr, LZR_FRAME_END);
    CHECK(in.pos == alen && out.pos == n);
    /* after the end, more calls keep reporting the end without reading */
    CHECK_RC(lzr_decompress_stream(ds, &out, &in), LZR_FRAME_END);
    CHECK(in.pos == alen);
    /* an error sticks until reset */
    lzr_dstream_reset(ds);
    in = (lzr_inbuf){"XYZW12345678", 12, 0};
    out = (lzr_outbuf){dec, 100, 0};
    CHECK_RC(lzr_decompress_stream(ds, &out, &in), LZR_ERR_CORRUPT);
    CHECK_RC(lzr_decompress_stream(ds, &out, &in), LZR_ERR_PARAM);
    lzr_dstream_reset(ds);

    /* memory limit */
    lzr_frame_info fi;
    CHECK_RC(lzr_get_frame_info(a, alen, &fi), LZR_OK);
    CHECK(fi.window_log > LZR_WINDOW_LOG_MIN);
    CHECK_RC(lzr_dstream_set(ds, LZR_D_MAX_WINDOW_LOG, LZR_WINDOW_LOG_MIN), LZR_OK);
    CHECK_RC(lzr_decompress_ds(ds, dec, 2 * n, &dlen, a, alen, &used), LZR_ERR_MEMLIMIT);
    in = (lzr_inbuf){a, alen, 0};
    out = (lzr_outbuf){dec, 2 * n, 0};
    CHECK_RC(lzr_decompress_stream(ds, &out, &in), LZR_ERR_MEMLIMIT);
    lzr_dstream_free(ds);

    free(cat);
    free(a);
    free(b);
    free(src);
    free(dec);
}

int main(void)
{
    test_basics();
    test_levels();
    test_stream_tiny_buffers();
    test_threads_segments();
    test_bound();
    test_errors_and_frames();
    if (failures) {
        fprintf(stderr, "test_api: %d check(s) failed\n", failures);
        return 1;
    }
    printf("test_api: all checks passed\n");
    return 0;
}
