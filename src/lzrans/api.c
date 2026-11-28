/*
 * api.c - liblzr public API (include/lzr.h) on top of the core (core.h).
 *
 * The core appends its output to growable buffers (lzr_buf) and consumes
 * whatever input it is given. The stream wrappers here add the zstd-style
 * in/out buffer interface: core output goes to a pending buffer that is
 * drained into the caller's output, and input is fed in bounded pieces only
 * while nothing is pending, so memory stays bounded by the window and one
 * piece's worth of output. The core's output does not depend on how its
 * input is split, so the bytes match a one-shot compress of the same data.
 */
#include "core.h"
#include <unistd.h>

#define CS_FEED ((size_t)1 << 20)   /* compressor input piece */
#define DS_FEED ((size_t)1 << 16)   /* decompressor input piece */

unsigned lzr_version_number(void) { return LZR_LIB_VERSION_NUMBER; }

#define LZR_STR2(x) #x
#define LZR_STR(x) LZR_STR2(x)
const char *lzr_version_string(void)
{
    return LZR_STR(LZR_LIB_VERSION_MAJOR) "." LZR_STR(LZR_LIB_VERSION_MINOR) "." LZR_STR(LZR_LIB_VERSION_PATCH);
}

const char *lzr_error_string(int code)
{
    switch (code) {
    case LZR_OK: return "no error";
    case LZR_ERR_NOMEM: return "out of memory";
    case LZR_ERR_CORRUPT: return "corrupt lzr data";
    case LZR_ERR_CHECKSUM: return "checksum mismatch: decompressed data is corrupt";
    case LZR_ERR_UNSUPPORTED: return "unsupported frame version, flags or block type (needs a newer decoder)";
    case LZR_ERR_PARAM: return "invalid parameter or API misuse";
    case LZR_ERR_MEMLIMIT: return "frame window exceeds the decoder memory limit (raise max window log, up to 31)";
    case LZR_ERR_DST_SIZE: return "destination buffer too small";
    case LZR_ERR_TRUNCATED: return "unexpected end of input (truncated lzr frame)";
    default: return "unknown error";
    }
}

/* 0 = one per online CPU, capped at LZR_THREADS_MAX. */
static int resolve_threads(int64_t t)
{
    if (t > 0) return (int)t;
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    return n < 1 ? 1 : (n > LZR_THREADS_MAX ? LZR_THREADS_MAX : (int)n);
}

/* Copy as much of the pending buffer into out as fits; 1 if some remains. */
static int drain(lzr_buf *pend, size_t *ppos, lzr_outbuf *out)
{
    size_t have = pend->len - *ppos, room = out->size - out->pos;
    size_t k = have < room ? have : room;
    if (k) {
        memcpy((uint8_t *)out->dst + out->pos, pend->p + *ppos, k);
        out->pos += k;
        *ppos += k;
    }
    if (*ppos < pend->len) return 1;
    pend->len = 0;
    *ppos = 0;
    return 0;
}

/* ======================= compression ======================= */

enum { ST_IDLE, ST_RUN, ST_DONE, ST_ERR };

struct lzr_cstream {
    int level, wlog, threads, filter;
    int64_t segment, csize;
    lzr_cctx *c;
    lzr_buf pend;
    size_t ppos;
    uint64_t fed;
    int state;
};

lzr_cstream *lzr_cstream_create(void)
{
    lzr_cstream *cs = (lzr_cstream *)calloc(1, sizeof(*cs));
    if (!cs) return NULL;
    cs->level = LZR_LEVEL_DEFAULT;
    cs->threads = 1;
    cs->segment = -1;
    cs->csize = -1;
    return cs;
}

void lzr_cstream_free(lzr_cstream *cs)
{
    if (!cs) return;
    lzr_cctx_free(cs->c);
    free(cs->pend.p);
    free(cs);
}

int lzr_cstream_reset(lzr_cstream *cs)
{
    if (!cs) return LZR_ERR_PARAM;
    lzr_cctx_free(cs->c);
    cs->c = NULL;
    cs->pend.len = 0;
    cs->ppos = 0;
    cs->fed = 0;
    cs->csize = -1;
    cs->state = ST_IDLE;
    return LZR_OK;
}

int lzr_cstream_set(lzr_cstream *cs, lzr_cparam param, int64_t v)
{
    if (!cs || cs->state != ST_IDLE) return LZR_ERR_PARAM;
    switch (param) {
    case LZR_C_LEVEL:
        if (v < LZR_LEVEL_MIN || v > LZR_LEVEL_MAX) return LZR_ERR_PARAM;
        cs->level = (int)v;
        return LZR_OK;
    case LZR_C_WINDOW_LOG:
        if (v && (v < LZR_WINDOW_LOG_MIN || v > LZR_WINDOW_LOG_ENC_MAX)) return LZR_ERR_PARAM;
        cs->wlog = (int)v;
        return LZR_OK;
    case LZR_C_THREADS:
        if (v < 0 || v > LZR_THREADS_MAX) return LZR_ERR_PARAM;
        cs->threads = (int)v;
        return LZR_OK;
    case LZR_C_SEGMENT:
        if (v < -1 || (v > 0 && (uint64_t)v < LZR_SEGMENT_MIN)) return LZR_ERR_PARAM;
        cs->segment = v;
        return LZR_OK;
    case LZR_C_FILTER:
        if (v != LZR_FILTER_NONE && v != LZR_FILTER_X86) return LZR_ERR_PARAM;
        cs->filter = (int)v;
        return LZR_OK;
    case LZR_C_CONTENT_SIZE:
        if (v < -1) return LZR_ERR_PARAM;
        cs->csize = v;
        return LZR_OK;
    }
    return LZR_ERR_PARAM;
}

static uint64_t cs_segment(const lzr_cstream *cs)
{
    return cs->segment < 0 ? lzr_default_segment(cs->level) : (uint64_t)cs->segment;
}

static int cs_start(lzr_cstream *cs)
{
    cs->c = lzr_cctx_new_seg(cs->level, cs->wlog, cs->csize, cs_segment(cs));
    if (!cs->c) return LZR_ERR_NOMEM;   /* parameters were validated by lzr_cstream_set */
    lzr_cctx_set_threads(cs->c, resolve_threads(cs->threads));
    int r = lzr_cctx_set_filter(cs->c, cs->filter);
    if (r) return r;
    cs->state = ST_RUN;
    return LZR_OK;
}

/* Whole-input feed: levels -3..6 with threads read a one-shot input in
 * place (core.h: jobs_borrow) when it arrives in a single write that
 * matches the pledged size. Elsewhere it would only buffer more output. */
static int cs_whole(const lzr_cstream *cs, size_t avail)
{
    return cs->fed == 0 && cs->csize >= 0 && (uint64_t)cs->csize == avail && avail < ((size_t)1 << 30) &&
           cs->level <= 6 && resolve_threads(cs->threads) > 1 && cs->filter == LZR_FILTER_NONE &&
           (cs_segment(cs) == 0 || (uint64_t)cs->csize <= cs_segment(cs));
}

int lzr_compress_stream(lzr_cstream *cs, lzr_outbuf *out, lzr_inbuf *in, lzr_end_op op)
{
    if (!cs || !out || !in || in->pos > in->size || out->pos > out->size || cs->state == ST_ERR)
        return LZR_ERR_PARAM;
    if (cs->state == ST_DONE) {
        if (drain(&cs->pend, &cs->ppos, out)) return 1;
        return in->pos < in->size ? LZR_ERR_PARAM : 0;
    }
    int r;
    if (cs->state == ST_IDLE && (r = cs_start(cs))) goto fail;
    for (;;) {
        if (drain(&cs->pend, &cs->ppos, out)) return 1;
        size_t avail = in->size - in->pos;
        if (avail) {
            size_t k = cs_whole(cs, avail) ? avail : (avail < CS_FEED ? avail : CS_FEED);
            if (cs->csize >= 0 && cs->fed + k > (uint64_t)cs->csize) { r = LZR_ERR_PARAM; goto fail; }
            r = lzr_cctx_write(cs->c, (const uint8_t *)in->src + in->pos, k, &cs->pend);
            if (r) goto fail;
            in->pos += k;
            cs->fed += k;
            continue;
        }
        if (op != LZR_FINISH) return 0;
        r = lzr_cctx_end(cs->c, &cs->pend);
        if (r) goto fail;
        cs->state = ST_DONE;
        lzr_cctx_free(cs->c);
        cs->c = NULL;
        return drain(&cs->pend, &cs->ppos, out);
    }
fail:
    cs->state = ST_ERR;
    return r;
}

size_t lzr_compress_bound(size_t n)
{
    /* stored blocks (9-byte header per <= 256 KiB, plus one per segment
     * end) + seek index (16 bytes per >= 1 MiB segment) + frame header,
     * end block and index trailer */
    return n + n / 8192 + 256;
}

int lzr_compress_cs(lzr_cstream *cs, void *dst, size_t dst_cap, size_t *dst_len, const void *src, size_t n)
{
    if (!cs || (!dst && dst_cap) || (!src && n)) return LZR_ERR_PARAM;
    lzr_cstream_reset(cs);
    cs->csize = (int64_t)n;
    lzr_buf out = {0};
    int r = cs_start(cs);
    if (!r) r = buf_reserve(&out, n / 2 + 64);
    if (!r) r = lzr_cctx_write(cs->c, (const uint8_t *)src, n, &out);
    if (!r) r = lzr_cctx_end(cs->c, &out);
    if (!r && out.len > dst_cap) r = LZR_ERR_DST_SIZE;
    if (!r) {
        if (out.len) memcpy(dst, out.p, out.len);
        if (dst_len) *dst_len = out.len;
    }
    free(out.p);
    lzr_cstream_reset(cs);
    return r;
}

int lzr_compress(void *dst, size_t dst_cap, size_t *dst_len, const void *src, size_t n, int level)
{
    lzr_cstream *cs = lzr_cstream_create();
    if (!cs) return LZR_ERR_NOMEM;
    int r = lzr_cstream_set(cs, LZR_C_LEVEL, level);
    if (!r) r = lzr_compress_cs(cs, dst, dst_cap, dst_len, src, n);
    lzr_cstream_free(cs);
    return r;
}

/* ======================= decompression ======================= */

struct lzr_dstream {
    int maxw, threads;
    lzr_dctx *d;
    lzr_buf pend;
    size_t ppos;
    int state;
};

lzr_dstream *lzr_dstream_create(void)
{
    lzr_dstream *ds = (lzr_dstream *)calloc(1, sizeof(*ds));
    if (!ds) return NULL;
    ds->maxw = LZR_MAX_WINDOW_LOG_DEFAULT;
    ds->threads = 1;
    return ds;
}

void lzr_dstream_free(lzr_dstream *ds)
{
    if (!ds) return;
    lzr_dctx_free(ds->d);
    free(ds->pend.p);
    free(ds);
}

int lzr_dstream_reset(lzr_dstream *ds)
{
    if (!ds) return LZR_ERR_PARAM;
    lzr_dctx_free(ds->d);
    ds->d = NULL;
    ds->pend.len = 0;
    ds->ppos = 0;
    ds->state = ST_IDLE;
    return LZR_OK;
}

int lzr_dstream_set(lzr_dstream *ds, lzr_dparam param, int64_t v)
{
    if (!ds || ds->state != ST_IDLE) return LZR_ERR_PARAM;
    switch (param) {
    case LZR_D_MAX_WINDOW_LOG:
        if (v < LZR_WINDOW_LOG_MIN || v > LZR_WINDOW_LOG_MAX) return LZR_ERR_PARAM;
        ds->maxw = (int)v;
        return LZR_OK;
    case LZR_D_THREADS:
        if (v < 0 || v > LZR_THREADS_MAX) return LZR_ERR_PARAM;
        ds->threads = (int)v;
        return LZR_OK;
    }
    return LZR_ERR_PARAM;
}

int lzr_decompress_stream(lzr_dstream *ds, lzr_outbuf *out, lzr_inbuf *in)
{
    if (!ds || !out || !in || in->pos > in->size || out->pos > out->size || ds->state == ST_ERR)
        return LZR_ERR_PARAM;
    int r;
    if (ds->state == ST_IDLE) {
        if (!(ds->d = lzr_dctx_new())) { r = LZR_ERR_NOMEM; goto fail; }
        if ((r = lzr_dctx_set_max_window_log(ds->d, ds->maxw))) goto fail;
        ds->state = ST_RUN;
    }
    for (;;) {
        if (drain(&ds->pend, &ds->ppos, out)) return 0;
        if (ds->state == ST_DONE) return LZR_FRAME_END;
        size_t avail = in->size - in->pos;
        if (!avail) return 0;
        size_t k = avail < DS_FEED ? avail : DS_FEED, used = 0;
        r = lzr_dctx_write(ds->d, (const uint8_t *)in->src + in->pos, k, &ds->pend, &used);
        if (r < 0) goto fail;
        if (r == 1) {
            in->pos += used;
            ds->state = ST_DONE;
            lzr_dctx_free(ds->d);
            ds->d = NULL;
        } else {
            in->pos += k;
        }
    }
fail:
    ds->state = ST_ERR;
    return r;
}

int lzr_decompress_ds(lzr_dstream *ds, void *dst, size_t dst_cap, size_t *dst_len,
                      const void *src, size_t n, size_t *src_used)
{
    if (!ds || (!dst && dst_cap) || (!src && n)) return LZR_ERR_PARAM;
    int64_t csize = -1;
    int hr = lzr_frame_header((const uint8_t *)src, n, NULL, NULL, &csize);
    if (hr == 0) return LZR_ERR_TRUNCATED;
    if (hr < 0) return hr;
    if (csize >= 0 && (uint64_t)csize > dst_cap) return LZR_ERR_DST_SIZE;
    lzr_buf out = {0};
    if (csize >= 0 && dst_cap - (size_t)csize >= LZR_DECODE_MARGIN) {
        /* in place: the decoder never writes past content_size + LZR_SLACK */
        out.p = (uint8_t *)dst;
        out.cap = dst_cap;
        out.fixed = 1;
    }
    size_t used = 0;
    int r = lzr_decompress_mt((const uint8_t *)src, n, &out, resolve_threads(ds->threads), ds->maxw, &used);
    if (r == 0) r = LZR_ERR_TRUNCATED;
    else if (r == 1) {
        r = LZR_OK;
        if (!out.fixed) {
            if (out.len > dst_cap) r = LZR_ERR_DST_SIZE;
            else if (out.len) memcpy(dst, out.p, out.len);
        }
        if (!r) {
            if (dst_len) *dst_len = out.len;
            if (src_used) *src_used = used;
        }
    }
    if (!out.fixed) free(out.p);
    return r;
}

int lzr_decompress(void *dst, size_t dst_cap, size_t *dst_len, const void *src, size_t n)
{
    lzr_dstream *ds = lzr_dstream_create();
    if (!ds) return LZR_ERR_NOMEM;
    size_t used = 0;
    int r = lzr_decompress_ds(ds, dst, dst_cap, dst_len, src, n, &used);
    if (!r && used != n) r = LZR_ERR_CORRUPT;
    lzr_dstream_free(ds);
    return r;
}

int lzr_get_frame_info(const void *src, size_t n, lzr_frame_info *info)
{
    if (!src || !info) return LZR_ERR_PARAM;
    int wlog = 0;
    uint8_t flags = 0;
    int64_t size = -1;
    int r = lzr_frame_header((const uint8_t *)src, n, &wlog, &flags, &size);
    if (r == 0) return LZR_ERR_TRUNCATED;
    if (r < 0) return r;
    info->content_size = size;
    info->window_log = wlog;
    info->version = ((const uint8_t *)src)[4];
    info->flags = flags;
    info->header_size = (size_t)r;
    return LZR_OK;
}
