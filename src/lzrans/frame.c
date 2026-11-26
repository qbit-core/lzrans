/*
 * frame.c - container assembly: header, per-block LZ77 + rANS, CRC32.
 *
 * Compression:  for each 1 MB block -> lz_parse -> split tokens into
 *               3 symbol streams -> normalise freqs -> rANS encode each.
 *               Falls back to a stored block if that is not smaller.
 * Decompression: build O(1) slot tables, decode tokens in order, copy
 *               literals / matches sequentially, verify CRC32 at the end.
 */
#include "frame.h"
#include "lz77.h"
#include "rans.h"
#include "crc32.h"
#include <stdlib.h>

/* ------------------------------------------------------------------ */
/* Compression                                                         */
/* ------------------------------------------------------------------ */

size_t lzr_compress_bound(size_t n)
{
    size_t nblocks = (n + LZR_BLOCK_SIZE - 1) / LZR_BLOCK_SIZE;
    return LZR_HEADER_SIZE + nblocks * 5 + n;
}

typedef struct {
    lz_state lz;
    uint16_t sym[LZR_BLOCK_SIZE];
    uint16_t dist[LZR_BLOCK_SIZE / LZ_MIN_MATCH + 1];
    uint16_t dhi[LZR_BLOCK_SIZE / LZ_MIN_MATCH + 1];
    uint16_t dlo[LZR_BLOCK_SIZE / LZ_MIN_MATCH + 1];
    uint8_t scratch[2 * LZR_BLOCK_SIZE + 8];     /* rans_bound(max nsym) */
    uint8_t blk[13 + 3 * (RANS_FREQ_MAX(SYM_LITLEN) + 4) + 2 * LZR_BLOCK_SIZE + 8
                + 2 * (2 * (LZR_BLOCK_SIZE / LZ_MIN_MATCH + 1) + 8)];
} enc_ctx;

/* Encode one stream (freq table + len + payload) at `o`. Returns bytes. */
static size_t put_stream(enc_ctx *c, uint8_t *o, const uint16_t *syms, size_t n, unsigned alphabet)
{
    if (n == 0) return 0;
    uint32_t freq[SYM_LITLEN];
    rans_build_freqs(syms, n, alphabet, freq);
    size_t k = rans_write_freqs(o, freq, alphabet);
    size_t len = rans_encode(syms, n, freq, alphabet, c->scratch, sizeof(c->scratch), o + k + 4);
    wr32le(o + k, (uint32_t)len);
    return k + 4 + len;
}

int64_t lzr_compress(const uint8_t *src, size_t n, uint8_t *dst, size_t cap, int level)
{
    if (cap < lzr_compress_bound(n)) return LZR_ERR_NOMEM;
    enc_ctx *c = (enc_ctx *)malloc(sizeof(enc_ctx));
    if (!c) return LZR_ERR_NOMEM;
    lz_init(&c->lz, level);

    uint8_t *o = dst;
    memcpy(o, LZR_MAGIC, 4);
    o[4] = LZR_VERSION;
    o[5] = 0;                       /* flags: plain (no encryption) */
    o[6] = (uint8_t)level;
    o[7] = 0;
    wr64le(o + 8, (uint64_t)n);
    wr32le(o + 16, lzr_crc32(0, src, n));
    o += LZR_HEADER_SIZE;

    for (size_t start = 0; start < n; start += LZR_BLOCK_SIZE) {
        size_t end = start + LZR_BLOCK_SIZE < n ? start + LZR_BLOCK_SIZE : n;
        size_t raw = end - start;
        size_t nm;
        size_t ns = lz_parse(&c->lz, src, n, start, end, c->sym, c->dist, &nm);
        for (size_t i = 0; i < nm; i++) {
            c->dhi[i] = c->dist[i] >> 8;
            c->dlo[i] = c->dist[i] & 0xFF;
        }

        uint8_t *b = c->blk;
        b[0] = 1;
        wr32le(b + 1, (uint32_t)raw);
        wr32le(b + 5, (uint32_t)ns);
        wr32le(b + 9, (uint32_t)nm);
        size_t k = 13;
        k += put_stream(c, b + k, c->sym, ns, SYM_LITLEN);
        k += put_stream(c, b + k, c->dhi, nm, SYM_BYTE);
        k += put_stream(c, b + k, c->dlo, nm, SYM_BYTE);

        if (k < 5 + raw) {
            memcpy(o, b, k);
            o += k;
        } else {                    /* incompressible: store */
            o[0] = 0;
            wr32le(o + 1, (uint32_t)raw);
            memcpy(o + 5, src + start, raw);
            o += 5 + raw;
        }
    }
    free(c);
    return (int64_t)(o - dst);
}

/* ------------------------------------------------------------------ */
/* Decompression                                                       */
/* ------------------------------------------------------------------ */

typedef struct {
    uint32_t x;
    const uint8_t *p, *end;
    int err;
} rdec;

static int rdec_init(rdec *d, const uint8_t *p, size_t len)
{
    if (len < 4) return 0;
    d->x = rd32le(p);
    d->p = p + 4;
    d->end = p + len;
    d->err = 0;
    return d->x >= RANS_L;
}

/* O(1) decode: one table lookup gives symbol, freq and slot bias. */
static inline unsigned rdec_get(rdec *d, const rans_slot *tbl)
{
    const rans_slot e = tbl[d->x & (RANS_M - 1)];
    uint32_t x = e.freq * (d->x >> RANS_SCALE_BITS) + e.bias;
    while (x < RANS_L) {
        if (d->p >= d->end) { d->err = 1; break; }
        x = (x << 8) | *d->p++;
    }
    d->x = x;
    return e.sym;
}

/* Parse "freq table + u32 len + payload". Returns bytes consumed or 0. */
static size_t get_stream(const uint8_t *in, size_t avail, unsigned alphabet,
                         rans_slot *tbl, rdec *d)
{
    size_t k = rans_read_freqs(in, avail, alphabet, tbl);
    if (!k || avail - k < 4) return 0;
    uint32_t len = rd32le(in + k);
    k += 4;
    if (len > avail - k) return 0;
    if (!rdec_init(d, in + k, len)) return 0;
    return k + len;
}

int lzr_frame_info(const uint8_t *src, size_t n, uint64_t *orig_size, uint32_t *crc, uint8_t *flags)
{
    if (n < LZR_HEADER_SIZE || memcmp(src, LZR_MAGIC, 4) != 0) return LZR_ERR_CORRUPT;
    if (src[4] != LZR_VERSION) return LZR_ERR_UNSUPPORTED;
    if (flags) *flags = src[5];
    if (src[5] & ~LZR_KNOWN_FLAGS) return LZR_ERR_UNSUPPORTED;
    uint64_t sz = rd64le(src + 8);
    /* Each block needs >= 5 header bytes and yields <= BLOCK_SIZE bytes. */
    uint64_t max_blocks = (n - LZR_HEADER_SIZE) / 5;
    if (sz > max_blocks * (uint64_t)LZR_BLOCK_SIZE) return LZR_ERR_CORRUPT;
    *orig_size = sz;
    if (crc) *crc = rd32le(src + 16);
    return LZR_OK;
}

typedef struct {
    rans_slot t1[RANS_M], t2[RANS_M], t3[RANS_M];
} dec_tables;

static int decode_block(const uint8_t *in, size_t avail, size_t *consumed,
                        uint8_t *base, uint8_t *op, uint8_t *bend, dec_tables *T)
{
    if (avail < 8) return LZR_ERR_CORRUPT;
    size_t ns = rd32le(in), nm = rd32le(in + 4);
    size_t k = 8;
    if (ns == 0 || ns > (size_t)(bend - op) || nm > ns) return LZR_ERR_CORRUPT;

    rdec d1, d2, d3;
    size_t c = get_stream(in + k, avail - k, SYM_LITLEN, T->t1, &d1);
    if (!c) return LZR_ERR_CORRUPT;
    k += c;
    if (nm) {
        c = get_stream(in + k, avail - k, SYM_BYTE, T->t2, &d2);
        if (!c) return LZR_ERR_CORRUPT;
        k += c;
        c = get_stream(in + k, avail - k, SYM_BYTE, T->t3, &d3);
        if (!c) return LZR_ERR_CORRUPT;
        k += c;
    }

    size_t matches = 0;
    for (size_t i = 0; i < ns; i++) {
        unsigned s = rdec_get(&d1, T->t1);
        if (s < SYM_BYTE) {
            if (op >= bend) return LZR_ERR_CORRUPT;
            *op++ = (uint8_t)s;
            continue;
        }
        if (++matches > nm) return LZR_ERR_CORRUPT;
        size_t len = s - SYM_BYTE + LZ_MIN_MATCH;
        size_t dist = ((size_t)rdec_get(&d2, T->t2) << 8 | rdec_get(&d3, T->t3)) + 1;
        if (dist > (size_t)(op - base) || len > (size_t)(bend - op)) return LZR_ERR_CORRUPT;

        const uint8_t *from = op - dist;
        if (dist >= len) {
            memcpy(op, from, len);           /* non-overlapping: one memcpy */
        } else {
            /* Overlapping run: copy the period, then keep doubling the
             * already-written span so each memcpy is non-overlapping. */
            uint8_t *dst = op;
            size_t left = len;
            while (left) {
                size_t chunk = (size_t)(dst - from);
                if (chunk > left) chunk = left;
                memcpy(dst, from, chunk);
                dst += chunk;
                left -= chunk;
            }
        }
        op += len;
    }

    /* Integrity: every stream must be fully consumed and back at state L. */
    if (op != bend || matches != nm) return LZR_ERR_CORRUPT;
    if (d1.err || d1.x != RANS_L || d1.p != d1.end) return LZR_ERR_CORRUPT;
    if (nm && (d2.err || d3.err || d2.x != RANS_L || d3.x != RANS_L ||
               d2.p != d2.end || d3.p != d3.end)) return LZR_ERR_CORRUPT;
    *consumed = k;
    return LZR_OK;
}

int lzr_decompress(const uint8_t *src, size_t n, uint8_t *dst, size_t cap)
{
    uint64_t orig;
    uint32_t crc;
    int r = lzr_frame_info(src, n, &orig, &crc, NULL);
    if (r) return r;
    if (orig != cap) return LZR_ERR_CORRUPT;

    dec_tables *T = NULL;
    const uint8_t *ip = src + LZR_HEADER_SIZE, *iend = src + n;
    uint8_t *op = dst, *oend = dst + cap;

    while (op < oend) {
        if (iend - ip < 5) { r = LZR_ERR_CORRUPT; goto done; }
        uint8_t type = ip[0];
        size_t raw = rd32le(ip + 1);
        ip += 5;
        if (raw == 0 || raw > LZR_BLOCK_SIZE || raw > (size_t)(oend - op)) { r = LZR_ERR_CORRUPT; goto done; }
        if (type == 0) {
            if ((size_t)(iend - ip) < raw) { r = LZR_ERR_CORRUPT; goto done; }
            memcpy(op, ip, raw);
            ip += raw;
        } else if (type == 1) {
            if (!T && !(T = (dec_tables *)malloc(sizeof(dec_tables)))) { r = LZR_ERR_NOMEM; goto done; }
            size_t used;
            r = decode_block(ip, (size_t)(iend - ip), &used, dst, op, op + raw, T);
            if (r) goto done;
            ip += used;
        } else {
            r = LZR_ERR_CORRUPT; goto done;
        }
        op += raw;
    }
    if (ip != iend) { r = LZR_ERR_CORRUPT; goto done; }
    r = lzr_crc32(0, dst, cap) == crc ? LZR_OK : LZR_ERR_CHECKSUM;
done:
    free(T);
    return r;
}
