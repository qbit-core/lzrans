/*
 * enc_internal.h - encoder state shared by the lazy parser (enc.c) and
 * the binary-tree / optimal parser (opt.c).
 */
#ifndef LZRANS_ENC_INTERNAL_H
#define LZRANS_ENC_INTERNAL_H
#include "lzr.h"
#include "codes.h"
#include "arc.h"

typedef struct {
    uint8_t wlog, hlog, lazy;
    uint16_t chain;   /* hash-chain steps (lazy) or tree depth (opt) per search */
    uint16_t nice;    /* stop searching once a match this long is found          */
    uint8_t opt;      /* 0 = lazy parser, N = optimal parser with N passes        */
} level_params;


struct lzr_cctx {
    level_params lp;
    size_t wsize, wmask;
    uint8_t *buf;          /* input history + pending data              */
    size_t cap;            /* 2*wsize + block                            */
    size_t pos;            /* first byte not yet compressed              */
    size_t end;            /* bytes filled                               */
    size_t ins;            /* next position to insert into the chains    */
    uint32_t *head, *prev; /* stored as buffer index + 1, 0 = empty      */
    uint32_t reps[LZR_NUM_REPS];
    uint32_t crc;
    uint64_t total_in;
    int64_t content_size;
    int header_done, finished;

    /* optimal parser: binary tree (son[2*wsize]), match cache, prices */
    uint32_t *son;
    void *opt;
    int have_stats;
    int threads;           /* match-finder threads (levels 7-9)           */
    size_t bt_maxdist;     /* farthest tree node followed (see opt_init)   */
    int adaptive;          /* emit block type 3 (adaptive range coding)    */
    arc_model *am;         /* [0] = live model, [1] = copy for rollback     */

    /* per-block scratch */
    uint32_t *seq_ll, *seq_ml, *seq_dist;  /* dist: 0..2 = rep idx, else dist + 2 */
    uint8_t *lits, *llc, *mlc, *ofc, *extra, *scratch, *blk;
};

#define MAX_SEQ (LZR_BLOCK_SIZE / LZR_ML_BASE + 1)
#define LZR_MAX_THREADS 16
/* Adaptive literal stream worst case: 8 decisions per literal, each
 * <= 7.05 bits (probabilities stay in [31, 4065] / 4096). */
#define ARC_CAP_LIT (8 * (size_t)LZR_BLOCK_SIZE + 16)
#ifndef LZR_ADAPTIVE_LEVEL
#define LZR_ADAPTIVE_LEVEL 6   /* levels >= this use adaptive blocks */
#endif
#ifndef LZR_HAVE_THREADS
#if defined(_WIN32)
#define LZR_HAVE_THREADS 0
#else
#define LZR_HAVE_THREADS 1
#endif
#endif
#define LZR_LOOKAHEAD 512   /* > max nice length; see opt.c gather() */

#ifndef LZR_HASH_LEN
#define LZR_HASH_LEN 5
#endif
/* Hash of the next LZR_HASH_LEN bytes (reads 8; the extra bytes are shifted
 * out, and the buffer has slack so the read stays in bounds). Hashing more
 * than MIN_MATCH bytes keeps chains free of short, useless candidates;
 * 4-byte matches are still found through the repeat offsets. */
static inline uint32_t hash4(const uint8_t *p, unsigned hlog)
{
    uint64_t v;
    memcpy(&v, p, 8);
    return (uint32_t)(((v << (64 - 8 * LZR_HASH_LEN)) * 0xCF1BBCDCB7A56463ull) >> (64 - hlog));
}

static inline size_t match_len(const uint8_t *a, const uint8_t *b, size_t maxlen)
{
    size_t len = 0;
    while (len + 8 <= maxlen) {
        uint64_t x, y;
        memcpy(&x, a + len, 8);
        memcpy(&y, b + len, 8);
        if (x != y) return len + ((size_t)__builtin_ctzll(x ^ y) >> 3);
        len += 8;
    }
    while (len < maxlen && a[len] == b[len]) len++;
    return len;
}

/* Map a chosen match to its offset code and update the repeat history. */
static inline uint32_t commit_offset(uint32_t *reps, uint32_t off)
{
    uint32_t dist = off < LZR_NUM_REPS ? reps[off] : off - 2;
    if (dist == reps[0]) return 0;
    if (dist == reps[1]) { reps[1] = reps[0]; reps[0] = dist; return 1; }
    if (dist == reps[2]) { reps[2] = reps[1]; reps[1] = reps[0]; reps[0] = dist; return 2; }
    reps[2] = reps[1]; reps[1] = reps[0]; reps[0] = dist;
    return dist + 2;
}


/* opt.c */
int opt_init(lzr_cctx *c);
void opt_free(lzr_cctx *c);
size_t opt_parse_block(lzr_cctx *c, size_t bs, size_t be, size_t *nlit_out);
void opt_slide(lzr_cctx *c, uint32_t delta);
void opt_sync(lzr_cctx *c);

#endif
