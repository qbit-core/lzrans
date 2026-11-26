/*
 * enc_internal.h - encoder state shared by the lazy parser (enc.c) and
 * the binary-tree / optimal parser (opt.c).
 */
#ifndef LZRANS_ENC_INTERNAL_H
#define LZRANS_ENC_INTERNAL_H
#include "lzr.h"
#include "codes.h"
#include "arc.h"
#include "x86.h"
#include "huf.h"
#include "fse.h"
#include "xxh64.h"

typedef struct {
    uint8_t wlog, hlog, lazy;
    uint16_t chain;   /* hash-chain steps (lazy) or tree depth (opt) per search */
    uint16_t nice;    /* stop searching once a match this long is found          */
    uint8_t opt;      /* 0 = lazy parser, N = optimal parser with N passes        */
    uint8_t row;      /* lazy levels: log2 of the row size (row match finder)    */
} level_params;


struct lzr_cctx {
    level_params lp;
    size_t wsize, wmask;
    uint8_t *buf;          /* input history + pending data              */
    size_t cap;            /* 2*wsize + block                            */
    size_t pos;            /* first byte not yet compressed              */
    size_t end;            /* bytes filled                               */
    size_t ins;            /* next position to insert into the chains    */
    uint32_t *head;        /* levels 7-9: buffer index + 1; levels 1-6 (row
                            * slots): base + buffer index + 1 (see slide()); 0 = empty */
    size_t base;           /* levels 1-6: bytes slid out since the last rebase */
    uint8_t *rtag, *rhead; /* row finder: 8-bit tags (head[] holds the
                            * positions) and each row's newest slot       */
    uint32_t reps[LZR_NUM_REPS];
    xxh64_state xs;        /* checksum of the original content */
    uint64_t total_in;
    int64_t content_size;
    int header_done, finished;

    /* optimal parser: binary tree (son[2*wsize]), match cache, prices */
    uint32_t *son;
    void *opt;
    int have_stats;
    int threads;           /* match-finder threads (levels 7-9)           */
    size_t bt_maxdist;     /* farthest tree node followed (see opt_init)   */
    int adaptive;          /* emit block type 4 (token range coding)        */
    tok_model *tm;         /* [0] = live model, [1] = copy for rollback     */
    /* block type 5: tables the decoder holds (repeat modes), and copies
     * for the stored-block rollback */
    huf_table *huf;        /* [0] live, [1] rollback                        */
    fse_table *fse;        /* [0..2] ll/ml/of live, [3..5] rollback         */
    lzr_x86 x86;           /* x86 filter state (if x86tmp)                 */
    uint8_t *x86tmp;       /* x86 filter work buffer, NULL = no filter     */

    /* per-block scratch */
    uint32_t *seq_ll, *seq_ml, *seq_dist;  /* dist: 0..2 = rep idx, else dist + 2 */
    uint8_t *lits, *llc, *mlc, *ofc, *blk;
    uint32_t *xlm, *xof;   /* type 5: extra bits (ll | ml << n, offset)     */
    uint8_t *nlm, *nof;    /* and their bit counts                          */
};

#define MAX_SEQ (LZR_BLOCK_SIZE / LZR_ML_BASE + 1)
#define LZR_MAX_THREADS 16
#ifndef LZR_ADAPTIVE_LEVEL
#define LZR_ADAPTIVE_LEVEL 9   /* levels >= this use block type 4 */
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
/* The lazy levels hash the next LZR_HASH_LEN bytes (enc.c row_hash): more
 * than MIN_MATCH keeps rows free of short, useless candidates; 4-byte
 * matches are still found through the repeat offsets. Hashes read 8
 * bytes; the buffer has slack so the read stays in bounds. */

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
