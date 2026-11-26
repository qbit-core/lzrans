#ifndef LZRANS_FRAME_H
#define LZRANS_FRAME_H
#include "common.h"

/*
 * Container ("frame") format, version 1. All integers little-endian.
 *
 *   offset size  field
 *   0      4     magic "LZRN"
 *   4      1     version (=1)
 *   5      1     flags   bit0 = LZR_FLAG_ENCRYPTED (reserved, not yet
 *                        implemented; decoders MUST reject unknown bits)
 *   6      1     level used by the encoder (informational)
 *   7      1     reserved (0)
 *   8      8     original (uncompressed) size
 *   16     4     CRC32 of the original data
 *   20     ...   blocks, until `original size` bytes are produced
 *
 * Block:
 *   u8  type     0 = stored (raw bytes follow), 1 = LZ77 + 3-stream rANS
 *   u32 raw_len  bytes this block produces (<= LZR_BLOCK_SIZE)
 *   type 0: raw_len bytes
 *   type 1: u32 nsym, u32 nmatch, then 3 streams
 *           (litlen[nsym], dist_hi[nmatch], dist_lo[nmatch]); each stream is
 *           freq table + u32 payload_len + rANS payload (omitted if count=0)
 *
 * Matches may reference earlier blocks (64 KB window spans blocks).
 */
#define LZR_MAGIC          "LZRN"
#define LZR_VERSION        1
#define LZR_HEADER_SIZE    20
#define LZR_BLOCK_SIZE     (1u << 20)
#define LZR_FLAG_ENCRYPTED 0x01u   /* hook for future AEAD envelope */
#define LZR_KNOWN_FLAGS    0x00u   /* flags this build can decode    */

size_t lzr_compress_bound(size_t n);

/* Returns compressed size (>0) or a negative LZR_ERR_* code. */
int64_t lzr_compress(const uint8_t *src, size_t n, uint8_t *dst, size_t cap, int level);

/* Validates the header and reports the original size. 0 or LZR_ERR_*. */
int lzr_frame_info(const uint8_t *src, size_t n, uint64_t *orig_size, uint32_t *crc, uint8_t *flags);

/* Decompresses into dst (cap must equal the original size). 0 or LZR_ERR_*. */
int lzr_decompress(const uint8_t *src, size_t n, uint8_t *dst, size_t cap);

#endif
