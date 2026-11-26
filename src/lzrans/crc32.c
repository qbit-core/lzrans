/*
 * crc32.c - CRC-32 (IEEE 802.3, reflected poly 0xEDB88320), slice-by-8.
 * Output is bit-identical to zlib.crc32 / binascii.crc32.
 */
#include "crc32.h"

static uint32_t crc_table[8][256];
static int crc_ready = 0;

static void crc32_init(void)
{
    for (uint32_t i = 0; i < 256; i++) {
        uint32_t c = i;
        for (int k = 0; k < 8; k++)
            c = (c & 1) ? (c >> 1) ^ 0xEDB88320u : c >> 1;
        crc_table[0][i] = c;
    }
    for (uint32_t i = 0; i < 256; i++)
        for (int t = 1; t < 8; t++)
            crc_table[t][i] = (crc_table[t - 1][i] >> 8) ^ crc_table[0][crc_table[t - 1][i] & 0xFF];
    crc_ready = 1;
}

/* Table init is idempotent, so a benign race between threads is harmless;
 * the module also calls lzr_crc32_init() once at import time. */
void lzr_crc32_init(void) { if (!crc_ready) crc32_init(); }

uint32_t lzr_crc32(uint32_t crc, const uint8_t *p, size_t n)
{
    lzr_crc32_init();
    crc = ~crc;
    /* Process 8 bytes per iteration (little-endian assumption on the word reads). */
    while (n >= 8) {
        uint32_t lo = rd32le(p) ^ crc, hi = rd32le(p + 4);
        crc = crc_table[7][lo & 0xFF] ^ crc_table[6][(lo >> 8) & 0xFF] ^
              crc_table[5][(lo >> 16) & 0xFF] ^ crc_table[4][lo >> 24] ^
              crc_table[3][hi & 0xFF] ^ crc_table[2][(hi >> 8) & 0xFF] ^
              crc_table[1][(hi >> 16) & 0xFF] ^ crc_table[0][hi >> 24];
        p += 8; n -= 8;
    }
    while (n--) crc = (crc >> 8) ^ crc_table[0][(crc ^ *p++) & 0xFF];
    return ~crc;
}
