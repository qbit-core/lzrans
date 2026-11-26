#ifndef LZRANS_CRC32_H
#define LZRANS_CRC32_H
#include "common.h"

void lzr_crc32_init(void);
/* Incremental: pass 0 as the initial crc, then feed the previous result. */
uint32_t lzr_crc32(uint32_t crc, const uint8_t *p, size_t n);

#endif
