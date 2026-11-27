"""Split a zstd frame into literal and sequence sections (RFC 8878).

Returns totals: literals regenerated, literal-section bytes, sequences,
sequence-section bytes (FSE tables + bit stream), raw/RLE block bytes and
frame overhead. Used by tools/entropy_report.py to compare against
lzrans block type 5 (LZR_STATS build).
"""


def breakdown(z: bytes) -> dict:
    st = dict(nlit=0, lit_bytes=0, nseq=0, seq_bytes=0, raw_bytes=0, overhead=0, blocks=0)
    p = 0
    while p < len(z):
        assert z[p:p + 4] == b'\x28\xb5\x2f\xfd', 'not a zstd frame'
        fhd = z[p + 4]
        q = p + 5
        single = fhd >> 5 & 1
        if not single:
            q += 1                                   # window descriptor
        q += (0, 1, 2, 4)[fhd & 3]                   # dictionary id
        q += (single, 2, 4, 8)[fhd >> 6]             # frame content size
        st['overhead'] += q - p
        while True:
            h = z[q] | z[q + 1] << 8 | z[q + 2] << 16
            last, btype, bsize = h & 1, h >> 1 & 3, h >> 3
            q += 3
            st['overhead'] += 3
            st['blocks'] += 1
            if btype == 0:
                st['raw_bytes'] += bsize; q += bsize
            elif btype == 1:
                st['raw_bytes'] += 1; q += 1
            else:
                b = z[q:q + bsize]
                lt, sf = b[0] & 3, b[0] >> 2 & 3
                if lt in (0, 1):                     # raw / RLE literals
                    if sf in (0, 2):
                        hl, regen = 1, b[0] >> 3
                    elif sf == 1:
                        hl, regen = 2, (b[0] >> 4) + (b[1] << 4)
                    else:
                        hl, regen = 3, (b[0] >> 4) + (b[1] << 4) + (b[2] << 12)
                    ls = hl + (regen if lt == 0 else 1)
                else:                                # Huffman (with or without tree)
                    hl, nb = {0: (3, 10), 1: (3, 10), 2: (4, 14), 3: (5, 18)}[sf]
                    v = int.from_bytes(b[:hl], 'little') >> 4
                    regen, comp = v & ((1 << nb) - 1), v >> nb & ((1 << nb) - 1)
                    ls = hl + comp
                st['nlit'] += regen
                st['lit_bytes'] += ls
                r = b[ls:]
                if r[0] == 0:
                    ns = 0
                elif r[0] < 128:
                    ns = r[0]
                elif r[0] < 255:
                    ns = ((r[0] - 128) << 8) + r[1]
                else:
                    ns = r[1] + (r[2] << 8) + 0x7F00
                st['nseq'] += ns
                st['seq_bytes'] += len(r)
                q += bsize
            if last:
                break
        if fhd >> 2 & 1:
            q += 4
            st['overhead'] += 4
        p = q
    return st
