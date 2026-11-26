# lzrans

LZ77 (64 KB hash-chain) + 3-stream rANS compressor. The core is C; Python
uses it through a CPython native extension.

```python
import lzrans
blob = lzrans.compress(data, level=5)   # level 1 (fast) .. 9 (best)
assert lzrans.decompress(blob) == data  # raises ChecksumError on CRC mismatch
```

## Build and test

```sh
pip install -e .[test]
pytest                   # Python round-trip and corruption tests
make -C tests/c run      # C core under ASan/UBSan, incl. mutated-input fuzzing
python bench/bench.py    # vs zlib
```

## Pipeline

| Stage | Implementation (`src/lzrans/`) |
|---|---|
| Match finding | `lz77.c`: 4-byte multiplicative hash (`HASH_BITS=16`), hash chains, 64 KB window, 32-bit word compares, one-step lazy matching at level >= 4. Matches are 4..259 bytes. |
| Tokens | Stream 1: literal/length symbols 0..511 (0-255 literal, 256-511 = match length 4..259). Stream 2: distance high byte. Stream 3: distance low byte. |
| Entropy coding | `rans.c`: byte-wise rANS, 32-bit state, 12-bit scale (M = 4096), one static frequency table per stream per 1 MB block. |
| Decoding | `frame.c`: 4096-entry slot table gives (symbol, freq, bias) in one lookup; literals and matches are copied sequentially with `memcpy` (overlapping matches are copied by doubling the period). |
| Integrity | `crc32.c`: slice-by-8 CRC-32, identical to `zlib.crc32`, verified after decompression. Each rANS stream must also end exactly at its initial state, which catches most corruption before the CRC. |

Incompressible blocks are stored raw, so output is at most `20 + 5 * blocks` bytes larger than the input.

## Frame format (v1, little-endian)

```
"LZRN" | version u8 | flags u8 | level u8 | reserved u8 | original_size u64 | crc32 u32
block*:  type u8 (0 = stored, 1 = lz+rans) | raw_len u32 | payload
  type 1: nsym u32 | nmatch u32 | 3 x [freq table | len u32 | rANS bytes]
```

Decoders reject unknown versions and unknown flag bits.

## Encryption hook (not implemented yet)

`flags` bit 0 (`lzrans.FLAG_ENCRYPTED`) is reserved. The plan is to wrap
the compressed frame in an AEAD envelope (libsodium XChaCha20-Poly1305 with
an Argon2id-derived key), so: compress, then encrypt; decrypt and verify the
tag, then decompress. This build raises `LzransError` on frames with that
flag set.
