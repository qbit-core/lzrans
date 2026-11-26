# lzrans

Streaming LZ77 + rANS compressor. The core is C; Python uses it through a
CPython native extension.

**Goal:** a better ratio than xz/LZMA with faster compression *and*
decompression, usable on streams. Status: v0.2 phase 2 (see Roadmap).

```python
import lzrans

blob = lzrans.compress(data, level=6)          # level 1 (fast) .. 9 (best)
assert lzrans.decompress(blob) == data         # ChecksumError on CRC mismatch

# streaming, bounded memory
c = lzrans.Compressor(level=6)
out = b"".join(c.compress(chunk) for chunk in chunks) + c.flush()

d = lzrans.Decompressor()
for chunk in compressed_chunks:
    sink.write(d.decompress(chunk))
assert d.eof
```

## Build and test

```sh
pip install -e .[test]
pytest                   # round-trip, streaming and corruption tests
make -C tests/c run      # C core under ASan/UBSan, incl. mutated-input fuzzing
python bench/bench.py    # vs zlib, xz and zstd (pip install zstandard)
```

## Design

| Stage | Implementation (`src/lzrans/`) |
|---|---|
| Match finding | Levels 1-6 (`enc.c`): hash chains over a 1-8 MB window, 5-byte hash, 3 repeat offsets probed first, lazy matching (depth 0-2). Levels 7-9 (`opt.c`): binary tree (BT4) over an 8-16 MB window collecting all match lengths per position. |
| Parsing | Levels 7-9: price-based optimal parsing (forward DP over 4096-position segments, repeat offsets tracked per path, prices from the previous pass/block statistics; 2 passes at levels 8-9). |
| Tokens | Sequences of (literal run, match length, offset). Offset code 0-2 = repeat offset, else a log-bucket code. Lengths/offsets use log-bucket codes (62 codes) + raw extra bits. |
| Entropy coding | `rans.c`: 4 rANS streams per 256 KB block (literals, LL codes, ML codes, offset codes), 12-bit static model per block, 2 interleaved states, raw/run fallbacks. |
| Decoding | `dec.c`: streams decoded into arrays with one table lookup per symbol, then a sequence loop with 16-byte wildcopy for literals and matches. |
| Streaming | Input is consumed and output produced block by block. Encoder memory is about 2x window + chain table; decoder memory is about 2x window. |
| Integrity | CRC-32 (`crc32.c`, identical to `zlib.crc32`) in the end block, plus structural checks: every rANS stream must end at its initial state. |

## Frame format v2 (little-endian)

```
"LZRN" | version u8 = 2 | flags u8 | window_log u8 | reserved u8 | [content_size u64]
block*: type u8 | raw_len u32 | comp_len u32 | payload
  type 0 stored      payload = raw bytes
  type 1 compressed  nseq u32 | nlit u32 | literals | ll | ml | of streams | extra_len u32 | extra bits
  type 2 end         payload = CRC32 of the whole content
```

`flags`: bit 0 = encrypted (reserved), bit 1 = content size present
(set by one-shot `compress`). Decoders reject unknown versions and flags.

## Encryption hook (not implemented yet)

Flag bit 0 (`lzrans.FLAG_ENCRYPTED`) is reserved for an AEAD envelope
(planned: libsodium XChaCha20-Poly1305, Argon2id key derivation): compress,
then encrypt per block; verify, then decompress. Frames with the flag set
are rejected for now.

## Roadmap

| Phase | Change | Status |
|---|---|---|
| P1 | Streaming frame, large window, sequence tokens, repeat offsets | done |
| P2 | Binary-tree match finder + price-based optimal parsing | done |
| P3 | Context-modelled entropy coding (literal contexts, match-byte context) | next |
| P4 | Decoder speed (more interleaved states, branchless refill), x86 filter | planned |
