# lzrans

Streaming LZ77 compressor with two entropy back ends: a fast one
(Huffman + tANS, zstd-class decode speed) and an LZMA-class one. The core is C; Python uses it through a
CPython native extension.

**Goal:** levels 1-8 match zstd -19 in ratio and speed; level 9 beats
xz -9 in ratio. Streaming throughout. Status: v0.7 phase 8 (see Roadmap).

```python
import lzrans

blob = lzrans.compress(data, level=8)          # 1 (fast) .. 8: fast decode; 9: best ratio
blob = lzrans.compress(data, level=9, threads=4)   # faster (levels 7-9: same bytes)
blob = lzrans.compress(exe_bytes, level=9, filter="x86")   # x86/x86-64 code
assert lzrans.decompress(blob) == data         # ChecksumError on checksum mismatch

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
| Match finding | Levels 1-6 (`enc.c`): row match finder over a 1-16 MB window: a 5-byte hash picks a row of the 16-64 most recent positions plus an 8-bit tag per slot; tags are compared 16 at a time (SSE2, scalar fallback) and only tag hits are checked, with rows prefetched ahead of insertion. 3 repeat offsets are probed first; lazy matching (depth 0-2). Levels 7-9 (`opt.c`): binary tree (BT4) over an 8-16 MB window collecting all match lengths per position. |
| Parsing | Levels 7-9: price-based optimal parsing (forward DP over segments of up to 4096 positions, repeat offsets and coder state tracked per path). Levels 7-8 price with the static code lengths of the previous block's statistics; level 9 prices from a shadow copy of the adaptive model, advanced over the chosen tokens and re-tabulated per segment; a segment ends early (after >= 512 positions) where all paths meet. Besides tree matches the DP tries 2-3 byte matches at short distances, 1-byte "short rep" matches, and match + literal + rep0 / literal + rep0 combinations. |
| Filter | `x86.h` (opt-in, `filter="x86"`): CALL (E8) rel32 operands with displacement < 2^25 are rewritten as absolute targets (mod 2^26), so repeated calls become repeated bytes. Bijective and streaming (holds back <= 4 bytes); the checksum covers the original content. ~4% smaller on an x86-64 shared object. |
| Threads | Levels 1-6: the input is cut into 4 MB jobs compressed independently (each first indexes the preceding window, or half of a 16 MB window); repeat offsets start unknown and the first block of a job builds new tables, so the format is unchanged. Threads take jobs as they finish; output is identical for any thread count > 1 and within ~0.1% of one thread. Levels 7-9: hash buckets are split between threads (each tree only links positions of its own bucket, so no locks), and the next block is gathered in the background while the current one is parsed. Output is byte-identical for any thread count. |
| Tokens | Sequences of (literal run, match length, offset); offsets 0-2 are the repeat offsets. |
| Entropy coding, levels 1-8 | Block type 5, `huf.c` + `fse.c` (zstd-like): literals as 4 canonical Huffman streams (codes <= 11 bits, table reused across blocks when cheaper); literal-length / match-length / offset codes by table-driven ANS (tANS, 512 states, tables new / repeated / single-symbol per block) sharing one backward bit stream with the raw extra bits. |
| Entropy coding, level 9 | Block type 4, `arc.h`: the coding layer of LZMA on our tokens. One adaptive binary range coder (12-bit probabilities); a 12-value state of the last token kinds selects literal/match/rep decisions; literals use the previous byte's top 3 bits as context and, right after a match, the byte at rep0 ("matched literal"); match lengths 2-273 by low/mid/high trees per position & 3; distances by 6-bit slot (context: length) plus reverse-tree or direct + 4 aligned bits. The model, state and repeat offsets carry over between blocks. |
| Decoding | `dec.c`: type 5: Huffman literals into an array (4 streams interleaved, one table lookup per symbol); then one loop decodes each sequence (one bit-stream read in the common case, branch-free repeat offsets) and copies with fixed 16-byte chunks. The hot loops are also compiled for BMI2 and picked at run time. Type 4: one token loop. One-shot `decompress` decodes straight into the output buffer (no history copy). |
| Streaming | Input is consumed and output produced block by block. Encoder memory is about 2x window + match tables (levels 1-6: 5 bytes per row slot, up to 40 MB at level 6); decoder memory is about 2x window. |
| Integrity | Low 32 bits of XXH64 (`xxh64.h`, same value as zstd's frame checksum) in the end block, plus structural checks: every bit stream must be consumed exactly. (`lzrans.crc32` stays available as a zlib-compatible helper.) |

## Frame format v5 (little-endian)

```
"LZRN" | version u8 = 5 | flags u8 | window_log u8 | reserved u8 | [content_size u64]
block*: type u8 | raw_len u32 | comp_len u32 | payload
  type 0 stored      payload = raw bytes
  type 2 end         payload = low 32 bits of XXH64 (seed 0) of the whole content
  type 4 tokens      one adaptive range-coded token stream (arc.h); level 9
  type 5 fast        nseq u32 | nlit u32 | Huffman literal section (huf.h)
                     | if nseq: ll/ml/of table modes | backward tANS bit stream (fse.h)
```

`flags`: bit 0 = encrypted (reserved), bit 1 = content size present
(set by one-shot `compress`), bit 2 = x86 filter applied to the content. Decoders reject unknown versions and flags.

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
| P3 | Adaptive literal coding (literal contexts, match-byte context), short repeat matches | done |
| P4 | x86 filter (binaries) | done |
| P5 | LZMA-style token coding (block type 4) + model-priced parsing: ratio at or above xz -9 | done |
| P6 | Block type 5 (Huffman + tANS) for levels 1-8: zstd -19 class ratio and speed | done |
| P7 | Row match finder for levels 1-6 (2-5x faster compression at about the same ratio) | done |
| P8 | Parallel jobs for levels 1-6 (threads > 1, format unchanged) + single-thread micro-optimizations | done |
| P9 | Type 4 decoder speed; compression speed of the priced parser | planned |
