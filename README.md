# lzrans

Streaming LZ77 + rANS compressor. The core is C; Python uses it through a
CPython native extension.

**Goal:** a better ratio than xz/LZMA with faster compression *and*
decompression, usable on streams. Status: v0.4 phase 5 (see Roadmap).

```python
import lzrans

blob = lzrans.compress(data, level=6)          # level 1 (fast) .. 9 (best)
blob = lzrans.compress(data, level=9, threads=4)   # same bytes, faster
blob = lzrans.compress(exe_bytes, level=9, filter="x86")   # x86/x86-64 code
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
| Parsing | Levels 7-9: price-based optimal parsing (forward DP over segments of up to 4096 positions, repeat offsets and coder state tracked per path). Prices come from a shadow copy of the adaptive model, advanced over the chosen tokens and re-tabulated per segment; a segment ends early (after >= 512 positions) where all paths meet. Besides tree matches the DP tries 2-3 byte matches at short distances, 1-byte "short rep" matches, and match + literal + rep0 / literal + rep0 combinations. |
| Filter | `x86.h` (opt-in, `filter="x86"`): CALL (E8) rel32 operands with displacement < 2^25 are rewritten as absolute targets (mod 2^26), so repeated calls become repeated bytes. Bijective and streaming (holds back <= 4 bytes); the CRC covers the original content. ~4% smaller on an x86-64 shared object. |
| Threads | Levels 7-9: hash buckets are split between threads (each tree only links positions of its own bucket, so no locks), and the next block is gathered in the background while the current one is parsed. Output is byte-identical for any thread count. |
| Tokens | Sequences of (literal run, match length, offset); offsets 0-2 are the repeat offsets. |
| Entropy coding, levels 1-5 | Block type 1, `rans.c`: literals and LL/ML/offset log-bucket codes (62 codes + raw extra bits) as rANS streams, 12-bit static model per 256 KB block, 2 interleaved states, raw/run fallbacks. Fast to decode. |
| Entropy coding, levels 6-9 | Block type 4, `arc.h`: the coding layer of LZMA on our tokens. One adaptive binary range coder (12-bit probabilities); a 12-value state of the last token kinds selects literal/match/rep decisions; literals use the previous byte's top 3 bits as context and, right after a match, the byte at rep0 ("matched literal"); match lengths 2-273 by low/mid/high trees per position & 3; distances by 6-bit slot (context: length) plus reverse-tree or direct + 4 aligned bits. The model, state and repeat offsets carry over between blocks. |
| Decoding | `dec.c`: type 1: code streams decoded into arrays first, then a sequence loop with 16-byte wildcopy. Type 4: one token loop decoding straight into the history. |
| Streaming | Input is consumed and output produced block by block. Encoder memory is about 2x window + chain table; decoder memory is about 2x window. |
| Integrity | CRC-32 (`crc32.c`, identical to `zlib.crc32`) in the end block, plus structural checks: every rANS stream must end at its initial state. |

## Frame format v4 (little-endian)

```
"LZRN" | version u8 = 4 | flags u8 | window_log u8 | reserved u8 | [content_size u64]
block*: type u8 | raw_len u32 | comp_len u32 | payload
  type 0 stored      payload = raw bytes
  type 1 compressed  nseq u32 | nlit u32 | literals | ll | ml | of streams | extra_len u32 | extra bits
  type 2 end         payload = CRC32 of the whole content
  type 4 compressed  one adaptive range-coded token stream (arc.h); levels 6-9
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
| P6 | Type 4 decoder speed; compression speed of the priced parser | planned |
