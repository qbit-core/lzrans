# lzrans

Streaming LZ77 compressor with two entropy back ends: a fast one
(Huffman + tANS, zstd-class decode speed) and an LZMA-class one. The core is C; Python uses it through a
CPython native extension.

**Goal:** levels 1-10 match zstd in ratio and speed (level 10: zstd -19);
level 11 beats xz -9 in ratio; level 0 is a zstd -3 class fast level,
levels -1..-3 are faster still (zstd -1 class ratio). Streaming throughout.
Status: v0.17 phase 21, format 1.0 (see Roadmap). Levels 7-11 were renumbered in v0.14:
the old levels 7, 8, 9 are now 9, 10, 11.

```python
import lzrans

blob = lzrans.compress(data, level=10)         # -3 (fastest) .. 10: fast decode; 11: best ratio
blob = lzrans.compress(data, level=11, threads=4)  # faster (levels 9-11: same bytes)
blob = lzrans.compress(exe_bytes, level=11, filter="x86")  # x86/x86-64 code
assert lzrans.decompress(blob) == data         # ChecksumError on checksum mismatch
assert lzrans.decompress(blob, threads=4) == data   # parallel: segments, or blocks (levels -3..10)
blob = lzrans.compress(data, level=6, segment=16 << 20)  # reset points every 16 MiB (any level)

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
| Match finding | Levels -1..-3 (`enc.c` `parse_block_fast1`): greedy single-hash finder like zstd's "fast": one table of the newest position per 6-byte hash (8K-16K entries, cache-resident) over a 0.5-1 MB window; the table slot and candidate bytes of the next position are loaded before the current one is resolved, a miss advances by 1-4 bytes plus one per 256 literals, and 2 positions per match are indexed. Level 0 (`enc.c` `parse_block_fast`): greedy double-hash finder over a 4 MB window, like zstd's "dfast": one table of the newest position per 8-byte hash and one per 5-byte hash (64K entries each); repeat offset 0 is tried one byte ahead, a short match is replaced by a long one starting a byte later when there is one, matches are extended backwards into the literals, only searched positions and 3 positions per match are indexed, and the search step grows by one every 256 literals. Levels 1-6 (`enc.c`): row match finder over a 1-16 MB window: a 5-byte hash picks a row of the 16-64 most recent positions plus an 8-bit tag per slot; tags are compared 16 at a time (SSE2, scalar fallback) and only tag hits are checked, with rows prefetched ahead of insertion. 3 repeat offsets are probed first; lazy matching (depth 0-2). Levels 7-8 (`enc.c` `row_gather`): the same rows (32-64 slots, 16 MB window) searched at every position (16-32 slots) for matches of increasing length; a match found one byte earlier is carried over one byte shorter, and inside a match of 16+ bytes only 2 slots are tried. Levels 9-11 (`opt.c`): binary tree (BT4) over a 32-64 MB window (level 9: 32 MB, 10-11: 64 MB; smaller when the input size is known and smaller) collecting all match lengths per position. |
| Parsing | Levels 7-11: price-based optimal parsing (forward DP over segments of up to 4096 positions, repeat offsets and coder state tracked per path). Levels 7-10 price with the static code lengths of the previous block's statistics; level 11 prices from a shadow copy of the adaptive model, advanced over the chosen tokens and re-tabulated per segment; a segment ends early (after >= 512 positions) where all paths meet. Besides tree matches the DP tries 2-3 byte matches at short distances, 1-byte "short rep" matches, and match + literal + rep0 / literal + rep0 combinations. Levels 7-8 skip the short matches, level 7 also the combinations, and both price every length only up to 12 (longer matches at full length). |
| Filter | `x86.h` (opt-in, `filter="x86"`): CALL (E8) rel32 operands with displacement < 2^25 are rewritten as absolute targets (mod 2^26), so repeated calls become repeated bytes. Bijective and streaming (holds back <= 4 bytes); the checksum covers the original content. ~4% smaller on an x86-64 shared object. |
| Threads | Levels -3..6: the input is cut into 4 MB jobs compressed independently (each first indexes the preceding window, or half of a 16 MB window; level 0 every third position of it); repeat offsets start unknown and the first block of a job builds new tables, so the format is unchanged. Threads take jobs as they finish; output is identical for any thread count > 1 and within ~0.1% of one thread. When `compress()` gets the whole input at once, the jobs read it in place (no copy) except the last ones, and the checksum is computed while the helpers compress; the jobs and the output are the same as when streamed. Levels 7-8: one thread. Levels 9-11: hash buckets are split between threads (each tree only links positions of its own bucket, so no locks), and the next block is gathered in the background while the current one is parsed. Output is byte-identical for any thread count. |
| Tokens | Sequences of (literal run, match length, offset); offsets 0-2 are the repeat offsets. |
| Entropy coding, levels -3..10 | Block type 5, `huf.c` + `fse.c` (zstd-like): literals as 4 canonical Huffman streams (codes <= 11 bits, table reused across blocks when cheaper); literal-length / match-length / offset codes by table-driven ANS (tANS, 512 states, tables new / repeated / single-symbol per block) sharing one backward bit stream with the raw extra bits. |
| Context literals, levels 9-10 | Block type 6 = type 5 with another literal section (`huf.c` `hufc_*`): each literal is coded with one of up to 8 Huffman tables, picked by the byte before it through a 256-entry map. The encoder groups the 256 previous-byte values by a few rounds of k-means on the block's (previous byte, literal) counts (2, 4 or 8 groups, seeded by the top bits or by character classes), gives every symbol of the block a code in every table so later blocks can reuse them, and uses type 6 only when the literal section is >= 2% smaller than type 5's. Levels 7-8 keep type 5: there the decode cost (-11..-16%) outweighed +0.7%. The decoder decodes the literals while executing the sequences: one table lookup per literal, whose entry also holds the group of the decoded byte for the next one. |
| Entropy coding, level 11 | Block type 4, `arc.h`: the coding layer of LZMA on our tokens. One adaptive binary range coder (12-bit probabilities); a 12-value state of the last token kinds selects literal/match/rep decisions; literals use the previous byte's top 3 bits as context and, right after a match, the byte at rep0 ("matched literal"); match lengths 2-273 by low/mid/high trees per position & 3; distances by 6-bit slot (context: length) plus reverse-tree or direct + 4 aligned bits. The model, state and repeat offsets carry over between blocks. |
| Segments | `segment=` (levels 9-11 default 48 MiB, others off): the content is cut into segments compressed independently, each starting at a reset point, with a seek index (skippable block 0x80, FORMAT.md 12.3) before the end block. Segments are compressed in parallel (as many as fit in about 1 GiB of encoder memory; spare threads go to the match finders) and the output is identical for any thread count. Inputs no longer than one segment give exactly the unsegmented output. `decompress(threads=N)` finds the segments from the block headers and decodes them straight into the output buffer in parallel; the last block of each segment is decoded without wildcopy past its end, and the calling thread hashes finished segments in order. With fewer segments than threads (or none), each segment goes through the block pipeline below with its share of the threads, the segments side by side. |
| Decoding | `dec.c`: type 5: Huffman literals into an array (4 streams interleaved, one table lookup per symbol); then one loop decodes each sequence (one bit-stream read in the common case, branch-free repeat offsets) and copies with fixed 16-byte chunks. The hot loops are also compiled for BMI2 and picked at run time. Type 4: one token loop; tree bits (literals, lengths, distance slots) are decoded without branches (conditional moves), and the next tree node's probabilities are loaded before the current bit is known. One-shot `decompress` decodes straight into the output buffer (no history copy). Block pipeline (`decompress(threads=N)`, types 5-6, format unchanged): a pass over the block headers finds where each block's Huffman and tANS tables are described; worker threads then run the entropy stage several blocks ahead (type 5 literals into an array, the sequence stream into (literal length, match length, offset value) triples, tables rebuilt per worker only when they change) and fault in the block's output pages, while the calling thread executes the blocks in order (repeat offsets, literal and match copies, type 6 literals, which need the byte before each, and type 0 / 4 blocks whole) and decodes a later block itself when the next one is not ready; a worker hashes executed blocks in order. Same output and the same errors as one thread (a malformed header sends the frame to the serial decoder). |
| Streaming | Input is consumed and output produced block by block. Encoder memory is about 2x window + match tables (level 0: 512 KB; levels 1-6: 5 bytes per row slot, up to 40 MB at level 6, 40-48 MB at levels 7-8; levels 9-11: 8 bytes per window byte for the tree, about 680 MB at 64 MB, close to xz -9); decoder memory is about 2x window. |
| Integrity | Low 32 bits of XXH64 (`xxh64.h`, same value as zstd's frame checksum) in the end block, plus structural checks: every bit stream must be consumed exactly. (`lzrans.crc32` stays available as a zlib-compatible helper.) |

## Frame format v7 = LZR 1.0 (little-endian)

The full specification, with golden test vectors and a pure-Python
reference decoder (`tools/refdec.py`), is [docs/FORMAT.md](docs/FORMAT.md);
this is a summary. Version 5 and 6 frames still decode.

```
"LZRN" | version u8 = 7 | flags u8 | window_log u8 | reserved u8 | [content_size u64]
block*: type u8 | raw_len u32 | comp_len u32 | payload
  type byte: bit 7 = skippable (raw_len 0, payload ignored); else bits 0-5 =
             type, bit 6 = reset point (state reset, no references before it)
  type 0 stored      payload = raw bytes
  type 2 end         payload = low 32 bits of XXH64 (seed 0) of the whole content
  type 4 tokens      one adaptive range-coded token stream (arc.h); level 11
  type 5 fast        nseq u32 | nlit u32 | Huffman literal section (huf.h)
                     | if nseq: ll/ml/of table modes | backward tANS bit stream (fse.h)
  type 6 ctx literals as type 5, literal section = mode u8 (0 new / 1 previous tables)
                     | mode 0: K u8, 256 nibbles (previous byte -> table), K tables
                     | stream size u32 | forward LSB-first stream (huf.h)
  types 7-63         newer formats: rejected as unsupported
```

`flags`: bit 0 = encrypted (reserved), bit 1 = content size present
(set by one-shot `compress`), bit 2 = x86 filter applied to the content,
bits 3-7 = reserved for required features. Decoders reject unknown versions
and flags.

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
| P9 | Type 4 decoder speed (branch-free tree decoding; format and output unchanged): 1.2-1.6x, about xz -9 decode speed or faster | done |
| P10 | Level 9 single-thread compression speed (parser and tree walk, output unchanged): 1.0-1.3x, about xz -9 except on source code | done |
| P11 | Level 0 (double-hash fast level, zstd -3 class: 0.84x its speed at 1 thread, 1.16x at 4 threads) + levels 1-6 1.0-1.1x faster with the same output | done |
| P12 | Same output, faster: job-mode threads overlap input copy/checksum with compression (L0 4 threads +17%, 1.22x zstd -3; L1-6 +5-8%); branch-free type 5 bit writers (L0 1 thread +8%) | done |
| P13 | Ultra-fast levels -1..-3 (single-hash finder; 4 threads 790-890 MB/s, ratio 3.2-3.8; about 0.65x zstd -1 speed at +3% ratio) + threaded one-shot compress() reads the input in place (L0 4 threads +23%, 1.42x zstd -3; L1-3 +9-16%; same output) | done |
| P15 | Redesign phase 1: levels 7-9 window 16 MB -> 32/64 MB (enwik8 L9 3.891 -> 4.032, above xz -9 4.022; Silesia L9 +0.4%, L8 +0.4%, L7 +0.7%; 1-thread compress 7-25% slower, memory about xz -9) | done |
| P16 | Redesign phase 2: cost-based block splitting for optimal-parse levels 7-8 (each block tried as up to 8 parts, kept only if smaller; format unchanged). Silesia L7 +0.4%, L8 +0.4% (xml +0.7%, mr +0.9%), speed unchanged; enwik8 unchanged | done |
| P17 | Redesign phase 3: levels 7-8 (new) = row finder feeding a light optimal parser; old levels 7-9 renumbered 9-11. Silesia 1 thread: L7 3.720 @ 5.5 MB/s (zstd -14 3.700 @ 5.4), L8 3.785 @ 3.6 MB/s (zstd -15 3.726 @ 3.7); enwik8 L7 +2.2% vs zstd -14 at 1.26x speed, L8 +3.4% vs zstd -15 at 1.28x. Speed gap between lazy level 6 and the tree levels 10x -> 3.3x | done |
| P18 | Frame format v6, block type 6: context literals for levels 9-10 (up to 8 Huffman tables chosen by the previous byte, picked per block, type 5 kept when the gain is under 2%). Silesia 1 thread vs v0.14: ratio +0.71% (L9) / +0.72% (L10, 4.027), compress -4% / 0%, decode -6% / -4%; enwik8 +0.02%. Levels 7-8 measured +0.7% at -11..-16% decode and were left on type 5 (output unchanged). v5 frames still decode | done |
| P19 | Format freeze: frame version 7 = LZR 1.0, specified in docs/FORMAT.md with golden test vectors (tests/vectors/) and a pure-Python reference decoder (tools/refdec.py). New in v7: reset points (parallel decoding / seeking later), skippable blocks, unsupported-vs-corrupt for newer block types, stored size of Huffman stream 3, complete Huffman codes only, comp_len <= raw_len. Compression output otherwise unchanged (+4 bytes per Huffman literal section). v5/v6 frames still decode | done |
| P20 | Parallel decompression: segment mode (reset points every 48 MiB by default at levels 9-11, `segment=` for any level) + seek index (skippable block 0x80, FORMAT.md 12.3) + `decompress(threads=N)`. Silesia L11: ratio -0.38% (still +0.6% vs xz -9), 4-thread decode 3.4x xz -9; L10 at zstd -19 ratio decodes 1.9x faster with 4 threads. enwik8 loses 1-2.6% (long-range repeats). Unsegmented output and 1-thread decode unchanged | done |
| P21 | Pipelined parallel decoding of one segment (format and output unchanged): `decompress(threads=N)` runs the type 5/6 entropy stage on worker threads ahead of the in-order block execution, and segments fewer than threads each get a pipeline. 4 threads vs 1 (v0.13): Silesia L8 2.5x, L10 2.3x (1325 MB/s, 1.6x zstd -19 at +0.4% ratio); enwik8 L10 unsegmented 2.8x (1140 MB/s, 1.8x zstd -19 at +5.1%), L10 default segments 2.1x. Single thread: sequences decoded 256 at a time before execution at windows >= 8 MB and for type 6 (enwik8 L10 1.37x, L6 1.16x, Silesia 1.03-1.05x). Level 11 (type 4) unchanged | done |
