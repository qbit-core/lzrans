# The LZR compressed data format

**Version 1.0** (frame version 7). Also defines how to read legacy frame versions 5 and 6 (Appendix C).

**Status:** frozen once this document is merged. From then on, the format changes only as section 12 allows: new block types and features are added behind a feature flag or a new version number, and nothing a 1.0 decoder accepts changes meaning. Appendix A records the decisions taken for the freeze.

A second decoder, `tools/refdec.py` (pure Python, about 730 lines), was written from this document alone. It decodes every golden vector in `tests/vectors/` (section 13) and everything the current encoder produces at levels -3 to 11, with the same result as the C decoder (`src/lzrans/dec.c`). Treat any difference between this text, `refdec.py` and `dec.c` as a bug.

---

## 1. Overview

An LZR **frame** holds one stream of bytes, the **content**. A frame has three parts:

```
frame  = header  block*  end-block
```

* The **header** (section 3.1) carries a magic number, the frame version, flags, and the window size. It can also carry the content size.
* **Blocks** (section 3.2) each decode to at most 256 KiB of content. Blocks are decoded in order. A block may copy bytes from earlier blocks, up to the window size back. A block may also reuse coding tables or model state left by earlier blocks. A block marked as a **reset point** (section 5.5) depends on nothing before it, so the stretch from one reset point to the next can be decoded on its own: in parallel, or after a seek. **Skippable blocks** (section 3.2) carry data that does not change the content.
* The **end block** carries a checksum of the content (section 9).

Most of the format is LZ77: sequences of *literal bytes* and *matches* that copy bytes from the decoded output. There are three ways to code a block:

| Type | Name | Coding | Written by |
|---|---|---|---|
| 0 | stored | raw bytes | any level, when coding does not help |
| 4 | tokens | LZMA-style adaptive binary range coding of literals and matches | level 11 |
| 5 | sequences | canonical Huffman literals + tANS-coded (literal length, match length, offset) sequences | levels -3 to 10 |
| 6 | context sequences | as type 5, but each literal uses one of up to 8 Huffman tables, picked by the byte before it | levels 9 and 10 |

An optional x86 filter (section 8) can be applied to the content before compression. The decoder then undoes it after decompression.

## 2. Conventions

* The key words **MUST**, **MUST NOT**, **SHOULD** and **MAY** are used as in RFC 2119.
* Multi-byte integers are unsigned and **little-endian**: `u8`, `u32` and `u64` are 1, 4 and 8 bytes.
* `x >> n`, `x << n`, `&`, `|` and `^` are the usual bit operations on unsigned integers. `floor_log2(x)` (for x ≥ 1) is the index of the highest set bit of x.
* "Corrupt" means the decoder MUST reject the frame (section 10).
* Pseudocode uses Python-like syntax. Integer arithmetic is exact unless a width is given.

Constants used throughout:

| Name | Value | Meaning |
|---|---|---|
| `BLOCK_MAX` | 2^18 = 262144 | largest `raw_len` of a block |
| `WLOG_MIN`, `WLOG_MAX` | 18, 26 | range of `window_log` |
| `ML_MIN` | 2 | smallest match length (type 5/6) |
| `HUF_MAXLEN` | 11 | longest Huffman code |
| `FSE_LOG` | 9 | tANS table log (512 states) |
| `COMP_MAX` | 2^19 + 2^16 = 589824 | largest `comp_len` of a block |

## 3. Frame

### 3.1 Frame header

```
offset size field
0      4    magic        "LZRN" (4C 5A 52 4E)
4      1    version      7 (5 and 6: legacy, Appendix C)
5      1    flags
6      1    window_log   18..26
7      1    reserved     0
8      8    content_size u64, only if flags bit 1 is set
```

The header is 8 bytes long, or 16 if it carries `content_size`.

**version.** This document defines version 7, the 1.0 format. A 1.0 decoder MUST also read versions 5 and 6, which are older lzrans formats with a few differences listed in Appendix C. It MUST reject any other version as *unsupported*. A higher version means a newer format, which this decoder cannot read (section 12).

**flags.**

| Bit | Mask | Name | Meaning |
|---|---|---|---|
| 0 | 0x01 | ENCRYPTED | reserved for a future encryption envelope; never set in valid frames |
| 1 | 0x02 | SIZE | the header carries `content_size` |
| 2 | 0x04 | X86 | the content went through the x86 filter (section 8) |
| 3-7 | 0xF8 | | reserved for required features (section 12) |

A decoder MUST reject a frame as *unsupported* if any flag other than SIZE or X86 is set. That includes ENCRYPTED.

**window_log.** The window size is `W = 2^window_log`. A match MUST NOT reach more than `W` bytes back (section 5.4). A decoder MUST reject `window_log` outside 18..26 as corrupt. Decoders need about `W + BLOCK_MAX` bytes of history (section 11). An encoder SHOULD use the smallest `window_log` that covers the content when the size is known.

**reserved.** MUST be 0. Otherwise the frame is corrupt.

**content_size.** The exact number of content bytes. It MUST be less than 2^63; otherwise the frame is corrupt. When it is present, the decoder MUST reject the frame as corrupt if the blocks decode to more bytes than this, or to fewer by the end block.

### 3.2 Blocks

Each block starts with a 9-byte header:

```
offset size field
0      1    type byte
1      4    raw_len   u32: content bytes this block decodes to
5      4    comp_len  u32: bytes of payload that follow
9      comp_len  payload
```

#### 3.2.1 Type byte

| Bits | Meaning |
|---|---|
| bit 7 = 1 | **skippable block**. The decoder skips its payload. Bits 0-6 name the kind of data (section 12.2). |
| bit 7 = 0 | a data block: bits 0-5 are the **block type**, bit 6 is the **reset** bit (section 5.5) |

Block types: 0 stored, 2 end, 4 tokens, 5 sequences, 6 context sequences. Types 1 and 3 belonged to pre-release formats and are invalid. Types 7-63 are reserved for newer formats.

#### 3.2.2 Checks

These checks use only the header, so a decoder can make them before the payload arrives. In this order, a decoder MUST reject the frame:

1. as corrupt if `raw_len > BLOCK_MAX` or `comp_len > COMP_MAX`;
2. skippable block: as corrupt if `raw_len != 0`; otherwise it skips `comp_len` bytes, and nothing else changes;
3. as *unsupported* if the block type is 7..63;
4. as corrupt if the block type is 1 or 3;
5. as corrupt if:
   * the type is 2 (end) with the reset bit set, or `raw_len != 0`, or `comp_len != 4`;
   * the type is 0 and `comp_len != raw_len`;
   * the type is 4, 5 or 6 and `raw_len == 0` or `comp_len > raw_len`.

A data block's payload MUST decode to exactly `raw_len` bytes. Sections 6, 7 and 10 give the other checks.

**Type 0, stored.** The payload is `raw_len` content bytes. A stored block may be empty (`raw_len == 0`). A stored block does not change any decoder state except the history (section 5.2). (A stored block with the reset bit still makes a reset point.)

**Type 2, end.** The payload is a `u32`: the low 32 bits of the XXH64 checksum of the content (section 9). The frame ends right after this block. Any bytes after it are not part of the frame (section 10.3).

**Types 4, 5 and 6** are defined in sections 6 and 7.

### 3.3 Content and history

The blocks decode to the **decoded stream**. If the X86 flag is clear, the decoded stream is the content. If the X86 flag is set, the decoded stream is the *filtered* content, and the content is found by running the inverse filter (section 8) over the whole decoded stream.

Everything that refers to earlier output works on the decoded stream: match copies, the "previous byte" contexts of types 4 and 6, and the "match byte" of type 4. The content size and the checksum are about the content.

## 4. Bit streams

Types 5 and 6 use two kinds of bit stream. Both number bits LSB-first: bit `i` of a byte string `b` is `(b[i >> 3] >> (i & 7)) & 1`.

### 4.1 Forward bit streams

Huffman-coded literals (sections 6.3 and 6.6) use **forward** streams. Bits are read from bit 0 upward. A field of `n` bits read at bit position `p` has the value `sum(bit(p + j) << j for j in 0..n-1)`.

A forward stream has a known byte length `S`. When decoding has consumed `P` bits, `ceil(P / 8)` MUST equal `S`. That means the stream ends inside its last byte. The unused high bits of the last byte SHOULD be zero. Decoders do not check them. Reading past the end of the stream makes the block corrupt.

Huffman codes are written so that the first bit of a code (its most significant bit, as in DEFLATE) is the first bit read. A decoder may look ahead 11 bits and use a table.

### 4.2 Backward bit streams

The sequence section (section 6.7) uses a **backward** stream. The encoder writes it forwards (LSB-first) and then appends one `1` bit, the **end marker**. It pads to a byte boundary with zeros.

Decoding a stream of `S` bytes:

1. The last byte MUST be non-zero; otherwise the block is corrupt. Let `P = 8*(S-1) + floor_log2(last byte)`. That is the position of the end marker, and the number of bits left to read.
2. `read(n)`: if `n > P`, the block is corrupt. Otherwise `P -= n` and the result is the `n`-bit field that starts at bit position `P`, read as in 4.1. So bits come out in the reverse of the order the encoder wrote them. Within a field, the value is the one the encoder wrote. `read(0) = 0`.
3. After the last sequence, `P` MUST be exactly 0.

## 5. Decoder state and matches

### 5.1 History

The decoder keeps the decoded stream produced so far in this frame, or at least its last `W` bytes. Matches copy from it.

### 5.2 State carried across blocks

Some state lives for the whole frame. It starts at the values below when the frame starts and changes only as described:

| State | Initial value | Changed by |
|---|---|---|
| repeat offsets `rep[0..2]` | 1, 4, 8 | types 4, 5, 6 (sections 5.3, 7.4) |
| Huffman literal table | none | type 5, literal mode 2 (6.3) |
| context literal tables | none | type 6, mode 0 (6.6) |
| tANS tables for LL, ML, OF | none | types 5 and 6, table mode 0 (6.4) |
| type 4 model and state | section 7.2 | type 4 |

Stored blocks and end blocks change none of this, and skippable blocks change nothing at all. A reset point (5.5) sets all of it back to the initial values. Block types share state: a type 5 block after a type 4 block sees the repeat offsets the type 4 block left. Types 5 and 6 share the tANS tables. The type 5 Huffman table survives type 6 blocks, and the type 6 tables survive type 5 blocks.

### 5.3 Repeat offsets

All match-coding block types share the three most-recently-used match distances. When a match with distance `d` is decoded:

* new distance: `rep = [d, rep[0], rep[1]]`
* `rep[0]`: unchanged
* `rep[1]`: `rep = [rep[1], rep[0], rep[2]]`
* `rep[2]`: `rep = [rep[2], rep[0], rep[1]]`

The update rule does not depend on the literal length. (In zstd, it does.)

### 5.4 Executing a match

A match `(length L, distance D)` at decoded-stream position `o` copies `L` bytes, one at a time:

```
for j in range(L): out[o + j] = out[o + j - D]
```

When `D < L` the copy overlaps itself, which repeats the last `D` bytes. A match is corrupt if `D > o - seg` (it reaches before the start of the current segment, 5.5), if `D > W`, or if it runs past the block's `raw_len`. Matches never span block boundaries, but they may copy from any earlier block of the segment within the window.

### 5.5 Reset points and segments

A data block whose type byte has the reset bit set is a **reset point**. Before decoding it, the decoder:

* sets every piece of state in 5.2 back to its initial value (repeat offsets 1, 4, 8; no Huffman, context or tANS tables; a fresh type 4 model, state 0);
* sets `seg`, the start of the current **segment**, to the current decoded-stream position. At the start of the frame, `seg` is 0.

Inside a segment, nothing before `seg` is used: a match may not reach before it (5.4), and every "previous byte" or "match byte" before `seg` reads as 0 or is corrupt as stated where it is used (6.6, 7.2, 7.4). The type 4 position context counts from `seg` (7.2). So a decoder can decode a segment knowing only where it starts in the frame. Encoders do not have to emit reset points. A future seek index (a skippable block, 12.2) will list them.

## 6. Block types 5 and 6

### 6.1 Payload layout

```
nseq     u32     number of sequences
nlit     u32     number of literal bytes in the block
literals         type 5: literal section (6.3), present only if nlit > 0
                 type 6: context literal section (6.6), always present
if nseq > 0:
  ll_mode, ml_mode, of_mode table descriptions (6.4.3), in that order
  sequence bit stream (6.7): all remaining payload bytes
```

### 6.2 Counts

The block is corrupt unless all of these hold:

* `nlit <= raw_len`
* `nseq <= BLOCK_MAX / 2 + 1`
* `2*nseq + nlit <= raw_len`
* type 6 only: `nlit > 0`

If `nseq == 0`, the block holds only literals. Then `nlit` MUST equal `raw_len`, and the literal section MUST end exactly at the end of the payload.

### 6.3 Huffman-coded literals (type 5)

#### 6.3.1 Canonical codes

A Huffman table is given by a code length `len[s]` in 0..11 for each byte value `s`. Zero means the symbol has no code. The table is corrupt if any length is greater than 11, if no symbol has a code, or if the Kraft sum `K = sum(2^(11 - len[s]) for len[s] > 0)` is greater than 2^11.

Codes are assigned canonically, as in DEFLATE (RFC 1951, 3.2.2): shorter codes first, and symbols of the same length in increasing order.

```
bl_count[l] = number of symbols with length l   (bl_count[0] = 0)
code = 0
for l in 1..11: code = (code + bl_count[l-1]) << 1; next_code[l] = code
for s in 0..255 with len[s] > 0: code[s] = next_code[len[s]]; next_code[len[s]] += 1
```

A code `c` of length `l` is stored with its most significant bit first in the bit stream (section 4.1).

The code MUST be **complete**: `K == 2^11`. Otherwise the table is corrupt. (So a table needs at least two symbols. A block whose literals are all one byte value uses the run mode, 6.3.3.)

#### 6.3.2 Table description

```
last      u8                 largest symbol with a code (or larger)
lengths   ceil((last+1)/2) bytes: len[0..last] as 4-bit values,
          len[2i] in the low nibble of byte i, len[2i+1] in the high nibble
```

The bytes take `(last + 2) // 2` bytes after `last`. Symbols above `last` have length 0. If `last` is even, the high nibble of the final byte is ignored.

#### 6.3.3 Literal section

```
mode u8
  0  raw:       nlit bytes follow, used as they are
  1  run:       one byte follows; all nlit literals are that byte
  2  new table: table description (6.3.2), then the streams below;
                the table becomes the frame's Huffman literal table
  3  repeat:    the streams below, coded with the frame's Huffman literal table
                (corrupt if no mode 2 section came earlier in the frame)
modes 2 and 3:
  size0 u32, size1 u32, size2 u32, size3 u32    byte lengths of streams 0-3
  stream0 | stream1 | stream2 | stream3
```

Any other mode is corrupt. Modes 0 and 1 leave the frame's Huffman table as it was. In mode 2, if the table description is corrupt, the frame is corrupt.

The `nlit` literals are split into four runs. With `q = ceil(nlit / 4)`, stream `k` (k = 0..3) holds literals `[min(k*q, nlit), min((k+1)*q, nlit))`. Some streams may be empty. Each stream is a forward bit stream (4.1) of the byte length given, and it MUST end exactly at that length (4.1). The streams follow each other with no gap, and the literal section ends after stream 3. The section is corrupt if the four sizes add up to more than the rest of the payload.

Because every size is stored, a decoder knows where the sequence section starts before decoding any literal. It can decode the four literal streams and the sequences in parallel.

### 6.4 tANS tables

The sequence section codes three symbol streams, each with its own table:

| Stream | Alphabet size | Symbols |
|---|---|---|
| LL (literal length) | 62 | codes 0..61 (6.5) |
| ML (match length) | 58 | codes 0..57 (6.5) |
| OF (offset) | 65 | codes 0..64 (6.5) |

#### 6.4.1 Normalized counts

A table is given by counts `norm[0..A-1]` with `sum(norm) == 512`. They are written as a run of varints:

```
varint: 1 to 3 bytes, 7 bits each, least significant group first;
        the high bit of a byte is set if another byte follows.
        A third byte with its high bit set is corrupt.
while s < A:
  v = varint
  if v > 0:  norm[s] = v; s += 1
  else:      r = varint; norm[s .. s+r] = 0; s += r + 1   (r+1 zero counts)
```

The description is corrupt if a zero run goes past `A`, if the running sum goes above 512, or if the sum is not exactly 512 at `s == A`. Encoders SHOULD write varints in their shortest form. Decoders accept longer forms.

#### 6.4.2 Decoding table

The 512 states are built from `norm`:

```
step = 512/2 + 512/8 + 3            # 323
pos = 0
for s in 0..A-1:
  repeat norm[s] times: spread[pos] = s; pos = (pos + step) & 511
next = copy of norm
for u in 0..511:
  s = spread[u]; x = next[s]; next[s] += 1
  nb = 9 - floor_log2(x)
  state[u] = (symbol s, nbits nb, base (x << nb) - 512)
```

Decoding a symbol in state `u`: the symbol is `state[u].symbol`, and the next state is `state[u].base + read(state[u].nbits)`. There is no special "less than one" probability as in zstd: every symbol with a non-zero count has at least one state.

#### 6.4.3 Table modes

For each of LL, ML and OF in turn:

```
mode u8
  0  new:     normalized counts (6.4.1) follow; the table replaces the
              frame's table for this stream
  1  repeat:  use the frame's table for this stream (corrupt if none)
  2  single:  one u8 symbol s follows (corrupt if s >= A); this block
              uses a one-state table: symbol s, nbits 0, base 0, initial
              state read with 0 bits. The frame's table is NOT changed.
```

Any other mode is corrupt. A new table is built and stored even if the rest of the block turns out to be corrupt. That does not matter, because decoding stops at the first error.

### 6.5 Value codes

Lengths and offsets are coded as a symbol (the **code**) followed by extra bits. A family of codes with `2^d` **direct codes** maps code `c` to `(base, nbits)`:

```
if c < 2^d: base = c, nbits = 0
else:
  k = c - 2^d; n = d + k // 2
  base = (1 << n) | ((k & 1) << (n - 1)); nbits = n - 1
value = base + read(nbits)
```

| Field | d | Value |
|---|---|---|
| LL | 4 | literal length = value |
| ML | 5 | match length = value + 2 |
| OF, code < 3 | | repeat offset index `code` (5.3) |
| OF, code ≥ 3 | 4 (applied to `code - 3`) | distance = value + 1 |

For example, distance 32 gives value 31, which is code 17 of the d=4 family (base 24, 3 extra bits = 7), so the OF code is 20. The alphabets are larger than valid blocks need. LL ≤ 2^18 needs codes up to 44, and distance ≤ 2^26 needs OF codes up to 62. Decoders MUST cope with the larger codes. A larger code always produces a value that fails a later check, so the block is corrupt.

### 6.6 Context literal section (type 6)

```
mode u8
  0  new:  K u8 (1..8)
           cmap: 128 bytes, 256 nibbles: cmap[p] for p = 0..255, low nibble
                 first; every entry MUST be < K
           K table descriptions (6.3.2), table 0 first
  1  repeat: use the frame's context tables (corrupt if none)
size u32                     byte length of the literal stream
stream                       one forward bit stream (4.1) of `size` bytes
```

Any other mode, or `K` outside 1..8, is corrupt. Each table follows 6.3.1. Mode 0 replaces the frame's context tables.

Literals are decoded *while the sequences run* (6.7), in content order. Each literal is decoded with table `cmap[prev]`, where `prev` is the byte just before it in the decoded stream. At the start of a segment (5.5), and so of the frame, `prev` is 0. After a match, `prev` is the last byte the match copied. After all `nlit` literals, the stream MUST end as in 4.1.

### 6.7 Sequence section

The sequence section has `nseq` sequences. Each sequence is (literal length `LL`, match length `ML`, offset `OF`) and is executed in turn: copy `LL` literals to the output, then the match. All payload bytes after the table descriptions form one backward bit stream (4.2).

```
s_ll = read(t_ll); s_ml = read(t_ml); s_of = read(t_of)   # initial states
                              # t = 9, or 0 for a single-symbol table
lit = 0                       # literals used so far
repeat nseq times:
  c_ll, c_ml, c_of = symbols of states s_ll, s_ml, s_of
  ofv = OF value  (c_of < 3: ofv = c_of; else ofv = value(c_of - 3) + 3)
  # extra bits, in this order:
  read the OF extra bits, then the ML extra bits, then the LL extra bits
  # state updates, in this order:
  s_ll = next(s_ll); s_ml = next(s_ml); s_of = next(s_of)
  # (each next() reads that state's nbits)
  LL, ML from their codes and extra bits (6.5)
  distance: ofv >= 3: ofv - 3; else rep[ofv]; then update rep (5.3)
  check: LL <= nlit - lit and LL + ML <= bytes left in the block
  output LL literals (type 5: the next LL decoded literals; type 6: decode
     them now, 6.6); lit += LL
  execute the match (5.4)
after the loop:
  P == 0 (4.2)
  nlit - lit == bytes left in the block; output those literals
```

The states are updated after every sequence, including the last one. So the encoder writes the initial states and one set of state bits per sequence.

**Bit order inside one sequence.** OF extra bits are read first, then ML extra, then LL extra, then the LL, ML and OF state bits. In the order the encoder writes (forwards), that is reversed: OF state bits, ML state bits, LL state bits, LL extra, ML extra, OF extra. The encoder writes sequences last to first, and the initial states last.

## 7. Block type 4

A type 4 payload is one stream from a binary adaptive range coder. It uses the token model of LZMA, applied to this format's own state, so it is *not* an LZMA stream. The model, the 12-value state and the repeat offsets carry over from block to block (5.2).

### 7.1 Range decoder

Probabilities are 12-bit (`P_ONE = 4096`) estimates that the next bit is 0. They adapt with shift 5.

```
init:   first byte MUST be 0 (else corrupt); range = 0xFFFFFFFF;
        code = next 4 bytes, big-endian. The payload MUST be >= 5 bytes.
normalize: if range < 2^24: range <<= 8 (mod 2^32); code = (code << 8 | next byte) mod 2^32
bit(p):
  bound = (range >> 12) * p
  if code < bound: range = bound;  p += (4096 - p) >> 5; b = 0
  else: range -= bound; code -= bound; p -= p >> 5;   b = 1
  normalize; return b
direct(n):                       # n raw bits, most significant first
  repeat n times: range >>= 1; b = code >= range; if b: code -= range
                  v = 2v + b; normalize
tree(probs, n):                  # n bits, most significant first
  m = 1; repeat n times: m = 2m + bit(probs[m]); return m - 2^n
rtree(probs, n):                 # n bits, least significant first
  m = 1; v = 0; for i in 0..n-1: b = bit(probs[m]); m = 2m + b; v |= b << i
```

One normalization step per bit is always enough, because adaptation keeps every probability in [31, 4065]. When the block's `raw_len` bytes have been decoded, the decoder MUST have read exactly `comp_len` bytes: the 5 initial bytes plus one byte per normalization. Reading past the end of the payload, or ending before it, makes the block corrupt.

### 7.2 Model

All probabilities start at 2048. They are reset only at the start of a frame.

| Array | Size | Index |
|---|---|---|
| `is_match` | 12 × 4 | state, pos |
| `is_rep`, `is_rep0`, `is_rep1` | 12 each | state |
| `is_rep0_long` | 12 × 4 | state, pos |
| `lit` | 8 × 768 | prev >> 5, node |
| `slot` | 4 × 64 | min(len - 2, 3), node |
| `spec` | 14 × 32 | slot, node (one tree per slot 4..13) |
| `align` | 16 | node |
| `mlen`, `rlen` | each: `choice`, `choice2`, `low[4][8]`, `mid[4][8]`, `high[256]` | |

`pos` is the decoded-stream position **relative to the start of the segment** (5.5), modulo 4; without reset points, that is the content position. `prev` is the byte before the current position in the decoded stream, or 0 at the start of a segment. The state starts at 0.

### 7.3 State transitions

```
after literal:    s < 4 ? 0 : s < 10 ? s - 3 : s - 6
after match:      s < 7 ? 7 : 10
after rep match:  s < 7 ? 8 : 11
after short rep:  s < 7 ? 9 : 11
```

### 7.4 Tokens

Tokens are decoded until the block holds `raw_len` bytes:

```
if bit(is_match[s][pos]) == 0:                     # literal
  probs = lit[prev >> 5]
  if s < 7: byte = tree(probs, 8)
  else:                                            # "matched" literal
    (corrupt if rep[0] > bytes decoded so far in the segment)
    mbyte = out[o - rep[0]]; v = 1; same = true
    for i = 7 down to 0:
      mb = (mbyte >> i) & 1
      b = bit(probs[0x100 + (mb << 8) + v] if same else probs[v])
      v = 2v + b; same = same and b == mb
    byte = v & 0xFF
  output byte; s = after literal
elif bit(is_rep[s]) == 0:                          # match, new distance
  len = length(mlen)                               # 7.5
  d0 = distance(len)                               # 7.6; corrupt if d0 >= W
  rep = [d0 + 1, rep[0], rep[1]]; s = after match; copy (len, d0 + 1)
else:                                              # repeat match
  if bit(is_rep0[s]) == 0:
    if bit(is_rep0_long[s][pos]) == 0:             # short rep: 1 byte
      copy (1, rep[0]); s = after short rep; continue
    D = rep[0]
  elif bit(is_rep1[s]) == 0: D = rep[1]; rep = [rep[1], rep[0], rep[2]]
  else:                      D = rep[2]; rep = [rep[2], rep[0], rep[1]]
  len = length(rlen); s = after rep match; copy (len, D)
```

The copies follow 5.4, with the same checks. A match of more than 273 bytes is coded as several tokens. The first is a new-distance or repeat match, and the rest are rep0 matches.

### 7.5 Lengths

```
length(l):
  if bit(l.choice) == 0:  return 2 + tree(l.low[pos], 3)      # 2..9
  if bit(l.choice2) == 0: return 10 + tree(l.mid[pos], 3)     # 10..17
  return 18 + tree(l.high, 8)                                 # 18..273
```

### 7.6 Distances

```
distance(len):
  slot = tree(slot[min(len - 2, 3)], 6)
  if slot < 4: return slot
  nb = (slot >> 1) - 1
  d0 = (2 | (slot & 1)) << nb
  if slot < 14: return d0 + rtree(spec[slot], nb)
  return d0 + (direct(nb - 4) << 4) + rtree(align, 4)
```

The result `d0` is the distance minus 1. The block is corrupt if `d0 >= W`.

## 8. x86 filter

When the X86 flag is set, the content was transformed before compression. The decoder applies the inverse to the whole decoded stream. Only the rel32 operand of x86 `CALL` (opcode E8) is changed, and only when it is small.

```
unfilter(b):                       # b: decoded stream, positions mod 2^32
  i = 0
  while i < len(b):
    if b[i] != 0xE8: i += 1; continue
    if len(b) - i < 5: break       # a truncated operand at the end is left as is
    v = u32 at b[i+1..i+4]
    if (v >> 25) in (0, 0x7F):     # v is a 26-bit signed value
      a = (v - (i + 5)) mod 2^26
      if a & 2^25: a |= 0xFC000000   # sign-extend from bit 25
      write a as u32 at b[i+1..i+4]
    i += 5                         # skip the operand, converted or not
```

The forward filter is the same, with `v + (i + 5)`. Both directions test the bytes they read, and the conversion maps small values to small values one-to-one. So the inverse is exact for any input. A streaming decoder may hold back up to 4 bytes after an E8 until it sees more data or the end of the frame.

## 9. Checksum

The end block holds `XXH64(content, seed 0) & 0xFFFFFFFF` as a `u32`. XXH64 is the 64-bit xxHash algorithm (https://github.com/Cyan4973/xxHash, `doc/xxhash_spec.md`), the same function zstd uses for its frame checksum. It is computed over the content, which means after the inverse x86 filter. A mismatch is a *checksum* error.

A frame has no other checksum: not for the header, not per block. Streaming decoders deliver content before they can verify it (section 10.2).

## 10. Errors and conformance

### 10.1 Error kinds

| Kind | When |
|---|---|
| **unsupported** | version not 5, 6 or 7; a reserved flag set; block type 7-63 (version 7) |
| **corrupt** | anything else this document calls corrupt, including a truncated frame |
| **checksum** | the end block's checksum does not match |

Errors are fatal for the frame. Partial output MUST NOT be presented as the complete content. A decoder MUST NOT read or write outside its buffers on any input, valid or not. The golden vectors include invalid frames (section 13).

### 10.2 Streaming

Every check in this document can be made with at most one block of payload buffered. A streaming decoder MAY deliver a block's output once the block has been checked. It MUST report a later error, including a checksum mismatch, when that error is found. Applications that need verified data before using it should buffer until the end block.

### 10.3 Trailing data and concatenation

A frame ends after its end block. What follows is outside the frame. The lzrans one-shot `decompress()` rejects trailing bytes, and the streaming `Decompressor` exposes them as `unused_data`. This document does not define a meaning for concatenated frames (Appendix A, item A12).

### 10.4 Conformance

* A **conforming decoder** decodes every valid frame to its content. It rejects every frame this document calls corrupt, unsupported or a checksum failure, with the kind given. Implementations may merge "corrupt" and "unsupported" into one error, as lzrans' Python API does today, but SHOULD keep them apart.
* A **conforming encoder** produces only valid frames. It need not use every feature.
* Encoders MUST NOT depend on padding bits being ignored.

## 11. Limits and resources

| Limit | Value |
|---|---|
| content size | < 2^63 bytes |
| block content | ≤ 2^18 bytes |
| block payload | ≤ `raw_len` (types 4-6); ≤ 589824 bytes (any block) |
| match distance | ≤ 2^window_log ≤ 64 MiB |
| sequences per block | ≤ 131073 |
| match length | type 5/6: ≤ 2^18; type 4 token: ≤ 273 |

A decoder needs about `2^window_log + 2^18` bytes of history, plus one block of input. It can decode straight into the output buffer when the whole output is in memory, unless the X86 flag is set. Tables take under 64 KiB, and the type 4 model about 15 KiB.

## 12. Versioning and extension

The rules below let the format grow without changing what a 1.0 decoder accepts.

**Version byte.** 7 is LZR 1.0. A change that old decoders must not misread (a different meaning for existing bytes) needs a new version number. Decoders reject versions they do not know as *unsupported*. Versions 5 and 6 stay readable (Appendix C).

**Required features.** A new block type that changes the content (for example the long-range dedup layer being prototyped) takes a number from 7-63 and a flag bit from 3-7 that frames using it MUST set. A 1.0 decoder then stops at the header with *unsupported*, before writing any output. If a decoder still meets an unknown type 7-63 (a frame without the flag), it reports *unsupported*, not *corrupt*, so users can tell "needs a newer decoder" from "damaged".

**Skippable blocks.** A type byte of 0x80 or more marks a skippable block (3.2.1). Its `raw_len` is 0, and decoders that do not know it skip `comp_len` bytes. Skippable blocks are for data that does not change the content: a seek index, a dedup manifest, encoder information. Kinds (bits 0-6) are assigned as they are defined:

| Type byte | Kind |
|---|---|
| 0x80 | seek index (reserved; format to be defined with the parallel decoder) |
| 0x81-0xFF | unassigned |

A decoder MUST NOT change the decoded content because of a skippable block it knows. It may use one, for example to find reset points.

## 13. Test vectors

`tests/vectors/` holds frozen vectors, listed in `tests/vectors/manifest.json` with the expected outcome. For each valid frame it gives the size and SHA-256 of the content and the format features it uses. For each invalid frame it gives the error kind.

* `enc-*`: written by the C encoder (all block types, every table mode, x86, window logs 18-26, streaming and one-shot headers).
* `reset-splice-*`, `skippable-blocks`: two independent encoder frames joined at a reset point (types 4, 5 and 6), and skippable blocks between data blocks.
* `hand-*`: built byte by byte by `tools/make_vectors.py`, for cases the encoder never writes. These are the initial repeat offsets and their reset, single-symbol tables with an empty bit stream, literal-only blocks, an empty stored block, a hand-coded Huffman section, a stored block between type 5 blocks that share tables, and a type 5 block that uses repeat offsets left by a type 4 block.
* `bad-*`: invalid frames, one fault each, with the error kind: header fields, type bytes (unknown types are *unsupported*), reset and skippable misuse, block limits, a wrong checksum or content size, truncation, trailing data, bit streams not consumed exactly, malformed tANS or Huffman tables (including incomplete codes and a wrong stream 3 size), oversized extra-bit reads, and malformed type 4 payloads.
* `v5-*`, `v6-*`: frames of the legacy versions, frozen when the format was version 6.

`tests/test_vectors.py` checks every vector with the C decoder (one-shot and streaming, including the error kind) and with `tools/refdec.py`. It also checks that the version 7 vectors cover every feature in this document. The vectors are frozen: `tools/make_vectors.py` adds new ones but keeps the committed bytes of existing ones.

## 14. Encoder notes (non-normative)

* A block is stored whenever its coded form would not be smaller than `raw_len`, which also keeps `comp_len <= raw_len` (3.2.2).
* The encoder does not emit reset points or skippable blocks yet.
* Blocks hold 256 KiB of input. The last block, and the halves made by cost-based block splitting (levels 7-10), are shorter.
* Level 11 uses only types 4 and 0. Levels -3 to 8 use types 5 and 0. Levels 9 and 10 use types 5, 6 and 0, choosing type 6 for a block only when its literal section is at least 2% smaller.
* With threads at levels -3 to 6, each 4 MiB job starts with fresh tables (mode 0) and does not count on the repeat offsets. The format is the same either way.

---

## Appendix A. Decisions taken for the freeze

Before 1.0, the draft of this document listed format problems found while writing it. The user chose to fix the recommended ones and freeze as version 7 (2026-11-28).

| # | Problem in version 6 | Decision in version 7 |
|---|---|---|
| A1 | No independent points: every block could depend on everything before it, which rules out parallel decoding and seeking. | **Fixed.** Reset bit in the type byte (5.5). Seek index reserved as skippable kind 0x80. |
| A2 | No extension mechanism: unknown block types were "corrupt", and no block could be skipped. | **Fixed.** Type byte layout, skippable blocks, *unsupported* for types 7-63, feature flags (3.2.1, 12). |
| A3 | Length of literal stream 3 not stored, so literals had to be decoded before the sequences could be found. | **Fixed.** `size3` stored (6.3.3). |
| A4 | Block size fixed at 256 KiB, with the ML alphabet tied to it. | **Kept** at 2^18. Larger blocks would need a new version. |
| A5 | Incomplete Huffman codes accepted ("symbol 0, 11 bits" for unused patterns). | **Fixed.** Codes must be complete (6.3.1). The encoder now always writes complete codes. |
| A6 | Only a 32-bit checksum, at the end. | **Kept**, same as zstd. A per-block checksum could come later as a feature. |
| A7 | Type 4 `pos` context relative to the block. | **Fixed.** Relative to the segment (7.2). Same result for everything the old encoder wrote. |
| A8 | `comp_len` limited only by an arbitrary constant. | **Fixed.** `comp_len <= raw_len` for types 4-6 (3.2.2). |
| A9 | Encryption flag reserved with no meaning. | **Kept** reserved. Encryption belongs in an outer envelope. |
| A10 | Version number of the frozen format. | **7**. Versions 5 and 6 stay readable. |
| A11 | LL and OF alphabets larger than needed. A crafted block could make the C decoder read 66 bits at once. | **Kept**; decoder fixed, vector `bad-extra-bits-66`. |
| A12 | Meaning of concatenated frames not defined. | **Open.** To be settled with the CLI (10.3). |

## Appendix B. Differences from zstd and LZMA (informative)

* Repeat offsets update the same way whatever the literal length (5.3). zstd shifts them when LL = 0.
* tANS tables have no "less than one" probability. Counts are written as varints rather than zstd's bit-packed form. The state updates after the last sequence are still read (6.7).
* Type 4 is LZMA's coding layer with these changes: lc = 3, lp = 0, pb = 2, but `pos` counts from the segment start; one reverse tree per distance slot (not a shared array); initial repeat offsets 1, 4, 8; and the model persists across blocks up to the next reset point.

## Appendix C. Legacy frame versions 5 and 6

A 1.0 decoder reads versions 5 and 6 with these differences from version 7:

* The type byte is the block type itself. Types 7-255 are *corrupt* (not unsupported, and not skippable), and there is no reset bit, so the whole frame is one segment.
* Block type 6 is valid only in version 6. Version 5 has no other difference from version 6.
* `comp_len` for types 4-6 may be larger than `raw_len` (up to `COMP_MAX`).
* Type 5 literal section, modes 2 and 3: only `size0`, `size1` and `size2` are stored (12 bytes). Stream 3 starts after stream 2 and ends at the byte where its last literal ends. The literal section ends there, and the section is corrupt if stream 3 would read past the payload.
* Huffman codes may be incomplete. A decoder MUST behave as if every 11-bit pattern that matches no code decodes to symbol 0 and uses 11 bits. (Valid streams never contain one.)
* Type 4 `pos` (7.2) is the position relative to the start of the **block**, modulo 4.
