# The LZR compressed data format

**Version 1.0-draft** (describes frame versions 5 and 6, as written and read by lzrans at commit `10f73c8`)

**Status:** draft for the 1.0 format freeze. Everything in sections 3 to 11 is what the current decoder (`src/lzrans/dec.c`) accepts, byte for byte. Two parts are not settled yet:

* **Appendix A** lists format problems worth fixing *before* the freeze. Some of the fixes change the bitstream, so they need a decision.
* **Section 12** describes the extension mechanism proposed for 1.0. Another branch of work (a long-range dedup layer) may add a block type. Section 12 reserves room for that.

Until this document drops the "-draft" label, a new lzrans release may still change the format. After the freeze, the format may only change as section 12 allows.

A second decoder, `tools/refdec.py` (pure Python, about 700 lines), was written from this document alone. It decodes every golden vector in `tests/vectors/` (section 13) and everything the current encoder produces at levels -3 to 11, with the same result as the C decoder. Treat any difference between this text, `refdec.py` and `dec.c` as a bug in this document.

---

## 1. Overview

An LZR **frame** holds one stream of bytes, the **content**. A frame has three parts:

```
frame  = header  block*  end-block
```

* The **header** (section 3.1) carries a magic number, the frame version, flags, and the window size. It can also carry the content size.
* **Blocks** (section 3.2) each decode to at most 256 KiB of content. Blocks are decoded in order. A block may copy bytes from earlier blocks, up to the window size back. A block may also reuse coding tables or model state left by earlier blocks.
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
4      1    version      5 or 6
5      1    flags
6      1    window_log   18..26
7      1    reserved     0
8      8    content_size u64, only if flags bit 1 is set
```

The header is 8 bytes long, or 16 if it carries `content_size`.

**version.** This document defines versions 5 and 6. They differ in one way: block type 6 is only valid in version 6 frames. A decoder MUST reject any other version as *unsupported*. A decoder for this document MUST accept both 5 and 6.

**flags.**

| Bit | Mask | Name | Meaning |
|---|---|---|---|
| 0 | 0x01 | ENCRYPTED | reserved for a future encryption envelope; never set in valid frames |
| 1 | 0x02 | SIZE | the header carries `content_size` |
| 2 | 0x04 | X86 | the content went through the x86 filter (section 8) |
| 3-7 | 0xF8 | | reserved |

A decoder MUST reject a frame as *unsupported* if any flag other than SIZE or X86 is set. That includes ENCRYPTED.

**window_log.** The window size is `W = 2^window_log`. A match MUST NOT reach more than `W` bytes back (section 5.4). A decoder MUST reject `window_log` outside 18..26 as corrupt. Decoders need about `W + BLOCK_MAX` bytes of history (section 11). An encoder SHOULD use the smallest `window_log` that covers the content when the size is known.

**reserved.** MUST be 0. Otherwise the frame is corrupt.

**content_size.** The exact number of content bytes. It MUST be less than 2^63; otherwise the frame is corrupt. When it is present, the decoder MUST reject the frame as corrupt if the blocks decode to more bytes than this, or to fewer by the end block.

### 3.2 Blocks

Each block starts with a 9-byte header:

```
offset size field
0      1    type      0, 2, 4, 5 or 6
1      4    raw_len   u32: content bytes this block decodes to
5      4    comp_len  u32: bytes of payload that follow
9      comp_len  payload
```

A decoder MUST reject the frame as corrupt if any of these hold:

* `type` is 1, 3, or greater than 6. Types 1 and 3 belonged to pre-release formats. For 7 and up, see section 12.
* `type` is 6 and the frame version is 5.
* `raw_len > BLOCK_MAX` or `comp_len > COMP_MAX`.
* `type` is 0 and `comp_len != raw_len`.
* `type` is 2 and (`raw_len != 0` or `comp_len != 4`).
* `type` is 4, 5 or 6 and `raw_len == 0`.

These checks use only the header, so a decoder can make them before the payload arrives. A block's payload MUST decode to exactly `raw_len` bytes. Section 10 gives the other per-type checks.

**Type 0, stored.** The payload is `raw_len` content bytes. A stored block may be empty (`raw_len == 0`). A stored block does not change any decoder state except the history (section 5.2).

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

Stored blocks and end blocks change none of this. Block types share state: a type 5 block after a type 4 block sees the repeat offsets the type 4 block left. Types 5 and 6 share the tANS tables. The type 5 Huffman table survives type 6 blocks, and the type 6 tables survive type 5 blocks.

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

When `D < L` the copy overlaps itself, which repeats the last `D` bytes. A match is corrupt if `D > o` (it reaches before the start of the frame's decoded stream), if `D > W`, or if it runs past the block's `raw_len`. Matches never span block boundaries, but they may copy from any earlier block within the window.

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

**Incomplete codes (K < 2^11).** These are currently allowed. A decoder MUST behave as if every 11-bit pattern with no matching code decodes to symbol 0 and uses 11 bits. A valid stream never contains such a pattern. Appendix A, item A5, proposes to forbid incomplete codes.

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
  size0 u32, size1 u32, size2 u32    byte lengths of streams 0-2
  stream0 | stream1 | stream2 | stream3
```

Any other mode is corrupt. Modes 0 and 1 leave the frame's Huffman table as it was. In mode 2, if the table description is corrupt, the frame is corrupt.

The `nlit` literals are split into four runs. With `q = ceil(nlit / 4)`, stream `k` (k = 0..3) holds literals `[min(k*q, nlit), min((k+1)*q, nlit))`. Some streams may be empty. Each stream is a forward bit stream (4.1). Streams 0-2 have the byte lengths given, and each MUST end exactly at its length (4.1). Stream 3 starts right after stream 2. **Its length is not stored:** it ends at the byte where its last literal ends, and the literal section ends there too. The section is corrupt if `size0 + size1 + size2` reaches past the payload, or if stream 3 would read past the payload. (Appendix A, item A3, proposes to store the length.)

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

Literals are decoded *while the sequences run* (6.7), in content order. Each literal is decoded with table `cmap[prev]`, where `prev` is the byte just before it in the decoded stream. At the start of the frame, `prev` is 0. After a match, `prev` is the last byte the match copied. After all `nlit` literals, the stream MUST end as in 4.1.

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

`pos` is the decoded-stream position **relative to the start of the block**, modulo 4. (Appendix A, item A7: LZMA uses the absolute position.) `prev` is the byte before the current position in the decoded stream, or 0 at the start of the frame. The state starts at 0.

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
    (corrupt if rep[0] > bytes decoded so far)
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
| **unsupported** | version not 5 or 6; a reserved flag set |
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
* Encoders MUST NOT depend on the incomplete-code rule in 6.3.1, or on padding bits being ignored.

## 11. Limits and resources

| Limit | Value |
|---|---|
| content size | < 2^63 bytes |
| block content | ≤ 2^18 bytes |
| block payload | ≤ 589824 bytes |
| match distance | ≤ 2^window_log ≤ 64 MiB |
| sequences per block | ≤ 131073 |
| match length | type 5/6: ≤ 2^18; type 4 token: ≤ 273 |

A decoder needs about `2^window_log + 2^18` bytes of history, plus one block of input. It can decode straight into the output buffer when the whole output is in memory, unless the X86 flag is set. Tables take under 64 KiB, and the type 4 model about 15 KiB.

## 12. Versioning and extension

### 12.1 Today (versions 5 and 6)

Today, any new block type or flag needs a new frame version. A decoder rejects unknown versions and flags as *unsupported*, and unknown block types as *corrupt*. There is no way to add optional data that old decoders can skip.

### 12.2 Proposed for 1.0 (not implemented)

This proposal is part of the freeze decision (Appendix A, items A1 and A2):

* **Frame version.** The frozen format gets one version byte value, called `V1` below (6, or 7 if bitstream fixes are made). Decoders keep reading 5 and 6. Versions above `V1` are *unsupported*.
* **Block type ranges.**

  | Types | Use |
  |---|---|
  | 0-6 | defined here (1 and 3 stay invalid) |
  | 7-63 | future **required** block types. A decoder that does not know one MUST stop with *unsupported*, not *corrupt*, so users can tell "needs a newer decoder" from "damaged". New required types must be announced by a feature flag or a version bump, so old decoders fail at the header, not halfway through the content. |
  | 64-127 | reserved |
  | 128-191 | **skippable** blocks: `raw_len` MUST be 0, and decoders that do not know the type skip `comp_len` bytes. For metadata that does not change the content: a seek index, a dedup manifest, the encoder name. |
  | 192-255 | reserved |

* **Feature flags.** Flag bits 3-7 stay reserved. A future required feature (for example the long-range dedup layer being prototyped now) sets a flag bit. Decoders without it reject the frame at the header.
* **Reset points** (A1): a way to mark a block as starting a new independent segment, so segments can be decoded in parallel and a seek index can point into them.

## 13. Test vectors

`tests/vectors/` holds 55 frozen vectors, listed in `tests/vectors/manifest.json` with the expected outcome. For each valid frame it gives the size and SHA-256 of the content and the format features it uses. For each invalid frame it gives the error kind.

* `enc-*`: written by the C encoder (all block types, every table mode, x86, window logs 18-26, streaming and one-shot headers, version 5).
* `hand-*`: built byte by byte by `tools/make_vectors.py`, for cases the encoder never writes. These are the initial repeat offsets, single-symbol tables with an empty bit stream, literal-only blocks, an empty stored block, a stored block between type 5 blocks that share tables, and a type 5 block that uses repeat offsets left by a type 4 block.
* `bad-*`: invalid frames, one fault each: header fields, block header limits, a wrong checksum or content size, truncation, trailing data, bit streams not consumed exactly, malformed tANS or Huffman tables, oversized extra-bit reads, and malformed type 4 payloads.

`tests/test_vectors.py` checks every vector with the C decoder (one-shot and streaming) and with `tools/refdec.py`. It also checks that together the vectors cover every feature listed in this document. The vectors are frozen: `tools/make_vectors.py` adds new ones but keeps the committed bytes of existing ones.

## 14. Encoder notes (non-normative)

* A block is stored whenever its coded form would not be smaller than `raw_len`. So in today's frames, `comp_len < raw_len` for types 4-6.
* Blocks hold 256 KiB of input. The last block, and the halves made by cost-based block splitting (levels 7-10), are shorter.
* Level 11 uses only types 4 and 0. Levels -3 to 8 use types 5 and 0. Levels 9 and 10 use types 5, 6 and 0, choosing type 6 for a block only when its literal section is at least 2% smaller.
* With threads at levels -3 to 6, each 4 MiB job starts with fresh tables (mode 0) and does not count on the repeat offsets. The format is the same either way.

---

## Appendix A. Format problems to settle before the freeze

Ordered by importance. "Breaking" means frames written with the fix can't be read by today's decoder, so the fix needs a new version byte (A10). Today's decoder can keep reading versions 5 and 6.

| # | Problem | Proposal | Breaking? |
|---|---|---|---|
| **A1** | **No independent points.** Every block may depend on all earlier blocks: history, Huffman/tANS tables, repeat offsets, type 4 model. So a frame can't be decoded in parallel or entered in the middle. That rules out the "xz-class ratio + parallel decompression" goal without a format change. | Add a *reset* marker: bit 7 of the block type byte (or a dedicated empty block type). It means "decoder state goes back to frame-start values, and no match reaches before this block." Add an optional seek index as a skippable block (12.2). It costs nothing when unused, and encoders can add parallel-decodable output later without another format change. | yes (new bit) |
| **A2** | **No extension mechanism.** Unknown block types are "corrupt". There are no skippable blocks. | Block-type ranges, feature flags, and an *unsupported* error for unknown required types (12.2). The dedup layer then fits as a new required type behind a flag. | rules only |
| **A3** | **The length of literal stream 3 is not stored** (6.3.3). The sequence section's start is only known after all literals are decoded. That forces literal and sequence decoding to run one after the other, and blocks the planned parallel entropy decoding. | Store `size3` (+4 bytes per block, about 0.002%). | yes |
| **A4** | **Block size is fixed at 256 KiB, and the ML alphabet (58 codes) is tied to it.** Phase 18 measured table cost per 256 KiB block as what holds back the context literal gain. | Decide now. Either keep 2^18 forever, or allow up to 2^20 to 2^22 and extend the ML alphabet. The second costs a bigger decoder buffer and higher streaming latency. Recommendation: keep 2^18 unless a measurement shows a clear gain. | if changed |
| **A5** | **Incomplete Huffman codes are accepted** and decode as "symbol 0, 11 bits" (6.3.1). That is an artifact of the table-based decoder that every other decoder must copy. | Require `K == 2^11` when the table is read: no cost in the hot loop. Every table the encoder wrote on the sampled data was complete (91 of 91), but a single-symbol table would be incomplete. That case is already coded as a run (mode 1). | yes (stricter) |
| A6 | Only a 32-bit checksum, only at the end. No header check, no per-block check. | Keep it (same as zstd), or add an optional per-block checksum flag. Low priority. | if added |
| A7 | **The type 4 `pos` context is relative to the block** (7.2). It matches the content position only if every earlier block's `raw_len` is a multiple of 4. | Redefine it as the content position mod 4. Frames from today's encoder decode the same, because level 11 blocks are always 256 KiB except the last. So this is a free cleanup. | no for existing frames |
| A8 | `COMP_MAX` (589824) is an arbitrary sanity limit. The encoder always keeps `comp_len < raw_len` for coded blocks. | Tighten to `comp_len <= raw_len` for types 4-6. That bounds input buffering by the block size. | yes (stricter) |
| A9 | Encryption flag (bit 0) is reserved with no meaning. | Keep it reserved ("MUST be 0"). Do encryption as an outer envelope, not inside the frame. | no |
| A10 | Version number for the frozen format. | No breaking fixes: freeze as version 6. With A1/A3/A5/A8: version 7 = "LZR 1.0". Decoders keep reading 5 and 6. | |
| A11 | The LL (62) and OF (65) alphabets are larger than valid blocks need. Before this draft, a crafted block could make the C decoder read 66 bits at once (undefined behavior in C). That is fixed in `dec.c` and covered by vector `bad-extra-bits-66`. | Keep the alphabets: changing them buys nothing. | no |
| A12 | The meaning of concatenated frames is not defined. | Define it like zstd and xz: a stream of frames decodes to the concatenation. Decide together with the CLI. | no |

## Appendix B. Differences from zstd and LZMA (informative)

* Repeat offsets update the same way whatever the literal length (5.3). zstd shifts them when LL = 0.
* tANS tables have no "less than one" probability. Counts are written as varints rather than zstd's bit-packed form. The state updates after the last sequence are still read (6.7).
* Type 4 is LZMA's coding layer with these changes: lc = 3, lp = 0, pb = 2, but `pos` is block-relative; one reverse tree per distance slot (not a shared array); initial repeat offsets 1, 4, 8; and the model persists across blocks within a frame.
