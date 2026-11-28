#!/usr/bin/env python3
"""make_vectors.py - (re)generate the golden test vectors in tests/vectors/.

The vectors are FROZEN once committed: tests decode the stored .lzr files
and compare with the stored hashes, so encoder changes cannot move them.
Run this script only to add vectors (it refuses to overwrite an existing
vector whose bytes would change, unless --force).

Three kinds of vectors:
  enc-*   produced by the C encoder from deterministic inputs
  hand-*  built byte by byte here, for format features the encoder never
          emits (initial repeat offsets, literal-only blocks, type 4 and
          type 5 blocks in one frame, ...)
  bad-*   invalid frames that a conforming decoder must reject

Every positive vector must decode identically with the C decoder and with
tools/refdec.py (written from docs/FORMAT.md); the script checks that.

    python tools/make_vectors.py [--force]
"""
import hashlib
import json
import os
import random
import struct
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import refdec  # noqa: E402
import lzrans  # noqa: E402

OUT = os.path.join(HERE, "..", "tests", "vectors")
BS = 1 << 18


# ---------------------------------------------------------------- inputs
def words(n, seed=0, vocab=300):
    rnd = random.Random(seed)
    ws = [bytes(rnd.choices(b"abcdefghijklmnopqrstuvwxyz", k=rnd.randint(2, 9))) for _ in range(vocab)]
    out = bytearray()
    while len(out) < n:
        out += rnd.choice(ws) + rnd.choice([b" ", b" ", b" ", b", ", b".\n"])
    return bytes(out[:n])


def mutated(n, seed=0, base=16384, every=400):
    """words() repeated with a random byte changed about every `every`
    bytes: long matches and few literals, so multi-block vectors stay small."""
    rnd = random.Random(seed)
    b = words(base, seed)
    out = bytearray()
    while len(out) < n:
        out += b
    for i in range(base, n - every, every):
        out[i + rnd.randrange(every // 2)] = rnd.randrange(32, 127)
    return bytes(out[:n])


def markov(n, seed=3):
    rnd = random.Random(seed)
    nxt = {a: rnd.sample(range(32, 127), 5) for a in range(256)}
    out, b = bytearray(), 32
    for _ in range(n):
        b = rnd.choice(nxt[b]) if rnd.random() < 0.93 else rnd.randrange(32, 127)
        out.append(b)
    return bytes(out)


def noise(n, seed):
    return random.Random(seed).randbytes(n)


def x86_like(n, seed=3):
    rnd = random.Random(seed)
    out = bytearray()
    targets = [rnd.randrange(0, n) for _ in range(200)]
    while len(out) < n:
        k = rnd.random()
        if k < 0.3:
            out += b"\xe8" + ((rnd.choice(targets) - len(out) - 5) & 0xFFFFFFFF).to_bytes(4, "little")
        elif k < 0.33:
            out += b"\xe8\xe8\xe8"
        else:
            out += bytes(rnd.choice(b"\x48\x89\x8b\xc3\x0f\x00\x01\xff") for _ in range(rnd.randint(1, 6)))
    return bytes(out[:n]) + b"\xe8\x01\x02"   # truncated operand at the end


def records(n):
    out = bytearray()
    i = 0
    while len(out) < n:
        out += b"id=%06d;name=item%d;flags=%s;\n" % (i * 7919 % 1000003, i % 97, b"ab" * (i % 5))
        i += 1
    return bytes(out[:n])


# ---------------------------------------------------------------- frames
def frame(blocks, version=7, flags=0, wlog=20, content_size=None, checksum_of=None, bad_checksum=False):
    """Assemble a frame from (type, raw_len, payload) blocks; the end block's
    checksum is computed over checksum_of (the content)."""
    h = bytearray(b"LZRN" + bytes([version, flags | (2 if content_size is not None else 0), wlog, 0]))
    if content_size is not None:
        h += struct.pack("<Q", content_size)
    for t, raw, p in blocks:
        h += bytes([t]) + struct.pack("<II", raw, len(p)) + p
    c = refdec.xxh64(checksum_of or b"") & 0xFFFFFFFF
    if bad_checksum:
        c ^= 1
    h += bytes([2]) + struct.pack("<II", 0, 4) + struct.pack("<I", c)
    return bytes(h)


def seq5(ns, lits, modes, bits, lit_mode=0):
    """Type 5 payload: counts, literal section (raw or run), table modes
    (list of byte strings), backward bit stream bytes."""
    p = struct.pack("<II", ns, len(lits))
    if lits:
        p += bytes([0]) + lits if lit_mode == 0 else bytes([1, lits[0]])
    return p + b"".join(modes) + bits


RLE = lambda sym: bytes([2, sym])   # noqa: E731  (table mode 2)


def split_blocks(blob):
    """[(type, raw, payload)] of an encoder frame, without the end block."""
    out = []
    for ev in refdec.blocks(blob):
        if ev[0] == "block" and ev[1] != 2:
            out.append(ev[1:4])
    return out


def header_len(blob):
    return refdec.parse_header(blob)[4]


# ---------------------------------------------------------------- vectors
def build():
    V = []   # (name, blob, expect, description, content or None)

    def enc(name, data, desc, **kw):
        V.append((name, lzrans.compress(data, **kw), "ok", desc, data))

    def stream(name, data, desc, level, **kw):
        c = lzrans.Compressor(level=level, **kw)
        blob = b"".join(c.compress(data[i:i + 100000]) for i in range(0, len(data), 100000)) + c.flush()
        V.append((name, blob, "ok", desc, data))

    # ---- encoder output (format version 7)
    enc("enc-empty", b"", "empty content: header + end block", level=6)
    enc("enc-stored", noise(3000, 1), "incompressible: one stored block", level=6)
    enc("enc-l1-text", words(40000, 1), "level 1, type 5: Huffman literals (mode 2), new tANS tables", level=1)
    enc("enc-l-3-multiblock", mutated(2 * BS + 999, 2), "level -3, three type 5 blocks: table repeat modes",
        level=-3)
    enc("enc-l6-records", records(BS // 2), "level 6, type 5, repeat offsets 1 and 2", level=6)
    enc("enc-l10-ctx", markov(60000), "level 10, type 6 context literals, new and repeated tables", level=10)
    enc("enc-l11-text", mutated(2 * BS + 77, 7, every=150), "level 11, type 4 over three blocks (model carries over)",
        level=11)
    enc("enc-l11-zeros", bytes(BS + 100), "level 11, long rep0 runs (matches split at 273)", level=11)
    enc("enc-l3-zeros", bytes(300000), "level 3, matches of 2^18 - 2 bytes", level=3)
    enc("enc-x86-l6", x86_like(20000), "x86 filter flag, type 5", level=6, filter="x86")
    enc("enc-x86-l11", x86_like(20000, 4), "x86 filter flag, type 4", level=11, filter="x86")
    enc("enc-wlog18", mutated(BS + BS // 3, 10), "window_log 18 (smallest window)", level=4, window_log=18)
    stream("enc-stream-l5", mutated(BS + 4321, 13), "streaming compressor: no content size in the header", 5)
    stream("enc-stream-l9", markov(30000, 14), "streaming compressor, level 9", 9)
    stream("enc-wlog26", words(30000, 11) * 2, "window_log 26 (largest window), level 11, streaming", 11)

    # ---- reset points (5.5): two independent encoder frames spliced into
    # one, the second's first block marked as a reset point. Valid only if
    # the second part decodes with no state from the first. The first part
    # has a length that is not a multiple of 4 (type 4 positions count from
    # the segment start).
    for name, lvl, a, b in (("reset-splice-type5", 1, words(5001, 30), mutated(BS + 999, 31)),
                            ("reset-splice-type6", 10, markov(5001, 32), markov(25000, 33)),
                            ("reset-splice-type4", 11, words(5001, 34), mutated(BS + 999, 35))):
        ba = split_blocks(lzrans.compress(a, level=lvl, window_log=20))
        bb = split_blocks(lzrans.compress(b, level=lvl, window_log=20))
        t0, r0, p0 = bb[0]
        bb[0] = (t0 | 0x40, r0, p0)
        V.append((name, frame(ba + bb, checksum_of=a + b), "ok",
                  "level %d frames spliced at a reset point" % lvl, a + b))
    # skippable blocks (12.2) between data blocks
    a = mutated(BS + 999, 36)
    bl = split_blocks(lzrans.compress(a, level=1))
    bl.insert(1, (0x80, 0, b"opaque metadata the decoder skips"))
    bl.insert(0, (0xBF, 0, b""))
    V.append(("skippable-blocks", frame(bl, checksum_of=a), "ok", "skippable blocks 0x80 and 0xBF (empty)", a))

    # ---- segment mode (14): reset points + seek index (12.3), and 0x80
    # blocks that are not a valid index, which decoders ignore
    a = mutated(2 * (1 << 20) + 5000, 40, every=2000)
    enc("seg-l11-type4", a, "segment mode, level 11: three segments (reset points) and a seek index",
        level=11, segment=1 << 20)
    enc("seg-l3-type5", a, "segment mode, level 3: three segments and a seek index", level=3, segment=1 << 20)
    good = lzrans.compress(a, level=3, segment=1 << 20)
    h = header_len(good)
    bl = [(ev[1], ev[2], ev[3]) if ev[0] == "block" else (ev[1], 0, ev[3])   # data blocks and the index block
          for ev in refdec.blocks(good) if ev[0] == "skip" or (ev[0] == "block" and ev[1] != 2)]
    bl = [(t | (0x40 if r else 0), raw, p) if t < 0x80 else (t, raw, p)
          for (t, raw, p), r in zip(bl, [ev[5] if ev[0] == "block" else False
                                         for ev in refdec.blocks(good) if ev[0] == "skip" or (ev[0] == "block" and ev[1] != 2)])]
    assert bl[-1][0] == 0x80
    idx = bytearray(bl[-1][2])
    idx[24] ^= 1                               # content_offset of entry 1
    V.append(("seg-index-bad-offset", frame(bl[:-1] + [(0x80, 0, bytes(idx))], content_size=len(a), checksum_of=a),
              "ok", "0x80 block whose entry 1 has a wrong content offset: not a valid index, ignored", a))
    V.append(("seg-index-bad-notlast", frame(bl + [(0, 0, b"")], content_size=len(a), checksum_of=a),
              "ok", "seek index followed by an empty stored block: not a valid index, ignored", a))
    assert h == 16

    # ---- encoder blocks with a hand-inserted stored block in between:
    # the decoder state (tables, repeat offsets) survives it. The content
    # differs from the encoder input (later matches now reach into the
    # stored bytes), so it is taken from the reference decoder. (Not
    # possible for type 4: its literal contexts depend on the history.)
    bl = split_blocks(lzrans.compress(mutated(BS + 30000, 21, every=100), level=-3))
    assert len(bl) == 2 and bl[0][0] == bl[1][0] == 5
    bl.insert(1, (0, 100, noise(100, 22)))
    cont = _decode_nocheck(*bl)
    V.append(("hand-stored-between-type5", frame(bl, checksum_of=cont), "ok",
              "type 5 tables repeated across a stored block", cont))

    # ---- hand-built. Version 7 needs comp_len <= raw_len for types 4-6,
    # so the blocks below are long enough to pay for their headers.
    # One sequence, every table single-symbol (mode 2): literal 'a', then a
    # match of 33 (ML code 31) at repeat offset 0, whose initial value is 1.
    c = b"a" * 34
    p = seq5(1, b"a", [RLE(1), RLE(31), RLE(0)], b"\x01")
    V.append(("hand-initial-rep0", frame([(5, 34, p)], checksum_of=c), "ok",
              "initial repeat offset 0 = 1; all tables mode 2; bit stream = end marker only", c))
    # initial rep2 = 8 (block 1), then rep1 (block 2):
    # reps [1,4,8] -> OF 2 -> [8,1,4] -> OF 1 -> [1,8,4]
    c = b"abcdefgh" * 3 + b"h" * 22
    p1 = seq5(1, b"abcdefgh", [RLE(8), RLE(14), RLE(2)], b"\x01")
    p2 = seq5(1, b"", [RLE(0), RLE(20), RLE(1)], b"\x01")
    V.append(("hand-initial-rep2-rep1", frame([(5, 24, p1), (5, 22, p2)], checksum_of=c), "ok",
              "initial repeat offset 2 = 8, then the repeat offset update rule", c))
    # reset point: block 2 uses repeat offset 0, which is 1 again only
    # because of the reset (it would be 8 without it)
    c2 = b"q" * 22
    p3 = seq5(1, b"q", [RLE(1), RLE(19), RLE(0)], b"\x01")
    V.append(("hand-reset-reps", frame([(5, 24, p1), (0x45, 22, p3)], checksum_of=b"abcdefgh" * 3 + c2), "ok",
              "a reset point restores the initial repeat offsets", b"abcdefgh" * 3 + c2))
    # literal-only blocks (nseq = 0): run literals, an empty stored block
    c = b"z" * 40 + b"!"
    V.append(("hand-literal-only", frame([(0, 0, b""), (5, 40, seq5(0, b"z" * 40, [], b"", lit_mode=1)),
                                          (0, 1, b"!")], checksum_of=c), "ok",
              "nseq = 0 block with run literals; an empty stored block", c))
    # a new offset with extra bits (6.5): literals "abcd" * 8, then distance
    # 32: OF value 31 -> code 3 + 17 = 20 with 3 extra bits (7); LL 32 ->
    # code 18 with 4 extra bits (0); ML 32 -> code 30. Bits, read from the
    # top: end marker, OF extra 111, LL extra 0000.
    c = b"abcd" * 16
    p = seq5(1, b"abcd" * 8, [RLE(18), RLE(30), RLE(20)], b"\xf0")
    V.append(("hand-new-offset", frame([(5, 64, p)], checksum_of=c), "ok",
              "a new offset and a literal length with extra bits", c))
    # Huffman literals with all 4 stream sizes (6.3.3): 64 literals 0/1,
    # lengths 1,1 (codes 0 and 1), 16 literals = 2 bytes per stream
    rnd = random.Random(40)
    c = bytes(rnd.randrange(2) for _ in range(64))
    streams = [huf_bits(c[16 * i:16 * i + 16], {0: "0", 1: "1"}) for i in range(4)]
    p = struct.pack("<II", 0, 64) + bytes([2, 1, 0x11]) + struct.pack("<4I", *map(len, streams)) + b"".join(streams)
    V.append(("hand-huffman", frame([(5, 64, p)], checksum_of=c), "ok", "Huffman literal section, 4 stream sizes", c))
    huf_ok = p
    # type 4 block from the encoder, then a hand-built type 5 block that
    # uses the repeat offsets the type 4 block left behind
    src = words(3000, 16)
    t4 = split_blocks(lzrans.compress(src, level=11))
    assert [t for t, _, _ in t4] == [4]
    p = seq5(1, b"", [RLE(0), RLE(18), RLE(0)], b"\x01")   # LL 0, ML 20, rep0
    fr_content = _decode_nocheck(t4[0], (5, 20, p))
    V.append(("hand-type4-then-type5", frame([t4[0], (5, 20, p)], checksum_of=fr_content), "ok",
              "repeat offsets carry over from a type 4 block to a type 5 block", fr_content))

    # ---- invalid frames
    good = lzrans.compress(words(5000, 17), level=6)
    hl = header_len(good)

    def bad(name, blob, expect, desc):
        V.append((name, bytes(blob), expect, desc, None))

    def patch(off, val):
        b = bytearray(good)
        b[off] = val
        return b

    bad("bad-magic", patch(0, ord("X")), "corrupt", "wrong magic")
    bad("bad-version4", patch(4, 4), "unsupported", "frame version 4")
    bad("bad-version8", patch(4, 8), "unsupported", "frame version 8 (newer than this spec)")
    bad("bad-flag-encrypted", patch(5, good[5] | 1), "unsupported", "flag bit 0 (encryption, reserved)")
    bad("bad-flag-bit3", patch(5, good[5] | 8), "unsupported", "unknown flag bit 3")
    bad("bad-wlog17", patch(6, 17), "corrupt", "window_log below 18")
    bad("bad-wlog27", patch(6, 27), "corrupt", "window_log above 26")
    bad("bad-reserved", patch(7, 1), "corrupt", "reserved header byte not zero")
    for t in (1, 3, 0x41, 0x43):
        bad("bad-block-type-0x%02x" % t, patch(hl, t), "corrupt", "block type byte 0x%02x" % t)
    for t in (7, 63, 0x47):
        bad("bad-block-type-0x%02x" % t, patch(hl, t), "unsupported",
            "block type byte 0x%02x: a block type of a newer format" % t)
    c = noise(100, 19)
    bad("bad-reset-on-end", frame([(0, 100, c)], checksum_of=c)[:-13] + bytes([0x42]) + struct.pack("<II", 0, 4)
        + struct.pack("<I", refdec.xxh64(c) & 0xFFFFFFFF), "corrupt", "reset bit on the end block")
    bad("bad-skippable-raw", frame([(0x80, 1, b"x"), (0, 100, c)], checksum_of=c), "corrupt",
        "skippable block with raw_len 1")
    bad("bad-match-before-reset", frame([(5, 24, p1), (0x45, 22, seq5(1, b"", [RLE(0), RLE(20), RLE(0)], b"\x01"))],
                                        checksum_of=b""), "corrupt", "a match reaching before a reset point")
    bad("bad-comp-gt-raw", frame([(5, 12, seq5(0, b"hello, world", [], b""))], checksum_of=b"hello, world"),
        "corrupt", "type 5 block with comp_len > raw_len (valid in version 6)")
    bad("bad-raw-too-big", frame([(0, BS + 1, b"")], checksum_of=b"")[:-13], "corrupt",
        "stored block header with raw_len 2^18 + 1 (rejected before its payload is needed)")
    bad("bad-stored-len", frame([(0, 100, c + b"x")], checksum_of=c), "corrupt", "stored block, comp_len != raw_len")
    bad("bad-checksum", frame([(0, 100, c)], checksum_of=c, bad_checksum=True), "checksum", "checksum mismatch")
    bad("bad-size-short", frame([(0, 100, c)], content_size=101, checksum_of=c), "corrupt",
        "content_size larger than the content")
    bad("bad-size-long", frame([(0, 100, c)], content_size=99, checksum_of=c), "corrupt",
        "content_size smaller than the content")
    bad("bad-truncated", good[:-9], "corrupt", "no end block")
    bad("bad-trailing", good + b"\0", "corrupt", "a byte after the end block (one-shot decoding)")
    bad("bad-end-len", frame([], checksum_of=b"")[:-13] + bytes([2]) + struct.pack("<II", 0, 5) + b"\0" * 5,
        "corrupt", "end block with comp_len 5")
    bad("bad-empty-type5", frame([(5, 0, b"")], checksum_of=b""), "corrupt", "type 5 block with raw_len 0")
    bad("bad-dist-before-start", frame([(5, 22, seq5(1, b"", [RLE(0), RLE(20), RLE(0)], b"\x01"))],
                                       checksum_of=b"\0" * 22), "corrupt", "match before the start of the content")
    a34 = b"a" * 34
    bad("bad-bits-left", frame([(5, 34, seq5(1, b"a", [RLE(1), RLE(31), RLE(0)], b"\x03"))], checksum_of=a34),
        "corrupt", "sequence bit stream not consumed exactly")
    bad("bad-no-end-marker", frame([(5, 34, seq5(1, b"a", [RLE(1), RLE(31), RLE(0)], b"\x01\x00"))],
                                   checksum_of=a34), "corrupt", "sequence bit stream: last byte zero")
    bad("bad-table-repeat", frame([(5, 34, seq5(1, b"a", [b"\x01", RLE(31), RLE(0)], b"\x01"))], checksum_of=a34),
        "corrupt", "tANS repeat mode without a table")
    bad("bad-table-sum", frame([(5, 34, seq5(1, b"a", [b"\x00\x80\x02" + b"\x00\x3c", RLE(31), RLE(0)],
                                             b"\x01\x01"))], checksum_of=a34), "corrupt", "tANS counts summing to 256")
    z40 = b"\0" * 40
    bad("bad-huf-oversub", frame([(5, 40, struct.pack("<II", 0, 40) + bytes([2, 2, 0x11, 0x01]) + b"\0" * 20)],
                                 checksum_of=z40), "corrupt", "Huffman lengths 1,1,1 (oversubscribed)")
    bad("bad-huf-repeat", frame([(5, 40, struct.pack("<II", 0, 40) + bytes([3]) + b"\0" * 20)],
                                checksum_of=z40), "corrupt", "Huffman repeat mode without a table")
    bad("bad-huf-len12", frame([(5, 40, struct.pack("<II", 0, 40) + bytes([2, 1, 0xC1]) + b"\0" * 20)],
                               checksum_of=z40), "corrupt", "Huffman code length 12")
    # incomplete code (lengths 1, 2): valid in version 6, corrupt in 7
    streams = [huf_bits(c64, {0: "0", 1: "10"}) for c64 in [bytes(16)] * 3 + [b"\x01" * 16]]
    p = struct.pack("<II", 0, 64) + bytes([2, 1, 0x21]) + struct.pack("<4I", *map(len, streams)) + b"".join(streams)
    bad("bad-huf-incomplete", frame([(5, 64, p)], checksum_of=bytes(48) + b"\x01" * 16), "corrupt",
        "incomplete Huffman code (lengths 1 and 2)")
    q = bytearray(huf_ok)
    q[8 + 3 + 12] += 1          # size of stream 3, one byte too large
    bad("bad-huf-size3", frame([(5, 64, bytes(q) + b"\0")], checksum_of=c), "corrupt",
        "stream 3 does not end at its stored size")
    bad("bad-extra-bits-66", frame([(5, 1000, seq5(1, b"", [RLE(61), RLE(57), RLE(64)], b"\xff" * 16))],
                                   checksum_of=b""), "corrupt",
        "codes with 25 + 16 + 25 extra bits (C decoder read > 64 bits at once before the fix)")
    t4b = bytearray(t4[0][2])
    t4b[0] = 1
    bad("bad-type4-first-byte", frame([(4, t4[0][1], bytes(t4b))], checksum_of=src), "corrupt",
        "type 4 payload whose first byte is not 0")
    bad("bad-type4-short", frame([(4, t4[0][1], t4[0][2][:-1])], checksum_of=src), "corrupt",
        "type 4 payload one byte short")
    return V


def huf_bits(syms, codes):
    """Forward LSB-first stream of the given codes (strings, first bit first)."""
    bits = "".join(codes[s] for s in syms)
    out = bytearray((len(bits) + 7) // 8)
    for i, b in enumerate(bits):
        out[i >> 3] |= (b == "1") << (i & 7)
    return bytes(out)


def _decode_nocheck(*blks):
    """Content of a version 7 frame made of blks (no reset points)."""
    fr = refdec.Frame(7, 20, 0, None)
    for t, raw, p in blks:
        if t == 0:
            fr.out += p
        elif t == 4:
            refdec.decode_tok_block(p, raw, fr)
        else:
            refdec.decode_seq_block(p, raw, fr, t == 6)
    return bytes(fr.out)


def c_decode(blob):
    try:
        return "ok", lzrans.decompress(blob)
    except lzrans.ChecksumError:
        return "checksum", None
    except lzrans.LzransError as e:
        return ("unsupported" if "unsupported" in str(e).lower() else "corrupt"), None


def ref_decode(blob, feats):
    try:
        return "ok", refdec.decode(blob, features=feats)
    except refdec.ChecksumError:
        return "checksum", None
    except refdec.UnsupportedError:
        return "unsupported", None
    except refdec.FormatError:
        return "corrupt", None


def check(name, blob, expect, content=None):
    """Both decoders must give `expect`; returns the manifest fields."""
    feats = set()
    rres, rout = ref_decode(blob, feats)
    cres, cout = c_decode(blob)
    assert rres == expect, (name, rres, expect)
    assert cres == expect, (name, cres, expect)
    entry = {}
    if expect == "ok":
        assert rout == cout, name
        if content is not None:
            assert rout == content, name
        entry.update(size=len(rout), sha256=hashlib.sha256(rout).hexdigest(), features=sorted(feats))
    return entry


def main():
    force = "--force" in sys.argv
    os.makedirs(OUT, exist_ok=True)
    mpath = os.path.join(OUT, "manifest.json")
    old = json.load(open(mpath))["vectors"] if os.path.exists(mpath) else []
    manifest, allfeat = [], set()
    # legacy vectors (frame versions 5 and 6): frozen, re-checked only
    for e in old:
        if e["file"].startswith(("v5-", "v6-")):
            blob = open(os.path.join(OUT, e["file"]), "rb").read()
            assert check(e["file"], blob, e["expect"]).get("sha256") == e.get("sha256"), e["file"]
            manifest.append(e)
    for name, blob, expect, desc, content in build():
        path = os.path.join(OUT, name + ".lzr")
        if os.path.exists(path) and open(path, "rb").read() != blob and not force:
            blob = open(path, "rb").read()      # frozen: keep the committed bytes
            content = None
            print("kept frozen", name)
        entry = {"file": name + ".lzr", "expect": expect, "description": desc}
        entry.update(check(name, blob, expect, content))
        allfeat |= set(entry.get("features", ()))
        open(path, "wb").write(blob)
        manifest.append(entry)
    with open(mpath, "w") as f:
        json.dump({"format": "LZR 1.0 (frame version 7; v5-/v6- vectors: legacy versions)",
                   "vectors": manifest}, f, indent=1)
        f.write("\n")
    print(len(manifest), "vectors;", sum(os.path.getsize(os.path.join(OUT, e["file"])) for e in manifest), "bytes")
    print("version 7 features covered:", " ".join(sorted(allfeat)))


if __name__ == "__main__":
    main()
