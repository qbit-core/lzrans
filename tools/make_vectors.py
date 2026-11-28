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
def frame(blocks, version=6, flags=0, wlog=20, content_size=None, checksum_of=None, bad_checksum=False):
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

    # ---- encoder output
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
    small = words(20000, 15)
    b = bytearray(lzrans.compress(small, level=1))
    assert 6 not in [t for t, _, _ in split_blocks(bytes(b))]
    b[4] = 5
    V.append(("enc-version5", bytes(b), "ok", "frame version 5 (still decoded)", small))

    # ---- encoder blocks with a hand-inserted stored block in between:
    # the decoder state (tables, repeat offsets, type 4 model) survives it.
    # The content differs from the encoder input (later matches now reach
    # into the stored bytes), so it is taken from the reference decoder.
    # (Not possible for type 4: its literal contexts depend on the history,
    # so changing the history changes how the rest decodes.)
    for name, lvl, desc in (("hand-stored-between-type5", -3, "type 5 tables repeated across a stored block"),):
        bl = split_blocks(lzrans.compress(mutated(BS + 30000, 21, every=100), level=lvl))
        assert len(bl) == 2 and bl[0][0] == bl[1][0] == (5 if lvl < 11 else 4)
        bl.insert(1, (0, 100, noise(100, 22)))
        cont = _decode_nocheck(*bl)
        V.append((name, frame(bl, checksum_of=cont), "ok", desc, cont))

    # ---- hand-built
    # one sequence, every table single-symbol (mode 2): literal 'a', then
    # a match of 10 at repeat offset 0, whose initial value is 1
    c = b"a" * 11
    p = seq5(1, b"a", [RLE(1), RLE(8), RLE(0)], b"\x01")
    V.append(("hand-initial-rep0", frame([(5, 11, p)], checksum_of=c), "ok",
              "initial repeat offset 0 = 1; all tables mode 2; bit stream = end marker only", c))
    # initial rep2 = 8 (block 1), then rep1 after the update (block 2):
    # reps [1,4,8] -> ov 2 -> [8,1,4] -> ov 1 -> [1,8,4]
    c = b"abcdefgh" * 2 + b"h" * 6
    p1 = seq5(1, b"abcdefgh", [RLE(8), RLE(6), RLE(2)], b"\x01")
    p2 = seq5(1, b"", [RLE(0), RLE(4), RLE(1)], b"\x01")
    V.append(("hand-initial-rep2-rep1", frame([(5, 16, p1), (5, 6, p2)], checksum_of=c), "ok",
              "initial repeat offset 2 = 8, then the repeat offset update rule", c))
    # literal-only blocks (nseq = 0): raw and run literal modes, an empty stored block
    c = b"hello, world" + b"z" * 40 + b"!"
    V.append(("hand-literal-only", frame([(5, 12, seq5(0, b"hello, world", [], b"")),
                                          (0, 0, b""),
                                          (5, 40, seq5(0, b"z" * 40, [], b"", lit_mode=1)),
                                          (0, 1, b"!")], checksum_of=c), "ok",
              "nseq = 0 blocks with raw and run literals; an empty stored block", c))
    # a new offset with extra bits (FORMAT.md 6.5): literals "abcd" * 8,
    # then distance 32: OF value 31 -> code 3 + 17 = 20 with 3 extra bits (7);
    # LL 32 -> code 18 with 4 extra bits (0); ML 32 -> code 30.
    # Bits, read from the top: end marker, OF extra 111, LL extra 0000.
    c = b"abcd" * 16
    p = seq5(1, b"abcd" * 8, [RLE(18), RLE(30), RLE(20)], b"\xf0")
    V.append(("hand-new-offset", frame([(5, 64, p)], checksum_of=c), "ok",
              "a new offset and a literal length with extra bits", c))
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
    bad("bad-version7", patch(4, 7), "unsupported", "frame version 7 (newer than this spec)")
    bad("bad-flag-encrypted", patch(5, good[5] | 1), "unsupported", "flag bit 0 (encryption, reserved)")
    bad("bad-flag-bit3", patch(5, good[5] | 8), "unsupported", "unknown flag bit 3")
    bad("bad-wlog17", patch(6, 17), "corrupt", "window_log below 18")
    bad("bad-wlog27", patch(6, 27), "corrupt", "window_log above 26")
    bad("bad-reserved", patch(7, 1), "corrupt", "reserved header byte not zero")
    for t in (1, 3, 7, 255):
        bad("bad-block-type%d" % t, patch(hl, t), "corrupt", "block type %d" % t)
    b6 = bytearray(lzrans.compress(markov(30000, 18), level=10))
    assert b6[header_len(b6)] == 6
    b6[4] = 5
    bad("bad-type6-in-v5", b6, "corrupt", "type 6 block in a version 5 frame")
    c = noise(100, 19)
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
    bad("bad-dist-before-start", frame([(5, 2, seq5(1, b"", [RLE(0), RLE(0), RLE(0)], b"\x01"))],
                                       checksum_of=b"\0\0"), "corrupt", "match before the start of the content")
    bad("bad-bits-left", frame([(5, 11, seq5(1, b"a", [RLE(1), RLE(8), RLE(0)], b"\x03"))],
                               checksum_of=b"a" * 11), "corrupt", "sequence bit stream not consumed exactly")
    bad("bad-no-end-marker", frame([(5, 11, seq5(1, b"a", [RLE(1), RLE(8), RLE(0)], b"\x01\x00"))],
                                   checksum_of=b"a" * 11), "corrupt", "sequence bit stream: last byte zero")
    bad("bad-table-repeat", frame([(5, 11, seq5(1, b"a", [b"\x01", RLE(8), RLE(0)], b"\x01"))],
                                  checksum_of=b"a" * 11), "corrupt", "tANS repeat mode without a table")
    bad("bad-table-sum", frame([(5, 11, seq5(1, b"a", [b"\x00\x80\x02" + b"\x00\x3c", RLE(8), RLE(0)],
                                             b"\x01\x01"))], checksum_of=b"a" * 11), "corrupt",
        "tANS counts summing to 256")
    bad("bad-huf-oversub", frame([(5, 4, struct.pack("<II", 0, 4) + bytes([2, 2, 0x11, 0x01]) + b"\0" * 13)],
                                 checksum_of=b"\0" * 4), "corrupt", "Huffman lengths 1,1,1 (oversubscribed)")
    bad("bad-huf-repeat", frame([(5, 4, struct.pack("<II", 0, 4) + bytes([3]) + b"\0" * 13)],
                                checksum_of=b"\0" * 4), "corrupt", "Huffman repeat mode without a table")
    bad("bad-huf-len12", frame([(5, 4, struct.pack("<II", 0, 4) + bytes([2, 1, 0xC1]) + b"\0" * 13)],
                               checksum_of=b"\0" * 4), "corrupt", "Huffman code length 12")
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


def _decode_nocheck(*blks):
    """Content of a frame made of blks (checksum not checked)."""
    fr = refdec.Frame(6, 20, 0, None)
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


def main():
    force = "--force" in sys.argv
    os.makedirs(OUT, exist_ok=True)
    manifest, allfeat = [], set()
    for name, blob, expect, desc, content in build():
        feats = set()
        rres, rout = ref_decode(blob, feats)
        cres, cout = c_decode(blob)
        assert rres == expect, (name, rres, expect)
        # the C decoder reports unsupported versions/flags through the same
        # exception class; only the "ok" / "checksum" / error split is compared
        assert (cres == "ok") == (expect == "ok") and (cres == "checksum") == (expect == "checksum"), (name, cres)
        entry = {"file": name + ".lzr", "expect": expect, "description": desc}
        if expect == "ok":
            assert rout == cout, name
            if content is not None:
                assert rout == content, name
            entry.update(size=len(rout), sha256=hashlib.sha256(rout).hexdigest(), features=sorted(feats))
            allfeat |= feats
        path = os.path.join(OUT, name + ".lzr")
        if os.path.exists(path) and open(path, "rb").read() != blob and not force:
            # frozen: keep the committed bytes (they still have to pass)
            blob = open(path, "rb").read()
            print("kept frozen", name)
        open(path, "wb").write(blob)
        manifest.append(entry)
    with open(os.path.join(OUT, "manifest.json"), "w") as f:
        json.dump({"format": "LZR 1.0-draft", "vectors": manifest}, f, indent=1)
        f.write("\n")
    print(len(manifest), "vectors;", sum(os.path.getsize(os.path.join(OUT, e["file"])) for e in manifest), "bytes")
    print("features covered:", " ".join(sorted(allfeat)))


if __name__ == "__main__":
    main()
