"""Golden test vectors (tests/vectors/, docs/FORMAT.md section 13).

Each vector is a frozen .lzr file with its expected outcome: "ok" (with the
size and SHA-256 of the content), "corrupt", "unsupported" or "checksum".
Both the C decoder and the pure-Python reference decoder (tools/refdec.py,
written from the specification) must give that outcome.
"""
import hashlib
import json
import os
import sys

import pytest

import lzrans

HERE = os.path.dirname(os.path.abspath(__file__))
VEC = os.path.join(HERE, "vectors")
sys.path.insert(0, os.path.join(HERE, "..", "tools"))
import refdec  # noqa: E402

MANIFEST = json.load(open(os.path.join(VEC, "manifest.json")))["vectors"]


def _load(e):
    with open(os.path.join(VEC, e["file"]), "rb") as f:
        return f.read()


def _check_content(e, out):
    assert len(out) == e["size"]
    assert hashlib.sha256(out).hexdigest() == e["sha256"]


@pytest.mark.parametrize("e", MANIFEST, ids=[e["file"] for e in MANIFEST])
def test_c_decoder_one_shot(e):
    blob = _load(e)
    if e["expect"] == "ok":
        _check_content(e, lzrans.decompress(blob))
    elif e["expect"] == "checksum":
        with pytest.raises(lzrans.ChecksumError):
            lzrans.decompress(blob)
    else:
        with pytest.raises(lzrans.LzransError) as ei:
            lzrans.decompress(blob)
        assert not isinstance(ei.value, lzrans.ChecksumError)


@pytest.mark.parametrize("e", [e for e in MANIFEST if e["expect"] == "ok"],
                         ids=[e["file"] for e in MANIFEST if e["expect"] == "ok"])
def test_c_decoder_streaming(e):
    blob = _load(e)
    d = lzrans.Decompressor()
    out = b"".join(d.decompress(blob[i:i + 777]) for i in range(0, len(blob), 777))
    assert d.eof
    _check_content(e, out)


@pytest.mark.parametrize("e", MANIFEST, ids=[e["file"] for e in MANIFEST])
def test_reference_decoder(e):
    blob = _load(e)
    exc = {"corrupt": refdec.FormatError, "unsupported": refdec.UnsupportedError,
           "checksum": refdec.ChecksumError}
    if e["expect"] == "ok":
        _check_content(e, refdec.decode(blob))
        return
    with pytest.raises(exc[e["expect"]]) as ei:
        refdec.decode(blob)
    if e["expect"] == "corrupt":   # not one of the more specific kinds
        assert not isinstance(ei.value, (refdec.UnsupportedError, refdec.ChecksumError))


def test_vectors_cover_the_format():
    feats = set()
    for e in MANIFEST:
        feats |= set(e.get("features", ()))
    need = {"version5", "version6", "flag.size", "flag.x86", "wlog18", "wlog26",
            "block0", "block4", "block5", "block6", "seq.none",
            "t5.lit.mode0", "t5.lit.mode1", "t5.lit.mode2", "t5.lit.mode3", "t6.ctx.mode0", "t6.ctx.mode1",
            "seq.of.new", "seq.of.rep0", "seq.of.rep1", "seq.of.rep2",
            "t4.lit.plain", "t4.lit.matched", "t4.shortrep", "t4.rep1", "t4.rep2",
            "t4.slot.direct", "t4.slot.tree", "t4.slot.align", "t4.len.low", "t4.len.mid", "t4.len.high"}
    need |= {"seq.%s.mode%d" % (s, m) for s in ("ll", "ml", "of") for m in (0, 1, 2)}
    assert need <= feats, need - feats
    assert {e["expect"] for e in MANIFEST} == {"ok", "corrupt", "unsupported", "checksum"}


def test_reference_decoder_matches_encoder():
    # beyond the frozen vectors: today's encoder output at every level
    import random
    rnd = random.Random(5)
    words = [bytes(rnd.choices(b"abcdefghij klmnop", k=rnd.randint(2, 8))) for _ in range(200)]
    data = b" ".join(rnd.choice(words) for _ in range(9000))
    for level in range(lzrans.LEVEL_MIN if hasattr(lzrans, "LEVEL_MIN") else -3, 12):
        blob = lzrans.compress(data, level=level)
        assert refdec.decode(blob) == data, level
