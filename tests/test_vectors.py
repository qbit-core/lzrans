"""Golden test vectors (tests/vectors/, docs/FORMAT.md section 13).
Names starting with v5- or v6- are frames of those legacy versions.

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
    maxw = e.get("max_window_log", 27)
    if e["expect"] == "ok":
        _check_content(e, lzrans.decompress(blob, max_window_log=maxw))
    elif e["expect"] == "checksum":
        with pytest.raises(lzrans.ChecksumError):
            lzrans.decompress(blob, max_window_log=maxw)
    else:
        with pytest.raises(lzrans.LzransError) as ei:
            lzrans.decompress(blob, max_window_log=maxw)
        assert not isinstance(ei.value, lzrans.ChecksumError)
        assert ("unsupported" in str(ei.value)) == (e["expect"] == "unsupported")


@pytest.mark.parametrize("e", [e for e in MANIFEST if e["expect"] == "ok"],
                         ids=[e["file"] for e in MANIFEST if e["expect"] == "ok"])
def test_c_decoder_streaming(e):
    blob = _load(e)
    d = lzrans.Decompressor(max_window_log=e.get("max_window_log", 27))
    out = b"".join(d.decompress(blob[i:i + 777]) for i in range(0, len(blob), 777))
    assert d.eof
    _check_content(e, out)


@pytest.mark.parametrize("e", MANIFEST, ids=[e["file"] for e in MANIFEST])
def test_c_decoder_parallel(e):
    # decompress(threads=...) must give the same result and error kind
    blob = _load(e)
    maxw = e.get("max_window_log", 27)
    if e["expect"] == "ok":
        _check_content(e, lzrans.decompress(blob, threads=3, max_window_log=maxw))
        return
    exc = lzrans.ChecksumError if e["expect"] == "checksum" else lzrans.LzransError
    with pytest.raises(exc) as ei:
        lzrans.decompress(blob, threads=3, max_window_log=maxw)
    if e["expect"] != "checksum":
        assert not isinstance(ei.value, lzrans.ChecksumError)
        assert ("unsupported" in str(ei.value)) == (e["expect"] == "unsupported")


def test_seek_index_vectors():
    feats = {e["file"]: set(e.get("features", ())) for e in MANIFEST}
    for name in ("seg-l11-type4.lzr", "seg-l3-type5.lzr"):
        assert {"reset", "seek_index"} <= feats[name], name
        idx = refdec.read_seek_index(_load({"file": name}))
        assert idx is not None and len(idx[0]) == 3, name
    for name in ("seg-index-bad-offset.lzr", "seg-index-bad-notlast.lzr", "skippable-blocks.lzr"):
        assert "seek_index" not in feats[name], name


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
    # version 7 vectors alone cover every feature; legacy v5-/v6- vectors
    # cover the older versions
    feats, legacy = set(), set()
    for e in MANIFEST:
        (legacy if e["file"].startswith(("v5-", "v6-")) else feats).update(e.get("features", ()))
    assert {"version5", "version6", "block6", "t5.lit.mode2"} <= legacy
    need = {"version7", "reset", "block.skippable", "seek_index", "flag.size", "flag.x86", "wlog18", "wlog26",
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
