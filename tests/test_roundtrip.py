"""Round-trip, format and corruption tests for lzrans."""
import os
import random
import zlib

import pytest

import lzrans

BS = lzrans.BLOCK_SIZE


def _text(n, seed=0):
    rnd = random.Random(seed)
    words = [bytes(rnd.choices(b"abcdefghijklmnopqrstuvwxyz", k=rnd.randint(2, 9))) for _ in range(400)]
    out = bytearray()
    while len(out) < n:
        out += rnd.choice(words) + b" "
    return bytes(out[:n])


CASES = {
    "empty": b"",
    "one_byte": b"x",
    "three_bytes": b"abc",          # shorter than the minimum match
    "four_bytes": b"abcd",
    "random_64k": os.urandom(65536),
    "zeros_1m": bytes(1 << 20),
    "period_3": b"abc" * 100_000,   # overlapping matches with dist < len
    "period_1_after_text": _text(1000) + b"z" * 50_000,
    "text_300k": _text(300_000),
    "binary_mix": os.urandom(5000) * 20 + _text(50_000),
    "block_edge_minus1": _text(BS - 1, 1),
    "block_edge": _text(BS, 2),
    "block_edge_plus1": _text(BS + 1, 3),
    "far_match_65536": (lambda r: r + os.urandom(65536 - len(r)) + r)(os.urandom(1000)),
    "far_match_3m": (lambda r: r + os.urandom(3 << 20) + r)(os.urandom(5000)),   # needs big window
    "rep_offsets": b"".join(b"key%05d=" % i + b"value" * 3 for i in range(20_000)),
}


@pytest.mark.parametrize("name", sorted(CASES))
@pytest.mark.parametrize("level", [1, 4, 7, 9])
def test_roundtrip(name, level):
    data = CASES[name]
    comp = lzrans.compress(data, level=level)
    assert lzrans.decompress(comp) == data
    size, wlog, flags = lzrans.frame_info(comp)
    assert size == len(data)
    assert 18 <= wlog <= 26
    assert flags == lzrans.FLAG_SIZE


def test_large_multi_block():
    data = _text(5 * BS + 12345, 7) + os.urandom(BS) + bytes(BS)
    comp = lzrans.compress(data)
    assert lzrans.decompress(comp) == data
    assert len(comp) < len(data)


def test_never_expands_much():
    data = os.urandom(3 * BS)
    comp = lzrans.compress(data)
    assert len(comp) <= len(data) + 16 + 9 * 3 + 13   # header + stored blocks + end block


def test_optimal_levels_beat_lazy():
    data = CASES["text_300k"] + CASES["binary_mix"]
    assert len(lzrans.compress(data, level=9)) < len(lzrans.compress(data, level=6))


def test_big_window_finds_far_match():
    r = os.urandom(200_000)
    data = r + os.urandom(3 << 20) + r
    assert len(lzrans.compress(data, level=6)) < len(data) - 150_000


def test_repetitive_compresses_well():
    assert len(lzrans.compress(bytes(1 << 20))) < 2000
    assert len(lzrans.compress(b"abc" * 100_000)) < 2000


def test_buffer_types():
    data = _text(10_000)
    for obj in (bytearray(data), memoryview(data)):
        assert lzrans.decompress(lzrans.compress(obj)) == data
    assert lzrans.decompress(bytearray(lzrans.compress(data))) == data


@pytest.mark.parametrize("data", [b"", b"a", b"hello world", os.urandom(1234), _text(100_000)])
def test_crc32_matches_zlib(data):
    assert lzrans.crc32(data) == zlib.crc32(data)
    half = len(data) // 2
    assert lzrans.crc32(data[half:], lzrans.crc32(data[:half])) == zlib.crc32(data)


def test_invalid_level():
    with pytest.raises(ValueError):
        lzrans.compress(b"x", level=0)
    with pytest.raises(ValueError):
        lzrans.compress(b"x", level=10)
    with pytest.raises(ValueError):
        lzrans.compress(b"x", window_log=30)


@pytest.mark.parametrize("bad", [b"", b"LZRN", b"NOPE" + bytes(30)])
def test_garbage_rejected(bad):
    with pytest.raises(lzrans.LzransError):
        lzrans.decompress(bad)


def test_truncation_rejected():
    comp = lzrans.compress(_text(200_000))
    for cut in (len(comp) - 1, len(comp) - 13, len(comp) // 2, 25, 7):
        with pytest.raises(lzrans.LzransError):
            lzrans.decompress(comp[:cut])


def test_trailing_garbage_rejected():
    comp = lzrans.compress(_text(1000))
    with pytest.raises(lzrans.LzransError):
        lzrans.decompress(comp + b"\x00")


def test_crc_mismatch_detected():
    data = _text(50_000)
    comp = bytearray(lzrans.compress(data))
    comp[-1] ^= 0x01                 # flip one bit of the CRC32 trailer
    with pytest.raises(lzrans.ChecksumError):
        lzrans.decompress(bytes(comp))


def test_stored_block_bitflip_detected_by_crc():
    data = os.urandom(10_000)        # incompressible -> stored block
    comp = bytearray(lzrans.compress(data))
    comp[-20] ^= 0x80                # inside the stored payload
    with pytest.raises(lzrans.ChecksumError):
        lzrans.decompress(bytes(comp))


def test_random_bitflips_never_return_wrong_data():
    data = _text(100_000, 9)
    comp = lzrans.compress(data)
    rnd = random.Random(1)
    for _ in range(500):
        mut = bytearray(comp)
        mut[rnd.randrange(len(mut))] ^= 1 << rnd.randrange(8)
        try:
            out = lzrans.decompress(bytes(mut))
        except (lzrans.LzransError, MemoryError):
            continue
        assert out == data


def test_encrypted_flag_reserved():
    comp = bytearray(lzrans.compress(b"secret"))
    comp[5] |= lzrans.FLAG_ENCRYPTED
    with pytest.raises(lzrans.LzransError, match="encrypted"):
        lzrans.decompress(bytes(comp))


def test_unknown_version_rejected():
    comp = bytearray(lzrans.compress(b"data"))
    comp[4] = 99
    with pytest.raises(lzrans.LzransError):
        lzrans.decompress(bytes(comp))


# ---------------- streaming ----------------

def _chunks(data, rnd, maxchunk):
    i = 0
    while i < len(data):
        n = rnd.randint(0, maxchunk)
        yield data[i:i + n]
        i += n


@pytest.mark.parametrize("name", ["empty", "text_300k", "binary_mix", "block_edge_plus1", "far_match_3m"])
@pytest.mark.parametrize("maxchunk", [1, 777, 300_000])
def test_streaming_roundtrip(name, maxchunk):
    data = CASES[name]
    if maxchunk == 1:
        data = data[:20_000]          # keep the byte-by-byte case fast
    rnd = random.Random(maxchunk)
    c = lzrans.Compressor(level=4)
    comp = b"".join(c.compress(ch) for ch in _chunks(data, rnd, maxchunk)) + c.flush()
    assert lzrans.frame_info(comp)[0] is None          # size unknown when streaming
    assert lzrans.decompress(comp) == data

    d = lzrans.Decompressor()
    out = b"".join(d.decompress(ch) for ch in _chunks(comp, rnd, maxchunk))
    assert out == data and d.eof and d.unused_data == b""


def test_streaming_output_is_incremental():
    c = lzrans.Compressor(level=1)
    first = c.compress(_text(3 * BS))
    assert len(first) > 0                               # full blocks emitted before flush
    d = lzrans.Decompressor()
    assert len(d.decompress(first)) >= 2 * BS           # and decodable before the end
    assert not d.eof


@pytest.mark.parametrize("level", [5, 9])
def test_streaming_matches_one_shot_chunking_independent(level):
    data = _text(700_000, 4)
    outs = set()
    for maxchunk in (1000, 65_536, 10**7):
        c = lzrans.Compressor(level=level)
        rnd = random.Random(0)
        outs.add(b"".join(c.compress(ch) for ch in _chunks(data, rnd, maxchunk)) + c.flush())
    assert len(outs) == 1                               # deterministic regardless of chunking


def test_decompressor_unused_data_and_eof():
    comp = lzrans.compress(b"hello" * 100)
    d = lzrans.Decompressor()
    assert d.decompress(comp + b"TAIL") == b"hello" * 100
    assert d.eof and d.unused_data == b"TAIL"
    with pytest.raises(EOFError):
        d.decompress(b"more")


def test_compressor_after_flush_rejected():
    c = lzrans.Compressor()
    c.compress(b"abc")
    c.flush()
    with pytest.raises(ValueError):
        c.compress(b"x")


def test_streaming_corruption_raises():
    comp = bytearray(lzrans.compress(_text(400_000)))
    comp[len(comp) // 2] ^= 0xFF
    d = lzrans.Decompressor()
    with pytest.raises(lzrans.LzransError):
        for i in range(0, len(comp), 4096):
            d.decompress(bytes(comp[i:i + 4096]))
