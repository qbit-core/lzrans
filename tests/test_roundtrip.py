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
}


@pytest.mark.parametrize("name", sorted(CASES))
@pytest.mark.parametrize("level", [1, 5, 9])
def test_roundtrip(name, level):
    data = CASES[name]
    comp = lzrans.compress(data, level=level)
    assert lzrans.decompress(comp) == data
    size, crc, flags = lzrans.frame_info(comp)
    assert size == len(data)
    assert crc == zlib.crc32(data)
    assert flags == 0


def test_large_multi_block():
    data = _text(5 * BS + 12345, 7) + os.urandom(BS) + bytes(BS)
    comp = lzrans.compress(data)
    assert lzrans.decompress(comp) == data
    assert len(comp) < len(data)


def test_never_expands_much():
    data = os.urandom(3 * BS)
    comp = lzrans.compress(data)
    assert len(comp) <= len(data) + 20 + 5 * 3   # header + stored-block headers


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


@pytest.mark.parametrize("bad", [b"", b"LZRN", b"NOPE" + bytes(30)])
def test_garbage_rejected(bad):
    with pytest.raises(lzrans.LzransError):
        lzrans.decompress(bad)


def test_truncation_rejected():
    comp = lzrans.compress(_text(200_000))
    for cut in (len(comp) - 1, len(comp) // 2, 25):
        with pytest.raises(lzrans.LzransError):
            lzrans.decompress(comp[:cut])


def test_trailing_garbage_rejected():
    comp = lzrans.compress(_text(1000))
    with pytest.raises(lzrans.LzransError):
        lzrans.decompress(comp + b"\x00")


def test_crc_mismatch_detected():
    data = _text(50_000)
    comp = bytearray(lzrans.compress(data))
    comp[16] ^= 0x01                 # flip one bit of the stored CRC32
    with pytest.raises(lzrans.ChecksumError):
        lzrans.decompress(bytes(comp))


def test_stored_block_bitflip_detected_by_crc():
    data = os.urandom(10_000)        # incompressible -> stored block
    comp = bytearray(lzrans.compress(data))
    comp[-1] ^= 0x80
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
        assert out == data           # only header bytes 6-7 are don't-care


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
