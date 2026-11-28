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
@pytest.mark.parametrize("level", [-3, -1, 0, 1, 4, 7, 8, 9, 11])
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
    assert len(lzrans.compress(data, level=11)) < len(lzrans.compress(data, level=6))


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
        lzrans.compress(b"x", level=-4)
    with pytest.raises(ValueError):
        lzrans.compress(b"x", level=12)
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


@pytest.mark.parametrize("level", [0, 5, 7, 11])
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


# ---------------- threads ----------------

@pytest.mark.parametrize("level", [7, 9, 11])
def test_threads_do_not_change_output(level):
    data = CASES["text_300k"] * 3 + CASES["binary_mix"] + CASES["far_match_3m"]
    ref = lzrans.compress(data, level=level, threads=1)
    for t in (2, 3, 4, 0):
        assert lzrans.compress(data, level=level, threads=t) == ref
    c = lzrans.Compressor(level=level, threads=4)
    rnd = random.Random(5)
    streamed = b"".join(c.compress(ch) for ch in _chunks(data, rnd, 500_000)) + c.flush()
    assert lzrans.decompress(streamed) == data
    assert lzrans.decompress(ref) == data


def test_threads_identical_on_window_sized_input():
    # Inputs larger than the window exercise tree-slot reuse across threads.
    rnd = random.Random(11)
    parts = [_text(200_000, k) for k in range(8)]
    data = b"".join(rnd.choice(parts)[: rnd.randint(1000, 200_000)] for _ in range(80))
    c1 = lzrans.compress(data, level=9, window_log=20, threads=1)
    assert lzrans.compress(data, level=9, window_log=20, threads=4) == c1
    assert lzrans.decompress(c1) == data


@pytest.mark.parametrize("level", [-3, -1, 0, 1, 3, 6])
def test_job_mode_levels_1_6(level):
    # Levels 1-6 with threads > 1 compress independent 4 MB jobs: the
    # output is the same for any thread count > 1 and any input chunking,
    # within ~1% of the single-thread size, and always decodes.
    rnd = random.Random(level)
    # seeded, not CASES["binary_mix"] (os.urandom): the 1% size bound is
    # then checked on the same bytes every run
    parts = [_text(200_000, k) for k in range(6)] + [random.Random(99).randbytes(5000) * 20 + _text(50_000)]
    data = b"".join(rnd.choice(parts)[: rnd.randint(1000, 200_000)] for _ in range(120))
    one = lzrans.compress(data, level=level, threads=1)
    ref = lzrans.compress(data, level=level, threads=2)
    for t in (3, 4, 0):
        assert lzrans.compress(data, level=level, threads=t) == ref
    streamed = []
    for t, maxchunk in ((2, 700_000), (4, 50_000)):   # (no content size in the header)
        c = lzrans.Compressor(level=level, threads=t)
        streamed.append(b"".join(c.compress(ch) for ch in _chunks(data, rnd, maxchunk)) + c.flush())
    assert streamed[0] == streamed[1]
    assert lzrans.decompress(ref) == data and lzrans.decompress(streamed[0]) == data
    assert len(ref) <= len(one) * 1.01


def test_job_mode_freed_mid_stream():
    # Dropping a threaded compressor with queued jobs (and helper threads
    # running) must not leak or crash; a fresh one still produces the
    # normal output.
    data = _text(200_000, 1) * 30
    c = lzrans.Compressor(level=0, threads=4)
    ref = c.compress(data) + c.flush()
    assert lzrans.decompress(ref) == data
    for cut in (1_000_000, 4_500_000, 5_000_000):
        c = lzrans.Compressor(level=0, threads=4)
        c.compress(data[:cut])
        del c
    c = lzrans.Compressor(level=0, threads=4)
    assert c.compress(data) + c.flush() == ref


@pytest.mark.parametrize("level", [-3, -1, 0, 2, 6])
def test_one_shot_jobs_match_streaming(level):
    # compress() with threads reads the input in place (no copy) for all
    # but the last jobs; the blocks must be those of the copying path
    # (streaming), whose frame differs only by the content size field.
    rnd = random.Random(level + 10)
    parts = [_text(300_000, k) for k in range(5)] + [CASES["binary_mix"]]
    for n in (4 * 2**20 + 63, 4 * 2**20 + 65, 9 * 2**20 + 12345, 13 * 2**20):
        data = b"".join(rnd.choice(parts)[: rnd.randint(1000, 300_000)] for _ in range(n // 100_000 + 1))[:n]
        for t in (2, 4):
            one = lzrans.compress(data, level=level, threads=t)
            c = lzrans.Compressor(level=level, threads=t)
            st = c.compress(data) + c.flush()
            assert one[16:] == st[8:], (n, t)
            assert lzrans.decompress(one) == data


def test_threads_validated():
    with pytest.raises(ValueError):
        lzrans.compress(b"x", threads=-1)
    with pytest.raises(ValueError):
        lzrans.Compressor(threads=99)


@pytest.mark.parametrize("level", [3, 8, 10, 11])
def test_adaptive_model_survives_stored_blocks(level):
    # Blocks carry state across blocks (type 5: Huffman / tANS tables for
    # the repeat modes; type 4: the adaptive token model); incompressible
    # blocks are stored and must leave it untouched on both sides.
    rnd = random.Random(7)
    parts = []
    for i in range(6):
        parts.append(_text(BS + 1234 * i, seed=i))
        parts.append(rnd.randbytes(BS // 2 + 17 * i))
    data = b"".join(parts)
    blob = lzrans.compress(data, level=level)
    assert lzrans.decompress(blob) == data
    d = lzrans.Decompressor()
    out = b"".join(d.decompress(blob[i:i + 1000]) for i in range(0, len(blob), 1000))
    assert out == data and d.eof


def _block_types(blob):
    pos = 16 if blob[5] & 0x02 else 8
    types = []
    while True:
        t, comp = blob[pos], int.from_bytes(blob[pos + 5:pos + 9], "little")
        types.append(t)
        pos += 9 + comp
        if t == 2:
            return types


@pytest.mark.parametrize("level,btype", [(0, 5), (1, 5), (6, 5), (7, 5), (8, 5), (10, 5), (11, 4)])
def test_block_type_by_level(level, btype):
    blob = lzrans.compress(_text(3 * BS), level=level)
    assert blob[4] == 7                              # format version
    assert set(_block_types(blob)) == {btype, 2}


def test_long_matches_split_into_tokens():
    # Matches longer than one token (273 bytes) and runs at block edges.
    data = (b"x" * 5000 + _text(2000) + b"ab" * 40000) * 3
    for level in (6, 8, 11):
        assert lzrans.decompress(lzrans.compress(data, level=level)) == data


def _x86_like(n, seed=3):
    # Code-like bytes: many E8 calls with small displacements, plus E8 at
    # chunk-hostile spots (stream end, operands that contain E8 bytes).
    rnd = random.Random(seed)
    out = bytearray()
    targets = [rnd.randrange(0, n) for _ in range(300)]
    while len(out) < n:
        k = rnd.random()
        if k < 0.3:
            out += b"\xe8" + ((rnd.choice(targets) - len(out) - 5) & 0xFFFFFFFF).to_bytes(4, "little")
        elif k < 0.35:
            out += b"\xe8\xe8\xe8"
        else:
            out += rnd.randbytes(rnd.randint(1, 6))
    return bytes(out[:n]) + b"\xe8\x01\x02"


@pytest.mark.parametrize("level", [0, 1, 6, 7, 11])
def test_x86_filter_roundtrip(level):
    data = _x86_like(700_000)
    blob = lzrans.compress(data, level=level, filter="x86")
    assert lzrans.frame_info(blob)[2] & lzrans.FLAG_X86
    assert lzrans.decompress(blob) == data
    # chunked streaming in both directions gives the same bytes
    c = lzrans.Compressor(level=level, filter="x86")
    rnd = random.Random(1)
    parts, i = [], 0
    while i < len(data):
        j = i + rnd.randint(1, 9000)
        parts.append(c.compress(data[i:j]))
        i = j
    blob2 = b"".join(parts) + c.flush()
    d = lzrans.Decompressor()
    out = b"".join(d.decompress(blob2[k:k + 777]) for k in range(0, len(blob2), 777))
    assert out == data and d.eof


@pytest.mark.parametrize("data", [b"", b"\xe8", b"\xe8\x00\x00\x00", b"\xe8\x00\x00\x00\x00", bytes(range(256)) * 50])
def test_x86_filter_edge_cases(data):
    assert lzrans.decompress(lzrans.compress(data, filter="x86")) == data


def test_x86_filter_helps_code():
    data = _x86_like(1 << 20)
    assert len(lzrans.compress(data, 11, filter="x86")) < len(lzrans.compress(data, 11))


def test_filter_validated():
    with pytest.raises(ValueError):
        lzrans.compress(b"x", filter="arm")
    with pytest.raises(ValueError):
        lzrans.Compressor(filter=1)


def _markov(n, seed=3):
    # bytes whose distribution depends on the byte before (block type 6)
    rnd = random.Random(seed)
    nxt = {a: rnd.sample(range(32, 127), 6) for a in range(256)}
    out, b = bytearray(), 32
    for _ in range(n):
        b = rnd.choice(nxt[b]) if rnd.random() < 0.9 else rnd.randrange(32, 127)
        out.append(b)
    return bytes(out)


@pytest.mark.parametrize("level", [9, 10])
def test_context_literals_block_type_6(level):
    data = _markov(3 * BS)
    blob = lzrans.compress(data, level=level)
    assert 6 in _block_types(blob)
    assert lzrans.decompress(blob) == data
    d = lzrans.Decompressor()
    out = b"".join(d.decompress(blob[i:i + 777]) for i in range(0, len(blob), 777))
    assert out == data and d.eof
    # corrupted type 6 blocks are rejected or decode to something, never crash
    rnd = random.Random(level)
    for _ in range(40):
        bad = bytearray(blob)
        for _ in range(3):
            bad[rnd.randrange(8, len(bad))] ^= 1 << rnd.randrange(8)
        try:
            lzrans.decompress(bytes(bad))
        except lzrans.LzransError:
            pass


@pytest.mark.parametrize("level", [1, 6, 7, 8, 11])
def test_no_type_6_outside_levels_9_10(level):
    # decode speed: levels other than 9-10 keep their literal coding
    assert 6 not in _block_types(lzrans.compress(_markov(3 * BS), level=level))


def test_writes_format_version_7():
    # Versions 5 and 6 still decode: see the frozen vectors in
    # tests/test_vectors.py (enc-version5, the other enc-* are version 6).
    assert lzrans.FORMAT_VERSION == 7
    assert lzrans.compress(b"abc" * 1000, level=1)[4] == 7


# ---- segment mode (reset points + seek index) and parallel decompression

MB = 1 << 20


def _seg_data(n=3 * MB + 12345, seed=8):
    # text with long-range repeats, so segments really lose matches
    t = _text(n // 2, seed)
    return (t + t[::-1])[:n]


@pytest.mark.parametrize("level", [-3, 0, 4, 7, 9, 10, 11])
def test_segment_mode_roundtrip_and_parallel_decode(level):
    data = _seg_data()
    blob = lzrans.compress(data, level=level, segment=MB)
    assert lzrans.decompress(blob) == data
    for t in (2, 3, 4, 8):
        assert lzrans.decompress(blob, threads=t) == data
    # the output does not depend on the thread count
    assert lzrans.compress(data, level=level, segment=MB, threads=3) == blob


def test_segment_mode_small_input_is_unchanged():
    data = _text(MB - 7)
    for level in (1, 9, 11):
        assert lzrans.compress(data, level=level, segment=MB) == lzrans.compress(data, level=level, segment=0)


def test_segment_mode_streaming_matches_one_shot():
    data = _seg_data(2 * MB + 999)
    one = lzrans.compress(data, level=6, segment=MB)
    rnd = random.Random(3)
    c = lzrans.Compressor(level=6, segment=MB, threads=2)
    parts, i = [], 0
    while i < len(data):
        k = rnd.randint(1, 400_000)
        parts.append(c.compress(data[i:i + k]))
        i += k
    blob = b"".join(parts) + c.flush()
    # same data blocks and checksum; the header has no content size when
    # streaming, so the seek index's frame offsets differ by 8
    ix = 9 + 16 * 3 + 16 + 13
    assert blob[8:-ix] == one[16:-ix] and blob[-13:] == one[-13:]
    assert lzrans.decompress(blob, threads=4) == data


def test_segment_mode_x86_filter():
    data = _seg_data(2 * MB + 5) + bytes([0xE8, 1, 2, 3]) * 1000
    blob = lzrans.compress(data, level=9, segment=MB, filter="x86", threads=2)
    assert lzrans.decompress(blob, threads=4) == data


def test_segment_validated():
    with pytest.raises(ValueError):
        lzrans.compress(b"x", segment=1000)
    with pytest.raises(ValueError):
        lzrans.Compressor(segment=5)


def _segment_offsets(blob):
    import struct
    pos, offs = 16, []
    while True:
        t, raw, comp = blob[pos], *struct.unpack_from("<II", blob, pos + 1)
        if t & 0xC0 == 0x40:
            offs.append(pos)
        if t == 2:
            return offs
        pos += 9 + comp


def test_parallel_decode_detects_corruption_in_each_segment():
    data = _seg_data()
    blob = lzrans.compress(data, level=9, segment=MB)
    offs = _segment_offsets(blob)
    assert len(offs) == 3
    for o in offs + [16]:
        bad = bytearray(blob)
        bad[o + 9 + 40] ^= 0x10   # inside the first block of that segment
        with pytest.raises(lzrans.LzransError):
            lzrans.decompress(bytes(bad), threads=4)
    bad = bytearray(blob)
    bad[-1] ^= 1
    with pytest.raises(lzrans.ChecksumError):
        lzrans.decompress(bytes(bad), threads=4)
    with pytest.raises(lzrans.LzransError):
        lzrans.decompress(blob[:-5], threads=4)
    with pytest.raises(lzrans.LzransError):
        lzrans.decompress(blob + b"x", threads=4)


def test_parallel_decode_without_reset_points():
    data = _seg_data(MB)
    for level in (3, 11):
        blob = lzrans.compress(data, level=level, segment=0)
        assert lzrans.decompress(blob, threads=4) == data
