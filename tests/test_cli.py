"""The lzr CLI (liblzr) must produce the same bytes as the Python package and
read what it writes. Skipped unless the CLI is built (`make`), or LZR_CLI
points to a binary."""
import os
import random
import subprocess

import pytest

import lzrans

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CLI = os.environ.get("LZR_CLI", os.path.join(ROOT, "build", "c", "lzr"))
pytestmark = pytest.mark.skipif(not os.access(CLI, os.X_OK), reason="lzr CLI not built (run make)")


def data(n=9 << 20, seed=7):
    r = random.Random(seed)
    words = [bytes(r.choice(b"abcdefghijklmnopqrstuvwxyz") for _ in range(r.randrange(2, 9))) for _ in range(3000)]
    out = bytearray()
    while len(out) < n:
        out += r.choice(words) + b" "
        if r.random() < 0.02:
            out += bytes(r.randrange(256) for _ in range(r.randrange(1, 40)))
    return bytes(out[:n])


@pytest.fixture(scope="module")
def sample(tmp_path_factory):
    d = data()
    p = tmp_path_factory.mktemp("cli") / "sample"
    p.write_bytes(d)
    return p, d


def run(args, stdin=None):
    return subprocess.run([CLI] + args, input=stdin, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True).stdout


@pytest.mark.parametrize("level,threads,segment", [
    (-3, 1, -1), (-1, 4, -1), (0, 1, -1), (0, 3, -1), (1, 4, -1), (3, 2, -1), (6, 1, -1), (6, 4, -1),
    (7, 2, -1), (8, 1, -1), (9, 4, 2 << 20), (10, 1, -1), (10, 2, 1 << 20), (11, 1, -1), (11, 4, 4 << 20),
])
def test_file_matches_python_oneshot(sample, level, threads, segment):
    p, d = sample
    args = ["-c", f"--level={level}", f"-T{threads}"] + ([f"--segment={segment}"] if segment >= 0 else [])
    out = run(args + [str(p)])
    assert out == lzrans.compress(d, level=level, threads=threads, segment=segment)
    assert lzrans.decompress(out, threads=2) == d


@pytest.mark.parametrize("level,threads", [(-2, 1), (0, 4), (5, 3), (8, 1), (10, 4)])
def test_pipe_matches_python_streaming(sample, level, threads):
    _, d = sample
    out = run(["-c", f"--level={level}", f"-T{threads}"], stdin=d)   # stdin is a pipe: size unknown
    c = lzrans.Compressor(level=level, threads=threads)
    ref = b"".join(c.compress(d[i:i + 65536]) for i in range(0, len(d), 65536)) + c.flush()
    assert out == ref
    assert lzrans.frame_info(out)[0] is None


@pytest.mark.parametrize("threads", [1, 4])
def test_cli_decodes_python_frames(sample, threads):
    _, d = sample
    for kw in ({"level": 11}, {"level": 10, "segment": 1 << 20}, {"level": 2, "filter": "x86"}):
        assert run(["-dc", f"-T{threads}"], stdin=lzrans.compress(d, **kw)) == d
    c = lzrans.Compressor(level=4)
    assert run(["-dc"], stdin=c.compress(d) + c.flush()) == d
