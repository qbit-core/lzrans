"""Benchmark: lzrans vs zlib, xz (lzma) and zstd (if installed).
Reports ratio, compress MB/s and decompress MB/s (single thread).

Usage:  python bench/bench.py [file ...]
Without arguments it uses built-in corpora (Python stdlib source, a binary,
synthetic text, zeros, random).
"""
import glob
import lzma
import os
import sys
import sysconfig
import time
import zlib

import lzrans


def corpora():
    stdlib = sysconfig.get_paths()["stdlib"]
    src = b"".join(open(p, "rb").read() for p in sorted(glob.glob(os.path.join(stdlib, "*.py"))))[: 8 << 20]
    yield "stdlib .py (text)", src
    yield "python binary", open(sys.executable, "rb").read()[: 8 << 20]
    yield "zeros 8MB", bytes(8 << 20)
    yield "random 4MB", os.urandom(4 << 20)


def best_time(fn, arg, min_total=0.3):
    best, total, out = float("inf"), 0.0, None
    while total < min_total:
        t = time.perf_counter()
        out = fn(arg)
        dt = time.perf_counter() - t
        best, total = min(best, dt), total + dt
    return best, out


def run(name, data):
    mb = len(data) / 1e6
    rows = []
    codecs = [(f"lzrans -{lv}", lambda d, lv=lv: lzrans.compress(d, lv), lzrans.decompress) for lv in (1, 6, 9)]
    codecs += [("zlib -9", lambda d: zlib.compress(d, 9), zlib.decompress)]
    codecs += [(f"xz -{lv}", lambda d, lv=lv: lzma.compress(d, preset=lv), lzma.decompress) for lv in (6, 9)]
    try:
        import zstandard  # optional: pip install zstandard
        dz = zstandard.ZstdDecompressor()
        codecs += [(f"zstd -{lv}", lambda d, lv=lv: zstandard.ZstdCompressor(level=lv).compress(d), dz.decompress)
                   for lv in (19, 22)]
    except ImportError:
        pass
    for label, comp, decomp in codecs:
        tc, c = best_time(comp, data)
        td, d = best_time(decomp, c)
        assert d == data
        rows.append((label, len(data) / len(c), mb / tc, mb / td))
    print(f"\n## {name}  ({len(data):,} bytes)\n")
    print("| codec | ratio | compress MB/s | decompress MB/s |")
    print("|---|---:|---:|---:|")
    for r in rows:
        print(f"| {r[0]} | {r[1]:.3f} | {r[2]:.0f} | {r[3]:.0f} |")


if __name__ == "__main__":
    items = [(p, open(p, "rb").read()) for p in sys.argv[1:]] or list(corpora())
    for name, data in items:
        run(name, data)
