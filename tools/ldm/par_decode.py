"""Feasibility of parallel decompression at level 9-11 ratios: split the
input into independent chunks (each its own frame), compress each, and
decode the chunks on N threads. Reports ratio loss vs one frame and the
decode speedup. Also xz -9 -T4 (multi-block, threaded decode) for reference.

usage: python par_decode.py FILE [chunk_mb ...]
"""
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import lzrans

path = sys.argv[1]
chunks_mb = [int(x) for x in sys.argv[2:]] or [64, 32, 16]
data = open(path, "rb").read()
n = len(data)


def best_of(f, k=3):
    return min(f() or 0 for _ in range(k)) if False else min(_t(f) for _ in range(k))


def _t(f):
    t = time.perf_counter(); f(); return time.perf_counter() - t


for lvl in (11, 10):
    whole = lzrans.compress(data, level=lvl, threads=4)
    d1 = min(_t(lambda: lzrans.decompress(whole)) for _ in range(3))
    print(f"L{lvl} whole: ratio {n/len(whole):.4f} dec1T {n/1e6/d1:.0f} MB/s", flush=True)
    for cm in chunks_mb:
        cs = cm << 20
        parts = [data[i:i + cs] for i in range(0, n, cs)]
        with ThreadPoolExecutor(4) as ex:
            comp = list(ex.map(lambda p: lzrans.compress(p, level=lvl, threads=1), parts))
        size = sum(map(len, comp))
        with ThreadPoolExecutor(4) as ex:
            dt = min(_t(lambda: list(ex.map(lzrans.decompress, comp))) for _ in range(3))
        print(f"L{lvl} {cm}MB x{len(parts)}: ratio {n/size:.4f} ({(size/len(whole)-1)*100:+.2f}% size) "
              f"dec4T {n/1e6/dt:.0f} MB/s ({d1/dt:.2f}x)", flush=True)

for t in ("1", "4"):
    out = "/home/user/data/tmp/x.xz"
    subprocess.run(f"xz -9 -T{t} -c {path} > {out}", shell=True, check=True)
    sz = os.path.getsize(out)
    dt = min(_t(lambda: subprocess.run(f"xz -d -T4 -c {out} > /dev/null", shell=True, check=True)) for _ in range(3))
    print(f"xz -9 -T{t}: ratio {n/sz:.4f}, decode with -T4 {n/1e6/dt:.0f} MB/s", flush=True)
