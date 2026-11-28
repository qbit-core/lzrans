"""Long-range benchmark: zstd --long / xz / lzrans with and without the LDM
prototype split (tools/ldm/ldm.c). Results append to a JSON-lines file, so
a rerun skips finished (set, method) pairs.

usage: python bench_ldm.py DATADIR OUT.jsonl LDM_BIN [sets...]
"""
import json
import os
import subprocess
import sys
import time

import lzrans

DATA, OUT, LDM = sys.argv[1], sys.argv[2], sys.argv[3]
SETS = sys.argv[4:] or ["src_versions.tar", "silesia_snap.bin", "db_snap.sql", "layers.tar"]
TMP = os.path.join(DATA, "tmp")
os.makedirs(TMP, exist_ok=True)
T = int(os.environ.get("BT", "4"))

done = set()
if os.path.exists(OUT):
    for line in open(OUT):
        r = json.loads(line)
        done.add((r["set"], r["method"]))


def sh(cmd):
    t = time.perf_counter()
    subprocess.run(cmd, shell=True, check=True)
    return time.perf_counter() - t


def cli(name, comp, decomp, src):
    out = os.path.join(TMP, "c.bin")
    ct = sh(f"{comp} < {src} > {out}")
    size = os.path.getsize(out)
    dt = sh(f"{decomp} < {out} > {TMP}/d.bin")
    ok = subprocess.run(f"cmp -s {TMP}/d.bin {src}", shell=True).returncode == 0
    return size, ct, dt, ok


def lz(src_bytes, level):
    t = time.perf_counter()
    c = lzrans.compress(src_bytes, level=level, threads=T)
    ct = time.perf_counter() - t
    t = time.perf_counter()
    d = lzrans.decompress(c)
    dt = time.perf_counter() - t
    return c, ct, dt, d == src_bytes


def ldm_split(src, minlen, min_dist):
    res, refs = os.path.join(TMP, "r.res"), os.path.join(TMP, "r.refs")
    t = sh(f"{LDM} e {src} {res} {refs} {minlen} {min_dist} 6 2>{TMP}/ldm.log")
    return res, refs, t, open(f"{TMP}/ldm.log").read().strip()


def record(r):
    with open(OUT, "a") as f:
        f.write(json.dumps(r) + "\n")
    print(json.dumps(r), flush=True)


ZT = f"-T{T}"
METHODS = [
    ("zstd-19-long27", f"zstd -q -19 --long=27 {ZT} -c", "zstd -q -d --long=27 -c"),
    ("zstd-19-long31", f"zstd -q -19 --long=31 {ZT} -c", "zstd -q -d --long=31 -c"),
    ("zstd-22-long31", f"zstd -q --ultra -22 --long=31 {ZT} -c", "zstd -q -d --long=31 -c"),
    ("xz-9", "xz -9 -T1 -c", "xz -d -c"),
]

for s in SETS:
    src = os.path.join(DATA, s)
    n = os.path.getsize(src)
    for name, c, d in (METHODS if not os.environ.get("NO_BASE") else []):
        if (s, name) in done:
            continue
        size, ct, dt, ok = cli(name, c, d, src)
        record(dict(set=s, method=name, n=n, size=size, ctime=ct, dtime=dt, ok=ok, threads=T))
    data = open(src, "rb").read()
    for lvl in ((11,) if not os.environ.get("NO_BASE") else ()):
        name = f"lzr-L{lvl}"
        if (s, name) in done:
            continue
        c, ct, dt, ok = lz(data, lvl)
        record(dict(set=s, method=name, n=n, size=len(c), ctime=ct, dtime=dt, ok=ok, threads=T))
    VARS = [tuple(int(x) for x in v.split(":")) for v in os.environ.get("LDM_VARS", "256:0").split(",")]
    for minlen, min_dist in VARS:
        tag = f"ldm{minlen}" + (f"d{min_dist}" if min_dist else "")
        todo = [m for m in [f"{tag}+" + x for x in os.environ.get("LDM_BACK", "L11,L10,zstd19").split(",")] if (s, m) not in done]
        if not todo:
            continue
        res, refs, lt, log = ldm_split(src, minlen, min_dist)
        rb, fb = open(res, "rb").read(), open(refs, "rb").read()
        # decode-side cost of the split: rebuild time
        ldt = sh(f"LDM_MINLEN={minlen} {LDM} d {res} {refs} {TMP}/rb.bin")
        for m in todo:
            if m.endswith("zstd19"):
                rs, rct, rdt, ok1 = cli(m, f"zstd -q -19 {ZT} -c", "zstd -q -d -c", res)
                fs, fct, fdt, ok2 = cli(m, f"zstd -q -19 -c", "zstd -q -d -c", refs)
            else:
                lvl = int(m.split("+L")[1])
                c1, rct, rdt, ok1 = lz(rb, lvl)
                c2, fct, fdt, ok2 = lz(fb, 11)
                rs, fs = len(c1), len(c2)
            record(dict(set=s, method=m, n=n, size=rs + fs, res_size=rs, refs_size=fs, res_raw=len(rb), refs_raw=len(fb),
                        ctime=lt + rct + fct, dtime=ldt + rdt + fdt, ldm_ctime=lt, ldm_dtime=ldt, ok=ok1 and ok2,
                        threads=T, ldm_log=log))
