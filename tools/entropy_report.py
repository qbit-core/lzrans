"""Compare where the bytes go: lzrans block type 5 vs zstd, same input.

Needs an lzrans build with -DLZR_STATS (emit_seqs5 prints one LZRSTAT5
line per block to stderr). Usage:
  python tools/entropy_report.py STATS_BUILD_DIR LZR_LEVEL ZSTD_LEVEL FILE...
Prints, per side: literals and bits per literal, sequences and bits per
sequence (tables + bit stream), stored/raw bytes, total.
"""
import json, os, subprocess, sys
import zstandard

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import zstd_breakdown

build, lv, zl, files = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4:]
CHILD = r'''
import sys; sys.path.insert(0, sys.argv[1]); import lzrans
d = open(sys.argv[3], 'rb').read()
z = lzrans.compress(d, int(sys.argv[2]), threads=1)
print(len(d), len(z))
'''
tot = {'lz': dict(orig=0, total=0, nlit=0, lit_bytes=0, nseq=0, seq_bytes=0, tab_bytes=0, extra_bits=0),
       'zs': dict(orig=0, total=0, nlit=0, lit_bytes=0, nseq=0, seq_bytes=0, raw_bytes=0, overhead=0)}
for f in files:
    r = subprocess.run([sys.executable, '-c', CHILD, build, str(lv), f], capture_output=True, text=True, check=True)
    orig, size = map(int, r.stdout.split())
    t = tot['lz']; t['orig'] += orig; t['total'] += size
    for line in r.stderr.splitlines():
        if not line.startswith('LZRSTAT5'):
            continue
        v = list(map(int, line.split()[1:]))
        t['nlit'] += v[0]; t['lit_bytes'] += v[1]; t['nseq'] += v[2]
        t['tab_bytes'] += v[3]; t['seq_bytes'] += v[3] + v[4]
        t['extra_bits'] += v[5] if len(v) > 5 else 0
    d = open(f, 'rb').read()
    z = zstandard.ZstdCompressor(level=zl).compress(d)
    s = zstd_breakdown.breakdown(z)
    t = tot['zs']; t['orig'] += len(d); t['total'] += len(z)
    for k in ('nlit', 'lit_bytes', 'nseq', 'seq_bytes', 'raw_bytes', 'overhead'):
        t[k] += s[k]
    print(os.path.basename(f), 'done', file=sys.stderr, flush=True)
out = {}
for side, t in tot.items():
    o = dict(t)
    o['ratio'] = t['orig'] / t['total']
    o['bits_per_lit'] = 8 * t['lit_bytes'] / max(1, t['nlit'])
    o['bits_per_seq'] = 8 * t['seq_bytes'] / max(1, t['nseq'])
    o['other_bytes'] = t['total'] - t['lit_bytes'] - t['seq_bytes']
    out[side] = o
print(json.dumps(dict(lzr_level=lv, zstd_level=zl, **out), indent=1))
