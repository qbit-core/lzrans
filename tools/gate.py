"""Check the level ladder against tools/levels_target.json.

Measures each built level and its reference codec on the same files in
the same process (VM speed drifts, so only ratios between codecs measured
together are compared), then prints pass/fail per level and the speed
step between adjacent built levels.

  python tools/gate.py CORPUS_DIR [--files a,b,c] [--levels 1,3,10] [--json out.json]

CORPUS_DIR holds the Silesia files (dickens ... xml). Use --files for a
quick subset while iterating; the gate itself is the full set.
"""
import argparse, json, lzma, os, sys, time
import zstandard
import lzrans

SILESIA = 'dickens mozilla mr nci ooffice osdb reymont samba sao webster x-ray xml'.split()
HERE = os.path.dirname(os.path.abspath(__file__))


def codec(name, threads):
    fam, lv = name.split()
    if fam == 'zstd':
        c = zstandard.ZstdCompressor(level=int(lv[1:]), write_checksum=True, threads=threads if threads > 1 else 0)
        return c.compress
    if fam == 'xz':
        return lambda d: lzma.compress(d, preset=int(lv[1:]))
    if fam == 'lzrans':
        return lambda d: lzrans.compress(d, int(lv), threads=threads)
    raise ValueError(name)


def measure(name, threads, data, reps):
    comp = codec(name, threads)
    size, t = 0, 0.0
    for d in data:
        best = None
        for _ in range(reps):
            s = time.perf_counter(); z = comp(d); e = time.perf_counter() - s
            best = e if best is None else min(best, e)
        size += len(z); t += best
    n = sum(map(len, data))
    return dict(ratio=n / size, speed=n / 1e6 / t)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('corpus')
    ap.add_argument('--files', default=','.join(SILESIA))
    ap.add_argument('--levels', default='')
    ap.add_argument('--json')
    a = ap.parse_args()
    cfg = json.load(open(os.path.join(HERE, 'levels_target.json')))
    data = [open(os.path.join(a.corpus, f), 'rb').read() for f in a.files.split(',')]
    want = {int(x) for x in a.levels.split(',') if x}
    rows, refs = [], {}
    for L in cfg['levels']:
        if L['api_level'] is None or (want and L['level'] not in want):
            continue
        t = L['threads']
        slow = L['level'] >= 10
        key = (L['ref'], 1 if L['ref'].startswith('xz') else t)
        if key not in refs:
            refs[key] = measure(L['ref'], key[1], data, 1 if slow else 2)
        m = measure(f"lzrans {L['api_level']}", t, data, 1 if slow else 2)
        r = refs[key]
        gain = m['ratio'] / r['ratio'] - 1
        rel = m['speed'] / r['speed']
        ok = gain >= L['min_ratio_gain'] - 1e-9 and rel >= L['min_speed_rel']
        rows.append(dict(level=L['level'], api=L['api_level'], threads=t, ref=L['ref'], ratio=m['ratio'], speed=m['speed'],
                         ref_ratio=r['ratio'], ref_speed=r['speed'], gain=gain, speed_rel=rel, ok=ok))
        print(f"L{L['level']:>2} (api {L['api_level']:>2}, {t}T) ratio {m['ratio']:.3f} {m['speed']:8.2f} MB/s | {L['ref']:8s} "
              f"{r['ratio']:.3f} {r['speed']:8.2f} | ratio {gain*100:+5.1f}% (need {L['min_ratio_gain']*100:+.1f}) "
              f"speed {rel:4.2f}x (need {L['min_speed_rel']:.2f}) {'PASS' if ok else 'FAIL'}", flush=True)
    one = sorted((r for r in rows), key=lambda r: r['level'])
    steps = []
    for p, q in zip(one, one[1:]):
        if p['threads'] == q['threads'] and q['level'] == p['level'] + 1:
            steps.append(dict(pair=f"L{p['level']}->L{q['level']}", step=p['speed'] / q['speed']))
    for s in steps:
        flag = 'OK' if s['step'] <= cfg['max_adjacent_speed_ratio'] else 'GAP'
        print(f"{s['pair']}: speed step {s['step']:.2f}x {flag}")
    missing = [L['level'] for L in cfg['levels'] if L['api_level'] is None]
    if missing:
        print('not built yet:', ', '.join(f'L{m}' for m in missing))
    if a.json:
        json.dump(dict(rows=rows, steps=steps, missing=missing), open(a.json, 'w'), indent=1)
    sys.exit(0 if all(r['ok'] for r in rows) else 1)


if __name__ == '__main__':
    main()
