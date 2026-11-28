#!/bin/sh
# CLI tests for lzr: file and pipe round trips, xz-style file handling,
# concatenated frames, error exits, and tar integration.
#   sh tests/cli/run.sh [path/to/lzr]        (make check runs this)
set -u
LZR=$(cd "$(dirname "${1:-build/c/lzr}")" && pwd)/$(basename "${1:-build/c/lzr}")
[ -x "$LZR" ] || { echo "lzr binary not found: $LZR" >&2; exit 1; }
T=$(mktemp -d "${TMPDIR:-/tmp}/lzr-cli.XXXXXX")
trap 'rm -rf "$T"' EXIT
cd "$T" || exit 1
ln -s "$LZR" unlzr
ln -s "$LZR" lzrcat
fail=0
pass=0
ok() { pass=$((pass + 1)); }
bad() { echo "FAIL: $*" >&2; fail=$((fail + 1)); }
check() { if "$@"; then ok; else bad "$*"; fi; }
# expect EXIT CMD...: run CMD, compare its exit status
expect() {
    want=$1; shift
    "$@" >/dev/null 2>err.txt; got=$?
    if [ "$got" -eq "$want" ]; then ok; else bad "exit $got (want $want): $* ($(head -c 200 err.txt))"; fi
}

# test data: text-like, ~3 MB, plus random and empty
i=0
: > text0
while [ $i -lt 12 ]; do cat "$LZR" >> text0; printf 'line %s of some text\n' $i >> text0; i=$((i + 1)); done
head -c 3000000 text0 > text; rm text0   # fixed size: 3 segments of 1 MiB
head -c 3000000 /dev/urandom > rand 2>/dev/null || dd if=/dev/urandom of=rand bs=1000 count=3000 2>/dev/null
: > empty

# 1. file round trip at every level; -k keeps the input; output identical
for lv in --fast=3 --fast=2 --fast=1 -0 -1 -2 -3 -4 -5 -6 -7 -8 -9 -10 -11; do
    cp text t
    "$LZR" -k $lv t || bad "compress $lv"
    check test -f t
    check test -f t.lzr
    rm t
    "$LZR" -d t.lzr || bad "decompress $lv"
    check test ! -f t.lzr
    check cmp -s t text
done

# 2. default: input removed, mode and mtime kept
cp text u; chmod 640 u; touch -t 202001020304 u
"$LZR" u
check test ! -f u
check test -f u.lzr
"$LZR" -d u.lzr
check cmp -s u text
check test "$(stat -c %a u 2>/dev/null || stat -f %Lp u)" = 640
check test "$(stat -c %Y u 2>/dev/null || stat -f %m u)" = "$(touch -t 202001020304 ref && (stat -c %Y ref 2>/dev/null || stat -f %m ref))"

# 3. refuse to overwrite without -f; skip names already .lzr / without .lzr
cp text v; cp text v2; "$LZR" -k v
expect 1 "$LZR" -k v
"$LZR" -kf v || bad "-f overwrite"
expect 2 "$LZR" v.lzr
expect 2 "$LZR" -d v2
check test -f v2

# 4. pipes (no content size in the header) and -c
"$LZR" -c text > c1.lzr
cat text | "$LZR" > c2.lzr
cat text | "$LZR" -6 -T3 > c3.lzr
check test "$(wc -c < c2.lzr)" -lt "$(wc -c < c1.lzr)"   # 8 bytes less: no size field
cat c2.lzr | "$LZR" -d > p2; check cmp -s p2 text
"$LZR" -dc c1.lzr > p1; check cmp -s p1 text
"$LZR" -dc < c3.lzr > p3; check cmp -s p3 text
./lzrcat c1.lzr > p4; check cmp -s p4 text
cp c1.lzr w.lzr; ./unlzr w.lzr; check cmp -s w text
cat rand | "$LZR" -1 | "$LZR" -dc > r1; check cmp -s r1 rand
"$LZR" -c empty | "$LZR" -dc > e1; check cmp -s e1 empty
check test ! -s e1
"$LZR" -c - < text > c5.lzr; "$LZR" -dc - < c5.lzr > p5; check cmp -s p5 text

# 5. same bytes from a file, a pipe with -T, and parallel decode
"$LZR" -c -10 --segment=1MiB text > s1.lzr
"$LZR" -c -10 --segment=1MiB -T4 text > s2.lzr
check cmp -s s1.lzr s2.lzr
"$LZR" -dc -T4 s1.lzr > q1; check cmp -s q1 text
"$LZR" -dc -T0 s1.lzr > q2; check cmp -s q2 text
"$LZR" -c --x86 -9 "$LZR" > x.lzr; "$LZR" -dc x.lzr > xb; check cmp -s xb "$LZR"
"$LZR" -c --window-log=18 -6 text > wl.lzr   # levels 9-11 use >= 2^20
"$LZR" -dc --max-window-log=18 wl.lzr > wl; check cmp -s wl text
expect 1 "$LZR" -dc --max-window-log=18 s1.lzr

# 6. concatenated frames, trailing garbage, -t, -l
cat c1.lzr c2.lzr > cat.lzr
"$LZR" -dc cat.lzr > two; cat text text > text2; check cmp -s two text2
cat c1.lzr c2.lzr | "$LZR" -dc > two; check cmp -s two text2
"$LZR" -dc -T4 cat.lzr > two; check cmp -s two text2
expect 0 "$LZR" -t cat.lzr
cat c1.lzr > g.lzr; printf 'garbage' >> g.lzr
expect 2 "$LZR" -t g.lzr
expect 0 "$LZR" -lv cat.lzr s1.lzr
check sh -c "'$LZR' -l s1.lzr | grep -q ' 3 '"     # 3 segments

# 7. corrupt and truncated input fail with status 1 and leave no output
head -c 1000 c1.lzr > tr.lzr
expect 1 "$LZR" -t tr.lzr
expect 1 "$LZR" -d tr.lzr
check test ! -f tr
check test -f tr.lzr
cp c1.lzr co.lzr
printf '\377\377\377\377' | dd of=co.lzr bs=1 seek=5000 conv=notrunc 2>/dev/null
expect 1 "$LZR" -t co.lzr
expect 1 "$LZR" -t -T4 co.lzr
expect 1 "$LZR" -t text
expect 1 "$LZR" -dc empty
expect 1 sh -c "cat tr.lzr | '$LZR' -dc"
expect 1 "$LZR" -t missing-file

# 8. options
expect 0 "$LZR" -h
expect 0 "$LZR" -V
expect 1 "$LZR" -12 text
expect 1 "$LZR" --fast=4 text
expect 1 "$LZR" --segment=1000 text
expect 1 "$LZR" -T 99 text
expect 1 "$LZR" --bogus
cp text o1; "$LZR" -9kT2 o1; check test -f o1.lzr
cp text o2; "$LZR" -S .z -k --level=3 o2; check test -f o2.z
"$LZR" -d -S .z -c o2.z > o2b; check cmp -s o2b text

# 9. tar integration: tar -I runs "lzr" to compress and "lzr -d" to extract
if tar --version 2>/dev/null | grep -q GNU; then
    mkdir -p src/sub; cp text src/a; cp rand src/sub/b; cp "$LZR" src/sub/c
    tar -I "$LZR" -cf x.tar.lzr src || bad "tar create"
    tar -I "$LZR -9 -T2" -cf y.tar.lzr src || bad "tar create with options"
    mkdir out1 out2
    tar -I "$LZR" -xf x.tar.lzr -C out1 || bad "tar extract"
    tar -I "$LZR -T4" -xf y.tar.lzr -C out2 || bad "tar extract with threads"
    check diff -r src out1/src
    check diff -r src out2/src
    tar -cf - src | "$LZR" -c > z.tar.lzr
    "$LZR" -dc z.tar.lzr | tar -tf - > list.txt
    check grep -q 'src/sub/c' list.txt
else
    echo "skip: GNU tar not found" >&2
fi

echo "cli tests: $pass passed, $fail failed"
[ "$fail" -eq 0 ]
