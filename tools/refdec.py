#!/usr/bin/env python3
"""refdec.py - reference decoder for the LZR frame format, written from
docs/FORMAT.md alone (no code shared with the C core).

It is deliberately simple and slow (pure Python, no dependencies): its job
is to show that the specification is complete, and to cross-check the C
decoder on the golden vectors in tests/vectors/. Section numbers in the
comments refer to docs/FORMAT.md.

    python tools/refdec.py FILE.lzr            # decode to stdout
    python tools/refdec.py --trace FILE.lzr    # list frame header and blocks
"""
import struct
import sys

__all__ = ["decode", "trace", "xxh64", "FormatError", "ChecksumError", "UnsupportedError"]


class FormatError(Exception):
    """Malformed frame (spec: "corrupt")."""


class UnsupportedError(FormatError):
    """Well-formed but uses a feature this decoder does not implement."""


class ChecksumError(FormatError):
    """Content checksum mismatch."""


def _need(cond, msg):
    if not cond:
        raise FormatError(msg)


# ---------------------------------------------------------------- 9. XXH64
_M64 = (1 << 64) - 1
_P1, _P2, _P3, _P4, _P5 = (0x9E3779B185EBCA87, 0xC2B2AE3D27D4EB4F, 0x165667B19E3779F9,
                          0x85EBCA77C2B2AE63, 0x27D4EB2F165667C5)


def _rotl(x, r):
    return ((x << r) | (x >> (64 - r))) & _M64


def _round(acc, v):
    acc = (acc + v * _P2) & _M64
    return (_rotl(acc, 31) * _P1) & _M64


def xxh64(data, seed=0):
    n, i = len(data), 0
    if n >= 32:
        v = [(seed + _P1 + _P2) & _M64, (seed + _P2) & _M64, seed, (seed - _P1) & _M64]
        while i + 32 <= n:
            for j in range(4):
                v[j] = _round(v[j], struct.unpack_from("<Q", data, i + 8 * j)[0])
            i += 32
        h = (_rotl(v[0], 1) + _rotl(v[1], 7) + _rotl(v[2], 12) + _rotl(v[3], 18)) & _M64
        for j in range(4):
            h = ((h ^ _round(0, v[j])) * _P1 + _P4) & _M64
    else:
        h = (seed + _P5) & _M64
    h = (h + n) & _M64
    while i + 8 <= n:
        h ^= _round(0, struct.unpack_from("<Q", data, i)[0])
        h = (_rotl(h, 27) * _P1 + _P4) & _M64
        i += 8
    if i + 4 <= n:
        h ^= (struct.unpack_from("<I", data, i)[0] * _P1) & _M64
        h = (_rotl(h, 23) * _P2 + _P3) & _M64
        i += 4
    while i < n:
        h ^= (data[i] * _P5) & _M64
        h = (_rotl(h, 11) * _P1) & _M64
        i += 1
    h ^= h >> 33
    h = (h * _P2) & _M64
    h ^= h >> 29
    h = (h * _P3) & _M64
    h ^= h >> 32
    return h


# ------------------------------------------------- 4. bit streams
class FwdBits:
    """Forward LSB-first bit stream (section 4.1) over data[start:end]."""

    def __init__(self, data, start, end):
        self.d, self.base, self.end, self.pos = data, start, end, 0   # pos in bits

    def peek(self, nb):
        # bits past `end` read as 0 (they can never be part of a valid code)
        byte, sh = self.base + (self.pos >> 3), self.pos & 7
        chunk = self.d[byte:min(byte + 3, self.end)]
        return (int.from_bytes(chunk, "little") >> sh) & ((1 << nb) - 1)

    def skip(self, nb):
        self.pos += nb
        _need(self.base + (self.pos + 7) // 8 <= self.end, "forward bit stream overrun")

    def used_bytes(self):
        return (self.pos + 7) // 8


class BackBits:
    """Backward bit stream (section 4.2) over data[start:end]."""

    def __init__(self, data, start, end):
        _need(end > start and data[end - 1] != 0, "backward stream: no end marker")
        self.d, self.base = data, start
        self.pos = (end - start - 1) * 8 + data[end - 1].bit_length() - 1   # bits left

    def read(self, nb):
        if nb == 0:
            return 0
        _need(self.pos >= nb, "backward bit stream overrun")
        self.pos -= nb
        lo, hi = self.pos >> 3, (self.pos + nb + 7) >> 3
        v = int.from_bytes(self.d[self.base + lo:self.base + hi], "little")
        return (v >> (self.pos & 7)) & ((1 << nb) - 1)


# ------------------------------------------------- 6.3 Huffman codes
HUF_MAXLEN = 11


def huf_table(lengths):
    """Decoding table (2^11 entries of (symbol, length)) from code lengths,
    canonical codes as in DEFLATE (section 6.3.1)."""
    _need(all(l <= HUF_MAXLEN for l in lengths), "Huffman code length > 11")
    kraft = sum(1 << (HUF_MAXLEN - l) for l in lengths if l)
    _need(0 < kraft <= 1 << HUF_MAXLEN, "Huffman lengths oversubscribed or empty")
    bl_count = [0] * (HUF_MAXLEN + 1)
    for l in lengths:
        if l:
            bl_count[l] += 1
    next_code, code = [0] * (HUF_MAXLEN + 2), 0
    for l in range(1, HUF_MAXLEN + 1):
        code = (code + bl_count[l - 1]) << 1 if l > 1 else 0
        next_code[l] = code
    # unassigned slots (incomplete code): symbol 0, 11 bits (section 6.3.1)
    table = [(0, HUF_MAXLEN)] * (1 << HUF_MAXLEN)
    for s, l in enumerate(lengths):
        if not l:
            continue
        c = next_code[l]
        next_code[l] += 1
        r = int(format(c, "0%db" % l)[::-1], 2)     # first code bit -> lowest stream bit
        for j in range(r, 1 << HUF_MAXLEN, 1 << l):
            table[j] = (s, l)
    return table


def read_lengths(p, k):
    """last u8 + nibbles (section 6.3.2). Returns (lengths[256], new k)."""
    _need(k < len(p), "truncated Huffman table")
    last = p[k]
    nb = (last + 2) // 2
    _need(k + 1 + nb <= len(p), "truncated Huffman table")
    lens = [0] * 256
    for s in range(last + 1):
        lens[s] = (p[k + 1 + s // 2] >> (4 * (s & 1))) & 15
    return lens, k + 1 + nb


def huf_sym(bits, table):
    s, l = table[bits.peek(HUF_MAXLEN)]
    bits.skip(l)
    return s


# ------------------------------------------------- 6.4 tANS
FSE_LOG = 9
FSE_L = 1 << FSE_LOG


def read_varint(p, k):
    v = 0
    for i in range(3):
        _need(k + i < len(p), "truncated varint")
        v |= (p[k + i] & 0x7F) << (7 * i)
        if not p[k + i] & 0x80:
            return v, k + i + 1
    raise FormatError("varint longer than 3 bytes")


def read_norm(p, k, alphabet):
    norm, s, cum = [0] * alphabet, 0, 0
    while s < alphabet:
        v, k = read_varint(p, k)
        if v == 0:
            r, k = read_varint(p, k)
            _need(s + r + 1 <= alphabet, "zero run past alphabet")
            s += r + 1
            continue
        _need(cum + v <= FSE_L, "counts exceed 512")
        norm[s] = v
        cum += v
        s += 1
    _need(cum == FSE_L, "counts do not sum to 512")
    return norm, k


def fse_table(norm):
    """[(symbol, nbits, base)] for the 512 states (section 6.4.2)."""
    spread, pos = [0] * FSE_L, 0
    step = (FSE_L >> 1) + (FSE_L >> 3) + 3
    for s, f in enumerate(norm):
        for _ in range(f):
            spread[pos] = s
            pos = (pos + step) & (FSE_L - 1)
    nxt = list(norm)
    out = []
    for u in range(FSE_L):
        s = spread[u]
        x = nxt[s]
        nxt[s] += 1
        nb = FSE_LOG - (x.bit_length() - 1)
        out.append((s, nb, (x << nb) - FSE_L))
    return out


# ------------------------------------------------- 6.5 value codes
def code_value(c, direct_log):
    """(base, extra bits) of code c with 2^direct_log direct codes."""
    d = 1 << direct_log
    if c < d:
        return c, 0
    k = c - d
    nb = direct_log + k // 2
    return (1 << nb) | ((k & 1) << (nb - 1)), nb - 1


LL_CODES, ML_CODES, OF_CODES = 62, 58, 65
ML_MIN = 2
BLOCK_MAX = 1 << 18
MAX_SEQ = BLOCK_MAX // 2 + 1


# ------------------------------------------------- decoder state
class Frame:
    def __init__(self, version, wlog, flags, content_size, features=None):
        self.feat = features if features is not None else set()
        self.version, self.wsize, self.flags, self.content_size = version, 1 << wlog, flags, content_size
        self.out = bytearray()            # decoded (filtered) content: the history
        self.reps = [1, 4, 8]             # 5.3
        self.huf = None                   # type 5 literal table (6.3.3 mode 2/3)
        self.ctx = None                   # type 6 tables: (cmap, [tables])
        self.fse = [None, None, None]     # sequence tables (6.4.3 mode 0/1)
        self.tok = TokModel()             # type 4 model (7)


# ------------------------------------------------- 6. block types 5 and 6
def huf_literals(p, k, n, fr):
    """Type 5 literal section (6.3.3): returns (literals, new k)."""
    _need(k < len(p), "truncated literal section")
    mode = p[k]
    k += 1
    fr.feat.add("t5.lit.mode%d" % mode)
    if mode == 0:
        _need(k + n <= len(p), "truncated raw literals")
        return bytes(p[k:k + n]), k + n
    if mode == 1:
        _need(k < len(p), "truncated run literal")
        return bytes([p[k]]) * n, k + 1
    if mode == 2:
        lens, k = read_lengths(p, k)
        fr.huf = None
        fr.huf = huf_table(lens)
    elif mode != 3 or fr.huf is None:
        raise FormatError("bad literal mode %d" % mode)
    _need(k + 12 <= len(p), "truncated stream sizes")
    sizes = list(struct.unpack_from("<3I", p, k))
    k += 12
    _need(sum(sizes) <= len(p) - k, "stream sizes past the payload")
    q = (n + 3) // 4
    out = bytearray()
    start = k
    for st in range(4):
        cnt = max(0, min((st + 1) * q, n) - min(st * q, n))
        end = start + sizes[st] if st < 3 else len(p)
        bits = FwdBits(p, start, end)
        for _ in range(cnt):
            out.append(huf_sym(bits, fr.huf))
        used = bits.used_bytes()
        if st < 3:
            _need(used == sizes[st], "Huffman stream %d does not end at its size" % st)
        start += used
    return bytes(out), start


def ctx_header(p, k, fr):
    """Type 6 context literal tables (6.6)."""
    _need(k < len(p), "truncated context literal section")
    mode = p[k]
    k += 1
    fr.feat.add("t6.ctx.mode%d" % mode)
    if mode == 1:
        _need(fr.ctx is not None, "no previous context tables")
        return k
    _need(mode == 0, "bad context literal mode %d" % mode)
    fr.ctx = None
    _need(k + 129 <= len(p), "truncated context map")
    K = p[k]
    _need(1 <= K <= 8, "bad table count %d" % K)
    cmap = [(p[k + 1 + i // 2] >> (4 * (i & 1))) & 15 for i in range(256)]
    _need(max(cmap) < K, "context map entry >= K")
    k += 129
    tables = []
    for _ in range(K):
        lens, k = read_lengths(p, k)
        tables.append(huf_table(lens))
    fr.ctx = (cmap, tables)
    return k


def decode_seq_block(p, raw, fr, ctx):
    _need(len(p) >= 8, "truncated block")
    ns, nlit = struct.unpack_from("<II", p, 0)
    _need(nlit <= raw and ns <= MAX_SEQ and ns * ML_MIN + nlit <= raw, "bad sequence/literal counts")
    k = 8
    lits, lbits = b"", None
    if ctx:
        _need(nlit > 0, "type 6 block without literals")
        k = ctx_header(p, k, fr)
        _need(k + 4 <= len(p), "truncated literal stream size")
        sz = struct.unpack_from("<I", p, k)[0]
        k += 4
        _need(sz <= len(p) - k, "literal stream past the payload")
        lbits = FwdBits(p, k, k + sz)
        lend = k + sz
        k += sz
        cmap, tables = fr.ctx
        prev = [fr.out[-1] if fr.out else 0]

        def take(n):
            o = bytearray()
            for _ in range(n):
                b = huf_sym(lbits, tables[cmap[prev[0]]])
                o.append(b)
                prev[0] = b
            return o
    elif nlit:
        lits, k = huf_literals(p, k, nlit, fr)
    out, start, li = fr.out, len(fr.out), 0

    if ns == 0:
        fr.feat.add("seq.none")
        _need(k == len(p) and nlit == raw, "literal-only block has trailing bytes")
        out += take(nlit) if ctx else lits
        if ctx:
            _need(lbits.used_bytes() == sz, "literal stream does not end at its size")
        return

    alpha = (LL_CODES, ML_CODES, OF_CODES)
    tabs = []
    for t in range(3):
        _need(k < len(p), "truncated table modes")
        mode = p[k]
        k += 1
        fr.feat.add("seq.%s.mode%d" % (("ll", "ml", "of")[t], mode))
        if mode == 0:
            fr.fse[t] = None
            norm, k = read_norm(p, k, alpha[t])
            fr.fse[t] = fse_table(norm)
            tabs.append((fr.fse[t], FSE_LOG))
        elif mode == 1:
            _need(fr.fse[t] is not None, "repeat mode without a table")
            tabs.append((fr.fse[t], FSE_LOG))
        elif mode == 2:
            _need(k < len(p) and p[k] < alpha[t], "bad single symbol")
            tabs.append(([(p[k], 0, 0)], 0))        # one state; mode 2 leaves fr.fse[t] alone
            k += 1
        else:
            raise FormatError("bad table mode %d" % mode)

    bits = BackBits(p, k, len(p))
    sll, sml, sof = (bits.read(tabs[t][1]) for t in range(3))
    for _ in range(ns):
        cl, nl, bl = tabs[0][0][sll]
        cm, nm, bm = tabs[1][0][sml]
        co, no, bo = tabs[2][0][sof]
        if co < 3:
            ob, ox = co, 0
        else:
            ob, ox = code_value(co - 3, 4)
            ob += 4                                  # value = distance + 3
        lb, lx = code_value(cl, 4)
        mb, mx = code_value(cm, 5)
        ofv = ob + bits.read(ox)
        ml = mb + ML_MIN + bits.read(mx)
        ll = lb + bits.read(lx)
        sll = bl + bits.read(nl)
        sml = bm + bits.read(nm)
        sof = bo + bits.read(no)
        r = fr.reps
        fr.feat.add("seq.of.%s" % ("new" if ofv >= 3 else "rep%d" % ofv))
        if ofv >= 3:
            dist = ofv - 3
            fr.reps = [dist, r[0], r[1]]
        elif ofv == 0:
            dist = r[0]
        elif ofv == 1:
            dist = r[1]
            fr.reps = [dist, r[0], r[2]]
        else:
            dist = r[2]
            fr.reps = [dist, r[0], r[1]]
        _need(ll <= nlit - li and ll + ml <= raw - (len(out) - start), "sequence past the block")
        if ctx:
            out += take(ll)
        else:
            out += lits[li:li + ll]
        li += ll
        _need(dist <= len(out) and dist <= fr.wsize, "distance too far")
        for _ in range(ml):
            out.append(out[-dist])
        if ctx:
            prev[0] = out[-1]
    _need(bits.pos == 0, "sequence bit stream not consumed exactly")
    _need(nlit - li == raw - (len(out) - start), "literal count mismatch")
    if ctx:
        out += take(nlit - li)
        _need(lbits.used_bytes() == sz, "literal stream does not end at its size")
    else:
        out += lits[li:]


# ------------------------------------------------- 7. block type 4
P_ONE, P_SHIFT = 1 << 12, 5


class TokModel:
    def __init__(self):
        h = P_ONE // 2
        self.is_match = [[h] * 4 for _ in range(12)]
        self.is_rep, self.is_rep0, self.is_rep1 = [h] * 12, [h] * 12, [h] * 12
        self.is_rep0_long = [[h] * 4 for _ in range(12)]
        self.lit = [[h] * 0x300 for _ in range(8)]
        self.slot = [[h] * 64 for _ in range(4)]
        self.spec = [[h] * 32 for _ in range(14)]
        self.align = [h] * 16
        self.mlen = LenModel()
        self.rlen = LenModel()
        self.state = 0


class LenModel:
    def __init__(self):
        h = P_ONE // 2
        self.choice, self.choice2 = [h], [h]
        self.low = [[h] * 8 for _ in range(4)]
        self.mid = [[h] * 8 for _ in range(4)]
        self.high = [h] * 256


class RangeDec:
    def __init__(self, p):
        _need(len(p) >= 5 and p[0] == 0, "range coder: bad first byte")
        self.p, self.i = p, 5
        self.range, self.code = 0xFFFFFFFF, int.from_bytes(p[1:5], "big")

    def _norm(self):
        if self.range < 1 << 24:
            self.range = (self.range << 8) & 0xFFFFFFFF
            b = self.p[self.i] if self.i < len(self.p) else 0
            self.i += 1
            self.code = ((self.code << 8) | b) & 0xFFFFFFFF

    def bit(self, probs, idx):
        pr = probs[idx]
        bound = (self.range >> 12) * pr
        if self.code < bound:
            self.range = bound
            probs[idx] = pr + ((P_ONE - pr) >> P_SHIFT)
            b = 0
        else:
            self.range -= bound
            self.code -= bound
            probs[idx] = pr - (pr >> P_SHIFT)
            b = 1
        self._norm()
        return b

    def direct(self, n):
        v = 0
        for _ in range(n):
            self.range >>= 1
            b = 1 if self.code >= self.range else 0
            if b:
                self.code -= self.range
            v = (v << 1) | b
            self._norm()
        return v

    def tree(self, probs, nbits):
        m = 1
        for _ in range(nbits):
            m = (m << 1) | self.bit(probs, m)
        return m - (1 << nbits)

    def rtree(self, probs, nbits):
        m, v = 1, 0
        for i in range(nbits):
            b = self.bit(probs, m)
            m = (m << 1) | b
            v |= b << i
        return v


def _len(rc, lm, ps):
    if not rc.bit(lm.choice, 0):
        return 2 + rc.tree(lm.low[ps], 3)
    if not rc.bit(lm.choice2, 0):
        return 10 + rc.tree(lm.mid[ps], 3)
    return 18 + rc.tree(lm.high, 8)


def decode_tok_block(p, raw, fr):
    rc, m, out = RangeDec(p), fr.tok, fr.out
    start = len(out)
    while len(out) - start < raw:
        _need(rc.i <= len(p), "range coder overrun")
        ps, s = (len(out) - start) & 3, m.state
        if not rc.bit(m.is_match[s], ps):                          # literal
            probs = m.lit[(out[-1] if out else 0) >> 5]
            fr.feat.add("t4.lit.%s" % ("matched" if s >= 7 else "plain"))
            if s >= 7:
                _need(fr.reps[0] <= len(out), "matched literal before content start")
                mbyte, v, same = out[-fr.reps[0]], 1, True
                for i in range(7, -1, -1):
                    mb = (mbyte >> i) & 1
                    b = rc.bit(probs, 0x100 + (mb << 8) + v if same else v)
                    v = (v << 1) | b
                    same = same and b == mb
                lit = v & 0xFF
            else:
                lit = rc.tree(probs, 8)
            out.append(lit)
            m.state = 0 if s < 4 else (s - 3 if s < 10 else s - 6)
            continue
        r = fr.reps
        if not rc.bit(m.is_rep, s):                                # new distance
            ln = _len(rc, m.mlen, ps)
            slot = rc.tree(m.slot[min(ln - 2, 3)], 6)
            d0 = slot
            fr.feat.add("t4.slot.%s" % ("direct" if slot < 4 else "tree" if slot < 14 else "align"))
            fr.feat.add("t4.len.%s" % ("low" if ln < 10 else "mid" if ln < 18 else "high"))
            if slot >= 4:
                nb = (slot >> 1) - 1
                d0 = (2 | (slot & 1)) << nb
                if slot < 14:
                    d0 += rc.rtree(m.spec[slot], nb)
                else:
                    d0 += rc.direct(nb - 4) << 4
                    d0 += rc.rtree(m.align, 4)
            _need(d0 < fr.wsize, "distance past the window")
            dist = d0 + 1
            fr.reps = [dist, r[0], r[1]]
            m.state = 7 if s < 7 else 10
        else:
            if not rc.bit(m.is_rep0, s):
                dist = r[0]
                if not rc.bit(m.is_rep0_long[s], ps):              # short rep
                    fr.feat.add("t4.shortrep")
                    _need(dist <= len(out), "distance before content start")
                    out.append(out[-dist])
                    m.state = 9 if s < 7 else 11
                    continue
            elif not rc.bit(m.is_rep1, s):
                fr.feat.add("t4.rep1")
                dist = r[1]
                fr.reps = [dist, r[0], r[2]]
            else:
                fr.feat.add("t4.rep2")
                dist = r[2]
                fr.reps = [dist, r[0], r[1]]
            ln = _len(rc, m.rlen, ps)
            fr.feat.add("t4.rep" if dist == r[0] else "t4.repN")
            m.state = 8 if s < 7 else 11
        _need(ln <= raw - (len(out) - start), "match past the block")
        _need(dist <= len(out) and dist <= fr.wsize, "distance too far")
        for _ in range(ln):
            out.append(out[-dist])
    _need(rc.i == len(p), "range coder did not end at the payload end")


# ------------------------------------------------- 8. x86 filter
def x86_unfilter(buf):
    b, i, n = bytearray(buf), 0, len(buf)
    while i < n:
        if b[i] != 0xE8:
            i += 1
            continue
        if n - i < 5:
            break
        v = struct.unpack_from("<I", b, i + 1)[0]
        top = v >> 25
        if top == 0 or top == 0x7F:
            a = (v - (i + 5)) & 0x3FFFFFF
            if a >> 25:
                a |= 0xFC000000
            struct.pack_into("<I", b, i + 1, a)
        i += 5
    return bytes(b)


# ------------------------------------------------- 3. frame
def parse_header(data):
    _need(len(data) >= 8 and data[:4] == b"LZRN", "bad magic or truncated header")
    ver, flags, wlog, res = data[4], data[5], data[6], data[7]
    if ver < 5 or ver > 6:
        raise UnsupportedError("frame version %d" % ver)
    if flags & ~0x06:
        raise UnsupportedError("flags 0x%02x" % flags)
    _need(18 <= wlog <= 26 and res == 0, "bad window_log or reserved byte")
    if flags & 2:
        _need(len(data) >= 16, "truncated header")
        size = struct.unpack_from("<Q", data, 8)[0]
        _need(size < 1 << 63, "content size too large")
        return ver, flags, wlog, size, 16
    return ver, flags, wlog, None, 8


def blocks(data):
    ver, flags, wlog, size, pos = parse_header(data)
    yield ("header", ver, flags, wlog, size)
    while True:
        _need(pos + 9 <= len(data), "truncated block header")
        t, raw, comp = data[pos], *struct.unpack_from("<II", data, pos + 1)
        _need(t in (0, 2, 4, 5, 6) and not (t == 6 and ver < 6), "bad block type %d" % t)
        _need(raw <= BLOCK_MAX and comp <= 2 * BLOCK_MAX + 65536, "block too large")
        _need(t != 0 or comp == raw, "stored block: comp_len != raw_len")
        _need(t != 2 or (raw == 0 and comp == 4), "bad end block")
        _need(t < 4 or raw > 0, "empty compressed block")
        _need(pos + 9 + comp <= len(data), "truncated block")
        pos += 9 + comp
        yield ("block", t, raw, data[pos - comp:pos], pos)
        if t == 2:
            return


def decode(data, allow_trailing=False, features=None):
    """Decode one frame. Returns the content; raises FormatError.
    features: optional set, filled with the format features met (for
    test-vector coverage)."""
    data = bytes(data)
    fr = None
    for ev in blocks(data):
        if ev[0] == "header":
            _, ver, flags, wlog, size = ev
            fr = Frame(ver, wlog, flags, size, features)
            fr.feat.update({"version%d" % ver, "wlog%d" % wlog})
            if flags & 2:
                fr.feat.add("flag.size")
            if flags & 4:
                fr.feat.add("flag.x86")
            continue
        _, t, raw, p, end = ev
        if t == 2:
            _need(allow_trailing or end == len(data), "trailing data after the frame")
            _need(fr.content_size is None or fr.content_size == len(fr.out), "content size mismatch")
            content = x86_unfilter(fr.out) if fr.flags & 4 else bytes(fr.out)
            if struct.unpack("<I", p)[0] != xxh64(content) & 0xFFFFFFFF:
                raise ChecksumError("checksum mismatch")
            return content
        _need(fr.content_size is None or len(fr.out) + raw <= fr.content_size, "more content than announced")
        before = len(fr.out)
        fr.feat.add("block%d" % t)
        if t == 0:
            fr.out += p
        elif t == 4:
            decode_tok_block(p, raw, fr)
        else:
            decode_seq_block(p, raw, fr, t == 6)
        assert len(fr.out) - before == raw
    raise FormatError("no end block")


def trace(data):
    """Human-readable block list (does not decode)."""
    lines = []
    for ev in blocks(bytes(data)):
        if ev[0] == "header":
            _, ver, flags, wlog, size = ev
            lines.append("header version=%d flags=0x%02x window_log=%d content_size=%s" % (ver, flags, wlog, size))
        elif ev[0] == "block":
            _, t, raw, p, _ = ev
            extra = ""
            if t in (5, 6) and len(p) >= 8:
                ns, nl = struct.unpack_from("<II", p, 0)
                extra = " nseq=%d nlit=%d lit_mode=%d" % (ns, nl, p[8] if nl else -1)
            lines.append("block type=%d raw=%d comp=%d%s" % (t, raw, len(p), extra))
    return "\n".join(lines)


if __name__ == "__main__":
    args = sys.argv[1:]
    if args and args[0] == "--trace":
        print(trace(open(args[1], "rb").read()))
    elif len(args) == 1:
        sys.stdout.buffer.write(decode(open(args[0], "rb").read()))
    else:
        print(__doc__)
        sys.exit(2)
