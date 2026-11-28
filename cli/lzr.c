/*
 * lzr - command-line front end of liblzr, used like xz / zstd.
 *
 *   lzr [OPTION]... [FILE]...       compress FILE to FILE.lzr (removes FILE unless -k)
 *   lzr -d FILE.lzr                 decompress to FILE
 *   lzr -c < in > out               stdin to stdout (so `tar -I lzr` works)
 *   unlzr = lzr -d, lzrcat = lzr -dc
 *
 * Structure: option parsing (main) -> per-file driver (do_file: names,
 * open/stat, output file, metadata, removal) -> one of the stream loops:
 * compress_fd, decompress_fd (streaming, bounded memory, or the parallel
 * one-shot path for mapped files with threads), list_fd (block headers
 * only). Exit status as xz: 0 ok, 1 error, 2 warning.
 */
#define _GNU_SOURCE
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <signal.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#include "lzr.h"

#define IO_SIZE ((size_t)1 << 20)

enum { M_COMPRESS, M_DECOMPRESS, M_TEST, M_LIST };

static struct {
    int mode, to_stdout, keep, force, verbose;   /* verbose: -1 quiet .. 2 */
    int level, threads, wlog, filter, maxw;
    int64_t segment;
    const char *suffix;
} opt = {M_COMPRESS, 0, 0, 0, 0, LZR_LEVEL_DEFAULT, 1, 0, LZR_FILTER_NONE, LZR_MAX_WINDOW_LOG_DEFAULT, -1, ".lzr"};

static const char *prog = "lzr";
static int status = 0;                        /* exit status */
static volatile sig_atomic_t sig_seen = 0;
static char *out_path_live = NULL;           /* output file to remove if interrupted */

static void msg(int level, const char *name, const char *fmt, ...)
{
    if (opt.verbose < 0 && level >= 1) return;   /* -q hides warnings */
    va_list ap;
    va_start(ap, fmt);
    fprintf(stderr, "%s: ", prog);
    if (name) fprintf(stderr, "%s: ", name);
    vfprintf(stderr, fmt, ap);
    fputc('\n', stderr);
    va_end(ap);
}
static void err(const char *name, const char *fmt, ...) __attribute__((format(printf, 2, 3)));
static void err(const char *name, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    fprintf(stderr, "%s: ", prog);
    if (name) fprintf(stderr, "%s: ", name);
    vfprintf(stderr, fmt, ap);
    fputc('\n', stderr);
    va_end(ap);
    status = 1;
}
#define warn(name, ...) do { msg(1, name, __VA_ARGS__); if (status == 0) status = 2; } while (0)

static void on_signal(int s)
{
    sig_seen = s;
    if (out_path_live) unlink(out_path_live);
    signal(s, SIG_DFL);
    raise(s);
}

static double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}

static int write_all(int fd, const void *p, size_t n)
{
    const char *c = (const char *)p;
    while (n) {
        ssize_t w = write(fd, c, n > ((size_t)1 << 30) ? ((size_t)1 << 30) : n);
        if (w < 0) {
            if (errno == EINTR) continue;
            return -1;
        }
        c += w;
        n -= (size_t)w;
    }
    return 0;
}

/* read up to n bytes; short only at EOF. -1 on error. */
static ssize_t read_full(int fd, void *p, size_t n)
{
    size_t got = 0;
    while (got < n) {
        ssize_t r = read(fd, (char *)p + got, n - got);
        if (r < 0) {
            if (errno == EINTR) continue;
            return -1;
        }
        if (r == 0) break;
        got += (size_t)r;
    }
    return (ssize_t)got;
}

/* Input: a mapped regular file, or read in pieces. */
typedef struct {
    int fd;
    const uint8_t *map;     /* whole file when mapped */
    size_t map_len;
    int64_t size;           /* -1 = unknown (pipe) */
    const char *name;
} input;

typedef struct {
    uint64_t in, out;
} counts;

/* ---------------- compression ---------------- */

static int compress_fd(const input *src, int ofd, counts *ct, const char *oname)
{
    lzr_cstream *cs = lzr_cstream_create();
    uint8_t *obuf = malloc(IO_SIZE), *ibuf = NULL;
    int ret = -1, r = 0;
    if (!cs || !obuf) { err(src->name, "out of memory"); goto done; }
    if ((r = lzr_cstream_set(cs, LZR_C_LEVEL, opt.level)) || (r = lzr_cstream_set(cs, LZR_C_THREADS, opt.threads)) ||
        (r = lzr_cstream_set(cs, LZR_C_WINDOW_LOG, opt.wlog)) || (r = lzr_cstream_set(cs, LZR_C_SEGMENT, opt.segment)) ||
        (r = lzr_cstream_set(cs, LZR_C_FILTER, opt.filter)) || (r = lzr_cstream_set(cs, LZR_C_CONTENT_SIZE, src->size))) {
        err(src->name, "%s", lzr_error_string(r));
        goto done;
    }
    lzr_outbuf out = {obuf, IO_SIZE, 0};
    lzr_inbuf in = {NULL, 0, 0};
    int eof = 0;
    if (src->map || src->size == 0) {
        in.src = src->map;
        in.size = src->map_len;
        eof = 1;
    } else if (!(ibuf = malloc(IO_SIZE))) {
        err(src->name, "out of memory");
        goto done;
    }
    for (;;) {
        if (!eof && in.pos == in.size) {
            ssize_t n = read_full(src->fd, ibuf, IO_SIZE);
            if (n < 0) { err(src->name, "read error: %s", strerror(errno)); goto done; }
            in.src = ibuf;
            in.size = (size_t)n;
            in.pos = 0;
            eof = n == 0;
            ct->in += (uint64_t)n;
        }
        r = lzr_compress_stream(cs, &out, &in, eof ? LZR_FINISH : LZR_CONTINUE);
        if (r < 0) { err(src->name, "%s", lzr_error_string(r)); goto done; }
        if (out.pos && (r == 1 || out.pos == out.size || (eof && r == 0))) {
            if (write_all(ofd, obuf, out.pos)) { err(oname, "write error: %s", strerror(errno)); goto done; }
            ct->out += out.pos;
            out.pos = 0;
        }
        if (eof && r == 0) break;
    }
    if (src->map) ct->in = src->map_len;
    ret = 0;
done:
    lzr_cstream_free(cs);
    free(obuf);
    free(ibuf);
    return ret;
}

/* ---------------- decompression ---------------- */

static uint64_t mem_budget(void)
{
    long pages = sysconf(_SC_PHYS_PAGES), psz = sysconf(_SC_PAGESIZE);
    if (pages <= 0 || psz <= 0) return (uint64_t)1 << 30;
    return (uint64_t)pages * (uint64_t)psz / 2;
}

static int is_magic(const uint8_t *p, size_t n)
{
    return n >= 4 && memcmp(p, "LZRN", 4) == 0;
}

/* ofd < 0: test only. Handles concatenated frames. */
static int decompress_fd(const input *src, int ofd, counts *ct, const char *oname)
{
    lzr_dstream *ds = lzr_dstream_create();
    uint8_t *obuf = malloc(IO_SIZE), *ibuf = NULL;
    int ret = -1, r;
    if (!ds || !obuf) { err(src->name, "out of memory"); goto done; }
    if ((r = lzr_dstream_set(ds, LZR_D_MAX_WINDOW_LOG, opt.maxw)) || (r = lzr_dstream_set(ds, LZR_D_THREADS, opt.threads))) {
        err(src->name, "%s", lzr_error_string(r));
        goto done;
    }
    lzr_inbuf in = {NULL, 0, 0};
    int eof = 0;
    if (src->map) {
        in.src = src->map;
        in.size = src->map_len;
        eof = 1;
    } else if (!(ibuf = malloc(IO_SIZE))) {
        err(src->name, "out of memory");
        goto done;
    }
    uint64_t frames = 0, consumed_before = 0;
    int in_frame = 0;   /* bytes of the current frame have been fed */
    for (;;) {
        if (!eof && in.pos == in.size) {
            consumed_before += in.size;
            ssize_t n = read_full(src->fd, ibuf, IO_SIZE);
            if (n < 0) { err(src->name, "read error: %s", strerror(errno)); goto done; }
            in.src = ibuf;
            in.size = (size_t)n;
            in.pos = 0;
            eof = n == 0;
        }
        if (!in_frame) {
            if (in.pos == in.size) {
                if (eof) break;
                continue;
            }
            /* Need 4 bytes to tell a frame from trailing garbage; a frame
             * split across reads is just fed (the decoder checks it). */
            const uint8_t *p = (const uint8_t *)in.src + in.pos;
            size_t avail = in.size - in.pos;
            if (avail >= 4 ? !is_magic(p, avail) : (eof && memcmp(p, "LZRN", avail) != 0)) {
                if (frames) warn(src->name, "trailing garbage after the last frame ignored");
                else err(src->name, "file format not recognized");
                if (!frames) goto done;
                break;
            }
            /* Parallel one-shot: mapped input, threads, size in header, fits in memory. */
            lzr_frame_info fi;
            if (src->map && opt.threads != 1 && lzr_get_frame_info(p, avail, &fi) == LZR_OK && fi.content_size >= 0 &&
                (uint64_t)fi.content_size + LZR_DECODE_MARGIN <= mem_budget()) {
                size_t cap = (size_t)fi.content_size + LZR_DECODE_MARGIN, len = 0, used = 0;
                uint8_t *buf = malloc(cap);
                if (buf) {
                    r = lzr_decompress_ds(ds, buf, cap, &len, p, avail, &used);
                    if (r) { free(buf); err(src->name, "%s", lzr_error_string(r)); goto done; }
                    if (ofd >= 0 && write_all(ofd, buf, len)) {
                        free(buf);
                        err(oname, "write error: %s", strerror(errno));
                        goto done;
                    }
                    free(buf);
                    ct->out += len;
                    in.pos += used;
                    frames++;
                    continue;
                }
            }
            lzr_dstream_reset(ds);
            in_frame = 1;
        }
        lzr_outbuf out = {obuf, IO_SIZE, 0};
        r = lzr_decompress_stream(ds, &out, &in);
        if (r < 0) { err(src->name, "%s", lzr_error_string(r)); goto done; }
        if (out.pos) {
            if (ofd >= 0 && write_all(ofd, obuf, out.pos)) { err(oname, "write error: %s", strerror(errno)); goto done; }
            ct->out += out.pos;
        }
        if (r == LZR_FRAME_END) {
            in_frame = 0;
            frames++;
        } else if (eof && in.pos == in.size && out.pos < out.size) {
            err(src->name, "%s", lzr_error_string(LZR_ERR_TRUNCATED));
            goto done;
        }
    }
    if (!frames) {   /* empty input */
        err(src->name, "%s", lzr_error_string(LZR_ERR_TRUNCATED));
        goto done;
    }
    ct->in = consumed_before + in.size;
    ret = 0;
done:
    lzr_dstream_free(ds);
    free(obuf);
    free(ibuf);
    return ret;
}

/* ---------------- list (-l): walks block headers, decodes nothing ---------------- */

static uint32_t rd32(const uint8_t *p) { return p[0] | p[1] << 8 | p[2] << 16 | (uint32_t)p[3] << 24; }

static void human(char *b, size_t bn, uint64_t v)
{
    const char *u[] = {"B", "KiB", "MiB", "GiB", "TiB"};
    double d = (double)v;
    int i = 0;
    while (d >= 1024 && i < 4) { d /= 1024; i++; }
    if (i == 0) snprintf(b, bn, "%llu B", (unsigned long long)v);
    else snprintf(b, bn, "%.1f %s", d, u[i]);
}

static int list_fd(const input *src, int header)
{
    if (!src->map) {
        err(src->name, "--list needs a regular file");
        return -1;
    }
    const uint8_t *p = src->map;
    size_t n = src->map_len, pos = 0;
    uint64_t frames = 0, segs = 0, raw = 0;
    int wlog = 0, ver = 0;
    unsigned flags = 0;
    while (pos < n) {
        lzr_frame_info fi;
        int r = lzr_get_frame_info(p + pos, n - pos, &fi);
        if (r) {
            if (frames) { warn(src->name, "trailing garbage after the last frame"); break; }
            err(src->name, "%s", r == LZR_ERR_TRUNCATED ? lzr_error_string(r) : "file format not recognized");
            return -1;
        }
        if (fi.window_log > wlog) wlog = fi.window_log;
        ver = fi.version;
        flags |= fi.flags;
        pos += fi.header_size;
        segs++;   /* a frame is one segment plus one per reset point after its first block */
        int end = 0, first = 1;
        while (!end) {
            if (n - pos < 9) { err(src->name, "%s", lzr_error_string(LZR_ERR_TRUNCATED)); return -1; }
            unsigned t = p[pos];
            uint32_t rl = rd32(p + pos + 1), cl = rd32(p + pos + 5);
            if (n - pos - 9 < cl) { err(src->name, "%s", lzr_error_string(LZR_ERR_TRUNCATED)); return -1; }
            if (fi.version >= 7) {
                if (!(t & 0x80)) {
                    if ((t & 0x40) && !first) segs++;
                    first = 0;
                    end = (t & 0x3F) == 2;
                    raw += rl;
                }
            } else {
                end = t == 2;
                raw += rl;
            }
            pos += 9 + (size_t)cl;
        }
        frames++;
    }
    if (header && opt.verbose >= 0)
        printf("%6s %8s %12s %12s %7s %6s  %s\n", "Frames", "Segments", "Compressed", "Uncompressed", "Ratio", "Window", "Filename");
    char cb[32], ub[32], wb[32];
    human(cb, sizeof cb, n);
    human(ub, sizeof ub, raw);
    human(wb, sizeof wb, (uint64_t)1 << wlog);
    printf("%6llu %8llu %12s %12s %7.3f %6s  %s%s%s\n", (unsigned long long)frames,
           (unsigned long long)segs, cb, ub, n ? (double)raw / (double)n : 0.0, wb, src->name,
           (flags & LZR_FRAME_FLAG_X86) ? " [x86]" : "", "");
    if (opt.verbose > 0) printf("    format version %d, checksum XXH64 (low 32 bits)\n", ver);
    return 0;
}

/* ---------------- per-file driver ---------------- */

static int has_suffix(const char *name)
{
    size_t n = strlen(name), s = strlen(opt.suffix);
    return n > s && strcmp(name + n - s, opt.suffix) == 0;
}

static void report(const char *name, const counts *ct, double secs)
{
    if (opt.verbose < 1) return;
    uint64_t raw = opt.mode == M_COMPRESS ? ct->in : ct->out;
    uint64_t comp = opt.mode == M_COMPRESS ? ct->out : ct->in;
    double mbs = secs > 0 ? raw / secs / 1e6 : 0;
    fprintf(stderr, "%s: %s: %llu -> %llu bytes (%.2f%%, ratio %.3f), %.1f MB/s\n", prog, name,
            (unsigned long long)ct->in, (unsigned long long)ct->out, raw ? 100.0 * comp / raw : 0.0,
            comp ? (double)raw / comp : 0.0, mbs);
}

/* name == NULL: stdin. */
static void do_file(const char *name, int first)
{
    input src = {0, NULL, 0, -1, name ? name : "(stdin)"};
    struct stat st;
    int use_stdout = opt.to_stdout || !name || opt.mode == M_TEST || opt.mode == M_LIST;
    char *oname = NULL;
    int ofd = -1;
    counts ct = {0, 0};

    if (name) {
        if (opt.mode == M_COMPRESS && !use_stdout && has_suffix(name) && !opt.force) {
            warn(name, "already has '%s' suffix, skipped", opt.suffix);
            return;
        }
        if (opt.mode == M_DECOMPRESS && !use_stdout && !has_suffix(name)) {
            warn(name, "unknown suffix, skipped (use -c or -S)");
            return;
        }
        src.fd = open(name, O_RDONLY | O_NOCTTY);
        if (src.fd < 0) { err(name, "%s", strerror(errno)); return; }
    } else {
        src.fd = STDIN_FILENO;
    }
    if (fstat(src.fd, &st)) { err(src.name, "%s", strerror(errno)); goto out; }
    if (name && !S_ISREG(st.st_mode) && !use_stdout) {
        warn(name, "not a regular file, skipped");
        goto out;
    }
    if (name && S_ISREG(st.st_mode) && st.st_nlink > 1 && !use_stdout && !opt.keep && !opt.force) {
        warn(name, "has %lu other hard links, skipped (use -k or -f)", (unsigned long)st.st_nlink - 1);
        goto out;
    }
    if (S_ISREG(st.st_mode)) {
        src.size = st.st_size;
        if (st.st_size > 0) {
            void *m = mmap(NULL, (size_t)st.st_size, PROT_READ, MAP_PRIVATE, src.fd, 0);
            if (m != MAP_FAILED) {
                src.map = m;
                src.map_len = (size_t)st.st_size;
                madvise(m, src.map_len, MADV_SEQUENTIAL);
            } else {
                src.size = opt.mode == M_COMPRESS ? st.st_size : -1;   /* fall back to reads */
            }
        }
    }
    if (opt.mode == M_LIST) {
        list_fd(&src, first);
        goto out;
    }
    if (opt.mode != M_COMPRESS && isatty(src.fd)) {
        err(src.name, "compressed data not read from a terminal (use -f to force)");
        goto out;
    }

    if (opt.mode == M_TEST) {
        ofd = -1;
    } else if (use_stdout) {
        ofd = STDOUT_FILENO;
        if (opt.mode == M_COMPRESS && isatty(ofd) && !opt.force) {
            err(NULL, "compressed data not written to a terminal (use -f to force)");
            goto out;
        }
    } else {
        size_t n = strlen(name), s = strlen(opt.suffix);
        oname = malloc(n + s + 1);
        if (!oname) { err(name, "out of memory"); goto out; }
        if (opt.mode == M_COMPRESS) {
            memcpy(oname, name, n);
            memcpy(oname + n, opt.suffix, s + 1);
        } else {
            memcpy(oname, name, n - s);
            oname[n - s] = 0;
        }
        if (opt.force) unlink(oname);
        ofd = open(oname, O_WRONLY | O_CREAT | O_EXCL | O_NOCTTY, 0600);
        if (ofd < 0) {
            err(oname, "%s", strerror(errno));
            free(oname);
            oname = NULL;
            goto out;
        }
        out_path_live = oname;
    }

    double t0 = now();
    int r = opt.mode == M_COMPRESS ? compress_fd(&src, ofd, &ct, oname ? oname : "(stdout)")
                                   : decompress_fd(&src, ofd, &ct, oname ? oname : "(stdout)");
    if (!r && oname) {
        /* metadata as xz: permissions, owner (best effort), times */
        if (fchown(ofd, (uid_t)-1, st.st_gid)) { /* not a member of the group: keep ours */ }
        if (fchmod(ofd, st.st_mode & 07777)) warn(oname, "cannot set permissions: %s", strerror(errno));
        struct timespec ts[2] = {st.st_atim, st.st_mtim};
        if (futimens(ofd, ts)) warn(oname, "cannot set timestamps: %s", strerror(errno));
    }
    if (oname && close(ofd)) {
        if (!r) err(oname, "write error: %s", strerror(errno));
        r = -1;
    }
    if (oname) ofd = -1;
    if (r) {
        if (oname) unlink(oname);
    } else {
        report(src.name, &ct, now() - t0);
        if (oname && !opt.keep && unlink(name)) warn(name, "cannot remove: %s", strerror(errno));
    }
    out_path_live = NULL;
    free(oname);
out:
    if (src.map) munmap((void *)src.map, src.map_len);
    if (name && src.fd >= 0) close(src.fd);
}

/* ---------------- options ---------------- */

static void usage(FILE *f)
{
    fprintf(f,
            "Usage: %s [OPTION]... [FILE]...\n"
            "Compress or decompress FILEs in the .lzr format (LZR 1.0). With no FILE, or\n"
            "when FILE is -, read standard input and write standard output.\n\n"
            "  -z, --compress        compress (default)\n"
            "  -d, --decompress      decompress (also: unlzr, lzrcat)\n"
            "  -t, --test            decompress and verify, write nothing\n"
            "  -l, --list            list frames, segments and sizes of .lzr files\n"
            "  -c, --stdout          write to standard output, keep input files\n"
            "  -k, --keep            keep input files\n"
            "  -f, --force           overwrite output, compress to a terminal, etc.\n"
            "  -0 ... -11            compression level (default %d; 11 = best ratio)\n"
            "      --fast[=N]        level -N, N = 1..3 (fastest; default 1)\n"
            "      --level=N         level -3..11\n"
            "  -T, --threads=N       worker threads (default 1, 0 = one per CPU); compression\n"
            "                        output does not depend on N (levels -3..6: any N > 1)\n"
            "      --segment=SIZE    reset points every SIZE bytes for parallel decompression\n"
            "                        (0 = off, default 48MiB at levels 9-11, else off; >= 1MiB)\n"
            "      --window-log=N    match window 2^N bytes, N = 18..26 (default by level)\n"
            "      --x86             x86/x86-64 CALL filter (executables, libraries)\n"
            "      --max-window-log=N  refuse frames with a window above 2^N (18..31, default %d)\n"
            "  -S, --suffix=.SUF     suffix of compressed files (default .lzr)\n"
            "  -q, --quiet           suppress warnings\n"
            "  -v, --verbose         report ratio and speed\n"
            "  -h, --help            this help\n"
            "  -V, --version         version\n\n"
            "Examples:  lzr -9 -T4 big.tar      tar -I 'lzr -9 -T0' -cf x.tar.lzr dir\n"
            "           lzr -dc x.lzr | less    tar -I lzr -xf x.tar.lzr\n",
            prog, LZR_LEVEL_DEFAULT, LZR_MAX_WINDOW_LOG_DEFAULT);
}

static int parse_int(const char *s, long lo, long hi, long *v)
{
    char *e;
    errno = 0;
    long x = strtol(s, &e, 10);
    if (!*s || *e || errno || x < lo || x > hi) return -1;
    *v = x;
    return 0;
}

/* 123, 4K, 16M, 1G, 2MiB, 1KB (all 1024-based) */
static int parse_size(const char *s, int64_t *v)
{
    char *e;
    errno = 0;
    unsigned long long x = strtoull(s, &e, 10);
    if (e == s || errno || *s == '-') return -1;
    int sh = 0;
    switch (*e) {
    case 0: break;
    case 'k': case 'K': sh = 10; e++; break;
    case 'm': case 'M': sh = 20; e++; break;
    case 'g': case 'G': sh = 30; e++; break;
    default: return -1;
    }
    if (sh && (!strcmp(e, "iB") || !strcmp(e, "B") || !strcmp(e, "i"))) e += strlen(e);
    if (*e) return -1;
    if (x > ((unsigned long long)INT64_MAX >> sh)) return -1;
    *v = (int64_t)(x << sh);
    return 0;
}

static void die_usage(const char *fmt, const char *arg)
{
    fprintf(stderr, "%s: ", prog);
    fprintf(stderr, fmt, arg);
    fprintf(stderr, "\nTry '%s --help' for more information.\n", prog);
    exit(1);
}

/* long option with value: --name=value or --name value */
static const char *long_val(const char *a, const char *name, int *i, int argc, char **argv)
{
    size_t n = strlen(name);
    if (strncmp(a, name, n)) return NULL;
    if (a[n] == '=') return a + n + 1;
    if (a[n] == 0) {
        if (*i + 1 >= argc) die_usage("option '%s' needs a value", a);
        return argv[++*i];
    }
    return NULL;
}

static void set_threads(const char *v)
{
    long x;
    if (parse_int(v, 0, LZR_THREADS_MAX, &x)) die_usage("invalid thread count '%s' (0..16)", v);
    opt.threads = (int)x;
}

int main(int argc, char **argv)
{
    const char *base = strrchr(argv[0], '/');
    base = base ? base + 1 : argv[0];
    if (!strncmp(base, "unlzr", 5)) opt.mode = M_DECOMPRESS;
    if (!strncmp(base, "lzrcat", 6)) { opt.mode = M_DECOMPRESS; opt.to_stdout = 1; }

    static char **files;   /* static: still reachable on the early exits below */
    files = calloc((size_t)argc + 1, sizeof(char *));
    int nfiles = 0, only_files = 0;
    if (!files) return 1;
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (only_files || a[0] != '-' || a[1] == 0) { files[nfiles++] = argv[i]; continue; }
        if (a[1] == '-') {
            const char *v;
            long x;
            if (!a[2]) only_files = 1;
            else if (!strcmp(a, "--compress")) opt.mode = M_COMPRESS;
            else if (!strcmp(a, "--decompress") || !strcmp(a, "--uncompress")) opt.mode = M_DECOMPRESS;
            else if (!strcmp(a, "--test")) opt.mode = M_TEST;
            else if (!strcmp(a, "--list")) opt.mode = M_LIST;
            else if (!strcmp(a, "--stdout") || !strcmp(a, "--to-stdout")) opt.to_stdout = 1;
            else if (!strcmp(a, "--keep")) opt.keep = 1;
            else if (!strcmp(a, "--force")) opt.force = 1;
            else if (!strcmp(a, "--quiet")) opt.verbose = -1;
            else if (!strcmp(a, "--verbose")) opt.verbose++;
            else if (!strcmp(a, "--x86")) opt.filter = LZR_FILTER_X86;
            else if (!strcmp(a, "--help")) { usage(stdout); return 0; }
            else if (!strcmp(a, "--version")) {
                printf("lzr %s (liblzr %s, format LZR 1.0)\n", lzr_version_string(), lzr_version_string());
                return 0;
            } else if (!strcmp(a, "--fast")) opt.level = -1;
            else if (!strncmp(a, "--fast=", 7)) {
                if (parse_int(a + 7, 1, -LZR_LEVEL_MIN, &x)) die_usage("invalid --fast value '%s' (1..3)", a + 7);
                opt.level = (int)-x;
            } else if ((v = long_val(a, "--level", &i, argc, argv))) {
                if (parse_int(v, LZR_LEVEL_MIN, LZR_LEVEL_MAX, &x)) die_usage("invalid level '%s' (-3..11)", v);
                opt.level = (int)x;
            } else if ((v = long_val(a, "--threads", &i, argc, argv))) set_threads(v);
            else if ((v = long_val(a, "--segment", &i, argc, argv))) {
                if (parse_size(v, &opt.segment) || (opt.segment && (uint64_t)opt.segment < LZR_SEGMENT_MIN))
                    die_usage("invalid segment size '%s' (0, or >= 1MiB)", v);
            } else if ((v = long_val(a, "--window-log", &i, argc, argv))) {
                if (parse_int(v, LZR_WINDOW_LOG_MIN, LZR_WINDOW_LOG_ENC_MAX, &x)) die_usage("invalid window log '%s' (18..26)", v);
                opt.wlog = (int)x;
            } else if ((v = long_val(a, "--max-window-log", &i, argc, argv))) {
                if (parse_int(v, LZR_WINDOW_LOG_MIN, LZR_WINDOW_LOG_MAX, &x)) die_usage("invalid max window log '%s' (18..31)", v);
                opt.maxw = (int)x;
            } else if ((v = long_val(a, "--suffix", &i, argc, argv))) {
                if (!*v) die_usage("empty suffix%s", "");
                opt.suffix = v;
            } else die_usage("unknown option '%s'", a);
            continue;
        }
        /* short options, clustered: -9kT4, -dc, -11 */
        for (const char *p = a + 1; *p; p++) {
            if (*p >= '0' && *p <= '9') {
                long lv = 0;
                while (*p >= '0' && *p <= '9') lv = lv * 10 + (*p++ - '0'), lv = lv > 99 ? 99 : lv;
                p--;
                if (lv > LZR_LEVEL_MAX) die_usage("invalid level in '%s' (0..11; --fast=N for faster)", a);
                opt.level = (int)lv;
                continue;
            }
            switch (*p) {
            case 'z': opt.mode = M_COMPRESS; break;
            case 'd': opt.mode = M_DECOMPRESS; break;
            case 't': opt.mode = M_TEST; break;
            case 'l': opt.mode = M_LIST; break;
            case 'c': opt.to_stdout = 1; break;
            case 'k': opt.keep = 1; break;
            case 'f': opt.force = 1; break;
            case 'q': opt.verbose = -1; break;
            case 'v': opt.verbose++; break;
            case 'h': usage(stdout); return 0;
            case 'V': printf("lzr %s (format LZR 1.0)\n", lzr_version_string()); return 0;
            case 'T':
            case 'S': {
                const char *v = p[1] ? p + 1 : (i + 1 < argc ? argv[++i] : NULL);
                if (!v) die_usage("option '-%s' needs a value", *p == 'T' ? "T" : "S");
                if (*p == 'T') set_threads(v);
                else if (!*v) die_usage("empty suffix%s", "");
                else opt.suffix = v;
                goto next_arg;   /* the value ends the cluster */
            }
            default: {
                char b[3] = {'-', *p, 0};
                die_usage("unknown option '%s'", b);
            }
            }
        }
    next_arg:;
    }

    signal(SIGINT, on_signal);
    signal(SIGTERM, on_signal);
    signal(SIGHUP, on_signal);
    signal(SIGPIPE, SIG_IGN);

    if (nfiles == 0) files[nfiles++] = "-";
    for (int i = 0; i < nfiles; i++) {
        if (!strcmp(files[i], "-")) {
            if (opt.mode == M_LIST) { err(NULL, "--list does not read standard input"); continue; }
            do_file(NULL, i == 0);
        } else {
            do_file(files[i], i == 0);
        }
    }
    free(files);
    if (fflush(stdout) || ferror(stdout)) { err(NULL, "write error on standard output"); }
    return status;
}
