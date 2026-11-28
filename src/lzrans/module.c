/*
 * module.c - CPython binding for the lzrans C core.
 *
 *   compress(data, level=6, window_log=0, threads=1, filter=None) -> bytes
 *   decompress(data) -> bytes
 *   Compressor(level=6, window_log=0, threads=1, filter=None): .compress(data), .flush()
 *   Decompressor(): .decompress(data) -> bytes, .eof, .unused_data
 *   crc32(data, value=0) -> int                     (zlib-compatible)
 *   frame_info(data) -> (content_size | None, window_log, flags)
 *
 * The GIL is released while the C core runs. A Compressor/Decompressor
 * object must not be used from two threads at once (a lock guards it).
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "lzr.h"
#include "crc32.h"

static PyObject *LzrError, *LzrChecksumError;

static PyObject *raise_code(int code)
{
    switch (code) {
    case LZR_ERR_NOMEM: return PyErr_NoMemory();
    case LZR_ERR_CHECKSUM:
        PyErr_SetString(LzrChecksumError, "checksum mismatch: decompressed data is corrupt");
        return NULL;
    case LZR_ERR_UNSUPPORTED:
        PyErr_SetString(LzrError, "unsupported frame version, flags or block type: needs a newer decoder (encrypted frames are not supported yet)");
        return NULL;
    case LZR_ERR_PARAM:
        PyErr_SetString(PyExc_ValueError, "invalid argument or stream already finished");
        return NULL;
    default:
        PyErr_SetString(LzrError, "corrupt or truncated lzrans frame");
        return NULL;
    }
}

/* Move a lzr_buf into a new bytes object and reset it. */
static PyObject *take_bytes(lzr_buf *b)
{
    PyObject *r = PyBytes_FromStringAndSize((const char *)b->p, (Py_ssize_t)b->len);
    b->len = 0;
    return r;
}

static int check_threads(int threads)
{
    if (threads < 0 || threads > 16) {
        PyErr_SetString(PyExc_ValueError, "threads must be between 0 (auto) and 16");
        return -1;
    }
    return 0;
}

/* 0 = auto: one per CPU, capped at 4. */
static int resolve_threads(int threads)
{
    if (threads) return threads;
    PyObject *os = PyImport_ImportModule("os");
    long n = 1;
    if (os) {
        PyObject *r = PyObject_CallMethod(os, "cpu_count", NULL);
        if (r && r != Py_None) n = PyLong_AsLong(r);
        Py_XDECREF(r);
        Py_DECREF(os);
    }
    PyErr_Clear();
    return n < 1 ? 1 : (n > 4 ? 4 : (int)n);
}

/* segment: -1 = the level's default (lzr_default_segment), 0 = off,
 * else the reset interval in bytes (>= 1 MiB). */
static int resolve_segment(long long seg, int level, uint64_t *out)
{
    if (seg < 0) { *out = lzr_default_segment(level); return 0; }
    if (seg > 0 && (uint64_t)seg < LZR_SEG_MIN) {
        PyErr_SetString(PyExc_ValueError, "segment must be -1 (level default), 0 (off) or >= 1 MiB");
        return -1;
    }
    *out = (uint64_t)seg;
    return 0;
}

static int check_level(int level, int wlog)
{
    if (level < LZR_LEVEL_MIN || level > LZR_LEVEL_MAX) {
        PyErr_SetString(PyExc_ValueError, "level must be between -3 and 11");
        return -1;
    }
    if (wlog && (wlog < LZR_WLOG_MIN || wlog > LZR_WLOG_MAX)) {
        PyErr_SetString(PyExc_ValueError, "window_log must be 0 (auto) or between 18 and 26");
        return -1;
    }
    return 0;
}

/* filter: None or "x86" -> LZR_FILTER_*, or -1 with an exception set. */
static int parse_filter(PyObject *f)
{
    if (!f || f == Py_None) return LZR_FILTER_NONE;
    if (PyUnicode_Check(f) && PyUnicode_CompareWithASCIIString(f, "x86") == 0) return LZR_FILTER_X86;
    PyErr_SetString(PyExc_ValueError, "filter must be None or \"x86\"");
    return -1;
}

/* ---------------- one-shot ---------------- */

static PyObject *py_compress(PyObject *self, PyObject *args, PyObject *kw)
{
    static char *kwlist[] = {"data", "level", "window_log", "threads", "filter", "segment", NULL};
    Py_buffer in;
    int level = 6, wlog = 0, threads = 1;
    long long segarg = -1;
    uint64_t seg = 0;
    PyObject *fobj = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kw, "y*|iiiOL:compress", kwlist, &in, &level, &wlog, &threads, &fobj, &segarg)) return NULL;
    int filter = parse_filter(fobj);
    if (filter < 0 || check_level(level, wlog) || check_threads(threads) || resolve_segment(segarg, level, &seg)) { PyBuffer_Release(&in); return NULL; }
    threads = resolve_threads(threads);

    lzr_buf out = {0};
    int r;
    Py_BEGIN_ALLOW_THREADS
    lzr_cctx *c = lzr_cctx_new_seg(level, wlog, (int64_t)in.len, seg);
    if (!c) r = LZR_ERR_NOMEM;
    else {
        lzr_cctx_set_threads(c, threads);
        r = lzr_cctx_set_filter(c, filter);
        if (!r) r = buf_reserve(&out, (size_t)in.len / 2 + 64);
        if (!r) r = lzr_cctx_write(c, (const uint8_t *)in.buf, (size_t)in.len, &out);
        if (!r) r = lzr_cctx_end(c, &out);
        lzr_cctx_free(c);
    }
    Py_END_ALLOW_THREADS
    PyBuffer_Release(&in);
    PyObject *res = r ? raise_code(r) : take_bytes(&out);
    free(out.p);
    return res;
}

static PyObject *py_decompress(PyObject *self, PyObject *args, PyObject *kw)
{
    static char *kwlist[] = {"data", "threads", NULL};
    Py_buffer in;
    int threads = 1;
    if (!PyArg_ParseTupleAndKeywords(args, kw, "y*|i:decompress", kwlist, &in, &threads)) return NULL;
    if (check_threads(threads)) { PyBuffer_Release(&in); return NULL; }
    threads = resolve_threads(threads);
    lzr_buf out = {0};
    int r;
    size_t used = 0;
    int64_t csize = -1;
    int hr = lzr_frame_header((const uint8_t *)in.buf, (size_t)in.len, NULL, NULL, &csize);
    /* When the header carries the content size (one-shot compress), decode
     * straight into the result bytes object: no growth, no final copy. The
     * decoder rejects frames that produce more than content_size bytes, so
     * the buffer (size + wildcopy slack) never needs to grow. */
    PyObject *bytes = NULL;
    const size_t extra = LZR_SLACK + 64;
    if (hr > 0 && csize > 0 && (uint64_t)csize <= (uint64_t)in.len * 4096u + (1u << 20) &&
        (uint64_t)csize <= (uint64_t)PY_SSIZE_T_MAX - extra) {
        bytes = PyBytes_FromStringAndSize(NULL, (Py_ssize_t)((size_t)csize + extra));
        if (!bytes) { PyBuffer_Release(&in); return NULL; }
        out.p = (uint8_t *)PyBytes_AS_STRING(bytes);
        out.cap = (size_t)csize + extra;
        out.fixed = 1;
    }
    Py_BEGIN_ALLOW_THREADS
    r = lzr_decompress_mt((const uint8_t *)in.buf, (size_t)in.len, &out, threads, &used);
    Py_END_ALLOW_THREADS
    Py_ssize_t inlen = in.len;
    PyBuffer_Release(&in);
    PyObject *res = NULL;
    if (r < 0) raise_code(r);
    else if (r != 1) PyErr_SetString(LzrError, "truncated lzrans frame (no end block)");
    else if ((Py_ssize_t)used != inlen) PyErr_SetString(LzrError, "trailing data after lzrans frame");
    else if (bytes) {
        if (_PyBytes_Resize(&bytes, (Py_ssize_t)out.len) == 0) { res = bytes; bytes = NULL; }
    } else res = take_bytes(&out);
    Py_XDECREF(bytes);
    if (!out.fixed) free(out.p);
    return res;
}

static PyObject *py_crc32(PyObject *self, PyObject *args)
{
    Py_buffer in;
    unsigned int value = 0;
    if (!PyArg_ParseTuple(args, "y*|I:crc32", &in, &value)) return NULL;
    uint32_t c;
    Py_BEGIN_ALLOW_THREADS
    c = lzr_crc32(value, (const uint8_t *)in.buf, (size_t)in.len);
    Py_END_ALLOW_THREADS
    PyBuffer_Release(&in);
    return PyLong_FromUnsignedLong(c);
}

static PyObject *py_frame_info(PyObject *self, PyObject *arg)
{
    Py_buffer in;
    if (PyObject_GetBuffer(arg, &in, PyBUF_SIMPLE) < 0) return NULL;
    int wlog = 0;
    uint8_t flags = 0;
    int64_t size = -1;
    int r = lzr_frame_header((const uint8_t *)in.buf, (size_t)in.len, &wlog, &flags, &size);
    PyBuffer_Release(&in);
    if (r == 0) r = LZR_ERR_CORRUPT;
    if (r < 0) return raise_code(r);
    if (size < 0) return Py_BuildValue("(OiB)", Py_None, wlog, flags);
    return Py_BuildValue("(LiB)", (long long)size, wlog, flags);
}

/* ---------------- Compressor ---------------- */

typedef struct {
    PyObject_HEAD
    lzr_cctx *c;
    lzr_buf out;
    PyThread_type_lock lock;
} CompressorObject;

static int Compressor_init(CompressorObject *self, PyObject *args, PyObject *kw)
{
    static char *kwlist[] = {"level", "window_log", "threads", "filter", "segment", NULL};
    int level = 6, wlog = 0, threads = 1;
    long long segarg = -1;
    uint64_t seg = 0;
    PyObject *fobj = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kw, "|iiiOL:Compressor", kwlist, &level, &wlog, &threads, &fobj, &segarg)) return -1;
    int filter = parse_filter(fobj);
    if (filter < 0 || check_level(level, wlog) || check_threads(threads) || resolve_segment(segarg, level, &seg)) return -1;
    threads = resolve_threads(threads);
    if (self->c) { lzr_cctx_free(self->c); self->c = NULL; }
    self->c = lzr_cctx_new_seg(level, wlog, -1, seg);
    if (!self->c) { PyErr_NoMemory(); return -1; }
    lzr_cctx_set_threads(self->c, threads);
    if (lzr_cctx_set_filter(self->c, filter)) { PyErr_NoMemory(); return -1; }
    if (!self->lock && !(self->lock = PyThread_allocate_lock())) { PyErr_NoMemory(); return -1; }
    return 0;
}

static void Compressor_dealloc(CompressorObject *self)
{
    lzr_cctx_free(self->c);
    free(self->out.p);
    if (self->lock) PyThread_free_lock(self->lock);
    Py_TYPE(self)->tp_free((PyObject *)self);
}

#define LOCK(o)   do { if (!PyThread_acquire_lock((o)->lock, 0)) { Py_BEGIN_ALLOW_THREADS PyThread_acquire_lock((o)->lock, 1); Py_END_ALLOW_THREADS } } while (0)
#define UNLOCK(o) PyThread_release_lock((o)->lock)

static PyObject *Compressor_compress(CompressorObject *self, PyObject *arg)
{
    if (!self->c) { PyErr_SetString(PyExc_ValueError, "Compressor not initialised"); return NULL; }
    Py_buffer in;
    if (PyObject_GetBuffer(arg, &in, PyBUF_SIMPLE) < 0) return NULL;
    LOCK(self);
    int r;
    Py_BEGIN_ALLOW_THREADS
    r = lzr_cctx_write(self->c, (const uint8_t *)in.buf, (size_t)in.len, &self->out);
    Py_END_ALLOW_THREADS
    PyBuffer_Release(&in);
    PyObject *res = r ? raise_code(r) : take_bytes(&self->out);
    UNLOCK(self);
    return res;
}

static PyObject *Compressor_flush(CompressorObject *self, PyObject *noargs)
{
    if (!self->c) { PyErr_SetString(PyExc_ValueError, "Compressor not initialised"); return NULL; }
    LOCK(self);
    int r;
    Py_BEGIN_ALLOW_THREADS
    r = lzr_cctx_end(self->c, &self->out);
    Py_END_ALLOW_THREADS
    PyObject *res = r ? raise_code(r) : take_bytes(&self->out);
    UNLOCK(self);
    return res;
}

static PyMethodDef Compressor_methods[] = {
    {"compress", (PyCFunction)Compressor_compress, METH_O,
     "compress(data) -> bytes. Feed data; returns any completed compressed output."},
    {"flush", (PyCFunction)Compressor_flush, METH_NOARGS,
     "flush() -> bytes. Finish the frame. The object cannot be used afterwards."},
    {NULL}};

static PyTypeObject CompressorType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "lzrans.Compressor",
    .tp_basicsize = sizeof(CompressorObject),
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_doc = "Compressor(level=6, window_log=0, threads=1, filter=None): streaming compressor.",
    .tp_new = PyType_GenericNew,
    .tp_init = (initproc)Compressor_init,
    .tp_dealloc = (destructor)Compressor_dealloc,
    .tp_methods = Compressor_methods,
};

/* ---------------- Decompressor ---------------- */

typedef struct {
    PyObject_HEAD
    lzr_dctx *d;
    lzr_buf out;
    int eof;
    PyObject *unused;
    PyThread_type_lock lock;
} DecompressorObject;

static int Decompressor_init(DecompressorObject *self, PyObject *args, PyObject *kw)
{
    static char *kwlist[] = {NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kw, ":Decompressor", kwlist)) return -1;
    if (self->d) lzr_dctx_free(self->d);
    self->d = lzr_dctx_new();
    if (!self->d) { PyErr_NoMemory(); return -1; }
    self->eof = 0;
    Py_XSETREF(self->unused, PyBytes_FromStringAndSize(NULL, 0));
    if (!self->lock && !(self->lock = PyThread_allocate_lock())) { PyErr_NoMemory(); return -1; }
    return self->unused ? 0 : -1;
}

static void Decompressor_dealloc(DecompressorObject *self)
{
    lzr_dctx_free(self->d);
    free(self->out.p);
    Py_XDECREF(self->unused);
    if (self->lock) PyThread_free_lock(self->lock);
    Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject *Decompressor_decompress(DecompressorObject *self, PyObject *arg)
{
    if (!self->d) { PyErr_SetString(PyExc_ValueError, "Decompressor not initialised"); return NULL; }
    if (self->eof) { PyErr_SetString(PyExc_EOFError, "end of lzrans frame already reached"); return NULL; }
    Py_buffer in;
    if (PyObject_GetBuffer(arg, &in, PyBUF_SIMPLE) < 0) return NULL;
    LOCK(self);
    int r;
    size_t used = 0;
    Py_BEGIN_ALLOW_THREADS
    r = lzr_dctx_write(self->d, (const uint8_t *)in.buf, (size_t)in.len, &self->out, &used);
    Py_END_ALLOW_THREADS
    PyObject *res = NULL;
    if (r < 0) {
        self->out.len = 0;
        raise_code(r);
    } else {
        if (r == 1) {
            self->eof = 1;
            Py_XSETREF(self->unused, PyBytes_FromStringAndSize((const char *)in.buf + used,
                                                               in.len - (Py_ssize_t)used));
        }
        res = self->unused ? take_bytes(&self->out) : NULL;
    }
    PyBuffer_Release(&in);
    UNLOCK(self);
    return res;
}

static PyObject *Decompressor_get_eof(DecompressorObject *self, void *c) { return PyBool_FromLong(self->eof); }
static PyObject *Decompressor_get_unused(DecompressorObject *self, void *c)
{
    return Py_NewRef(self->unused ? self->unused : Py_None);
}

static PyMethodDef Decompressor_methods[] = {
    {"decompress", (PyCFunction)Decompressor_decompress, METH_O,
     "decompress(data) -> bytes. Feed compressed data; returns decoded output so far."},
    {NULL}};

static PyGetSetDef Decompressor_getset[] = {
    {"eof", (getter)Decompressor_get_eof, NULL, "True once the end block was verified.", NULL},
    {"unused_data", (getter)Decompressor_get_unused, NULL, "Bytes found after the end of the frame.", NULL},
    {NULL}};

static PyTypeObject DecompressorType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "lzrans.Decompressor",
    .tp_basicsize = sizeof(DecompressorObject),
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_doc = "Decompressor(): streaming decompressor.",
    .tp_new = PyType_GenericNew,
    .tp_init = (initproc)Decompressor_init,
    .tp_dealloc = (destructor)Decompressor_dealloc,
    .tp_methods = Decompressor_methods,
    .tp_getset = Decompressor_getset,
};

/* ---------------- module ---------------- */

static PyMethodDef methods[] = {
    {"compress", (PyCFunction)(void (*)(void))py_compress, METH_VARARGS | METH_KEYWORDS,
     "compress(data, level=6, window_log=0, threads=1, filter=None, segment=-1) -> bytes\n"
     "threads: worker threads (0 = auto, up to 4). Levels 7-11 and segment mode: output is identical for any value.\n"
     "segment: reset interval in bytes for parallel decompression (-1 = level default, 0 = none, >= 1 MiB).\n"
     "filter: None, or \"x86\" for x86/x86-64 machine code (executables, libraries)."},
    {"decompress", (PyCFunction)(void (*)(void))py_decompress, METH_VARARGS | METH_KEYWORDS,
     "decompress(data, threads=1) -> bytes\n"
     "threads: segments decoded in parallel (0 = auto, up to 4). Frames written with segments\n"
     "(compress(..., segment=...)) have reset points; others decode on one thread."},
    {"crc32", py_crc32, METH_VARARGS, "crc32(data, value=0) -> int (zlib-compatible)"},
    {"frame_info", py_frame_info, METH_O, "frame_info(data) -> (content_size or None, window_log, flags)"},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef moddef = {PyModuleDef_HEAD_INIT, "_lzrans",
                                    "LZ77 + rANS streaming compressor (C core).", -1, methods};

PyMODINIT_FUNC PyInit__lzrans(void)
{
    lzr_crc32_init();
    if (PyType_Ready(&CompressorType) < 0 || PyType_Ready(&DecompressorType) < 0) return NULL;
    PyObject *m = PyModule_Create(&moddef);
    if (!m) return NULL;
    LzrError = PyErr_NewException("lzrans.LzransError", PyExc_ValueError, NULL);
    LzrChecksumError = PyErr_NewException("lzrans.ChecksumError", LzrError, NULL);
    if (PyModule_AddObjectRef(m, "LzransError", LzrError) < 0 ||
        PyModule_AddObjectRef(m, "ChecksumError", LzrChecksumError) < 0 ||
        PyModule_AddObjectRef(m, "Compressor", (PyObject *)&CompressorType) < 0 ||
        PyModule_AddObjectRef(m, "Decompressor", (PyObject *)&DecompressorType) < 0 ||
        PyModule_AddIntConstant(m, "FLAG_ENCRYPTED", LZR_FLAG_ENCRYPTED) < 0 ||
        PyModule_AddIntConstant(m, "FLAG_SIZE", LZR_FLAG_SIZE) < 0 ||
        PyModule_AddIntConstant(m, "FLAG_X86", LZR_FLAG_X86) < 0 ||
        PyModule_AddIntConstant(m, "BLOCK_SIZE", LZR_BLOCK_SIZE) < 0 ||
        PyModule_AddIntConstant(m, "FORMAT_VERSION", LZR_VERSION) < 0) {
        Py_DECREF(m);
        return NULL;
    }
    return m;
}
