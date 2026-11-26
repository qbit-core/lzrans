/*
 * module.c - CPython binding for the lzrans C core.
 *
 *   compress(data, level=5) -> bytes
 *   decompress(data) -> bytes
 *   crc32(data, value=0) -> int        (zlib-compatible)
 *   frame_info(data) -> (original_size, crc32, flags)
 *
 * Accepts any C-contiguous buffer (bytes, bytearray, memoryview, ...).
 * The GIL is released while the C core runs, so threads can compress
 * in parallel.
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "frame.h"
#include "crc32.h"

static PyObject *LzrError;   /* lzrans.LzransError (subclass of ValueError) */
static PyObject *LzrChecksumError;

static PyObject *raise_code(int code)
{
    switch (code) {
    case LZR_ERR_NOMEM: return PyErr_NoMemory();
    case LZR_ERR_CHECKSUM:
        PyErr_SetString(LzrChecksumError, "CRC32 mismatch: decompressed data is corrupt");
        return NULL;
    case LZR_ERR_UNSUPPORTED:
        PyErr_SetString(LzrError, "unsupported frame version or flags (encrypted frames are not supported yet)");
        return NULL;
    default:
        PyErr_SetString(LzrError, "corrupt or truncated lzrans frame");
        return NULL;
    }
}

static PyObject *py_compress(PyObject *self, PyObject *args, PyObject *kw)
{
    static char *kwlist[] = {"data", "level", NULL};
    Py_buffer in;
    int level = 5;
    if (!PyArg_ParseTupleAndKeywords(args, kw, "y*|i:compress", kwlist, &in, &level)) return NULL;
    if (level < 1 || level > 9) {
        PyBuffer_Release(&in);
        PyErr_SetString(PyExc_ValueError, "level must be between 1 and 9");
        return NULL;
    }
    size_t cap = lzr_compress_bound((size_t)in.len);
    PyObject *out = PyBytes_FromStringAndSize(NULL, (Py_ssize_t)cap);
    if (!out) { PyBuffer_Release(&in); return NULL; }

    int64_t r;
    Py_BEGIN_ALLOW_THREADS
    r = lzr_compress((const uint8_t *)in.buf, (size_t)in.len,
                     (uint8_t *)PyBytes_AS_STRING(out), cap, level);
    Py_END_ALLOW_THREADS
    PyBuffer_Release(&in);
    if (r < 0) { Py_DECREF(out); return raise_code((int)r); }
    if (_PyBytes_Resize(&out, (Py_ssize_t)r) < 0) return NULL;
    return out;
}

static PyObject *py_decompress(PyObject *self, PyObject *arg)
{
    Py_buffer in;
    if (PyObject_GetBuffer(arg, &in, PyBUF_SIMPLE) < 0) return NULL;
    uint64_t orig;
    int r = lzr_frame_info((const uint8_t *)in.buf, (size_t)in.len, &orig, NULL, NULL);
    if (r) { PyBuffer_Release(&in); return raise_code(r); }
    if (orig > (uint64_t)PY_SSIZE_T_MAX) { PyBuffer_Release(&in); return PyErr_NoMemory(); }

    PyObject *out = PyBytes_FromStringAndSize(NULL, (Py_ssize_t)orig);
    if (!out) { PyBuffer_Release(&in); return NULL; }
    Py_BEGIN_ALLOW_THREADS
    r = lzr_decompress((const uint8_t *)in.buf, (size_t)in.len,
                       (uint8_t *)PyBytes_AS_STRING(out), (size_t)orig);
    Py_END_ALLOW_THREADS
    PyBuffer_Release(&in);
    if (r) { Py_DECREF(out); return raise_code(r); }
    return out;
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
    uint64_t orig = 0;
    uint32_t crc = 0;
    uint8_t flags = 0;
    int r = lzr_frame_info((const uint8_t *)in.buf, (size_t)in.len, &orig, &crc, &flags);
    PyBuffer_Release(&in);
    if (r) return raise_code(r);
    return Py_BuildValue("(KkB)", (unsigned long long)orig, (unsigned long)crc, flags);
}

static PyMethodDef methods[] = {
    {"compress", (PyCFunction)(void (*)(void))py_compress, METH_VARARGS | METH_KEYWORDS,
     "compress(data, level=5) -> bytes"},
    {"decompress", py_decompress, METH_O, "decompress(data) -> bytes"},
    {"crc32", py_crc32, METH_VARARGS, "crc32(data, value=0) -> int (zlib-compatible)"},
    {"frame_info", py_frame_info, METH_O, "frame_info(data) -> (original_size, crc32, flags)"},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef moddef = {PyModuleDef_HEAD_INIT, "_lzrans",
                                    "LZ77 + 3-stream rANS compressor (C core).", -1, methods};

PyMODINIT_FUNC PyInit__lzrans(void)
{
    lzr_crc32_init();
    PyObject *m = PyModule_Create(&moddef);
    if (!m) return NULL;
    LzrError = PyErr_NewException("lzrans.LzransError", PyExc_ValueError, NULL);
    LzrChecksumError = PyErr_NewException("lzrans.ChecksumError", LzrError, NULL);
    if (PyModule_AddObjectRef(m, "LzransError", LzrError) < 0 ||
        PyModule_AddObjectRef(m, "ChecksumError", LzrChecksumError) < 0 ||
        PyModule_AddIntConstant(m, "FLAG_ENCRYPTED", LZR_FLAG_ENCRYPTED) < 0 ||
        PyModule_AddIntConstant(m, "BLOCK_SIZE", LZR_BLOCK_SIZE) < 0 ||
        PyModule_AddIntConstant(m, "FORMAT_VERSION", LZR_VERSION) < 0) {
        Py_DECREF(m);
        return NULL;
    }
    return m;
}
