"""Build script for the lzrans C extension (setuptools ext_modules).

The extension compiles the same core sources as liblzr (Makefile / CMakeLists.txt);
module.c binds the core API (src/lzrans/core.h) directly, which lets one-shot
decompress write straight into the result bytes object."""
import glob
import sys
from setuptools import Extension, setup

extra = ["/O2"] if sys.platform == "win32" else ["-O3", "-std=gnu11", "-Wall", "-Wextra", "-Wno-unused-parameter", "-Wno-missing-field-initializers"]

setup(
    ext_modules=[
        Extension(
            "lzrans._lzrans",
            sources=[
                "src/lzrans/module.c",
                "src/lzrans/enc.c",
                "src/lzrans/dec.c",
                "src/lzrans/opt.c",
                "src/lzrans/crc32.c",
                "src/lzrans/huf.c",
                "src/lzrans/fse.c",
            ],
            depends=sorted(glob.glob("src/lzrans/*.h")) + ["include/lzr.h"],   # rebuild when a header changes
            include_dirs=["src/lzrans"],
            extra_compile_args=extra,
            libraries=[] if sys.platform == "win32" else ["m", "pthread"],
        )
    ],
)
