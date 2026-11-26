"""Build script for the lzrans C extension (setuptools ext_modules)."""
import sys
from setuptools import Extension, setup

extra = ["/O2"] if sys.platform == "win32" else ["-O3", "-std=c11", "-Wall", "-Wextra", "-Wno-unused-parameter", "-Wno-missing-field-initializers"]

setup(
    ext_modules=[
        Extension(
            "lzrans._lzrans",
            sources=[
                "src/lzrans/module.c",
                "src/lzrans/enc.c",
                "src/lzrans/dec.c",
                "src/lzrans/rans.c",
                "src/lzrans/crc32.c",
            ],
            include_dirs=["src/lzrans"],
            extra_compile_args=extra,
            libraries=[] if sys.platform == "win32" else ["m"],
        )
    ],
)
