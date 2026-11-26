"""lzrans - LZ77 (64 KB window) + 3-stream rANS compression with a C core.

Public API:
    compress(data, level=6, window_log=0) -> bytes
    decompress(data) -> bytes
    Compressor(level=6, window_log=0)   .compress(chunk) / .flush()
    Decompressor()                      .decompress(chunk), .eof, .unused_data
    crc32(data, value=0) -> int
    frame_info(data) -> (content_size or None, window_log, flags)

Encryption hook: frames carry a flags byte. FLAG_ENCRYPTED is reserved for a
future AEAD envelope (planned: libsodium XChaCha20-Poly1305 + Argon2id).
This build rejects encrypted frames with LzransError.
"""
from ._lzrans import (
    BLOCK_SIZE,
    FLAG_ENCRYPTED,
    FLAG_SIZE,
    FLAG_X86,
    Compressor,
    Decompressor,
    FORMAT_VERSION,
    ChecksumError,
    LzransError,
    compress,
    crc32,
    decompress,
    frame_info,
)

__all__ = [
    "compress", "decompress", "Compressor", "Decompressor", "crc32", "frame_info",
    "LzransError", "ChecksumError",
    "BLOCK_SIZE", "FLAG_ENCRYPTED", "FLAG_SIZE", "FLAG_X86", "FORMAT_VERSION",
]
__version__ = "0.3.0.dev0"
