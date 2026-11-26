"""lzrans - LZ77 (64 KB window) + 3-stream rANS compression with a C core.

Public API:
    compress(data, level=5) -> bytes
    decompress(data) -> bytes
    crc32(data, value=0) -> int
    frame_info(data) -> (original_size, crc32, flags)

Encryption hook: frames carry a flags byte. FLAG_ENCRYPTED is reserved for a
future AEAD envelope (planned: libsodium XChaCha20-Poly1305 + Argon2id).
This build rejects encrypted frames with LzransError.
"""
from ._lzrans import (
    BLOCK_SIZE,
    FLAG_ENCRYPTED,
    FORMAT_VERSION,
    ChecksumError,
    LzransError,
    compress,
    crc32,
    decompress,
    frame_info,
)

__all__ = [
    "compress", "decompress", "crc32", "frame_info",
    "LzransError", "ChecksumError",
    "BLOCK_SIZE", "FLAG_ENCRYPTED", "FORMAT_VERSION",
]
__version__ = "0.1.0"
