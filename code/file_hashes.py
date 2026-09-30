"""sha256 helpers shared by the snapshot builder and the archive tooling.

Kept apart from ``scripts/archive/build_archive_subsets.py`` so that
``build_snapshot.py`` -- part of running the pipeline -- does not import the
archive tooling, which a release tree may not carry.
"""

from __future__ import annotations

import gzip
import hashlib
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def open_text(path: Path):
    """Read a source file, transparently un-gzipping the committed copies."""
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def content_sha256(path: Path) -> str:
    """sha256 of the *uncompressed* content, so the archive file and its gzipped copy agree."""
    if path.suffix != ".gz":
        return sha256(path)
    h = hashlib.sha256()
    with gzip.open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()
