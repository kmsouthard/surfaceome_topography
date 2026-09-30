#!/usr/bin/env python
"""Download STRING protein links and aliases for one organism.

    python code/database/inputs/download_string.py --out-dir DIR [--version 12.0] [--taxon 9606]

The interaction set the pipeline reads is derived from STRING **v11.0**, pulled in 2020. This
fetches a chosen release so the derivation can be repeated against a newer one; v12.0 is
current.

Three files, all gzipped as served:

    <taxon>.protein.links.detailed.<v>.txt.gz   the channel scores and combined_score
    <taxon>.protein.aliases.<v>.txt.gz          STRING id -> UniProt accession, among others
    <taxon>.protein.info.<v>.txt.gz             preferred names

The links file is 140 MB for human and is **not** vendored: it is a public, versioned download,
and what the pipeline reads is the derived pair table. Existing files are skipped.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import requests

BASE = "https://stringdb-downloads.org/download"
FILES = ("protein.links.detailed", "protein.aliases", "protein.info")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--version", default="12.0")
    ap.add_argument("--taxon", default="9606")
    a = ap.parse_args()

    a.out_dir.mkdir(parents=True, exist_ok=True)
    for kind in FILES:
        name = f"{a.taxon}.{kind}.v{a.version}.txt.gz"
        dest = a.out_dir / name
        if dest.exists() and dest.stat().st_size > 0:
            print(f"  cached  {name}")
            continue
        url = f"{BASE}/{kind}.v{a.version}/{name}"
        with requests.get(url, stream=True, timeout=900) as r:
            if r.status_code != 200:
                print(f"  FAILED  {url} -> HTTP {r.status_code}", file=sys.stderr)
                return 1
            with open(dest, "wb") as fh:
                shutil.copyfileobj(r.raw, fh)
        print(f"  {dest.stat().st_size / 1e6:>7.1f} MB  {name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
