#!/usr/bin/env python
"""Download the mmCIF files a chain list needs, from RCSB.

    python code/figures/validation/download_pdb_structures.py --chains CSV --out-dir DIR [--have DIR]

`--chains` is the table `figures/notebooks/S2_methods_vs_structures.ipynb` writes at cell 17
(`data/human_surfaceome_pdbs.csv`): the solved-structure chains that span an ectodomain.
`--have` names a directory of structures already on disk -- the archive's `code/tidy_code/pdb/`
holds the set downloaded in 2021 -- so only what is missing is fetched.

Files are named `<pdbid>.cif`, lowercase, which is what `measure_pdb_chains.py` expects.
Existing files are skipped, so an interrupted run can simply be repeated.

These are not vendored. At 4.9 GB for the archived set alone they are far too large, and
unlike GlycoMine they are a public, versioned, scripted download -- the same standing as the
AlphaFold models. What the pipeline reads is the *measurement* output, which is vendored.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import pandas as pd
import requests

BASE = "https://files.rcsb.org/download"


def fetch(session, pdb: str, out_dir: Path, retries: int = 3) -> str:
    dest = out_dir / f"{pdb.lower()}.cif"
    if dest.exists() and dest.stat().st_size > 0:
        return "cached"
    url = f"{BASE}/{pdb.lower()}.cif"
    for attempt in range(retries):
        try:
            r = session.get(url, timeout=120)
            if r.status_code == 200:
                dest.write_bytes(r.content)
                return "ok"
            if r.status_code == 404:
                return "missing"
        except requests.RequestException:
            pass
        time.sleep(2 ** attempt)
    return "failed"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--chains", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--have", type=Path, action="append", default=[],
                    help="directory of structures already on disk; repeatable")
    a = ap.parse_args()

    a.out_dir.mkdir(parents=True, exist_ok=True)
    want = sorted({str(p).lower() for p in pd.read_csv(a.chains)["PDB"].dropna().unique()})
    have = {p.stem.lower() for d in list(a.have) + [a.out_dir] if d.is_dir()
            for p in d.glob("*.cif")}
    todo = [p for p in want if p not in have]
    print(f"  {len(want):,} structures wanted, {len(want) - len(todo):,} on disk, "
          f"{len(todo):,} to fetch")

    counts: dict[str, int] = {}
    missing = []
    with requests.Session() as s:
        for i, pdb in enumerate(todo, 1):
            state = fetch(s, pdb, a.out_dir)
            counts[state] = counts.get(state, 0) + 1
            if state in ("missing", "failed"):
                missing.append(f"{pdb}\t{state}")
            if i % 200 == 0 or i == len(todo):
                print(f"    {i:,}/{len(todo):,}  " +
                      "  ".join(f"{k}={v}" for k, v in sorted(counts.items())), flush=True)

    if missing:
        (a.out_dir / "missing.txt").write_text("\n".join(missing) + "\n")
        print(f"  {len(missing)} not retrieved; listed in {a.out_dir / 'missing.txt'}")
    return 1 if counts.get("failed") else 0


if __name__ == "__main__":
    sys.exit(main())
