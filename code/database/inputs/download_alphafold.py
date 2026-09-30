#!/usr/bin/env python
"""Download AlphaFold DB models for a list of UniProt accessions.

    python code/database/inputs/download_alphafold.py --accessions FILE --out-dir DIR [--version 6]

Skips files already present, so it is safe to re-run after an interruption. Writes
`missing.txt` listing accessions AlphaFold has no model for.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import requests

BASE = "https://alphafold.ebi.ac.uk/files"


def fetch(session, acc: str, out_dir: Path, version: int, retries: int = 3) -> str:
    dest = out_dir / f"AF-{acc}-F1-model_v{version}.pdb"
    if dest.exists() and dest.stat().st_size > 0:
        return "cached"
    url = f"{BASE}/AF-{acc}-F1-model_v{version}.pdb"
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
    ap.add_argument("--accessions", required=True, type=Path,
                    help="one UniProt accession per line")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--version", type=int, default=6)
    a = ap.parse_args()

    accs = [l.strip() for l in a.accessions.read_text().splitlines() if l.strip()]
    a.out_dir.mkdir(parents=True, exist_ok=True)
    s = requests.Session()
    tally = {"ok": 0, "cached": 0, "missing": 0, "failed": 0}
    missing = []

    for i, acc in enumerate(accs, 1):
        res = fetch(s, acc, a.out_dir, a.version)
        tally[res] += 1
        if res in ("missing", "failed"):
            missing.append(f"{acc}\t{res}")
        if i % 100 == 0 or i == len(accs):
            print(f"  {i}/{len(accs)}  " + "  ".join(f"{k}={v}" for k, v in tally.items()), flush=True)

    (a.out_dir / "missing.txt").write_text("\n".join(missing) + ("\n" if missing else ""))
    print("\n" + "  ".join(f"{k}={v}" for k, v in tally.items()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
