#!/usr/bin/env python
"""Fetch the SIFTS PDB→UniProt residue mapping from the EBI.

    python code/figures/validation/download_sifts.py --out DIR

`pdb_chain_uniprot.tsv` is what decides, for every solved structure, which residue range of
the protein a chain actually covers. `figures/notebooks/S2_methods_vs_structures.ipynb` uses it to
pick the chains that span an ectodomain, and `measure_pdb_chains.py` measures exactly the
mapped range -- which is how the archived measurements were reproduced to 0.0003 Å once the
convention was recovered.

The file is a weekly rolling release: its first line is a vintage stamp, and this prints it,
because that stamp is the only version this database has.

    # 2021/11/14 - 17:19 | PDB: 45.21 | UniProt: 2021.04     the archived copy
    # 2026/09/08 - 13:56 | PDB: 36.26 | UniProt: 2026.04     current

**A refresh is not a drop-in.** The measurements in `structure_sizes_211220.txt` were taken
over the ranges the *old* mapping gave. Re-run `measure_pdb_chains.py` against the refreshed
chain list, or S2 will merge new ranges onto old measurements and say nothing about it.
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import sys
from pathlib import Path

import requests

URL = ("https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/"
       "pdb_chain_uniprot.tsv.gz")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, required=True, help="directory to write into")
    ap.add_argument("--url", default=URL)
    ap.add_argument("--keep-gz", action="store_true",
                    help="also keep the downloaded .gz beside the .tsv")
    a = ap.parse_args()

    a.out.mkdir(parents=True, exist_ok=True)
    packed = a.out / "pdb_chain_uniprot.tsv.gz"
    dest = a.out / "pdb_chain_uniprot.tsv"

    with requests.get(a.url, stream=True, timeout=300) as r:
        r.raise_for_status()
        with open(packed, "wb") as fh:
            shutil.copyfileobj(r.raw, fh)

    with gzip.open(packed, "rb") as src, open(dest, "wb") as out:
        shutil.copyfileobj(src, out, 1 << 20)
    if not a.keep_gz:
        packed.unlink()

    with open(dest) as fh:
        stamp = fh.readline().strip()
        rows = sum(1 for _ in fh) - 1        # minus the column header
    print(f"  {stamp}")
    print(f"  {rows:,} chain mappings -> {dest}")
    print("\n  Re-measure before using this: the chain ranges have moved, and")
    print("  structure_sizes_211220.txt was measured over the old ones.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
