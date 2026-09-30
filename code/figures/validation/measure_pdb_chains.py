#!/usr/bin/env python
"""Measure solved-structure chains with the inertia-axis-aligned bounding box.

    python code/figures/validation/measure_pdb_chains.py --structures DIR --chains CSV --out TSV

Regenerates `data/measurements/pdb_structure_sizes.txt`, the solved-structure measurements
`figures/notebooks/S2_methods_vs_structures.ipynb` compares the predicted heights with. It is a
measurement output, not source data: it is recomputed from the mmCIF files of the chains SIFTS
lists (`download_pdb_structures.py`), the chain list the S2 notebook writes
(`database/surfaceome_pdbs.csv`), and the bounding-box routine in `measure_structures.py`,
validated to 0.05 A against the published AlphaFold measurements.

Output matches what the notebook reads -- headerless, tab separated:

    CHAIN  PDB  PDB_len  x  y  z

`PDB_len` is a residue count and `x >= y >= z` are the box dimensions in Angstrom.

What is measured had to be recovered from the archived output, the same way the AlphaFold
trim rule was. Three candidates were tested against
`unique_extracellular_pdb_structures_and_measurements.csv`, comparing `max_dim`:

    whole chain                     3/20 within 0.05 A   median |diff| 3.79 A
    CA atoms only                   0/20                 median |diff| 2.83 A
    the SIFTS-mapped range          19/20 within 0.05 A  median |diff| 0.000 A

So the measured selection is `polymer and chain <C> and resi <PDB_BEG>-<PDB_END>` -- the
mapped range, not the chain as deposited. `polymer` keeps waters, ions and ligands out of
the box.

**`PDB_len` counts what was measured, in state 1.** Counting without a state returns atoms
the geometry pass never sees -- 453 against 84 for 7bg7 chain B -- and since `S2` sorts on
`PDB_len` to choose each protein's representative structure, an inflated count wins the sort
and the chosen structure then measures small. That is how ICAM1 came out at 45 A against the
190 A of the structure the archived table used.

**One caveat, stated rather than papered over.** `PDB_len` does *not* reproduce exactly.
Every residue-count definition tried -- CA atoms, unique `resi`, PyMOL `guide` atoms, the
range span -- agrees with the others and gives 11 of 25 exact against the archived column,
median difference 1 residue, 7 of 25 off by more than 5. The archived measurements were
evidently taken against a different generation of the PDB->UniProt mapping than the
`human_surfaceome_pdbs.csv` archived beside them. The dimensions are unaffected; `PDB_len`
feeds only a filter discarding structures longer than the ectodomain, so a regenerated S2
will keep a slightly different set of structures from the submitted one. Worth knowing
before comparing a regenerated Figure S2 against the published one.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # code/, for surfaceome_config and the package
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "database" / "inputs"))  # measure_structures lives with the input build
from measure_structures import _setup, iabb  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--structures", required=True, type=Path, action="append",
                    help="directory of mmCIF/PDB files named <pdbid>.cif; repeatable, "
                         "searched in order -- the archived 2021 pull and a newer download "
                         "of what a refreshed SIFTS adds are two directories")
    ap.add_argument("--chains", required=True, type=Path,
                    help="CSV with PDB, CHAIN, PDB_BEG, PDB_END (human_surfaceome_pdbs.csv)")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()

    want = (pd.read_csv(a.chains)[["PDB", "CHAIN", "PDB_BEG", "PDB_END"]]
            .dropna().drop_duplicates(subset=["PDB", "CHAIN"]))
    want["PDB"] = want.PDB.astype(str)
    want["CHAIN"] = want.CHAIN.astype(str)
    if a.limit:
        want = want.head(a.limit)
    print(f"{len(want):,} chains over {want.PDB.nunique():,} entries")

    cmd = _setup()
    rows, missing, empty = [], set(), 0
    for n, (pdb, grp) in enumerate(want.groupby("PDB"), 1):
        path = next((d / f"{pdb.lower()}.cif" for d in a.structures
                     if (d / f"{pdb.lower()}.cif").exists()), None)
        if path is None:
            missing.add(pdb)
            continue
        cmd.delete("all")
        try:
            cmd.load(str(path), "s")
        except Exception:
            missing.add(pdb)
            continue
        for row in grp.itertuples():
            chain = row.CHAIN
            sel = (f"s and polymer and chain {chain} "
                   f"and resi {int(row.PDB_BEG)}-{int(row.PDB_END)}")
            try:
                # state=1 matters: without it PyMOL counts atoms the geometry pass does
                # not see, so PDB_len stops describing what was actually boxed.  7bg7 chain
                # B counted 453 residues against the 84 whose coordinates were measured, and
                # PDB_len is what S2 sorts on to pick a protein's representative structure --
                # so an inflated count wins the sort and then measures small.
                n_res = cmd.count_atoms(f"{sel} and name CA and alt ''+A", state=1)
                if n_res < 3:
                    empty += 1
                    continue
                dims = iabb(cmd, sel)
            except Exception:
                empty += 1
                cmd.delete("all")
                cmd.load(str(path), "s")
                continue
            if dims is None:
                empty += 1
            else:
                rows.append((chain, pdb, n_res, *dims))
            #iabb() transforms the object in place, so reload before the next chain
            cmd.delete("all")
            cmd.load(str(path), "s")
        if n % 200 == 0:
            print(f"  {n:,}/{want.PDB.nunique():,} entries, {len(rows):,} chains measured",
                  flush=True)

    out = pd.DataFrame(rows, columns=["CHAIN", "PDB", "PDB_len", "x", "y", "z"])
    a.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, sep="\t", header=False, index=False)
    print(f"\nmeasured {len(out):,} chains")
    print(f"  entries not found:      {len(missing):,}")
    print(f"  chains skipped (empty): {empty:,}")
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
