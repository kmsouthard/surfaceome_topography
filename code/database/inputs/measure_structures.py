#!/usr/bin/env python
"""Measure ectodomain dimensions with the inertia-axis-aligned bounding box, in PyMOL.

    python code/database/inputs/measure_structures.py --models DIR --regions FILE --out FILE [--pattern ...]

This is `Draw_Protein_Dimensions.py` (Pablo Guardado Calvo) reimplemented for Python 3 and
driven in batch. The original is a Python 2 PyMOL plugin that writes distances to a file as a
side effect; the maths here is identical, calling the same PyMOL API:

  1. translate the selection's centre of mass to the origin
  2. build the mass-weighted inertia tensor about the origin
  3. take its eigenvectors, sorted by ascending eigenvalue, rounded to 3 decimals, as rows
     (rows reversed if the determinant is negative, to keep a proper rotation)
  4. apply that rotation and take the coordinate extents

`--regions` is a CSV with columns `ID`, `ecd_start`, `ecd_end`; each row is measured over that
residue range. Output columns are `ID`, `n_atoms`, `dim_1` >= `dim_2` >= `dim_3` (Angstroms),
`mean_plddt` (the B-factor column, which AlphaFold uses for per-residue confidence).

Which residues to pass
----------------------
To reproduce the published AlphaFold measurements, measure the **confidence-trimmed** ECD, not
the whole annotated one. In `af_confidence.csv`, `n_trim` and `c_trim` are 0-based offsets
*within* the ECD rather than absolute residue positions -- `c_trim <= ecd_length` holds for
every row, and `ecd_length == ecd_end - ecd_start + 1`. So the measured range is::

    ecd_start + n_trim  ..  ecd_start + c_trim

Measuring that range reproduces the published `max_dim_first` and `min_dim_first` exactly:
100% of 559 proteins agree to within 0.05 A, median error 0.0003 A. Measuring the untrimmed
`ecd_start..ecd_end` instead matches only ~48%, because the untrimmed ends are the
low-confidence tails the trimming exists to remove.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def _setup():
    import pymol
    pymol.finish_launching(["pymol", "-qc"])
    from pymol import cmd
    return cmd


def iabb(cmd, selection: str):
    """Return the IABB dimensions of a selection, following Draw_Protein_Dimensions."""
    model = cmd.get_model(selection)
    if len(model.atom) < 3:
        return None

    # 1. centre of mass to the origin
    tot = sum(a.get_mass() for a in model.atom)
    cx = sum(a.coord[0] * a.get_mass() for a in model.atom) / tot
    cy = sum(a.coord[1] * a.get_mass() for a in model.atom) / tot
    cz = sum(a.coord[2] * a.get_mass() for a in model.atom) / tot
    cmd.transform_selection(selection, [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, -cx, -cy, -cz, 1])

    # 2. mass-weighted inertia tensor about the origin
    model = cmd.get_model(selection)
    I = np.zeros((3, 3))
    for a in model.atom:
        m = a.get_mass()
        x, y, z = a.coord
        I[0, 0] += m * (y * y + z * z)
        I[1, 1] += m * (x * x + z * z)
        I[2, 2] += m * (x * x + y * y)
        I[0, 1] -= m * x * y
        I[0, 2] -= m * x * z
        I[1, 2] -= m * y * z
    I[1, 0], I[2, 0], I[2, 1] = I[0, 1], I[0, 2], I[1, 2]

    # 3. eigenvectors as rows, ascending eigenvalue, rounded as the original does
    val, vec = np.linalg.eig(I)
    R = np.around(vec[:, np.argsort(val)].T, 3)
    if np.linalg.det(R) < 0:
        R = R[::-1]

    # 4. rotate and take extents
    transf = np.transpose(R)
    cmd.transform_selection(
        selection,
        [transf[0][0], transf[0][1], transf[0][2], 0,
         transf[1][0], transf[1][1], transf[1][2], 0,
         transf[2][0], transf[2][1], transf[2][2], 0,
         0, 0, 0, 1],
        homogenous=0, transpose=1)
    (mn, mx) = cmd.get_extent(selection)
    return sorted((mx[0] - mn[0], mx[1] - mn[1], mx[2] - mn[2]), reverse=True)


def measure(models: Path, regions: pd.DataFrame,
            pattern: str = "AF-{acc}-F1-model_v6.pdb") -> pd.DataFrame:
    """Measure each `ID, ecd_start, ecd_end` row and return the dimensions frame.

    Shared with ``build_af_dimensions.py``, which derives the regions from the confidence
    table rather than reading them from a file -- so both go through one implementation of
    the measurement and cannot drift apart.
    """
    cmd = _setup()
    reg = regions.dropna(subset=["ecd_start", "ecd_end"])

    rows = []
    for i, r in enumerate(reg.itertuples(), 1):
        acc = r.ID
        # allow 'ACC|tag' ids so one model can be measured over several regions
        path = models / pattern.format(acc=acc.split('|')[0])
        if not path.exists():
            continue
        lo, hi = int(r.ecd_start), int(r.ecd_end)
        if hi <= lo:
            continue
        cmd.reinitialize()
        cmd.load(str(path), "m")
        cmd.remove("solvent or hydrogens")
        sel = f"m and resi {lo}-{hi}"
        n = cmd.count_atoms(sel)
        if n < 3:
            continue
        plddt = float(np.mean([at.b for at in cmd.get_model(sel).atom]))
        dims = iabb(cmd, sel)
        if dims is None:
            continue
        rows.append({"ID": acc, "n_atoms": n,
                     "dim_1": dims[0], "dim_2": dims[1], "dim_3": dims[2],
                     "mean_plddt": plddt})
        if i % 200 == 0:
            print(f"  {i}/{len(reg)} measured {len(rows)}", flush=True)

    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", required=True, type=Path)
    ap.add_argument("--regions", required=True, type=Path,
                    help="CSV with ID, ecd_start, ecd_end")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--pattern", default="AF-{acc}-F1-model_v6.pdb")
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()

    reg = pd.read_csv(a.regions)
    if a.limit:
        reg = reg.head(a.limit)
    out = measure(a.models, reg, a.pattern)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, index=False)
    print(f"\nmeasured {len(out):,} structures -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
