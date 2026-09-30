#!/usr/bin/env python
"""Measure assembled interaction models with the same bounding box as every other measurement.

    python code/figures/validation/measure_interaction_models.py [--models DIR] [--out FILE]

An assembled model is two or more structures superposed in PyMOL into one bound complex, or one
protein built from fragments in a conformation no single structure shows -- the extended-open
integrin, for one. Each `.pdb.gz` in `--models` is measured whole, its protein atoms only, with
`measure_structures.iabb`, and written as one row: `model`, `n_atoms`, `dim_1` >= `dim_2` >= `dim_3`
(Angstroms). Every model reproduces the measurement stored with it in 2018 to within 0.02 A;
`interaction_models.csv` records where each came from and the partners it joins.

The models were assembled in 2018 in PyMOL sessions kept outside the surfaceome archive
(`macrophage_surface_analysis/bmdm/ligand receptor models/`). The committed files are those
sessions' molecule objects saved as PDB; see `data/README.md` for which objects.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "database" / "inputs"))  # measure_structures lives with the input build
from measure_structures import _setup, iabb

REPO = Path(__file__).resolve().parents[3]  # this file is code/figures/validation/measure_interaction_models.py
MODELS = REPO / "data/measurements/models_2018"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--models", type=Path, default=MODELS)
    parser.add_argument("--out", type=Path, default=MODELS / "interaction_model_dimensions.csv")
    args = parser.parse_args()

    cmd = _setup()
    rows = []
    for path in sorted(args.models.glob("*.pdb.gz")):
        cmd.reinitialize("everything")
        cmd.load(str(path))
        cmd.remove("not polymer")
        dims = iabb(cmd, "all")
        rows.append({"model": path.name.removesuffix(".pdb.gz"), "n_atoms": cmd.count_atoms("all"),
                     "dim_1": round(dims[0], 2), "dim_2": round(dims[1], 2), "dim_3": round(dims[2], 2)})
    table = pd.DataFrame(rows)
    table.to_csv(args.out, index=False)
    print(table.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
