#!/usr/bin/env python
"""AlphaFold models → the dimensions table the pipeline reads, in one step.

    python code/database/inputs/build_af_dimensions.py --models DIR --af-confidence CSV --out FILE
    python code/database/inputs/build_af_dimensions.py ... --check data/inputs/alphafold/alphafold_dimensions.txt

Two steps used to sit between `measure_structures.py` and the file the combined-estimates
notebook reads, and neither existed as code -- they were done by hand and the reasoning was
recoverable only by inspecting the outputs:

1. **Which residues to measure.** Not the annotated ectodomain but the *confidence-trimmed*
   one: `ecd_start + n_trim … ecd_start + c_trim`, from `af_confidence.csv`. Measuring the
   untrimmed span reproduces only ~48% of the published values; the trimmed span reproduces
   100% to within 0.05 Å. See `measure_structures.py` for how that was established.

2. **The format.** The notebook reads `210923_alphafold_dist_UP000000589.txt` with
   `sep='\\t', names=['ID','x','y','z']` -- headerless, tab separated, four columns. The
   measurement writes six, with a header. The reduction is `dim_1..3` -> `x,y,z`, dropping
   `n_atoms` and `mean_plddt`.

Both are now here, so the chain from a downloaded model to the table a height is computed
from is scripted end to end. `--check` compares the result against an existing copy and
reports the largest disagreement, which is how this was verified against the vendored file.

Rows with no model, or with a trimmed span of nothing, drop out -- 2,742 of 2,786 in the 2026
run. `--regions-out` keeps the intermediate for inspection.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from measure_structures import measure


def regions(af_confidence: pd.DataFrame) -> pd.DataFrame:
    """The confidence-trimmed span to measure, per ectodomain."""
    d = af_confidence.dropna(subset=["ecd_start", "n_trim", "c_trim"])
    out = pd.DataFrame({
        "ID": d["ID"],
        "ecd_start": (d["ecd_start"] + d["n_trim"]).astype(int),
        "ecd_end": (d["ecd_start"] + d["c_trim"]).astype(int),
    })
    return out[out.ecd_end > out.ecd_start].reset_index(drop=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", required=True, type=Path)
    ap.add_argument("--af-confidence", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--pattern", default="AF-{acc}-F1-model_v6.pdb")
    ap.add_argument("--regions-out", type=Path, help="keep the trimmed spans")
    ap.add_argument("--check", type=Path, help="compare against an existing dimensions file")
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()

    conf = pd.read_csv(a.af_confidence)
    reg = regions(conf)
    if a.limit:
        reg = reg.head(a.limit)
    print(f"  {len(conf):,} ectodomains in the confidence table, {len(reg):,} with a "
          f"trimmed span to measure")
    if a.regions_out:
        a.regions_out.parent.mkdir(parents=True, exist_ok=True)
        reg.to_csv(a.regions_out, index=False)

    dims = measure(a.models, reg, a.pattern)
    print(f"  measured {len(dims):,}")

    a.out.parent.mkdir(parents=True, exist_ok=True)
    dims[["ID", "dim_1", "dim_2", "dim_3"]].to_csv(a.out, sep="\t", header=False, index=False)
    print(f"  wrote {a.out}")

    if a.check:
        want = pd.read_csv(a.check, sep="\t", header=None, names=["ID", "x", "y", "z"])
        got = dims.rename(columns={"dim_1": "x", "dim_2": "y", "dim_3": "z"})
        m = want.merge(got[["ID", "x", "y", "z"]], on="ID", how="outer",
                       suffixes=("_want", "_got"), indicator=True)
        counts = m._merge.value_counts().to_dict()
        both = m[m._merge == "both"]
        worst = max(float(np.nanmax((both[f"{c}_got"] - both[f"{c}_want"]).abs()))
                    for c in "xyz") if len(both) else float("nan")
        print(f"\n  against {a.check.name}: {counts.get('both', 0):,} shared, "
              f"{counts.get('right_only', 0)} new, {counts.get('left_only', 0)} lost")
        print(f"  largest disagreement on any axis: {worst:.6g} A")
        return 0 if (worst == 0 or worst < 5e-2) and not counts.get("left_only") else 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
