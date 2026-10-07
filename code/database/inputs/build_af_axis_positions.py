#!/usr/bin/env python
"""AlphaFold models → where each residue sits along its model's box axes, for the epitope heights.

    python code/database/inputs/build_af_axis_positions.py --models DIR --af-confidence CSV --out FILE

An epitope's height needs its residues' positions inside the antigen's model
(`surfaceomeTopography.epitope`), and a run reads no model files. This measures them once, for
the antigens of ``data/curated/antibody_epitopes.csv``: per residue of the model's confident
span, its C-alpha along each inertia axis of the span, nm from the low end of the box, with the
box's three dimensions. The span is `build_af_dimensions.regions`'s, so the box is the one the
model's height was measured on.

``--accessions`` measures other proteins instead, one accession per line.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "code"))
from surfaceomeTopography.epitope import model_axes  # noqa: E402

from build_af_dimensions import regions  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", required=True, type=Path)
    ap.add_argument("--af-confidence", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--epitopes", type=Path, default=REPO / "data/curated/antibody_epitopes.csv")
    ap.add_argument("--accessions", type=Path, help="measure these instead of the epitopes' antigens")
    ap.add_argument("--pattern", default="AF-{acc}-F1-model_v6.pdb")
    a = ap.parse_args()

    wanted = (a.accessions.read_text().split() if a.accessions
              else sorted(set(pd.read_csv(a.epitopes)["accession"])))
    spans = regions(pd.read_csv(a.af_confidence))
    spans = spans.assign(accession=spans["ID"].str.split("[").str[0])
    out = []
    for acc in wanted:
        model = a.models / a.pattern.format(acc=acc)
        for _, r in spans[spans["accession"] == acc].iterrows():
            if not model.is_file():
                print(f"  {acc}: no model in {a.models}")
                continue
            axes = model_axes(model, r["ecd_start"], r["ecd_end"]).round(3)
            out.append(axes.assign(accession=acc, start=r["ecd_start"], end=r["ecd_end"]))
    table = pd.concat(out)[["accession", "start", "end", "residue", "axis_1", "axis_2", "axis_3",
                            "box_1", "box_2", "box_3"]]
    a.out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(a.out, index=False)
    print(f"  {table['accession'].nunique()} models, {len(table):,} residues")
    print(f"  wrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
