#!/usr/bin/env python
"""Per-residue AlphaFold confidence for each ectodomain, and where to trim it.

    python code/database/inputs/build_af_confidence.py --models DIR --surfaceome CSV --out CSV

Reproduces the columns of the archived `af_confidence.csv`, which the pipeline reads to decide
which part of an AlphaFold model to measure and which stretches to call disordered:

    ID, length, ecd_start, ecd_end, ecd_length, mean_confidence,
    n_trim, c_trim, ecd_length_trimmed, mean_confidence_trimmed, confidence_string

`confidence_string` is one digit per ECD residue, pLDDT // 10, indexed from `ecd_start`.

The trimming rule was not documented anywhere and had to be recovered from the archived table:
**n_trim is the first ECD residue with pLDDT >= 70 and c_trim the last**, both as 0-based
offsets within the ECD. That reproduces the archived values for 97.7% of 3,000 entries checked,
and measuring `ecd_start + n_trim .. ecd_start + c_trim` reproduces the published bounding-box
dimensions exactly (100% within 0.05 A). pLDDT 70 is also AlphaFold's own boundary between
"confident" and "low".

pLDDT is read from the B-factor column, which is where AlphaFold stores it, taking one value per
residue from the CA atom.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

TRIM_PLDDT = 70


def residue_plddt(pdb_path: Path) -> dict[int, float]:
    """{residue number: pLDDT} from a model's CA atoms."""
    out = {}
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    out[int(line[22:26])] = float(line[60:66])
                except ValueError:
                    continue
    return out


def confidence_row(acc: str, plddt: dict[int, float], ecd_start: int, ecd_end: int) -> dict | None:
    span = [plddt.get(r) for r in range(int(ecd_start), int(ecd_end) + 1)]
    vals = [v for v in span if v is not None]
    if not vals:
        return None
    digits = "".join(str(min(int(v // 10), 9)) if v is not None else "0" for v in span)
    idx = [i for i, v in enumerate(span) if v is not None and v >= TRIM_PLDDT]
    n_trim, c_trim = (idx[0], idx[-1]) if idx else (np.nan, np.nan)
    trimmed = [v for v in span[int(n_trim):int(c_trim) + 1] if v is not None] if idx else []
    return {
        "ID": acc,
        "length": max(plddt) if plddt else np.nan,
        "ecd_start": ecd_start,
        "ecd_end": ecd_end,
        "ecd_length": len(span),
        "mean_confidence": float(np.mean(vals)),
        "n_trim": n_trim,
        "c_trim": c_trim,
        "ecd_length_trimmed": (c_trim - n_trim + 1) if idx else np.nan,
        "mean_confidence_trimmed": float(np.mean(trimmed)) if trimmed else np.nan,
        "confidence_string": digits,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", required=True, type=Path)
    ap.add_argument("--surfaceome", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--pattern", default="AF-{acc}-F1-model_v6.pdb")
    ap.add_argument("--sources", default="ecd,gpi,long",
                    help="comma-separated `source` values to include (default: modelled ones)")
    a = ap.parse_args()

    surf = pd.read_csv(a.surfaceome, low_memory=False)
    surf = surf[surf["source"].isin(a.sources.split(","))].dropna(subset=["start", "end"])
    # one row per (protein, ECD); take the longest span per protein to match the archived table
    surf = surf.sort_values("seq_len", ascending=False).drop_duplicates("ID link")

    rows, missing = [], []
    for i, r in enumerate(surf.itertuples(), 1):
        acc = getattr(r, "_1")          # 'ID link'
        path = a.models / a.pattern.format(acc=acc)
        if not path.is_file():
            missing.append(acc)
            continue
        row = confidence_row(acc, residue_plddt(path), r.start, r.end)
        if row:
            rows.append(row)
        if i % 500 == 0:
            print(f"  {i}/{len(surf)}  built {len(rows)}", flush=True)

    out = pd.DataFrame(rows)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, index=False)
    print(f"\n{len(out):,} ectodomains -> {a.out}")
    if missing:
        print(f"{len(missing):,} without a model, e.g. {missing[:5]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
