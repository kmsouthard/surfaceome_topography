#!/usr/bin/env python
"""Check the expression-weighted immune panels against a primary-cell surfaceome (nanoMAPS).

    python code/figures/validation/check_nanomaps_surfaceome.py --run OUT_DIR \\
        --source-data "Source Data Main Figures.xlsx" [--supplementary 41467_2026_77876_MOESM3_ESM.xlsx] \\
        [--out CSV]

Lorentzian, Bartleson et al., Nature Communications 2026 (doi:10.1038/s41467-026-77876-4;
raw data MassIVE MSV000097303) profile the surface proteins of sorted primary human monocytes,
B cells and T cells from 500 cells each: surface biotinylation, streptavidin capture in a
microlitre droplet, label-free LC-MS, and a call of "cell-surface protein" for anything in the
Surfy list that is enriched over an unlabelled control. That is a *presence* measurement at
the membrane of a primary cell, which the whole-cell proteome our panels weight by (Expression
Atlas E-PROT-1) cannot give. This script asks, per panel:

* how many of the proteins the panel weights were captured at the surface, and what share
  of the panel's abundance weight they carry;
* which abundant proteins the surface capture never saw;
* which captured surface proteins the panel has no expression value for.

With ``--supplementary`` it also compares the Surfy list (their Supplementary Data 1) with
our surfaceome. The paper's files are CC BY-NC-ND and are not redistributed here: download
"Source Data" (zip, Main Figures workbook) and "Supplementary Data" from the article page.
Intensities from an enriched fraction fold in capture efficiency, so this is a check on
presence, not a replacement for the abundance weights.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from surfaceome_config import resolve  # noqa: E402

#: our per-cell-type expression tables (written by F3) -> the paper's Figure 4B block
PANELS = {"monocyte": "Monocytes", "bcell": "B-Cells", "cd4tcell": "T-Cells", "cd8tcell": "T-Cells"}


def acc(s: pd.Series) -> pd.Series:
    """First accession of a protein group, isoform suffix dropped."""
    return s.astype(str).str.split(";").str[0].str.split("-").str[0]


def figure4b_blocks(xlsx: Path) -> dict[str, pd.DataFrame]:
    """The three cell-type blocks stacked in the Figure 4B sheet, keyed by cell type."""
    raw = pd.read_excel(xlsx, sheet_name="Figure4B", header=None)
    starts = [i for i in range(len(raw) - 1)
              if isinstance(raw.iat[i, 0], str) and raw.iat[i + 1, 0] == "Protein.ID"]
    blocks = {}
    for k, i in enumerate(starts):
        j = starts[k + 1] if k + 1 < len(starts) else len(raw)
        blk = raw.iloc[i + 2:j].copy()
        blk.columns = list(raw.iloc[i + 1])
        blk = blk.dropna(subset=["Protein.ID"])
        blk["acc"] = acc(blk["Protein.ID"])
        blk["surface"] = blk["CS_protein"].astype(str) == "True"
        blk["enriched"] = blk["Expression"].astype(str) == "Biotin_Enriched"
        blocks[raw.iat[i, 0]] = blk
    return blocks


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", type=Path, required=True, help="a run's output directory (tables/expression_eprot1_*.csv)")
    ap.add_argument("--source-data", type=Path, required=True, help='the paper\'s "Source Data Main Figures.xlsx"')
    ap.add_argument("--supplementary", type=Path, help="the paper's Supplementary Data workbook (Surfy list)")
    ap.add_argument("--out", type=Path, help="per-protein table (default: RUN/tables/nanomaps_validation.csv)")
    ap.add_argument("--top", type=int, default=8, help="abundant proteins to name per panel")
    a = ap.parse_args()

    blocks = figure4b_blocks(a.source_data)
    heights = pd.read_csv(resolve("database/height_estimates.csv"), usecols=["ID link_first", "Entry name_first"])
    surfaceome = set(acc(heights["ID link_first"]))
    rows = []
    print(f"{'panel':<10} {'weighted':>8} {'captured':>9} {'weight share':>12} {'enriched':>9} {'surface, no value':>18}")
    for panel, cell in PANELS.items():
        table = a.run / "tables" / f"expression_eprot1_{panel}.csv"
        if not table.is_file():
            print(f"{panel:<10} (no {table.name} in the run)")
            continue
        ours = pd.read_csv(table)
        ours["acc"] = acc(ours["ID link"])
        blk = blocks[cell]
        cs = blk[blk["surface"]].set_index("acc")
        ours["captured"] = ours["acc"].isin(cs.index)
        ours["enriched"] = ours["acc"].map(cs["enriched"]).fillna(False).astype(bool)
        weight = ours["percent_expression"].sum()
        share = ours.loc[ours["captured"], "percent_expression"].sum() / weight if weight else float("nan")
        unseen = set(cs.index) & surfaceome - set(ours["acc"])
        print(f"{panel:<10} {len(ours):>8} {int(ours['captured'].sum()):>9} {share:>11.0%} "
              f"{int(ours['enriched'].sum()):>9} {len(unseen):>18}")
        top = ours[~ours["captured"]].nlargest(a.top, "percent_expression")
        print("           abundant, never captured: " + ", ".join(f"{g} ({p:.1f}%)" for g, p in
                                                                 zip(top["Gene Name"], top["percent_expression"])))
        names = heights.set_index(acc(heights["ID link_first"]))["Entry name_first"]
        print("           captured, no expression value: " + ", ".join(
            sorted(str(names.get(x, x)).replace("_HUMAN", "") for x in unseen)[:16])
              + (" ..." if len(unseen) > 16 else ""))
        for _, r in ours.iterrows():
            rows.append({"panel": panel, "nanomaps_cell_type": cell, "accession": r["acc"], "gene": r["Gene Name"],
                         "percent_expression": r["percent_expression"], "total_height": r["total_height"],
                         "captured_at_surface": bool(r["captured"]), "biotin_enriched": bool(r["enriched"])})
        for x in sorted(unseen):
            rows.append({"panel": panel, "nanomaps_cell_type": cell, "accession": x, "gene": str(names.get(x, "")).replace("_HUMAN", ""),
                         "percent_expression": float("nan"), "total_height": float("nan"),
                         "captured_at_surface": True, "biotin_enriched": bool(cs.loc[x, "enriched"]) if x in cs.index else False})

    if a.supplementary and a.supplementary.is_file():
        surfy = pd.read_excel(a.supplementary, sheet_name="Supplementary Data 1")
        s = set(acc(surfy["UniProt accession"]))
        print(f"\nSurfy list: {len(s)} proteins; in our surfaceome {len(s & surfaceome)}; "
              f"Surfy only {len(s - surfaceome)}; ours only {len(surfaceome - s)} of {len(surfaceome)}")

    out = a.out or (a.run / "tables" / "nanomaps_validation.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
