#!/usr/bin/env python
"""An epitope's height, the gap an antibody on it holds, and what that gap excludes.

    python code/predict_epitope.py ERBB2 579-625 --models DIR
    python code/predict_epitope.py --antibody trastuzumab --models DIR
    python code/predict_epitope.py --curated --models DIR --out antibody_gaps.csv

An epitope is typed, as a gene name or UniProt accession of the height table and its residues
(UniProt positions, as ranges and single residues), or named by its antibody in
``data/curated/antibody_epitopes.csv``; ``--curated`` runs every antibody of that table. Typed
residues are uncurated: nothing checks that an antibody binds them, and the output says so.

For each Fc receptor the table gives the gap, whether the epitope is within 10 nm of its
membrane, and whether each protein of ``data/curated/segregation_probes.csv`` (CD45, its short
isoform CD45RO and CD148) fits in the gap; ``--protein`` adds others from the height table.

The antigens of the curated antibodies have their models' residue positions committed
(``data/inputs/alphafold/alphafold_axis_positions.csv``). Any other protein needs ``--models``, a
directory of AlphaFold v6 models (``code/database/inputs/download_alphafold.py`` fetches them);
without its model a residue inside it is placed by its share of the model's residues, and the
output says so. The tables are a run's, read from the working
directory if it is one and from the committed snapshot otherwise.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from surfaceome_config import resolve
from surfaceomeTopography.epitope import antigen_epitope, predict_bridge, probe_heights


def residue_list(text: str) -> list[int]:
    """``579-625,630`` as residues."""
    out = []
    for part in text.split(","):
        first, _, last = part.partition("-")
        out.extend(range(int(first), int(last or first) + 1))
    return out


def cited(row) -> str:
    """Where a curated epitope comes from, for the output."""
    how = f"PDB {row['pdb']}" if row["basis"] == "structure" else row["basis"]
    return f"curated: {row['antibody']}, {how}, PMID {row['pmids']}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("antigen", nargs="?", help="gene name or UniProt accession")
    ap.add_argument("residues", nargs="?", help="UniProt positions, like 579-625,630")
    ap.add_argument("--antibody", help="an antibody of antibody_epitopes.csv, in place of antigen and residues")
    ap.add_argument("--curated", action="store_true", help="every antibody of antibody_epitopes.csv")
    ap.add_argument("--models", type=Path, help="directory of AlphaFold v6 models")
    ap.add_argument("--protein", action="append", default=[], help="another protein to place in the gap, by gene")
    ap.add_argument("--margin", type=float, default=5.0, help="nm taller than the gap at which a protein is excluded")
    ap.add_argument("--out", type=Path, help="write the table here as well")
    a = ap.parse_args()

    estimates = pd.read_csv(resolve("database/height_estimates.csv"))
    assignments = pd.read_csv(resolve("database/domain_disorder_assignments_alphafold.csv"), index_col=0, low_memory=False)
    ecds = pd.read_csv(resolve("database/UP000005640_surfaceome_largest_ecds_classified.csv"))
    axes = pd.read_csv(resolve("inputs/alphafold/alphafold_axis_positions.csv"))
    named = pd.read_csv(resolve("database/protein_heights_for_interactions.csv"))
    gene = named["Gene names"].fillna("").str.split().str[0]
    height = pd.Series(named["total_height"].to_numpy(), index=gene).groupby(level=0).first()
    accession = dict(zip(gene, named["ID link"]))
    by_accession = dict(zip(named["ID link"], named["total_height"]))

    curated = pd.read_csv(resolve("curated/antibody_epitopes.csv"), dtype=str, keep_default_na=False)
    curated = curated[curated["included"] == "yes"]
    if a.curated:
        asked = [(r["accession"], r["residues"], cited(r)) for _, r in curated.iterrows()]
    elif a.antibody:
        match = curated[curated["antibody"].str.contains(a.antibody, case=False, regex=False)]
        if len(match) != 1:
            ap.error(f"{a.antibody} names {len(match)} antibodies with an epitope; there are: {', '.join(curated['antibody'])}")
        r = match.iloc[0]
        asked = [(r["accession"], r["residues"], cited(r))]
    elif a.antigen and a.residues:
        asked = [(accession.get(a.antigen, a.antigen), a.residues, "uncurated: residues as typed")]
    else:
        ap.error("give an antigen and its residues, --antibody NAME or --curated")

    bridges = pd.read_csv(resolve("curated/antibody_bridges.csv"))
    receptors = bridges[(bridges["role"] == "Fc receptor") & (bridges["included"] == "yes")]
    receptors = receptors.assign(height=receptors["accession"].map(by_accession)).dropna(subset=["height"])
    missing = [g for g in a.protein if g not in height]
    if missing:
        ap.error(f"not in the height table: {', '.join(missing)}")
    probes = pd.concat([probe_heights(pd.read_csv(resolve("curated/segregation_probes.csv")), estimates, assignments),
                        height[a.protein]])
    print("placed in each gap: " + ", ".join(f"{name} {h:.1f} nm" for name, h in probes.items()) + "\n")

    tables = []
    for antigen, residues, source in asked:
        try:
            epitope = antigen_epitope(antigen, residue_list(residues), estimates, assignments, ecds, a.models, axes)
        except (KeyError, ValueError) as problem:
            ap.error(str(problem).strip("'\""))
        table = predict_bridge(epitope, receptors, probes, a.margin)
        print(f"{antigen} {epitope['ecd']}, {epitope['antigen_height']:.1f} nm tall. {source}")
        print(f"epitope at {epitope['height']:.1f} nm ({epitope['lowest']:.1f} to {epitope['highest']:.1f}), in "
              f"{', '.join(epitope['kinds'])}; {'within' if table['within_phagocytosis_range'].all() else 'beyond'} "
              f"10 nm of its membrane")
        for note in epitope["notes"]:
            print(f"  note: {note}")
        if epitope["outside"]:
            print(f"  outside the ectodomain, left out: {len(epitope['outside'])} residues")
        shown = table.drop(columns=["antigen_height", "epitope_height", "within_phagocytosis_range"])
        print(shown.round(1).to_string(index=False) + "\n")
        tables.append(table.assign(antigen=antigen, residues=residues, source=source,
                                   epitope_lowest=epitope["lowest"], epitope_highest=epitope["highest"],
                                   notes="; ".join(epitope["notes"])))
    if a.out:
        pd.concat(tables, ignore_index=True).to_csv(a.out, index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
