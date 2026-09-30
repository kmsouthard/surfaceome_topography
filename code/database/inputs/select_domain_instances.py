#!/usr/bin/env python
"""Choose the solved-structure chains to measure each Pfam family on.

    python code/database/inputs/select_domain_instances.py --sifts-pfam pdb_chain_pfam.tsv.gz \\
        --sifts-taxonomy pdb_chain_taxonomy.tsv.gz --assignments CSV [--clans CSV] \\
        [--per-family 20] --out CSV

The 2020 domain heights were measured on PDB chains listed by an RCSB advanced search over
human and mouse entries, one row per chain and Pfam family. SIFTS now publishes that mapping
(``pdb_chain_pfam.tsv``: PDB, CHAIN, SP_PRIMARY, PFAM_ID, COVERAGE) weekly, with
``pdb_chain_taxonomy.tsv`` for the species, so the selection is a download plus this script.

Which families: every family the surfaceome's Pfam assignments use (``--assignments``), plus,
with ``--clans``, every family in the clan-height table, so families that inherit a clan height
are measured too. Which chains: human first, then mouse, keeping at most ``--per-family``
distinct PDB entries per family and preferring chains SIFTS marks as fully covering the
family (``COVERAGE = 1``). Within that, entries are taken in id order, which is arbitrary but
reproducible; resolution is not in SIFTS and is not used, as it was not in 2020.

The output lists one row per (PDB, CHAIN, PFAM_ID) with the species and the UniProt accession,
and is what ``download_pdb_structures.py`` fetches and ``measure_domain_instances.py`` measures.
Families with no human or mouse chain are written beside it as ``<out>.unmeasurable.txt``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

HUMAN, MOUSE = 9606, 10090


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sifts-pfam", type=Path, required=True)
    ap.add_argument("--sifts-taxonomy", type=Path, required=True)
    ap.add_argument("--assignments", type=Path, required=True,
                    help="the pipeline's Pfam assignment table (a `Pfam` column, PFxxxxx.n)")
    ap.add_argument("--clans", type=Path, help="the clan-height table, to add the families it measured (height_pfam set)")
    ap.add_argument("--per-family", type=int, default=20, help="max distinct PDB entries per family")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()

    fams = set(pd.read_csv(a.assignments, usecols=["Pfam"], low_memory=False)["Pfam"]
               .dropna().str.split(".").str[0])
    if a.clans:
        # only families with a measurement of their own; the rest of the table inherits a clan height
        clans = pd.read_csv(a.clans, usecols=["Pfam", "height_pfam"])
        fams |= set(clans.dropna(subset=["height_pfam"])["Pfam"].dropna().str.split(".").str[0])
    print(f"families to measure: {len(fams)}")

    pf = pd.read_csv(a.sifts_pfam, sep="\t", comment="#", low_memory=False,
                     names=["PDB", "CHAIN", "SP_PRIMARY", "PFAM_ID", "COVERAGE"])
    tx = pd.read_csv(a.sifts_taxonomy, sep="\t", comment="#", low_memory=False,
                     names=["PDB", "CHAIN", "TAX_ID", "SCIENTIFIC_NAME"])
    tx["TAX_ID"] = pd.to_numeric(tx["TAX_ID"], errors="coerce")   # the header line reads as a row
    tx = tx[tx["TAX_ID"].isin([HUMAN, MOUSE])].drop_duplicates(["PDB", "CHAIN"])
    pf = pf[pf["PFAM_ID"].isin(fams)].merge(tx[["PDB", "CHAIN", "TAX_ID"]], on=["PDB", "CHAIN"])
    pf["species"] = pf["TAX_ID"].map({HUMAN: "human", MOUSE: "mouse"})
    pf["full"] = pd.to_numeric(pf["COVERAGE"], errors="coerce").eq(1)
    print(f"chains in SIFTS with one of them, human or mouse: {len(pf):,} over {pf['PDB'].nunique():,} entries")

    # human before mouse, full coverage before partial, then id order; cap distinct entries
    pf = pf.sort_values(["PFAM_ID", "species", "full", "PDB", "CHAIN"],
                        ascending=[True, True, False, True, True])
    keep = []
    for fam, g in pf.groupby("PFAM_ID", sort=False):
        entries = g["PDB"].drop_duplicates().head(a.per_family)
        keep.append(g[g["PDB"].isin(entries)].drop_duplicates(["PDB", "CHAIN"]))
    sel = pd.concat(keep, ignore_index=True)[["PFAM_ID", "PDB", "CHAIN", "SP_PRIMARY", "species", "full"]]
    sel.to_csv(a.out, index=False)

    missing = sorted(fams - set(sel["PFAM_ID"]))
    Path(str(a.out) + ".unmeasurable.txt").write_text("\n".join(missing) + "\n")
    by = sel.groupby("PFAM_ID")["PDB"].nunique()
    print(f"selected: {len(sel):,} chains, {sel['PDB'].nunique():,} entries, "
          f"{sel['PFAM_ID'].nunique()} families ({(sel['species'] == 'mouse').sum():,} mouse chains)")
    print(f"entries per family: median {by.median():.0f}, min {by.min()}, families at the cap: {(by >= a.per_family).sum()}")
    print(f"families with no human or mouse chain: {len(missing)} -> {a.out}.unmeasurable.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
