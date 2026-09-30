#!/usr/bin/env python
"""Map each surfaceome accession to the AlphaFold model used for it.

    python code/database/inputs/build_af_mapping.py --surfaceome CSV --af-confidence CSV --out TSV

Rebuilds `data/inputs/proteome/surfaceome_current.tsv`, which the pipeline reads to join
AlphaFold measurements onto surfaceome proteins:

    current<TAB>af_id       'current' is a comma-separated list of accessions
                            sharing the model 'af_id'

The archived copy is a **manual UniProt ID-mapping export from 2022** and was never
regenerated. That is a problem once the proteome moves, because it fixes the accession
universe at its 2022 state: UniProt 2026_02 merged the HLA allele entries, and this file
puts them straight back. Its four comma-bearing rows carry 79 accessions between them --
P01889 x35, P04439 x21, P10321 x14, P01911 x9 -- so 75 accessions UniProt has retired
re-enter the run, each inheriting a model it does not have of its own.

What that cost, measured on the 2026 run before this script existed:

* the estimates table reported **3,285 ectodomains against a 3,210-protein surfaceome**;
* the same model's row was duplicated once per alias ahead of the per-ECD `groupby.agg`,
  so P01889 summed `alphafold_len` to 9,940 residues over a 285-residue ectodomain --
  the 3,488% coverage outlier previously read as a coverage-metric defect;
* the retired alleles took **12% of the B cell expression weight**, moving that panel's
  weighted mean by 1.19 nm;
* `F2` cell 37 and `03_interaction_heights` cell 56 both failed outright, because the
  outer merge against a stale accession list produces all-NaN rows.

Heights themselves escaped: all four affected ectodomains carry `methods == 'structure'`,
so the solved structure overrode the inflated sum.

Construction
------------
One row per surfaceome accession that has an AlphaFold model. Models are downloaded per
accession (`download_alphafold.py` fetches `AF-<accession>-F1`), so `af_id` equals the
accession and no aliasing arises -- which is the point: an alias only ever appeared because
the 2022 export folded UniProt's then-separate allele entries onto one model.

Validation against the archived 2022 file: it holds 3,079 rows, 3,077 of them identity, and
its exploded `current` column reproduces the 2022 surfaceome exactly (3,153 accessions).
`--check` reports that agreement rather than asserting it, since the 77 alias rows come from
UniProt's ID-mapping service at a particular release and cannot be re-derived from the model
files alone.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def build(surfaceome: pd.DataFrame, af_ids: set[str]) -> pd.DataFrame:
    accs = sorted(set(surfaceome["ID link"].dropna().unique()) & af_ids)
    return pd.DataFrame({"current": accs, "af_id": accs})


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--surfaceome", required=True,
                    help="UP000005640_surfaceome_largest_ecds_classified.csv")
    ap.add_argument("--af-confidence", required=True,
                    help="af_confidence.csv -- its ID column is the set of modelled accessions")
    ap.add_argument("--out", required=True)
    ap.add_argument("--check", help="an existing surfaceome_current.tsv to compare against")
    args = ap.parse_args()

    sf = pd.read_csv(args.surfaceome, low_memory=False)
    afc = pd.read_csv(args.af_confidence)
    out = build(sf, set(afc["ID"].dropna().unique()))

    modelled = set(afc["ID"].dropna().unique())
    surf = set(sf["ID link"].dropna().unique())
    print(f"surfaceome proteins        {len(surf):,}")
    print(f"accessions with a model    {len(modelled):,}")
    print(f"  modelled but not in the surfaceome: {len(modelled - surf)}")
    print(f"mapping rows written       {len(out):,}")

    if args.check:
        old = pd.read_csv(args.check, sep="\t", header=None, names=["current", "af_id"])
        exp = old.assign(current=old.current.astype(str).str.split(",")).explode("current")
        alias = exp[exp.current != exp.af_id]
        print(f"\nagainst {Path(args.check).name}:")
        print(f"  rows {len(old):,}, exploded {len(exp):,}, identity {(exp.current == exp.af_id).sum():,}")
        print(f"  alias rows (accession borrowing another's model): {len(alias):,}")
        if len(alias):
            print(f"  models shared: {sorted(alias.af_id.unique())}")
        stale = set(exp.current) - surf
        print(f"  accessions it names that are NOT in this surfaceome: {len(stale)}")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", header=False, index=False)
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
