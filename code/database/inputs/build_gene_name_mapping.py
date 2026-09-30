#!/usr/bin/env python
"""Rebuild the surfaceome accession -> names mapping from UniProt.

    python code/database/inputs/build_gene_name_mapping.py --surfaceome CSV --out TSV [--check TSV]

Regenerates `data/inputs/mapping/uniprot_surfaceome_gene_name.tab`, which the pipeline
joins onto the surfaceome to attach `Entry name`, `Protein names` and `Gene names`. Those
columns are not cosmetic: `Entry name` labels the interaction heatmaps, and
`Gene names  (primary )` is the column the expression datasets merge on, so a protein
missing from this table carries no expression weight at all.

The archived copy is a 2022 UniProt export of the 2022 surfaceome, and it was not refreshed
with the rest. Against the 2026 surfaceome that leaves **169 proteins with no names** --
exactly the 169 UniProt 2026_02 added -- which surfaces as
`ValueError: Cannot mask with non-boolean array containing NA / NaN values` the first time
anything calls `.str.contains` on the column, and silently as missing expression weight
everywhere else.

The 2026 proteome download cannot supply this on its own: `TAB_FIELDS` omits `Entry name`,
`Status`, `Organism` and the primary gene name, so those are fetched here.

Column order matches the archived file exactly, because the notebooks index these by name
and one of them (`Gene names  (primary )`) has two spaces before the parenthesis.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # code/, for surfaceome_config and the package
from surfaceomeTopography.download_proteome import download_proteome  # noqa: E402

#: UniProt field ids for the archived file's columns, in its order.
FIELDS = ("accession", "id", "reviewed", "protein_name", "gene_names",
          "organism_name", "length", "gene_primary")

#: The archived header, reproduced verbatim -- note the double space in the last column.
HEADER = ["ID link", "Entry", "Entry name", "Status", "Protein names", "Gene names",
          "Organism", "Length", "Gene names  (primary )"]


def fetch(proteome: str, cache: Path) -> pd.DataFrame:
    """Whole-proteome pull with the name fields.

    Fetching per accession is not an option: a batched ``accession:X OR ...`` query long
    enough to be worth batching is rejected by UniProt with HTTP 400.  The proteome query
    is the one ``download_proteome`` already uses and is known to page correctly.
    """
    cache.parent.mkdir(parents=True, exist_ok=True)
    if not cache.is_file():
        written = download_proteome(proteome, "tsv", cache.parent, fields=FIELDS,
                                    legacy_headers=False, filename=cache.name)
        print(f"  cached at {written}")
    return pd.read_csv(cache, sep="\t", dtype=str)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--surfaceome", required=True)
    ap.add_argument("--proteome", default="UP000005640")
    ap.add_argument("--cache", required=True,
                    help="where to keep the whole-proteome name pull (reused if present)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--check", help="the archived mapping, to report coverage against")
    args = ap.parse_args()

    sf = pd.read_csv(args.surfaceome, low_memory=False)
    accs = set(sf["ID link"].dropna().unique())
    print(f"surfaceome proteins: {len(accs):,}")

    full = fetch(args.proteome, Path(args.cache))
    print(f"proteome rows fetched: {len(full):,}")
    acc_col = full.columns[0]
    df = full[full[acc_col].isin(accs)].drop_duplicates(subset=[acc_col]).copy()
    df.insert(0, "ID link", df[acc_col])
    df = df.iloc[:, :len(HEADER)]
    df.columns = HEADER

    got = set(df["ID link"])
    print(f"fetched {len(df):,} rows for {len(got):,} accessions")
    missing = accs - got
    if missing:
        print(f"  NOT returned by UniProt: {len(missing)}  {sorted(missing)[:10]}")

    if args.check:
        old = pd.read_csv(args.check, sep="\t", dtype=str)
        old_accs = set(old["ID link"])
        print(f"\nagainst {Path(args.check).name}: {len(old):,} rows")
        print(f"  surfaceome accessions it covers:      {len(accs & old_accs):,} of {len(accs):,}")
        print(f"  surfaceome accessions it MISSES:      {len(accs - old_accs):,}")
        print(f"  accessions it names but are now gone: {len(old_accs - accs):,}")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
