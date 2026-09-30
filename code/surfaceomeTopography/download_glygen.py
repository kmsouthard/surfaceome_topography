#!/usr/bin/env python3
"""Experimental glycosylation sites for the human proteome, from GlyGen.

The pipeline's glycan density has always come from GlycoMine (Li et al., *Bioinformatics*
2015;31:1411-19), a random forest whose archived output is dated 2014-03-16. It was trained on
68 C-linked, 416 N-linked and 649 O-linked verified sites, treating every unverified N/S/T/W as
a negative -- which in 2014 meant most genuine glycosites were training negatives.

GlyGen publishes per-site human glycosylation curated from experiment, mapped to UniProt
accessions with ECO evidence codes. Benchmarked against it, GlycoMine has 51.4% recall and a
per-segment density correlation of 0.705.

Neither source is used alone. On heavily glycosylated mucins the per-site databases are badly
incomplete, because a VNTR's glycans are numerous and never individually mapped -- for CD43,
GlycoMine gives 0.404 and GlyGen 0.123 against Cyster's measured 0.333. Elsewhere GlycoMine
misses sites GlyGen has, notably PSGL-1 (0.012 vs 0.039). Taking the **union** of the two
recovers the second class without deflating the first; it raises total disorder height by ~5%.

O-GlcNAcylation is excluded. It is a nucleocytoplasmic modification, so it cannot contribute to
an ectodomain's glycan density, and GlyGen carries it in the same files.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import requests

BASE = "https://data.glygen.org/ln2data/releases/data/current/reviewed"

#: Experimental and curated sources. Deliberately excludes the `oglcnac_*` files
#: (nucleocytoplasmic) and the `predicted_*` files, which would defeat the purpose of using
#: this as an independent check on a predictor.
SOURCES = (
    "uniprotkb", "unicarbkb", "glyconnect", "literature",
    "literature_mining_manually_verified", "gptwiki", "harvard",
    "pdb", "rcsb_pdb", "platelet", "c_man",
)

#: Subtypes to drop even when they appear in an included file.
EXCLUDE_SUBTYPE = ("GlcNAcylation",)


def _fetch(source: str, cache: Path | None) -> pd.DataFrame:
    name = f"human_proteoform_glycosylation_sites_{source}.csv"
    if cache is not None:
        cache.mkdir(parents=True, exist_ok=True)
        local = cache / name
        if not local.is_file():
            r = requests.get(f"{BASE}/{name}", timeout=300)
            r.raise_for_status()
            local.write_bytes(r.content)
        return pd.read_csv(local, low_memory=False)
    r = requests.get(f"{BASE}/{name}", timeout=300)
    r.raise_for_status()
    from io import StringIO
    return pd.read_csv(StringIO(r.text), low_memory=False)


def download_glygen_sites(sources=SOURCES, cache: Path | None = None) -> pd.DataFrame:
    """Return a frame of (acc, site, amino_acid, glycosylation_type, src), deduplicated."""
    keep = ["uniprotkb_canonical_ac", "glycosylation_site_uniprotkb",
            "amino_acid", "glycosylation_type", "glycosylation_subtype"]
    frames = []
    for s in sources:
        d = _fetch(s, cache)
        cols = [c for c in keep if c in d.columns]
        if "uniprotkb_canonical_ac" not in cols or "glycosylation_site_uniprotkb" not in cols:
            continue
        d = d[cols].copy()
        d["src"] = s
        frames.append(d)
    g = pd.concat(frames, ignore_index=True)
    if "glycosylation_subtype" in g.columns:
        pat = "|".join(EXCLUDE_SUBTYPE)
        g = g[~g["glycosylation_subtype"].fillna("").str.contains(pat, case=False)]
    g["acc"] = g["uniprotkb_canonical_ac"].astype(str).str.split("-").str[0]
    g["site"] = pd.to_numeric(g["glycosylation_site_uniprotkb"], errors="coerce")
    g = g.dropna(subset=["site"])
    g["site"] = g["site"].astype(int)
    return g[["acc", "site", "glycosylation_type", "src"]].drop_duplicates()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, type=Path, help="CSV to write")
    ap.add_argument("--cache", type=Path, default=None, help="directory to cache raw downloads")
    a = ap.parse_args()
    g = download_glygen_sites(cache=a.cache)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    g.to_csv(a.out, index=False)
    print(f"{len(g):,} unique sites across {g.acc.nunique():,} proteins -> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
