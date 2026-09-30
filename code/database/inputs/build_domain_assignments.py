#!/usr/bin/env python
"""Pfam domain assignments for ectodomains, from an hmmsearch table.

    python code/database/inputs/build_domain_assignments.py --domtbl FILE --clans CSV --out CSV

Turns raw HMMER output into the schema the pipeline reads, which had no generating script --
only the 2021 output file existed. Recovered from that file's own columns:

    ID link   accession parsed from the ECD id, e.g. 'P21589[27:549]' -> 'P21589'
    start,end ECD bounds from the same id
    seq_len   end - start          (the pipeline's convention throughout; not +1)
    dom_start start + align_from   verified against the archived table
    dom_end   start + align_to
    height    the measured clan height for that Pfam family

Search direction matters. The parser in `hmmer_assignment.hmmer_import` names column 1
`target_name` and column 4 `query_name`, and the archived table has sequence ids in the first
and family names in the fourth -- that is **hmmsearch** layout (target = sequence, query = HMM),
not hmmscan. `tlen` is therefore the ECD length and `qlen` the model length, which is what makes
the Methods' "90% coverage" filter `tlen > 0.9 * qlen` mean what it says. Run as::

    hmmsearch --domtblout OUT --cpu 8 -o /dev/null Pfam-A.hmm ecds.fasta

Families with no measured clan height are **dropped**, because the archived table is
post-join: all 346 of its families have heights, and common domains absent from the height
table -- PF00084 Sushi among them -- appear in neither. Keeping them is actively harmful, not
merely inert. The overlap resolution selects the highest-scoring non-overlapping set, and
Pfam 38.2 adds long composite families that outscore the individual repeats they span:
PF25024 (EGF_TEN, median span 178 aa) displaces NOTCH1's 30-residue EGF repeats, and since it
carries no height the protein's domain contribution collapses from 88.4 nm to 19.8 nm. Twelve
of 83 domain-counted proteins lost more than 30% of their height that way, every one of them
because a heightless family had been kept.

Pass ``--keep-heightless`` to retain them for inspection.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # code/, for surfaceome_config and the package
from surfaceomeTopography.hmmer_assignment import hmmer_import

ECD_ID = re.compile(r"^([^\[]+)\[(\d+):(\d+)\]$")


def build(domtbl: Path, clans: Path, evalue=0.0001, i_evalue=2,
          keep_heightless: bool = False) -> pd.DataFrame:
    d = hmmer_import(str(domtbl), filter=True, Evalue=evalue, i_Evalue=i_evalue).copy()
    d = d.loc[:, ~d.columns.str.startswith("Unnamed")]

    parsed = d["target_name"].astype(str).str.extract(ECD_ID)
    d["ID link"] = parsed[0]
    d["start"] = pd.to_numeric(parsed[1], errors="coerce")
    d["end"] = pd.to_numeric(parsed[2], errors="coerce")
    d = d.dropna(subset=["start", "end"])
    d["seq_len"] = d["end"] - d["start"]
    d["dom_start"] = d["start"] + d["align_from"]
    d["dom_end"] = d["start"] + d["align_to"]

    heights = (pd.read_csv(clans, low_memory=False)[["Pfam", "height"]]
               .dropna().drop_duplicates("Pfam"))
    d = d.merge(heights, on="Pfam", how="left")
    if not keep_heightless:
        d = d[d["height"].notna()]

    cols = ["target_name", "tlen", "query_name", "Pfam", "qlen", "E-value_full", "score_full",
            "bias_full", "dom_num", "type_count", "c-Evalue_dom", "i-Evalue_dom", "score_dom",
            "bias_dom", "hmm_from", "hmm_to", "align_from", "align_to", "env_from", "env_to",
            "accuracy", "height", "ID link", "start", "end", "dom_start", "dom_end", "seq_len"]
    return d[[c for c in cols if c in d.columns]]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--domtbl", required=True, type=Path)
    ap.add_argument("--clans", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--keep-heightless", action="store_true",
                    help="retain families with no measured clan height (default: drop)")
    a = ap.parse_args()

    d = build(a.domtbl, a.clans, keep_heightless=a.keep_heightless)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(a.out, index=True)          # the pipeline drops the leading 'Unnamed: 0'

    print(f"{len(d):,} assignments over {d['target_name'].nunique():,} ectodomains, "
          f"{d['Pfam'].nunique()} families -> {a.out}")
    n_missing = d["height"].isna().sum()
    if n_missing:
        fams = d.loc[d['height'].isna(), 'Pfam'].nunique()
        print(f"WARNING: {n_missing:,} assignments across {fams} families have no measured "
              f"clan height. They contribute nothing and can displace domains that do.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
