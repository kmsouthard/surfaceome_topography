#!/usr/bin/env python
"""The STRING-derived interaction table: every STRING pair between surfaceome proteins, labelled.

    python code/database/inputs/build_string_interactions.py --string-dir DIR --surfaceome CSV \\
        --types data/curated/interaction_types.csv [--review data/curated/string_v12_pair_review.csv] \\
        --out inputs/string/string_surfaceome_pairs.csv [--version 12.0] [--min-score 400]

The interaction database starts from STRING's protein links (``download_string.py`` fetches
the release): each link whose two proteins both map to surfaceome accessions is a candidate
pair, carrying its channel scores and ``combined_score``. Whether a pair is *trans* (binds
across two cells), *cis*, *secreted* or a *pathway* association is not something STRING says;
it is a call made here, per pair, and kept in ``data/curated/interaction_types.csv`` (the 2020
calls, made on STRING v11) with new calls in ``string_v12_pair_review.csv``. This script joins
the two: the STRING release supplies the pairs and scores, the decision tables the labels, so
refreshing STRING is a re-join, not a re-annotation.

Output columns follow the table the interaction notebook always read: ``Human ID link_prot1``,
``string_id_prot1``, ``Human ID link_prot2``, ``string_id_prot2``, the seven channel scores and
``combined_score_human``, ``Interaction Type (cis, trans, or secreted)``, ``Human interaction_id``
(the two accessions sorted), ``_merge`` (``string``), ``Human gene name_prot1/2`` (STRING's
preferred names). Every pair appears in both orders, as STRING lists it and as the notebooks
expect. A pair with several labels (the 2020 calls were made mostly per protein, so a pair can
be cis one way and trans the other) keeps one row per label and order; a pair with no call has
an empty label and is excluded from the trans set until it is called.

``--min-score`` keeps pairs at or above a combined score (default 400, STRING's medium-confidence
cutoff, which the interaction notebook applies as well; at 150, the 2020 table's floor, the table
is five times larger for pairs nothing reads). Changing the cutoff means rebuilding this table.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import pandas as pd

TYPE = "Interaction Type (cis, trans, or secreted)"
CHANNELS = ("neighborhood", "fusion", "cooccurence", "coexpression", "experimental", "database",
            "textmining", "combined_score")
ALIAS_SOURCES = ("UniProt_AC", "Ensembl_UniProt")


def string_to_uniprot(aliases: Path, universe: set[str]) -> dict[str, str]:
    """STRING id -> surfaceome accession. Restricted to the universe before choosing one per id,
    because an id carries several accessions (isoforms, secondary) and an arbitrary pick loses
    the join."""
    al = pd.read_csv(aliases, sep="\t", comment=None)
    al.columns = ["sid", "alias", "source"]
    al = al[al.source.isin(ALIAS_SOURCES) & al.alias.isin(universe)]
    al = al.sort_values(["sid", "source"]).drop_duplicates("sid")
    return dict(zip(al.sid, al.alias))


def surfaceome_links(links: Path, ids: dict[str, str], min_score: int) -> pd.DataFrame:
    rows = []
    with gzip.open(links, "rt") as fh:
        header = fh.readline().split()
        assert header[:2] == ["protein1", "protein2"] and header[2:] == list(CHANNELS), header
        for line in fh:
            p = line.split()
            if p[0] in ids and p[1] in ids and int(p[-1]) >= min_score:
                rows.append((p[0], p[1], *(int(x) for x in p[2:])))
    d = pd.DataFrame(rows, columns=["string_id_prot1", "string_id_prot2", *CHANNELS])
    d["Human ID link_prot1"] = d.string_id_prot1.map(ids)
    d["Human ID link_prot2"] = d.string_id_prot2.map(ids)
    # two STRING ids on one accession: keep the best-scored link per ordered accession pair
    d = (d.sort_values("combined_score", ascending=False)
         .drop_duplicates(["Human ID link_prot1", "Human ID link_prot2"]))
    return d[d["Human ID link_prot1"] != d["Human ID link_prot2"]]


def labels(types: Path, review: Path | None) -> pd.DataFrame:
    """(accession_1, accession_2, type) per call; a reviewed call overrides the 2020 one."""
    t = pd.read_csv(types, dtype=str)[["accession_1", "accession_2", "type"]]
    if review is not None and review.is_file():
        r = pd.read_csv(review, dtype=str, keep_default_na=False)
        r = r[r["call"].str.strip() != ""]
        if len(r):
            called = set(zip(r.accession_1, r.accession_2))
            t = t[[p not in called for p in zip(t.accession_1, t.accession_2)]]
            t = pd.concat([t, r[["accession_1", "accession_2"]].assign(type=r["call"].str.strip())])
    return t.drop_duplicates()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--string-dir", type=Path, required=True,
                    help="directory holding <taxon>.protein.{links.detailed,aliases,info}.v<version>.txt.gz")
    ap.add_argument("--surfaceome", type=Path, required=True,
                    help="a surfaceome table with an 'ID link' accession column (stage 0's <proteome>.csv)")
    ap.add_argument("--types", type=Path, required=True, help="data/curated/interaction_types.csv")
    ap.add_argument("--review", type=Path, default=None,
                    help="data/curated/string_v12_pair_review.csv; rows with a call override --types")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--version", default="12.0")
    ap.add_argument("--taxon", default="9606")
    ap.add_argument("--min-score", type=int, default=400)
    a = ap.parse_args()

    f = lambda kind: a.string_dir / f"{a.taxon}.{kind}.v{a.version}.txt.gz"
    universe = set(pd.read_csv(a.surfaceome, usecols=["ID link"])["ID link"].dropna().astype(str))
    ids = string_to_uniprot(f("protein.aliases"), universe)
    print(f"{len(ids)} STRING ids map to {len(set(ids.values()))} of {len(universe)} surfaceome accessions")

    d = surfaceome_links(f("protein.links.detailed"), ids, a.min_score)
    info = pd.read_csv(f("protein.info"), sep="\t", usecols=[0, 1])
    info.columns = ["sid", "name"]
    names = dict(zip(info.sid, info.name))
    d["Human gene name_prot1"] = d.string_id_prot1.map(names)
    d["Human gene name_prot2"] = d.string_id_prot2.map(names)
    d = d.rename(columns={c: f"{c}_human" for c in CHANNELS})

    lab = labels(a.types, a.review)
    key = pd.DataFrame({"accession_1": d[["Human ID link_prot1", "Human ID link_prot2"]].min(axis=1),
                        "accession_2": d[["Human ID link_prot1", "Human ID link_prot2"]].max(axis=1)},
                       index=d.index)
    d = pd.concat([d, key], axis=1).merge(lab, on=["accession_1", "accession_2"], how="left")
    d[TYPE] = d.pop("type").fillna("")
    d["Human interaction_id"] = d.accession_1 + "," + d.accession_2
    d["_merge"] = "string"
    d = d.drop(columns=["accession_1", "accession_2"])

    cols = ["Human ID link_prot1", "string_id_prot1", "Human ID link_prot2", "string_id_prot2",
            *(f"{c}_human" for c in CHANNELS), TYPE, "Human interaction_id", "_merge",
            "Human gene name_prot1", "Human gene name_prot2"]
    d = d[cols].sort_values(["Human interaction_id", "Human ID link_prot1", TYPE]).reset_index(drop=True)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(a.out, index=False)

    pairs = d.drop_duplicates("Human interaction_id")
    hc = d[d.combined_score_human >= 400]
    print(f"{len(pairs)} pairs at combined score >= {a.min_score} ({len(d)} rows, both orders) -> {a.out}")
    print(f"  labelled: {(pairs[TYPE] != '').sum()} pairs; unlabelled: {(pairs[TYPE] == '').sum()}")
    print("  at >= 400, by label:", hc.drop_duplicates(['Human interaction_id', TYPE])[TYPE]
          .replace('', '(no call)').value_counts().to_dict())
    return 0


if __name__ == "__main__":
    sys.exit(main())
