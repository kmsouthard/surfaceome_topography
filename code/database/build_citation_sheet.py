#!/usr/bin/env python
"""Gather every citation behind a pipeline decision into one sheet for review.

    python code/database/build_citation_sheet.py [--offline]

Reads the tables that record decisions, each with its PMIDs and what they are cited for, and writes
`data/curated/citation_sheet.csv`: one row per decision and PMID, with the assertion, a PubMed link, and the
paper's first author, year, journal and title. A decision with no PMID gets one row with an empty
PMID and its status. The `reviewed` and `review_notes` columns are for the reader; whatever is in
them is kept when the sheet is rebuilt, matched on source, decision, subject and PMID.

Titles come from `data/curated/citation_pubmed.csv`, a cache filled from PubMed E-utilities for any PMID not
yet in it (skipped with --offline). `data/curated/CITATIONS.md` is rewritten from the same rows.

Run it after changing any decision table, as part of the checks.
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]  # this file is code/database/build_citation_sheet.py
TABLES = REPO / "data/curated"
MODELS = REPO / "data/measurements/models_2018"
SHEET = TABLES / "citation_sheet.csv"
PUBMED = TABLES / "citation_pubmed.csv"
REGISTER = TABLES / "CITATIONS.md"
KEY = ["source", "decision", "subject", "pmid"]
REVIEW = ["reviewed", "review_notes"]
COLUMNS = ["source", "decision", "subject", "status", "assertion", "pmid", "link", "first_author", "year",
           "journal", "title", *REVIEW]


def _read(name: str) -> pd.DataFrame:
    return pd.read_csv(TABLES / name, dtype=str, keep_default_na=False)


def decisions() -> pd.DataFrame:
    """One row per recorded decision: source, decision, subject, included, assertion, pmids, status."""
    out = []

    t = _read("interaction_type_corrections.csv")
    out.append(pd.DataFrame({
        "source": "interaction_type_corrections.csv",
        "decision": ["pair removed" if c == "removed" else f"relabelled {l} -> {c}" for l, c in zip(t.label_2020, t.corrected)],
        "subject": t.gene_prot1 + "–" + t.gene_prot2, "included": ["no" if c == "removed" else "yes" for c in t.corrected],
        "assertion": t.evidence, "pmids": t.pmids}))

    t = _read("curated_trans_interactions.csv")
    out.append(pd.DataFrame({
        "source": "curated_trans_interactions.csv",
        "decision": ["trans pair added" if i == "yes" else "trans pair checked and left out" for i in t.included],
        "subject": t["Human gene name_prot1"] + "–" + t["Human gene name_prot2"], "included": t.included,
        "assertion": t.evidence, "pmids": t.pmids}))

    t = _read("antibody_bridges.csv")
    out.append(pd.DataFrame({
        "source": "antibody_bridges.csv",
        "decision": ["antibody bridge " + r + ("" if i == "yes" else ", left out") for r, i in zip(t.role, t.included)],
        "subject": [g + (f" ({a})" if a else "") for g, a in zip(t.gene, t.antibody)], "included": t.included,
        "assertion": t.evidence, "pmids": t.pmids}))

    t = _read("antibody_epitopes.csv")
    out.append(pd.DataFrame({
        "source": "antibody_epitopes.csv",
        "decision": [f"antibody epitope ({b})" if i == "yes" else "antibody epitope, none recorded" for i, b in zip(t.included, t.basis)],
        "subject": t.gene + " (" + t.antibody + ")", "included": t.included,
        "assertion": [e or n for e, n in zip(t.evidence, t.note)], "pmids": t.pmids,
        "status": ["" if i == "yes" else "no epitope known" for i in t.included]}))

    t = _read("curated_complexes.csv")
    out.append(pd.DataFrame({
        "source": "curated_complexes.csv",
        "decision": ["complex added" if i == "yes" else "complex considered and left out" for i in t.included],
        "subject": t.gene_1 + "–" + t.gene_2, "included": t.included,
        "assertion": [e or n for e, n in zip(t.evidence, t.note)], "pmids": t.pmids}))

    t = _read("segregation_probes.csv")
    out.append(pd.DataFrame({
        "source": "segregation_probes.csv", "decision": "protein placed in an antibody's gap",
        "subject": t.name, "included": t.included, "assertion": t.evidence, "pmids": t.pmids}))

    t = _read("interaction_height_rules.csv")
    out.append(pd.DataFrame({
        "source": "interaction_height_rules.csv", "decision": "height rule: " + t.rule,
        "subject": [g if p in ("", "any") else f"{g}–{p}" for g, p in zip(t.gene, t.partner_genes)], "included": "yes",
        "assertion": t.evidence + [f" ({h} nm)" if h else "" for h in t.height_nm], "pmids": t.pmids}))

    t = _read("cis_pairs.csv")
    t = t[t.pmids != ""]
    out.append(pd.DataFrame({
        "source": "cis_pairs.csv", "decision": "cis pair (" + t.subtype + ")",
        "subject": t.gene_prot1 + "–" + t.gene_prot2, "included": "yes",
        "assertion": [n or s for n, s in zip(t.note, t.source)], "pmids": t.pmids}))

    t = pd.read_csv(TABLES / "model_templates.csv", dtype=str, keep_default_na=False)
    genes = pd.read_csv(MODELS / "interaction_models.csv", dtype=str, keep_default_na=False).set_index("model").human_genes
    out.append(pd.DataFrame({
        "source": "model_templates.csv", "decision": "assembled model template PDB " + t.pdb_id,
        "subject": [f"{m} ({genes.get(m, '').replace(';', '–') or 'no human pair'})" for m in t.model], "included": "yes",
        "assertion": "PDB " + t.pdb_id + ": " + t.pdb_title + "; " + t.check + [f"; {n}" if n else "" for n in t.note],
        "pmids": t.pmid, "status": ["" if p else "PDB entry, no publication" for p in t.pmid]}))

    t = _read("method_citations.csv")
    out.append(pd.DataFrame({
        "source": "method_citations.csv", "decision": t.decision, "subject": t.subject, "included": "yes",
        "assertion": t.assertion + " [" + t.applied_in + "]", "pmids": t.pmids, "status": t.status}))

    d = pd.concat(out, ignore_index=True)
    d["status"] = d["status"].fillna("")
    cited = d.pmids.str.strip() != ""
    d.loc[(d.status == "") & cited, "status"] = "cited"
    d.loc[(d.status == "") & ~cited, "status"] = "needs citation"
    d.loc[(d.included == "no") & (d.status == "cited"), "status"] = "cited, left out"
    return d


def pubmed_details(pmids: set, offline: bool) -> pd.DataFrame:
    """Author, year, journal and title per PMID, from the cache, filling it from PubMed when online."""
    columns = ["pmid", "first_author", "year", "journal", "title"]
    cache = pd.read_csv(PUBMED, dtype=str, keep_default_na=False) if PUBMED.exists() else pd.DataFrame(columns=columns)
    missing = sorted(pmids - set(cache.pmid))
    if missing and not offline:
        rows = []
        for i in range(0, len(missing), 150):
            batch = missing[i:i + 150]
            url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id="
                   + ",".join(batch))
            result = json.load(urllib.request.urlopen(url))["result"]
            for p in batch:
                r = result.get(p, {})
                if "title" not in r:
                    raise SystemExit(f"PMID {p} is not in PubMed")
                rows.append({"pmid": p, "first_author": r["authors"][0]["name"] if r.get("authors") else "",
                             "year": r.get("pubdate", "")[:4], "journal": r.get("source", ""), "title": r["title"]})
        cache = pd.concat([cache, pd.DataFrame(rows)], ignore_index=True).sort_values("pmid", key=lambda s: s.astype(int))
        cache.to_csv(PUBMED, index=False)
    return cache


def write_register(sheet: pd.DataFrame) -> None:
    """`data/curated/CITATIONS.md`: counts per table and status, and every decision still without a citation."""
    per_decision = sheet.drop_duplicates(["source", "decision", "subject"])
    counts = per_decision.groupby(["source", "status"]).size().unstack(fill_value=0)
    open_items = per_decision[per_decision.status.isin(["needs citation", "database only", "PDB entry, no publication"])]
    lines = [
        "# Citation register",
        "",
        "Generated by `code/database/build_citation_sheet.py`; do not edit by hand. The full list, one row per",
        "decision and citation with what each citation is used to assert, is `data/curated/citation_sheet.csv`,",
        "with `reviewed` and `review_notes` columns for review. Each PMID's title was taken from PubMed.",
        "",
        f"{sheet.pmid.replace('', pd.NA).dropna().nunique()} distinct papers cited across "
        f"{len(per_decision)} decisions.",
        "",
        "| decision table | " + " | ".join(counts.columns) + " |",
        "|---|" + "---|" * len(counts.columns),
        *[f"| `{src}` | " + " | ".join(str(v) for v in row) + " |" for src, row in counts.iterrows()],
        "",
        "## Still needing a citation",
        "",
        "| decision | subject | status | what it rests on |",
        "|---|---|---|---|",
        *[f"| {r.decision} | {r.subject} | {r.status} | {r.assertion} |" for r in open_items.itertuples()],
        "",
    ]
    REGISTER.write_text("\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--offline", action="store_true", help="do not query PubMed for new PMIDs")
    args = parser.parse_args()

    d = decisions()
    d["pmid"] = d.pmids.str.split(";")
    d = d.explode("pmid")
    d["pmid"] = d.pmid.str.strip()
    papers = pubmed_details(set(d.pmid) - {""}, args.offline)
    d = d.merge(papers, on="pmid", how="left").fillna("")
    d["link"] = ["" if p == "" else f"https://pubmed.ncbi.nlm.nih.gov/{p}/" for p in d.pmid]

    if SHEET.exists():
        old = pd.read_csv(SHEET, dtype=str, keep_default_na=False)[KEY + REVIEW]
        d = d.merge(old, on=KEY, how="left")
    d[REVIEW] = d.reindex(columns=REVIEW).fillna("")
    sheet = d[COLUMNS]
    sheet.to_csv(SHEET, index=False)
    write_register(sheet)

    untitled = sheet[(sheet.pmid != "") & (sheet.title == "")]
    print(f"{len(sheet)} rows, {sheet.pmid.replace('', pd.NA).dropna().nunique()} distinct PMIDs; "
          f"status: {sheet.drop_duplicates(['source', 'decision', 'subject']).status.value_counts().to_dict()}")
    if len(untitled):
        print(f"{len(untitled)} rows have a PMID with no PubMed title (run without --offline)")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())