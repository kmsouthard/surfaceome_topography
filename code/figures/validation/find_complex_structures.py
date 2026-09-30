#!/usr/bin/env python
"""Find PDB entries that contain both partners of an interaction, by UniProt accession.

    python code/figures/validation/find_complex_structures.py --pairs CSV --out CSV [--any-species]

``--pairs`` has ``uniprot_1`` and ``uniprot_2`` columns (a ``pair`` label is carried through).
For each pair the RCSB Search API is asked for entries with a polymer entity mapped to each
accession; an entry that has both is a candidate template for the bound complex. The output
lists, per pair, every such entry with its title, experimental method, resolution, release
date and the chain ids carrying each partner, so a template can be chosen and cited.

By default the accessions are human, so only structures of the human proteins match. Pass
``--any-species`` to also search the mouse orthologs' accessions given in ``uniprot_1_mouse``
and ``uniprot_2_mouse`` columns, when present.

This replaces the hand search behind the 2018 assembled models
(``data/README.md``), whose templates were mostly mouse complexes.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import pandas as pd
import requests

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
DATA = "https://data.rcsb.org/graphql"


def _acc_node(acc: str) -> dict:
    return {"type": "group", "logical_operator": "and", "nodes": [
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
            "operator": "in", "value": [acc]}},
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name",
            "operator": "exact_match", "value": "UniProt"}},
    ]}


def entries_with_both(a: str, b: str, session: requests.Session) -> list[str]:
    """PDB ids with a polymer entity for each accession (two entities, or one if a == b)."""
    if a == b:
        query = _acc_node(a)
    else:
        query = {"type": "group", "logical_operator": "and", "nodes": [_acc_node(a), _acc_node(b)]}
    body = {"query": query, "return_type": "entry",
            "request_options": {"paginate": {"start": 0, "rows": 500}}}
    r = session.post(SEARCH, json=body, timeout=120)
    if r.status_code == 204:
        return []
    r.raise_for_status()
    return [x["identifier"] for x in r.json().get("result_set", [])]


GQL = """
query($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    struct { title }
    exptl { method }
    rcsb_entry_info { resolution_combined }
    rcsb_accession_info { initial_release_date }
    rcsb_primary_citation { pdbx_database_id_PubMed }
    polymer_entities {
      rcsb_polymer_entity_container_identifiers { auth_asym_ids reference_sequence_identifiers { database_accession database_name } }
      rcsb_entity_source_organism { ncbi_taxonomy_id }
    }
  }
}
"""


def describe(ids: list[str], session: requests.Session) -> dict[str, dict]:
    out = {}
    for i in range(0, len(ids), 50):
        r = session.post(DATA, json={"query": GQL, "variables": {"ids": ids[i:i + 50]}}, timeout=120)
        r.raise_for_status()
        for e in r.json()["data"]["entries"] or []:
            chains = {}
            for pe in e["polymer_entities"] or []:
                ids_ = pe["rcsb_polymer_entity_container_identifiers"]
                for ref in ids_.get("reference_sequence_identifiers") or []:
                    if ref["database_name"] == "UniProt":
                        chains.setdefault(ref["database_accession"], []).extend(ids_["auth_asym_ids"] or [])
            res = e["rcsb_entry_info"]["resolution_combined"]
            out[e["rcsb_id"]] = {
                "title": e["struct"]["title"],
                "method": ";".join(m["method"] for m in e["exptl"]),
                "resolution": res[0] if res else None,
                "released": (e["rcsb_accession_info"] or {}).get("initial_release_date", "")[:10],
                "pmid": (e["rcsb_primary_citation"] or {}).get("pdbx_database_id_PubMed"),
                "chains": chains,
            }
        time.sleep(0.2)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--any-species", action="store_true")
    a = ap.parse_args()

    pairs = pd.read_csv(a.pairs)
    s = requests.Session()
    rows = []
    for p in pairs.itertuples():
        accs = [(p.uniprot_1, p.uniprot_2, "human")]
        if a.any_species and {"uniprot_1_mouse", "uniprot_2_mouse"} <= set(pairs.columns) \
                and isinstance(p.uniprot_1_mouse, str) and isinstance(p.uniprot_2_mouse, str):
            accs.append((p.uniprot_1_mouse, p.uniprot_2_mouse, "mouse"))
        for u1, u2, species in accs:
            try:
                ids = entries_with_both(u1, u2, s)
            except requests.RequestException as exc:
                print(f"  {getattr(p, 'pair', u1 + '-' + u2)}: search failed ({exc})", file=sys.stderr)
                continue
            info = describe(ids, s) if ids else {}
            for pdb, d in info.items():
                rows.append({"pair": getattr(p, "pair", f"{u1}-{u2}"), "species": species,
                             "uniprot_1": u1, "uniprot_2": u2, "pdb": pdb,
                             "chains_1": ",".join(d["chains"].get(u1, [])),
                             "chains_2": ",".join(d["chains"].get(u2, [])),
                             "method": d["method"], "resolution": d["resolution"],
                             "released": d["released"], "pmid": d["pmid"], "title": d["title"]})
            print(f"  {getattr(p, 'pair', u1 + '-' + u2)} ({species}): {len(ids)} entries", flush=True)
            time.sleep(0.3)
    out = pd.DataFrame(rows).sort_values(["pair", "species", "resolution"], na_position="last")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, index=False)
    print(f"wrote {a.out}: {len(out)} entries over {out['pair'].nunique() if len(out) else 0} pairs")
    return 0


if __name__ == "__main__":
    sys.exit(main())
