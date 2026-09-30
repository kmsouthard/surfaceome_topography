#!/usr/bin/env python
"""Flatten a CellphoneDB v5 database bundle into the four tables the pipeline reads.

    python code/database/inputs/build_cellphonedb_tables.py --bundle cellphonedb.zip --out-dir DIR

``02_cellphonedb_interactions`` was written against the 2022-02-25 download, four flat files:
``interaction_curated.csv`` (one row per pair, partners named by UniProt accession or complex
name), ``protein_curated.csv`` (one row per protein with its membrane flags),
``complex_curated.csv`` (one row per complex with its subunits in ``uniprot_1`` to
``uniprot_4``) and ``hla_curated.csv`` (unused; the 2022 copy is a single malformed line).

CellphoneDB v5 ships one bundle, ``cellphonedb.zip`` in the ``cellphonedb-data`` repository,
holding eight normalised tables instead: ``multidata_table`` (every protein and complex with
its flags), ``protein_table``, ``gene_table``, ``complex_table``, ``complex_composition_table``
and ``interaction_table``. This script joins them back into the 2022 shape, so the notebook
runs unchanged against either snapshot. Three things differ in v5 and are handled here:

* **Non-protein ligands.** v5 lists hormones and small molecules as partners. They have no
  height and no accession, so only interactions with ``is_ppi = True`` are kept; the count
  dropped is reported.
* **Provenance.** v5 records ``curator``, ``directionality`` and ``classification`` per
  interaction. They are carried as extra columns, which the notebook ignores.
* **Protein names.** v5 names proteins by UniProt entry name (``ESR1_HUMAN``) in
  ``protein_table``; the 2022 file did the same, so ``protein_name`` is taken from there.

Every output row is traceable to the bundle: ``id_cp_interaction`` is CellphoneDB's own
stable id, and the bundle's version tag is written to ``VERSION.txt`` beside the tables.
"""

from __future__ import annotations

import argparse
import io
import sys
import zipfile
from pathlib import Path

import pandas as pd

FLAGS = ["transmembrane", "peripheral", "secreted", "secreted_desc", "secreted_highlight",
         "receptor", "receptor_desc", "integrin", "other", "other_desc"]


def _read(z: zipfile.ZipFile, name: str) -> pd.DataFrame:
    return pd.read_csv(io.BytesIO(z.read(name)), low_memory=False)


def flatten(bundle: Path) -> dict[str, pd.DataFrame]:
    with zipfile.ZipFile(bundle) as z:
        md = _read(z, "multidata_table.csv")
        prot = _read(z, "protein_table.csv")
        comp = _read(z, "complex_table.csv")
        comp_comp = _read(z, "complex_composition_table.csv")
        inter = _read(z, "interaction_table.csv")

    md = md.set_index("id_multidata")
    is_complex = md["is_complex"].astype(str).str.lower().eq("true")

    # --- proteins: one row per non-complex multidata, named by accession -----------------
    p = md[~is_complex].copy()
    p = p.join(prot.set_index("protein_multidata_id")[["protein_name", "tags", "tags_reason",
                                                        "tags_description"]])
    proteins = pd.DataFrame({"uniprot": p["name"], "protein_name": p["protein_name"]})
    for c in FLAGS:
        proteins[c] = p[c]
    for c in ("tags", "tags_description", "tags_reason"):
        proteins[c] = p[c]
    proteins["pfam"] = ""
    proteins = proteins.reset_index(drop=True)

    # --- complexes: subunits from the composition table, in the bundle's order ----------
    subunits = (comp_comp.merge(md[["name"]], left_on="protein_multidata_id", right_index=True)
                .groupby("complex_multidata_id")["name"].apply(list))
    c = md[is_complex].copy()
    c = c.join(comp.set_index("complex_multidata_id")[["pdb_id", "pdb_structure",
                                                        "stoichiometry", "comments_complex"]])
    complexes = pd.DataFrame({"complex_name": c["name"]})
    subs = c.index.map(lambda i: subunits.get(i, []))
    for k in range(4):
        complexes[f"uniprot_{k + 1}"] = [s[k] if len(s) > k else None for s in subs]
    over = [len(s) for s in subs if len(s) > 4]
    for col in FLAGS + ["pdb_id", "pdb_structure", "stoichiometry", "comments_complex"]:
        complexes[col] = c[col]
    complexes = complexes.reset_index(drop=True)

    # --- interactions: partners by multidata name; protein names for proteins only ------
    name_of = md["name"]
    pname_of = prot.set_index("protein_multidata_id")["protein_name"]
    inter = inter.copy()
    ppi = inter["is_ppi"].astype(str).str.lower().eq("true")
    dropped = int((~ppi).sum())
    inter = inter[ppi]
    interactions = pd.DataFrame({
        "id_cp_interaction": inter["id_cp_interaction"],
        "partner_a": inter["multidata_1_id"].map(name_of),
        "partner_b": inter["multidata_2_id"].map(name_of),
        "protein_name_a": inter["multidata_1_id"].map(pname_of),
        "protein_name_b": inter["multidata_2_id"].map(pname_of),
        "annotation_strategy": inter["annotation_strategy"],
        "source": inter["source"],
        "curator": inter["curator"],
        "directionality": inter["directionality"],
        "classification": inter["classification"],
    }).reset_index(drop=True)

    hla = pd.DataFrame(columns=["ensembl", "gene_name", "hgnc_symbol", "uniprot"])
    return {"protein_curated.csv": proteins, "complex_curated.csv": complexes,
            "interaction_curated.csv": interactions, "hla_curated.csv": hla,
            "_dropped_non_protein": dropped, "_complexes_over_four_subunits": over}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bundle", type=Path, required=True, help="cellphonedb.zip from cellphonedb-data")
    ap.add_argument("--version", default="v5.0.0", help="the bundle's release tag, recorded beside the tables")
    ap.add_argument("--out-dir", type=Path, required=True)
    a = ap.parse_args()

    tables = flatten(a.bundle)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    for name, df in tables.items():
        if name.startswith("_"):
            continue
        df.to_csv(a.out_dir / name, index=False)
        print(f"  {name:<26} {len(df):>6} rows")
    (a.out_dir / "VERSION.txt").write_text(
        f"CellphoneDB {a.version}, from cellphonedb-data/cellphonedb.zip; flattened by "
        f"build_cellphonedb_tables.py\n")
    print(f"  interactions with a non-protein partner dropped: {tables['_dropped_non_protein']}")
    if tables["_complexes_over_four_subunits"]:
        print(f"  complexes with more than four subunits (truncated to four, as the 2022 schema): "
              f"{len(tables['_complexes_over_four_subunits'])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
