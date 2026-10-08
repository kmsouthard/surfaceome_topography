#!/usr/bin/env python
"""A run's tables → the data files the web explorer reads.

    python code/explorer/build_explorer.py RUN_DIR [--models DIR] [--dest explorer/data]

The explorer (``explorer/index.html``) is a static page: it computes nothing the pipeline
computes, and reads what this writes from a finished run (``python code/run_notebooks.py RUN_DIR``).

    meta.json        versions, the constants of the bridge and exclusion rules, the Fc receptors,
                     the antibodies and their epitopes, the proteins placed in a gap
    proteins.json    every protein of the height table: names, height, classification, its
                     ectodomain's stack of segments and which end the membrane holds
    pairs.json       the trans pairs and their gaps
    calls.json       every pair's call (trans, cis, secreted, pathway), read when a user uploads an
                     interaction table of their own
    residues/NN.json every ectodomain residue's height (`surfaceomeTopography.epitope`), in 64
                     files by a hash of the accession; the page fetches the one it needs
    surfaces.json    each cell type's surface: its proteins and their shares of the abundance (the
                     blood cells of E-PROT-1, not its tissues)
    contacts/*.json  each contact between two cell types: the units on both sides, and the pairs
                     with their gaps and shares, without an antibody and with each one that can
                     bridge it (`surfaceomeTopography.contact.contact_pairs`)

``--models`` is a directory of AlphaFold v6 models, read for the residue heights. Without it
only the curated antigens' residues are placed from their models (the committed
``alphafold_axis_positions.csv``); the rest take their share of the model's residues, flagged.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "code"))
from surfaceomeTopography import (cell_surface, contact_pairs, expressed_complexes, percent_expression,  # noqa: E402
                                  surface_units)
from surfaceomeTopography.epitope import (ANTIBODY_LENGTH, PHAGOCYTOSIS_RANGE, ectodomain_segments,  # noqa: E402
                                          probe_heights, residue_heights)
from surfaceomeTopography.interactions import with_curated_complexes  # noqa: E402
from surfaceomeTopography.interface import BULLSEYE_BANDS, _unit_heights  # noqa: E402
from surfaceomeTopography.topography import HEIGHT_BINS  # noqa: E402

DATA = REPO / "data"
#: The cell types a contact is offered for: one side from the first list, the other from both.
BLOOD = ["adult, B cell", "adult, CD4-positive T cell", "adult, CD8-positive T cell", "adult, monocyte",
         "adult, natural killer cell", "adult, platelet"]
#: E-PROT-27 is tumour tissue (40 breast tumours, Tyanova et al. 2016), and is named as such.
TUMOURS = {"HER2 Positive Breast Carcinoma": "HER2-positive breast tumour", "breast tumor luminal": "luminal breast tumour",
           "triple-negative breast cancer": "triple-negative breast tumour"}
MONOCYTE_SUBSETS = ["Classical monocytes", "Intermediate monocytes", "Non-classical monocytes"]
SHARDS = 64


def shard(accession: str) -> int:
    """The residue file an accession is in; `app.js` computes the same."""
    return sum(ord(c) for c in accession) % SHARDS


def write(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, separators=(",", ":"), allow_nan=False))


def clean(value):
    """A cell as JSON: NaN as None, numbers rounded to what the page shows."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    if isinstance(value, (float, np.floating)):
        return round(float(value), 3)
    if isinstance(value, (int, np.integer)):
        return int(value)
    return value


def slug(sample: str) -> str:
    return "".join(c if c.isalnum() else "-" for c in sample.lower()).strip("-")


def build_proteins(run: Path, models: Path | None, dest: Path) -> pd.DataFrame:
    estimates = pd.read_csv(run / "database/height_estimates.csv")
    assignments = pd.read_csv(run / "database/domain_disorder_assignments_alphafold.csv", index_col=0, low_memory=False)
    ecds = pd.read_csv(run / "database/UP000005640_surfaceome_largest_ecds_classified.csv")
    ecds["ID"] = ecds["ID link"] + "[" + ecds["start"].astype(int).astype(str) + ":" + ecds["end"].astype(int).astype(str) + "]"
    anchor = ecds.drop_duplicates("ID").set_index("ID")["ecd_domain_order"]
    axes = pd.read_csv(DATA / "inputs/alphafold/alphafold_axis_positions.csv")
    names = pd.concat([pd.read_csv(DATA / "inputs/mapping/uniprot_surfaceome_gene_name.tab", sep="\t"),
                       pd.read_csv(DATA / "inputs/mapping/uniprot_mapping_cleanup.tab", sep="\t")])
    names = names.drop_duplicates("ID link").set_index("ID link")
    glyco = pd.read_csv(run / "database/glycosites.csv").groupby("ID link")["glycosite"].apply(
        lambda s: sorted(set(int(x) for x in s.dropna())))

    proteins, shards = [], [dict() for _ in range(SHARDS)]
    for _, r in estimates.iterrows():
        acc, ecd = r["ID link_first"], r["ID"]
        start, end = (int(x) for x in ecd.split("[")[1].rstrip("]").split(":"))
        held = anchor.get(ecd)
        stack = ectodomain_segments(assignments, ecd, r["total_height"])
        notes = []
        if isinstance(held, str) and end > start:
            model = None if models is None else models / f"AF-{acc}-F1-model_v6.pdb"
            heights = residue_heights(ecd, assignments, r["total_height"], held,
                                      model if model is not None and model.exists() else None,
                                      axes[axes["accession"] == acc])
            notes = sorted(n for n in heights["note"].unique() if n)
            if len(heights):
                first = int(heights["residue"].min())
                by_residue = heights.set_index("residue")["height"].reindex(range(first, int(heights["residue"].max()) + 1))
                shards[shard(acc)][acc] = [first, [None if np.isnan(h) else int(round(h * 100)) for h in by_residue]]
        gene = names["Gene names  (primary )"].get(acc)
        proteins.append({
            "acc": acc, "gene": gene if isinstance(gene, str) else r["Entry name_first"].split("_")[0],
            "entry": r["Entry name_first"], "name": clean(names["Protein names"].get(acc)),
            "height": clean(r["total_height"]), "methods": clean(r["methods"]), "location": clean(r["uniprot_location"]),
            "tm": clean(r["tm_count"]), "surfy": clean(r["surfy_label"]), "cd": clean(r["cd_number"]),
            "family": clean(r["almen_main"]), "ecd": [start, end], "anchor": clean(held),
            "stack": [[k, int(a), int(b), clean(float(h))] for k, a, b, h in stack.itertuples(index=False)],
            "glycosites": [g for g in glyco.get(acc, []) if start <= g <= end], "notes": notes})
    write(dest / "proteins.json", proteins)
    for i, content in enumerate(shards):
        write(dest / "residues" / f"{i:02d}.json", content)
    print(f"  {len(proteins):,} proteins, {sum(len(s) for s in shards):,} with residue heights")
    return estimates


def interaction_classes(run: Path) -> dict:
    """CellphoneDB's class of each of its interactions (``Signaling by Notch``, ``Adhesion by ICAM``), by its id."""
    cp = pd.read_csv(run / "database/cellphonedb_tm_interaction_sizes.csv").dropna(subset=["classification"])
    return dict(zip(cp["id_cp_interaction"], cp["classification"]))


def build_pairs(run: Path, dest: Path) -> pd.DataFrame:
    interactions = pd.read_csv(run / "database/interaction_heights.csv")
    classes = interaction_classes(run)
    trans = interactions[interactions["_merge"] != "FcR"].drop_duplicates("Human interaction_id")
    write(dest / "pairs.json", [[r["Human ID link_prot1"], r["Human ID link_prot2"], clean(r["Human gene name_prot1"]),
                                 clean(r["Human gene name_prot2"]), round(float(r["interaction_dim"]), 2),
                                 "CellphoneDB" if isinstance(r["id_cp_interaction"], str) and not isinstance(r["_merge"], str)
                                 else {"string": "STRING", "curated": "curated", "complex": "complex"}.get(r["_merge"], "CellphoneDB"),
                                 classes.get(r["id_cp_interaction"])]
                                for _, r in trans.iterrows() if pd.notna(r["interaction_dim"])])
    print(f"  {len(trans):,} trans pairs")
    return interactions


def build_calls(interactions: pd.DataFrame, dest: Path) -> None:
    """Every pair's call, for annotating an uploaded interaction table: ``{"ACC1|ACC2": "t" | "c" | "s" | "p"}``.

    The calls of ``interaction_types.csv`` with the corrections of ``interaction_type_corrections.csv``;
    then every pair of the run's trans table is trans and every pair of ``cis_pairs.csv`` cis, as the
    pipeline has them. Accessions are in sorted order.
    """
    key = lambda a, b: "|".join(sorted((a, b)))
    types = pd.read_csv(DATA / "curated/interaction_types.csv")
    calls = {key(a, b): t[0] for a, b, t in zip(types["accession_1"], types["accession_2"], types["type"])}
    corrections = pd.read_csv(DATA / "curated/interaction_type_corrections.csv")
    for a, b, to in zip(corrections["Human ID link_prot1"], corrections["Human ID link_prot2"], corrections["corrected"]):
        if to == "removed":
            calls.pop(key(a, b), None)
        else:
            calls[key(a, b)] = to[0]
    cis = pd.read_csv(DATA / "curated/cis_pairs.csv")
    for a, b in zip(cis["Human ID link_prot1"], cis["Human ID link_prot2"]):
        calls[key(a, b)] = "c"
    trans = interactions[interactions["_merge"] != "FcR"]
    for a, b in zip(trans["Human ID link_prot1"], trans["Human ID link_prot2"]):
        calls[key(a, b)] = "t"
    write(dest / "calls.json", calls)
    counts = pd.Series(list(calls.values())).value_counts().to_dict()
    print(f"  {len(calls):,} pair calls: {counts}")


def expression_tables(run: Path):
    """The two Expression Atlas proteomes as long surface tables, as the Figure 6 notebook builds them."""
    names = pd.concat([pd.read_csv(DATA / "inputs/mapping/uniprot_surfaceome_gene_name.tab", sep="\t"),
                       pd.read_csv(DATA / "inputs/mapping/uniprot_mapping_cleanup.tab", sep="\t")])
    surfaceome = names.drop_duplicates(subset=["ID link"])
    heights = pd.read_csv(run / "database/height_estimates.csv").rename(
        columns={"ID link_first": "ID link", "Entry name_first": "Entry name"})
    heights = heights.drop(columns=["max_dim_first", "min_dim_first", "seq_len_first", "source_last"])
    ids = ["ID link", "Entry", "Entry name", "Status", "Protein names", "Gene names", "Organism", "Length",
           "Gene names  (primary )", "Gene ID", "Gene Name", "ID", "methods_<lambda>", "total_height"]
    complexes = with_curated_complexes(pd.read_csv(DATA / "inputs/cellphonedb/complex_curated.csv"),
                                       pd.read_csv(DATA / "curated/curated_complexes.csv"))
    out = {}
    for key, file in (("immune", "E-PROT-1-query-results.tsv.gz"), ("cancer", "E-PROT-27-query-results.tsv")):
        expression = pd.read_csv(DATA / "inputs/expression" / file, sep="\t", header=4)
        surf = pd.merge(surfaceome, expression, left_on="Gene names  (primary )", right_on="Gene Name", how="inner")
        surf = pd.merge(surf, heights, on=["ID link", "Entry name"])
        long = pd.melt(surf, id_vars=ids, value_vars=list(expression.columns[2:]), var_name="sample", value_name="expression")
        long["percent_expression"] = percent_expression(long)
        out[key] = (long, expressed_complexes(long, complexes))
    return out


def build_surfaces(run: Path, tables: dict, dest: Path) -> list:
    surfaces, listing = {}, []
    for key, dataset in (("immune", "Expression Atlas E-PROT-1 (Kim et al. 2014)"),
                         ("cancer", "Expression Atlas E-PROT-27 (breast tumours, Tyanova et al. 2016)")):
        long = tables[key][0].dropna(subset=["percent_expression"]).drop_duplicates(["sample", "ID link"])
        for sample, rows in long.groupby("sample", sort=False):
            if key == "immune" and sample not in BLOOD:      # the explorer shows cells, not E-PROT-1's tissues
                continue
            surfaces[slug(sample)] = [[a, round(float(p), 4)] for a, p in zip(rows["ID link"], rows["percent_expression"])
                                      if p > 0]
            listing.append({"id": slug(sample), "label": TUMOURS.get(sample, sample.replace("adult, ", "")), "dataset": dataset,
                            "proteins": len(surfaces[slug(sample)])})
    monocytes = pd.read_csv(run / "tables/ravenhill_2020_monocyte_expression_dataset.csv")
    for subset in MONOCYTE_SUBSETS:
        donors = [f"{subset}, donation {i}" for i in (1, 2, 3)]
        surface = cell_surface(monocytes, donors, columns=donors + ["total_height", "ID link"])
        surface = surface[surface["percent_expression"] > 0].drop_duplicates("ID link")
        surfaces[slug(subset)] = [[a, round(float(p), 4)] for a, p in zip(surface["ID link"], surface["percent_expression"])]
        listing.append({"id": slug(subset), "label": subset.lower(), "dataset": "Ravenhill et al. 2020 (surface proteomics)",
                        "proteins": len(surfaces[slug(subset)])})
    write(dest / "surfaces.json", surfaces)
    print(f"  {len(surfaces)} cell surfaces")
    return listing


def unit_table(surface: pd.DataFrame) -> pd.DataFrame:
    units = _unit_heights(surface)
    return units[["unit", "unit_name", "Entry", "unit_height", "percent_expression"]]


def build_contacts(run: Path, tables: dict, interactions: pd.DataFrame, bridges: pd.DataFrame, dest: Path) -> list:
    interactions = interactions.drop(columns=["ID link_prot1", "ID_prot1", "ID link_prot2", "ID_prot2"])
    antigens = bridges[(bridges["role"] == "antigen") & (bridges["included"] == "yes")]
    classes = interaction_classes(run)
    sides = [("immune", s) for s in BLOOD] + [("cancer", s) for s in tables["cancer"][0]["sample"].unique()]
    listing = []
    for sample_a in BLOOD:
        long_a, cplx_a = tables["immune"]
        a = surface_units(long_a, cplx_a, sample_a)
        for key_b, sample_b in sides:
            long_b, cplx_b = tables[key_b]
            expressed = set(long_b.loc[(long_b["sample"] == sample_b) & long_b["percent_expression"].notna(), "ID link"])
            variants, units = {}, {"a": {}, "b": {}}

            def add(side: str, surface: pd.DataFrame) -> None:
                for u in unit_table(surface).itertuples(index=False):
                    units[side].setdefault(u.unit, [u.unit_name, u.Entry, round(float(u.unit_height), 2),
                                                    round(float(u.percent_expression), 4)])

            def variant(b: pd.DataFrame, target: str | None) -> list | None:
                pairs = contact_pairs(a, b, interactions, antibody_target=target)
                if target is not None and not (pairs["_merge"] == "FcR").any():
                    return None
                add("a", a)
                add("b", b)
                unit_a = pairs["complex_name_immune"].fillna(pairs["Entry_immune"])
                unit_b = pairs["complex_name_cancer"].fillna(pairs["Entry_cancer"])
                # a bridge row keeps its sum and its antigen's height: the page sets the gap from the epitope
                bridge = (pairs["_merge"] == "FcR").to_numpy()
                raw = contact_pairs(a, b, interactions, antibody_target=target, epitope_heights=1e6)["interaction_dim"] - 1e6
                return [[ua, ub, round(float(r if f else g), 3), round(float(s), 4), int(f),
                         round(float(h), 3) if f else 0, classes.get(c)]
                        for ua, ub, g, r, s, f, h, c in zip(unit_a, unit_b, pairs["interaction_dim"], raw + pairs["total_height_prot1"],
                                                            pairs["interface_shared"], bridge, pairs["total_height_prot1"],
                                                            pairs["id_cp_interaction"])]

            variants[""] = variant(surface_units(long_b, cplx_b, sample_b), None)
            for g in antigens.itertuples():
                if g.accession in expressed:
                    found = variant(surface_units(long_b, cplx_b, sample_b, include=[g.entry_name]), g.entry_name)
                    if found:
                        variants[g.accession] = found
            index = {side: {u: i for i, u in enumerate(units[side])} for side in units}
            name = f"{slug(sample_a)}__{slug(sample_b)}"
            write(dest / "contacts" / f"{name}.json", {
                "a": list(units["a"].values()), "b": list(units["b"].values()),
                "variants": {k: [[index["a"][p[0]], index["b"][p[1]], *p[2:]] for p in v] for k, v in variants.items()}})
            listing.append({"id": name, "a": slug(sample_a), "b": slug(sample_b), "antigens": [k for k in variants if k]})
        print(f"  contacts of {sample_a}")
    return listing


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", type=Path, help="a finished run's output directory")
    ap.add_argument("--models", type=Path, help="directory of AlphaFold v6 models")
    ap.add_argument("--dest", type=Path, default=REPO / "explorer" / "data")
    a = ap.parse_args()

    estimates = build_proteins(a.run, a.models, a.dest)
    interactions = build_pairs(a.run, a.dest)
    build_calls(interactions, a.dest)
    tables = expression_tables(a.run)
    surfaces = build_surfaces(a.run, tables, a.dest)
    bridges = pd.read_csv(DATA / "curated/antibody_bridges.csv")
    contacts = build_contacts(a.run, tables, interactions, bridges, a.dest)

    assignments = pd.read_csv(a.run / "database/domain_disorder_assignments_alphafold.csv", index_col=0, low_memory=False)
    height = estimates.set_index("ID link_first")["total_height"]
    receptors = bridges[(bridges["role"] == "Fc receptor") & (bridges["included"] == "yes")]
    epitopes = pd.read_csv(a.run / "database/antibody_epitope_heights.csv")
    cited = pd.read_csv(DATA / "curated/antibody_epitopes.csv", dtype=str, keep_default_na=False).set_index("accession")
    manifest = json.loads((DATA / "snapshot/MANIFEST.json").read_text())["provenance"]
    write(a.dest / "meta.json", {
        "built_from": manifest.get("code_commit", "")[:12], "databases": manifest.get("databases", {}),
        "models_measured": a.models is not None,
        "antibody_length": ANTIBODY_LENGTH, "phagocytosis_range": PHAGOCYTOSIS_RANGE, "exclusion_margin": 5.0,
        "bands": list(BULLSEYE_BANDS), "height_bins": [round(float(b), 4) for b in HEIGHT_BINS], "shards": SHARDS,
        "receptors": [{"gene": r.gene, "acc": r.accession, "height": round(float(height[r.accession]), 3)}
                      for r in receptors.itertuples() if r.accession in height],
        "probes": [{"name": n, "height": round(float(h), 2)} for n, h in
                   probe_heights(pd.read_csv(DATA / "curated/segregation_probes.csv"), estimates, assignments).items()],
        "antibodies": [{"acc": r.accession, "gene": r.gene, "antibody": r.antibody, "pdb": clean(r.pdb),
                        "pmid": cited.loc[r.accession, "pmids"], "residues": clean(r.residues),
                        "height": clean(r.epitope_height), "lowest": clean(r.epitope_lowest),
                        "highest": clean(r.epitope_highest), "basis": r.basis, "notes": clean(r.notes)}
                       for r in epitopes.itertuples()],
        "surfaces": surfaces, "contacts": contacts})
    size = sum(f.stat().st_size for f in a.dest.rglob("*.json"))
    print(f"  wrote {a.dest} ({size / 1e6:.1f} MB in {sum(1 for _ in a.dest.rglob('*.json'))} files)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
