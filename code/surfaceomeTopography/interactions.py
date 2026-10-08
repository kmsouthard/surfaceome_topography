"""Tables that size interactions from the height estimates.

The CellphoneDB mapping notebook sizes every curated interaction and complex from a per-protein
table of heights, and 03_interaction_heights then joins those sizes to STRING's. Both start
from the same view of the height table -- each protein keyed by its UniProt accession and
carrying the id of its AlphaFold model -- which is built here once rather than inline in each.

Before this, the per-protein table was built inside 03_interaction_heights and saved for the
mapping notebook, which runs *first*. So every run sized its CellphoneDB interactions from the
stored 2022 copy of that table rather than from the heights the same run had just computed.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

TYPE = "Interaction Type (cis, trans, or secreted)"


def correct_interaction_types(annotation: pd.DataFrame, corrections: pd.DataFrame) -> pd.DataFrame:
    """The 2020 interaction annotation with the pairs in ``corrections`` relabelled or removed.

    The 2020 labels were set mostly per protein, so a pair listed in both orders can carry a
    different label each way round: TIGIT-PVR is cis one way and trans the other. A row whose
    label a correction contradicts is dropped where the pair also has a row with the right
    label, and relabelled where it has none, so each corrected pair ends with one label and no
    extra rows. A correction of ``removed`` drops every row of a pair with no evidence of binding.

    ``corrections`` has ``Human ID link_prot1``, ``Human ID link_prot2`` and ``corrected``,
    plus the evidence for each.
    """
    def pairs(df):
        return pd.Series([",".join(sorted(p)) for p in zip(df["Human ID link_prot1"], df["Human ID link_prot2"])],
                         index=df.index, dtype=object)

    label = dict(zip(pairs(corrections), corrections["corrected"]))
    pair = pairs(annotation)
    missing = set(label) - set(pair)
    if missing:
        # a corrected pair the current STRING release no longer lists: the correction is moot,
        # and the CellphoneDB side is still covered by ``drop_corrected_pairs``
        print(f"{len(missing)} corrected pair(s) not in this STRING release: {sorted(missing)}")
    corrected = pair.map(label)
    removed = corrected == "removed"
    wrong = corrected.notna() & ~removed & (annotation[TYPE] != corrected)
    has_right = pair.isin(pair[corrected.notna() & ~removed & ~wrong])
    out = annotation.copy()
    out.loc[wrong, TYPE] = corrected[wrong]
    return out[~(wrong & has_right) & ~removed]


def drop_corrected_pairs(interactions: pd.DataFrame, corrections: pd.DataFrame) -> pd.DataFrame:
    """Trans interactions without the pairs ``corrections`` says are not trans.

    `correct_interaction_types` fixes the 2020 STRING labels; CellphoneDB carries no labels, so a pair
    corrected to cis, or removed, is dropped from its rows here. Pairs are matched on their two
    accessions in either order.
    """
    not_trans = corrections[corrections["corrected"] != "trans"]
    drop = {frozenset(p) for p in zip(not_trans["Human ID link_prot1"], not_trans["Human ID link_prot2"])}
    pair = [frozenset(p) for p in zip(interactions["Human ID link_prot1"], interactions["Human ID link_prot2"])]
    return interactions[[p not in drop for p in pair]]


def with_curated_complexes(complexes: pd.DataFrame, curated: pd.DataFrame) -> pd.DataFrame:
    """CellphoneDB's complexes with the ones curated here added.

    ``curated`` is `curated_complexes.csv`: ``complex_name``, its subunits' accessions in
    ``uniprot_1`` and ``uniprot_2``, ``included`` and the citations. Rows with ``included`` =
    ``no`` are left out. A curated complex is a transmembrane complex like any other: it is one
    unit on a cell's surface, expressed when both subunits are, standing at its taller subunit.
    """
    used = curated[curated["included"] == "yes"]
    added = pd.DataFrame({"complex_name": used["complex_name"], "uniprot_1": used["uniprot_1"],
                          "uniprot_2": used["uniprot_2"], "transmembrane": True, "peripheral": False,
                          "secreted": False, "comments_complex": used["evidence"]})
    clash = set(added["complex_name"]) & set(complexes["complex_name"])
    if clash:
        raise ValueError(f"curated complexes already in CellphoneDB: {sorted(clash)}")
    return pd.concat([complexes, added], ignore_index=True)


def inherit_through_complexes(trans: pd.DataFrame, curated: pd.DataFrame) -> pd.DataFrame:
    """Trans pairs extended to the other chain of a curated complex.

    A pair recorded with one chain of a complex in `curated_complexes.csv` holds for the complex,
    so the same pair is added with each sibling chain: CD4 binds HLA-DRA in the table, and so
    binds HLA-DRB1. A proteome often detects one chain of a heterodimer and not the other, and
    the complex then never forms as a unit; without this the detected chain would have no pair.
    Added rows carry ``_merge`` = ``complex``; a pair already in the table is not added again.
    """
    used = curated[curated["included"] == "yes"]
    siblings = {}
    for r in used.itertuples():
        siblings.setdefault(r.uniprot_1, set()).add((r.uniprot_2, r.gene_2))
        siblings.setdefault(r.uniprot_2, set()).add((r.uniprot_1, r.gene_1))
    have = {frozenset(p) for p in zip(trans["Human ID link_prot1"], trans["Human ID link_prot2"])}
    added = []
    for side, other in (("prot1", "prot2"), ("prot2", "prot1")):
        for _, row in trans[trans[f"Human ID link_{side}"].isin(siblings)].iterrows():
            for accession, gene in sorted(siblings[row[f"Human ID link_{side}"]]):
                pair = frozenset((accession, row[f"Human ID link_{other}"]))
                if pair in have or len(pair) == 1:
                    continue
                have.add(pair)
                new = row.copy()
                new[f"Human ID link_{side}"], new[f"Human gene name_{side}"], new["_merge"] = accession, gene, "complex"
                if f"Entry name_{side}" in new:
                    new[f"Entry name_{side}"] = np.nan
                added.append(new)
    return pd.concat([trans, pd.DataFrame(added)], ignore_index=True) if added else trans


def antibody_bridge_pairs(bridges: pd.DataFrame) -> pd.DataFrame:
    """Antibody bridges: every antigen of a therapeutic antibody paired with every Fc-gamma receptor.

    ``bridges`` is `antibody_bridges.csv`: one row per antigen or receptor (``role``), with
    ``accession``, ``gene``, ``entry_name``, ``included`` and the citations. Rows with ``included`` =
    ``no`` are left out. The antigen is ``prot1`` and the receptor ``prot2``, as
    `interface.antibody_bridge_heights` expects; ``_merge`` is ``FcR``.
    """
    used = bridges[bridges["included"] == "yes"]
    antigens, receptors = used[used["role"] == "antigen"], used[used["role"] == "Fc receptor"]
    pairs = antigens.merge(receptors, how="cross", suffixes=("_prot1", "_prot2"))
    return pd.DataFrame({"Human ID link_prot1": pairs["accession_prot1"], "Human ID link_prot2": pairs["accession_prot2"],
                         "Human gene name_prot1": pairs["gene_prot1"], "Human gene name_prot2": pairs["gene_prot2"],
                         "Entry name_prot1": pairs["entry_name_prot1"], "Entry name_prot2": pairs["entry_name_prot2"],
                         "Interaction Type (cis, trans, or secreted)": "trans", "_merge": "FcR"})


def _swap_sides(df: pd.DataFrame) -> pd.DataFrame:
    """The same rows with every ``_prot1`` column exchanged for its ``_prot2`` column."""
    swap = {c: c[:-1] + ("2" if c.endswith("1") else "1") for c in df.columns if c.endswith(("_prot1", "_prot2"))}
    return df.rename(columns=swap)[df.columns]


def combine_interaction_sources(string_rows: pd.DataFrame, cellphonedb_rows: pd.DataFrame) -> pd.DataFrame:
    """STRING and CellphoneDB trans interactions in one table, a pair found in both on one row.

    A pair is matched on its two accessions in either order. The two sources often list a pair in
    opposite orders, and before this they were merged on their order and on their computed
    heights, so 70 of the 170 pairs in both appeared twice. A CellphoneDB row keeps its order,
    because ``partner_a`` and ``partner_b`` name the complex on each side; the STRING row's
    ``_prot1``/``_prot2`` columns are swapped to match. CellphoneDB can list one pair several times,
    once per complex, and each of those rows carries the STRING annotation.

    Both sources are sized from the same height table, so a matched pair must have the same
    heights; a disagreement raises ``ValueError``. Returns the STRING columns followed by
    CellphoneDB's own.
    """
    heights = ["total_height_prot1", "total_height_prot2", "interaction_dim"]
    ids = ["Human ID link_prot1", "Human ID link_prot2"]
    extra = [c for c in cellphonedb_rows.columns if c not in string_rows.columns]
    both_orders = pd.concat([string_rows, _swap_sides(string_rows)], ignore_index=True).drop_duplicates(ids)
    matched = pd.merge(cellphonedb_rows, both_orders, on=ids, how="left", suffixes=("", "_string"), indicator="_source")
    found = matched["_source"] == "both"
    disagree = found & ~np.isclose(matched[heights].to_numpy(dtype=float),
                                   matched[[h + "_string" for h in heights]].to_numpy(dtype=float)).all(axis=1)
    if disagree.any():
        raise ValueError(f"heights differ between STRING and CellphoneDB for {sorted(matched.loc[disagree, 'Human interaction_id'])}")
    string_only = string_rows[~string_rows["Human interaction_id"].isin(matched.loc[found, "Human interaction_id"])]
    return pd.concat([string_only, matched[string_rows.columns.tolist() + extra]], ignore_index=True)


#: Rules that set an interaction's height other than as the sum of its two partners.
HEIGHT_RULES = {
    "extended": "the protein stands at ``height_nm`` in every interaction, its height when bound, "
                "in place of its predicted height",
    "side by side": "the protein binds its partner's membrane-proximal domains, so the two stand side by side "
                    "and the gap is the taller of the two",
    "model length": "the pair's height is ``height_nm``, the longest axis of an assembled model of the bound pair",
}


def apply_height_rules(interactions: pd.DataFrame, rules: pd.DataFrame) -> pd.DataFrame:
    """Interaction heights with the recorded rules applied (``HEIGHT_RULES``).

    ``rules`` has ``rule``, ``protein`` (an accession), ``partners`` (accessions separated by ``;``)
    and ``height_nm``, plus the evidence for each. *Extended* sets the protein's height, and moves
    the interaction's height by the change, so anything added to the sum -- an antibody's length --
    is kept. *Side by side*, then *model length*, apply to an interaction between the protein and
    one of its partners, in either order; a measured model outranks every other rule. A rule not in
    ``HEIGHT_RULES`` raises ``ValueError``.
    """
    unknown = set(rules["rule"]) - set(HEIGHT_RULES)
    if unknown:
        raise ValueError(f"unknown height rules: {sorted(unknown)}")
    out = interactions.copy()
    extended = rules[rules["rule"] == "extended"]
    bound_height = dict(zip(extended["protein"], extended["height_nm"].astype(float)))
    for side in ("prot1", "prot2"):
        new = out[f"Human ID link_{side}"].map(bound_height).fillna(out[f"total_height_{side}"])
        out["interaction_dim"] += new - out[f"total_height_{side}"]
        out[f"total_height_{side}"] = new

    pair = pd.Series(list(zip(out["Human ID link_prot1"], out["Human ID link_prot2"])), index=out.index)
    beside = pair.isin(_rule_pairs(rules[rules["rule"] == "side by side"]))
    out["interaction_dim"] = out["interaction_dim"].where(~beside, out[["total_height_prot1", "total_height_prot2"]].max(axis=1))

    model_length = _rule_pairs(rules[rules["rule"] == "model length"])
    measured = pair.map(model_length)
    out["interaction_dim"] = measured.fillna(out["interaction_dim"])
    return out


def _rule_pairs(rules: pd.DataFrame) -> dict:
    """Each (protein, partner) pair a set of rules names, in both orders, with the rule's ``height_nm``."""
    pairs = {(r.protein, partner): r.height_nm for r in rules.itertuples() for partner in r.partners.split(";")}
    pairs.update({(b, a): h for (a, b), h in pairs.items()})
    return {p: pd.to_numeric(h, errors="coerce") for p, h in pairs.items()}


def estimates_with_af_ids(estimates: pd.DataFrame, af_mappings: pd.DataFrame) -> pd.DataFrame:
    """The height table keyed by accession, with each protein's AlphaFold model id attached.

    ``estimates`` is the height table as written by the height notebook. ``af_mappings`` is
    ``surfaceome_current.tsv`` as read: ``current``, a comma-separated list of accessions that
    share one model, and ``af_id``. The join keeps the height table's rows: a model with no
    height row would only add an empty row, and break anything that expects a height.
    """
    heights = estimates.rename(columns={"ID link": "ID link structure",
                                        "ID link_first": "ID link"})
    af = af_mappings.copy()
    af["current"] = af.apply(lambda row: row.current.split(","), axis=1)
    af = af.explode("current")
    return pd.merge(heights, af, left_on="ID link", right_on="current", how="left")


def protein_heights_for_interactions(estimates: pd.DataFrame, af_mappings: pd.DataFrame,
                                     names: pd.DataFrame) -> pd.DataFrame:
    """Per-protein heights for sizing interactions, with names attached.

    ``names`` is the accession-to-name mapping -- ``uniprot_surfaceome_gene_name.tab`` and its
    cleanup file, concatenated. One row per (ectodomain, height).

    Two things the original inline version did are deliberately not reproduced:

    * it filled a missing ``Entry name`` with "". The mapping notebook read the table back from
      CSV, which turns "" into NaN again, so NaN is what it has always seen, and is kept.
    * it tried to set every integrin to its extended 20 nm conformation, but never assigned the
      result, so no integrin was changed -- its own check printed ``modification check: False``.
      Applying it would change 26 proteins' interaction heights, which is a decision to make on
      purpose rather than as a side effect of moving code.
    """
    heights = estimates_with_af_ids(estimates, af_mappings)
    named = pd.merge(heights, names[["ID link", "Entry name", "Entry", "Gene names"]],
                     on="ID link", how="left", indicator=True)
    return named.drop_duplicates(subset=["ID", "total_height"])
