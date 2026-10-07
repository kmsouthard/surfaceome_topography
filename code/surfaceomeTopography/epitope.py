"""Where a residue stands above the membrane: the height of an epitope.

A protein's height is its ectodomain's segments laid end to end along the membrane normal
(01_height_estimates): an AlphaFold model, disordered chains, counted domains, and 0.04 nm for
each residue left over. A residue's height is the same stack read part-way up: every segment
between it and the membrane, plus how far it sits into its own segment.

    inside a model            its C-alpha's position along the inertia axis the model's height is
                              measured on, from the membrane end of the box
    inside a disordered chain  its share of the chain's residues, times the chain's height
    inside a counted domain    its share of the domain's residues, times the domain's height
    anywhere else              its share of the residues left over, times their height

The membrane end comes from the ectodomain's place in the topology: the first ectodomain of a
chain is held at its C-terminal end, the last at its N-terminal end. A loop between two
transmembrane helices is held at both; it is placed only inside a model, with the membrane at the
end of the box its two ends lie nearer, and is flagged.

The top residue of an ectodomain stands at the protein's height, by construction. The share rule
for disordered chains and domains is an assumption of this work, not a measurement.

`predict_bridge` carries an epitope to the contact an antibody on it would hold: the gap to each
Fc receptor (`interface.antibody_bridge_heights`), whether the epitope is within 10 nm of the
target membrane, and whether each probe protein (`probe_heights`: CD45, its short isoform and
CD148) fits in that gap or is excluded from it (the rule of `interface.contact_exclusion`). Heights are the
extended estimates, so gaps are upper bounds; a protein's own height is too, so its exclusion is
overstated where it bends.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .disorder_height import disorder_height
from .interface import antibody_bridge_heights

__all__ = ["ANTIBODY_LENGTH", "PHAGOCYTOSIS_RANGE", "ectodomain_segments", "model_axes",
           "residue_heights", "epitope_height", "antigen_epitope", "curated_epitope_heights",
           "probe_heights", "predict_bridge"]

#: The length an antibody adds between its antigen and an Fc receptor, nm (PDB 5DK3 and 4X4M;
#: 03_interaction_heights).
ANTIBODY_LENGTH = 3.24

#: How far above the target membrane an antibody can be held before phagocytosis falls, nm
#: (Bakalar et al. 2018).
PHAGOCYTOSIS_RANGE = 10.0

_MASS = {"C": 12.011, "N": 14.007, "O": 15.999, "S": 32.06, "SE": 78.971}
_SPANS = (("alphafold", "alphafold_start", "alphafold_end"),
          ("disorder", "disorder_start", "disorder_end"),
          ("domain", "domain_start", "domain_end"))


def ectodomain_segments(assignments: pd.DataFrame, ecd_id: str, total_height: float) -> pd.DataFrame:
    """One ectodomain's stack in sequence order: ``kind``, ``start``, ``end``, ``height``.

    ``assignments`` is ``domain_disorder_assignments_alphafold.csv``; ``ecd_id`` its ``ID``, like
    ``P04626[23:652]``. Residues no segment covers become ``sequence`` rows sharing whatever height
    the total holds beyond the segments' sum: 0.04 nm a residue where the pipeline imputed it,
    nothing where it did not.
    """
    start, end = (int(x) for x in ecd_id.split("[")[1].rstrip("]").split(":"))
    rows = []
    for _, r in assignments[assignments["ID"] == ecd_id].iterrows():
        for kind, a, b in _SPANS:
            if r["methods"] == kind and pd.notna(r[a]):
                rows.append((kind, int(r[a]), int(r[b]), float(r["height"])))
    rows.sort(key=lambda s: (s[1], s[2]))

    stack, at = [], start
    for kind, a, b, height in rows:
        a = max(a, at)                      # a segment starting inside the last one begins where it ends
        if b <= a:
            continue
        if a > at:
            stack.append(("sequence", at, a, np.nan))
        stack.append((kind, a, b, height))
        at = b
    if at < end:
        stack.append(("sequence", at, end, np.nan))
    stack = pd.DataFrame(stack, columns=["kind", "start", "end", "height"])

    left = stack["kind"] == "sequence"
    spare = max(total_height - stack.loc[~left, "height"].sum(), 0.0)
    n_left = (stack.loc[left, "end"] - stack.loc[left, "start"]).sum()
    stack.loc[left, "height"] = (stack["end"] - stack["start"]) * (spare / n_left if n_left else 0.0)
    return stack


def model_axes(model: Path, start: int, end: int) -> pd.DataFrame:
    """Each residue's C-alpha along a model's three box axes, nm from each one's low end.

    The axes and the box are those of `measure_structures.iabb`: the inertia axes of ``start``..
    ``end``, heavy atoms, mass weighted. One row per residue: ``residue``, ``axis_1`` to
    ``axis_3``, longest box dimension first, and the dimensions themselves, ``box_1`` to ``box_3``.
    """
    xyz, mass, ca = [], [], {}
    for line in Path(model).read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        resi, element = int(line[22:26]), line[76:78].strip().upper()
        if not start <= resi <= end or element == "H":
            continue
        if line[12:16].strip() == "CA":
            ca[resi] = len(xyz)
        xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        mass.append(_MASS.get(element, 12.011))
    xyz, mass = np.array(xyz), np.array(mass)
    xyz = xyz - np.average(xyz, axis=0, weights=mass)
    inertia = np.einsum("i,ij,ij->", mass, xyz, xyz) * np.eye(3) - np.einsum("i,ij,ik->jk", mass, xyz, xyz)
    with np.errstate(all="ignore"):         # Accelerate raises spurious matmul warnings on macOS
        boxed = xyz @ np.linalg.eigh(inertia)[1]
    boxed = (boxed - boxed.min(axis=0)) / 10
    boxed = boxed[:, np.argsort(-boxed.max(axis=0))]
    out = pd.DataFrame(boxed[list(ca.values())], columns=["axis_1", "axis_2", "axis_3"])
    out.insert(0, "residue", list(ca))
    for i in range(3):
        out[f"box_{i + 1}"] = boxed[:, i].max()
    return out.sort_values("residue", ignore_index=True)


def _height_axis(axes: pd.DataFrame, height: float) -> pd.Series:
    """The axis a model's height is measured on: the one whose box dimension is ``height``.

    A model's height is the box's longest dimension, except the plexins', which take its shortest
    (01_height_estimates). ``attrs["box"]`` is the dimension.
    """
    box = axes[["box_1", "box_2", "box_3"]].iloc[0].to_numpy()
    i = int(np.argmin(np.abs(box - height)))
    along = pd.Series(axes[f"axis_{i + 1}"].to_numpy(), index=axes["residue"].to_numpy())
    along.attrs["box"] = box[i]
    return along


def residue_heights(ecd_id: str, assignments: pd.DataFrame, total_height: float, anchor: str,
                    model: Path | None = None, axes: pd.DataFrame | None = None) -> pd.DataFrame:
    """The height of every residue of an ectodomain: ``residue``, ``kind``, ``height``, ``note``.

    ``anchor`` is the ectodomain's ``ecd_domain_order``: ``first`` (held at its C-terminal end),
    ``last`` (at its N-terminal end) or ``loop`` (both). The model's residues are placed from
    ``model``, the protein's AlphaFold file, or from ``axes``, `model_axes` of it stored with the
    span it was measured over (``start``, ``end``); with neither, or a stored span that is not
    this model's, a residue inside the model takes its share of the model's residues, and says so.
    """
    stack = ectodomain_segments(assignments, ecd_id, total_height)
    if anchor == "first":                   # membrane after the last residue: stack from the C-terminal end
        stack = stack.iloc[::-1]
    below, out = 0.0, []
    for kind, a, b, height in stack.itertuples(index=False):
        residues = np.arange(a, b + 1) if kind == "alphafold" else np.arange(a + 1, b + 1)
        share = (residues - a) / (b - a)
        if anchor == "first":
            share = 1 - share
        within, note = share * height, ""
        measured = None
        if kind == "alphafold" and model is not None:
            measured = model_axes(model, a, b)
        elif kind == "alphafold" and axes is not None and ((axes["start"] == a) & (axes["end"] == b)).all():
            measured = axes if len(axes) else None
        if measured is not None:
            along = _height_axis(measured, height)
            ends = along.reindex([a, b]).to_numpy()
            if anchor == "loop":
                flipped = ends.mean() > along.attrs["box"] / 2
                note = "loop: membrane at the end of the model its two ends lie nearer"
            else:
                membrane_end = ends[1] if anchor == "first" else ends[0]
                flipped = membrane_end > along.attrs["box"] / 2
                if 0.25 < membrane_end / along.attrs["box"] < 0.75:
                    note = "the model's membrane end lies mid-box: which way up is uncertain"
            placed = along.attrs["box"] - along if flipped else along
            # the stored height is the same box measured in PyMOL; the two agree to rounding
            within = placed.reindex(residues).to_numpy() * (height / along.attrs["box"])
        elif kind == "alphafold":
            note = "no model: placed by its share of the model's residues"
        elif anchor == "loop":
            within, note = np.full(len(residues), np.nan), "loop outside a model: not placed"
        out.append(pd.DataFrame({"residue": residues, "kind": kind, "height": below + within, "note": note}))
        below += height
    if not out:                             # an ectodomain of no residues
        return pd.DataFrame(columns=["residue", "kind", "height", "note"])
    return pd.concat(out).drop_duplicates("residue").sort_values("residue", ignore_index=True)


def epitope_height(residues, heights: pd.DataFrame) -> dict:
    """An epitope's height from `residue_heights`: the mean of its residues, with their range.

    Residues outside the ectodomain are returned in ``outside`` and left out of the mean.
    """
    residues = sorted(set(int(r) for r in residues))
    hit = heights[heights["residue"].isin(residues)]
    return {"height": hit["height"].mean(), "lowest": hit["height"].min(), "highest": hit["height"].max(),
            "kinds": sorted(hit["kind"].unique()), "notes": sorted(n for n in hit["note"].unique() if n),
            "outside": sorted(set(residues) - set(hit["residue"]))}


def antigen_epitope(accession: str, residues, estimates: pd.DataFrame, assignments: pd.DataFrame,
                    ecds: pd.DataFrame, model_dir: Path | None = None, axes: pd.DataFrame | None = None,
                    pattern: str = "AF-{acc}-F1-model_v6.pdb") -> dict:
    """An epitope's height on one protein of the height table, from its accession and residues.

    ``estimates`` is ``height_estimates.csv``, ``assignments``
    ``domain_disorder_assignments_alphafold.csv`` and ``ecds``
    ``UP000005640_surfaceome_largest_ecds_classified.csv``. The model is read from ``model_dir``,
    a directory of AlphaFold models, or failing that from ``axes``,
    ``alphafold_axis_positions.csv``. Returns `epitope_height`'s values with the
    ``antigen_height`` and the ectodomain (``ecd``) they were read on.
    """
    row = estimates[estimates["ID link_first"] == accession]
    if row.empty:
        raise KeyError(f"{accession} has no row in the height table")
    row = row.iloc[0]
    here = ecds[ecds["ID link"] == accession]
    span = here["start"].astype(int).astype(str) + ":" + here["end"].astype(int).astype(str)
    anchor = here.loc[accession + "[" + span + "]" == row["ID"], "ecd_domain_order"]
    if anchor.empty:
        raise KeyError(f"{row['ID']} has no topology row, so no membrane end")
    model = None if model_dir is None else Path(model_dir) / pattern.format(acc=accession)
    heights = residue_heights(row["ID"], assignments, row["total_height"], anchor.iloc[0],
                              model if model is not None and model.exists() else None,
                              None if axes is None else axes[axes["accession"] == accession])
    out = epitope_height(residues, heights)
    if not len(set(int(r) for r in residues)) > len(out["outside"]):
        raise ValueError(f"no residue given lies in the ectodomain {row['ID']}")
    return {**out, "antigen_height": float(row["total_height"]), "ecd": row["ID"]}


def curated_epitope_heights(epitopes: pd.DataFrame, estimates: pd.DataFrame, assignments: pd.DataFrame,
                            ecds: pd.DataFrame, axes: pd.DataFrame | None = None,
                            model_dir: Path | None = None) -> pd.DataFrame:
    """The height of each antibody's epitope in ``antibody_epitopes.csv``, one row per antibody.

    ``basis`` is the table's: a ``structure`` of the antibody on its antigen, ``peptide mapping``,
    or the ``domain`` it is reported to bind, the mean then taken over the whole stretch. An
    antibody with no epitope recorded is kept at ``epitope_height`` 0, a membrane-proximal
    epitope, the least gap its bridge can hold.
    """
    rows = []
    for _, r in epitopes.iterrows():
        row = {"accession": r["accession"], "gene": r["gene"], "antibody": r["antibody"], "pdb": r["pdb"],
               "residues": r["residues"], "epitope_height": 0.0, "basis": "none: taken at the membrane"}
        if r["included"] == "yes":
            residues = [n for part in r["residues"].split(",") for n in
                        range(int(part.split("-")[0]), int(part.split("-")[-1]) + 1)]
            found = antigen_epitope(r["accession"], residues, estimates, assignments, ecds, model_dir, axes)
            row.update(epitope_height=found["height"], epitope_lowest=found["lowest"], epitope_highest=found["highest"],
                       antigen_height=found["antigen_height"], basis=r["basis"], notes="; ".join(found["notes"]))
        rows.append(row)
    return pd.DataFrame(rows, columns=["accession", "gene", "antibody", "pdb", "residues", "epitope_height",
                                       "epitope_lowest", "epitope_highest", "antigen_height", "basis", "notes"])


def probe_heights(probes: pd.DataFrame, estimates: pd.DataFrame, assignments: pd.DataFrame) -> pd.Series:
    """The height of each protein a prediction places in the gap, by name.

    ``probes`` is ``segregation_probes.csv``. A row takes its protein's height from the height
    table. A row naming an ``isoform`` stands on the canonical protein's folded segments, with a
    disordered chain of its own ``disorder_residues`` and ``glycosites`` in place of the
    canonical one (`disorder_height.disorder_height`).
    """
    out = {}
    for _, r in probes[probes["included"] == "yes"].iterrows():
        row = estimates[estimates["ID link_first"] == r["accession"]].iloc[0]
        height = row["total_height"]
        if pd.notna(r["isoform"]) and r["isoform"] != "":
            stack = ectodomain_segments(assignments, row["ID"], height)
            n = float(r["disorder_residues"])
            height = (stack.loc[stack["kind"] != "disorder", "height"].sum()
                      + float(disorder_height(n, float(r["glycosites"]) / n)))
        out[r["name"]] = float(height)
    return pd.Series(out)


def _fit(height: float, gap: pd.Series, margin: float) -> np.ndarray:
    return np.select([height <= gap, height >= gap + margin], ["fits", "excluded"], "partly excluded")


def predict_bridge(epitope: dict, receptors: pd.DataFrame, proteins: pd.Series,
                   margin: float = 5.0) -> pd.DataFrame:
    """The gap an antibody on an epitope holds to each Fc receptor, and which proteins it excludes.

    ``epitope`` is `antigen_epitope`'s; ``receptors`` has a ``gene`` and a ``height`` per Fc
    receptor; ``proteins`` is the height of each protein to place, by name. One row per receptor:
    the ``gap`` at the epitope's mean height, ``gap_low`` and ``gap_high`` at its lowest and
    highest residue, ``held_by`` (``epitope``, or ``antigen`` where the antigen is taller than the
    bridge and sets the gap itself), ``within_phagocytosis_range`` (the epitope stands no more
    than `PHAGOCYTOSIS_RANGE` above its membrane), and for each protein ``fits``, ``partly
    excluded`` (up to ``margin`` nm taller than the gap) or ``excluded``.
    """
    antigen = epitope["antigen_height"]
    pairs = pd.DataFrame({"gene": receptors["gene"].to_numpy(), "receptor_height": receptors["height"].to_numpy(),
                          "_merge": "FcR", "total_height_prot1": antigen})
    pairs["interaction_dim"] = antigen + pairs["receptor_height"] + ANTIBODY_LENGTH
    out = pairs[["gene", "receptor_height"]].assign(epitope_height=epitope["height"], antigen_height=antigen)
    for column, at in (("gap", "height"), ("gap_low", "lowest"), ("gap_high", "highest")):
        out[column] = antibody_bridge_heights(pairs, epitope[at])["interaction_dim"]
    out["held_by"] = np.where(out["gap"] > antigen, "epitope", "antigen")
    out["within_phagocytosis_range"] = epitope["height"] <= PHAGOCYTOSIS_RANGE
    for name, height in proteins.items():
        out[name] = _fit(height, out["gap"], margin)
    return out
