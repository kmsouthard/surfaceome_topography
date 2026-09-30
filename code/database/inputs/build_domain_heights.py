#!/usr/bin/env python
"""Per-family domain heights from instance measurements, in the 2020 table's shape.

    python code/database/inputs/build_domain_heights.py --measurements CSV --orientations CSV \\
        --hmm-dat Pfam-A.hmm.dat --out CSV [--min-instances 1] [--max-extent-per-residue 4.0]

Rebuilds ``human_domain_sizes_oriented_clans_20200123.csv``, the table ``build_clan_heights.py``
derives every clan height from, using 2026 measurements (``measure_domain_instances.py``).
Each family gets the mean, standard deviation, standard error and count of its instances'
longest and shortest box dimensions, and an ``orientation`` call saying which one is the
domain's height when it stacks in an ectodomain: ``max`` if it stands along its longest axis,
``min`` if it lies along it. ``height`` is the chosen mean.

**Orientation** is the one judgement in the chain (METHODS.md §2), taken in this order:

1. ``--orientations``: a table of hand calls, one row per family (``Pfam``, ``orientation``,
   ``basis``). The 2020 calls are seeded into it; new families are added as they are called.
   Rows whose basis is a stacking call from an earlier run are dropped and recomputed.
2. Stacking, for a family whose copies occur in tandem on one chain (each covering at least
   80% of the model, so split HMMER fragments of one domain do not count): the vector between
   the centroids of consecutive copies is compared with each copy's longest and shortest
   inertia axes. If it aligns with the longest axis in at least 60% of pairs the family stacks end
   to end (``max``); in at most 40%, it lies across the stack (``min``); between, no call.
   The majority rule reproduced 79% of the 2020 hand calls where it could be tested.
3. Clan majority, for a family in a Pfam clan whose called members (rules 1 and 2) agree on
   an orientation in at least 60% of cases: members of one clan share a fold and so the
   geometry of chain entry and exit, and the rule reproduces 95% of the 2020 hand calls when
   each is predicted from its clan-mates alone.
4. Otherwise no call: the family carries no height of its own and inherits its clan's, which
   is how the 2020 table treated its 88 uncalled families. These are listed in
   ``<out>.needs_call.csv`` with their dimensions, for a call by hand.

Calls from rules 2 and 3 are written back to ``--orientations`` with their basis and
recomputed on every run; only hand calls are carried over.

**Stack pitch.** For a family with runs of three or more whole copies on one chain whose
centroids lie on a line, ``stack_pitch`` is the median advance per copy along that line: the
height one copy adds to a string, tilt and overlap included. It is recorded beside the box
height as a check (it runs 10-30% under the longest box dimension for bead-like folds) and
becomes the height only with ``--height-from-stack-pitch``, which is off so every family's
height is the same quantity, a box dimension.

**Neighbour pitch** (``--height-from-neighbour-pitch``, the adopted rule; see
``neighbour_pitch``): every domain on every measured chain, whatever its family, gives the
domains beside it a stacking axis, and a family's height is the median advance of its
instances along that axis. The box dimension along the called axis remains only for families
never seen beside a neighbour. **Composite** models spanning several domains of other
families carry no height at all (``composite_families``), so the domains they cover are the
ones counted.

**Sanity filter.** An instance whose longest dimension exceeds ``--max-extent-per-residue``
angstrom per residue (a fully extended chain is 3.8) is not a folded domain and is dropped;
the count dropped is reported. This is the filter the 2020 run lacked when whole chains were
measured as Sushi domains.

Columns match the 2020 table so ``build_clan_heights.py`` runs unchanged: ``Pfam``,
``max_dim_mean/std/sem/len``, ``min_dim_mean/std/sem/len``, ``query_name_x`` (the family
name), ``orientation``, ``Clan_x`` (the clan), ``height``, plus ``orientation_basis``,
``n_structures``, ``median_n_res``, ``height_basis``, ``stack_pitch`` and ``stack_runs`` for
the record.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def pfam_metadata(hmm_dat: Path) -> pd.DataFrame:
    """Pfam, name, clan for every family in Pfam-A.hmm.dat."""
    rows, cur = [], {}
    with open(hmm_dat) as fh:
        for line in fh:
            if line.startswith("# STOCKHOLM"):
                cur = {}
            elif line.startswith("#=GF ID"):
                cur["query_name_x"] = line.split(None, 2)[2].strip()
            elif line.startswith("#=GF AC"):
                cur["Pfam"] = line.split(None, 2)[2].strip().split(".")[0]
            elif line.startswith("#=GF CL"):
                cur["Clan_x"] = line.split(None, 2)[2].strip()
            elif line.startswith("//"):
                rows.append(cur)
    return pd.DataFrame(rows)


def stacking_calls(m: pd.DataFrame, min_coverage: float = 0.8) -> pd.DataFrame:
    """Orientation from consecutive same-family copies on one chain (METHODS.md §2).

    Only copies covering at least ``min_coverage`` of the family's model count: HMMER splits
    some large domains into fragments that would otherwise pass as tandem copies.
    """
    calls = []
    m = m[m["n_res"] >= min_coverage * m["model_len"]]
    m = m.sort_values(["PFAM_ID", "PDB", "CHAIN", "seq_from"])
    for fam, g in m.groupby("PFAM_ID"):
        along_long = along_short = 0
        for (_, _), c in g.groupby(["PDB", "CHAIN"]):
            if len(c) < 2:
                continue
            c = c.reset_index(drop=True)
            for i in range(len(c) - 1):
                v = c.loc[i + 1, ["cx", "cy", "cz"]].values - c.loc[i, ["cx", "cy", "cz"]].values
                n = np.linalg.norm(v)
                if n == 0:
                    continue
                v = v / n
                for j in (i, i + 1):
                    lo = abs(np.dot(v, c.loc[j, ["long_x", "long_y", "long_z"]].values))
                    sh = abs(np.dot(v, c.loc[j, ["short_x", "short_y", "short_z"]].values))
                    mid = np.sqrt(max(0.0, 1 - lo * lo - sh * sh))
                    if lo >= max(sh, mid):
                        along_long += 1
                    else:
                        along_short += 1
        n = along_long + along_short
        if n:
            frac = along_long / n
            # a clear majority either way; an even split is left for a call by hand
            call = "max" if frac >= 0.6 else ("min" if frac <= 0.4 else None)
            if call:
                calls.append((fam, call, n, frac))
    return pd.DataFrame(calls, columns=["Pfam", "orientation", "n_pairs", "frac_along_long"])


def stack_pitch(m: pd.DataFrame, min_coverage: float = 0.8, gap: int = 40,
                min_copies: int = 3, min_linearity: float = 0.9) -> pd.DataFrame:
    """Median advance per copy along runs of stacked whole copies, per family.

    A run is consecutive copies on one chain separated by at most ``gap`` residues; it counts
    when it has ``min_copies`` copies, their centroids fall on a line (first principal
    component carrying ``min_linearity`` of their spread) and advance one way along it.
    """
    m = m[m["n_res"] >= min_coverage * m["model_len"]]
    m = m.sort_values(["PFAM_ID", "PDB", "CHAIN", "seq_from"])
    rows = []
    for (fam, _, _), c in m.groupby(["PFAM_ID", "PDB", "CHAIN"]):
        if len(c) < min_copies:
            continue
        c = c.reset_index(drop=True)
        runs = [[0]]
        for i in range(1, len(c)):
            if c.loc[i, "seq_from"] - c.loc[i - 1, "seq_to"] > gap:
                runs.append([i])
            else:
                runs[-1].append(i)
        for r in runs:
            if len(r) < min_copies:
                continue
            C = c.loc[r, ["cx", "cy", "cz"]].values.astype(float)
            _, sv, vt = np.linalg.svd(C - C.mean(0))
            if sv[0] ** 2 / (sv ** 2).sum() < min_linearity:
                continue
            steps = np.diff(C @ vt[0])
            if not (np.all(steps > 0) or np.all(steps < 0)):
                continue
            rows.append((fam, np.abs(steps).mean(), len(r)))
    t = pd.DataFrame(rows, columns=["Pfam", "pitch", "copies"])
    if t.empty:
        return pd.DataFrame(columns=["Pfam", "stack_pitch", "stack_runs"])
    return t.groupby("Pfam").agg(stack_pitch=("pitch", "median"),
                                 stack_runs=("pitch", "size")).reset_index()


def neighbour_pitch(m: pd.DataFrame, min_coverage: float = 0.8, gap: int = 40,
                    overlap: int = 4) -> pd.DataFrame:
    """Per-family median of each domain's centre-to-centre spacing to its sequence neighbours.

    Each chain is tiled with its best-scoring non-overlapping whole domains (at most
    ``overlap`` shared residues, copies covering ``min_coverage`` of their model), the same
    resolution the pipeline applies to an ectodomain's assignments. A tile's neighbours are
    the adjacent tiles within ``gap`` residues, of any family. The axis through its
    neighbours (or to its one neighbour) is the local stacking direction. A tile between two
    neighbours advances the chain by half its spacing to each; a tile with one neighbour by
    half that spacing plus half its own extent along the axis (the box projected as an
    ellipsoid), so a small domain beside a large one is not credited with the large one's
    size. A domain with no neighbour within reach has no pitch.
    """
    m = m[m["n_res"] >= min_coverage * m["model_len"]]
    rows = []
    for (pdb, chain), c in m.groupby(["PDB", "CHAIN"]):
        c = c.sort_values("score", ascending=False)
        tiles = []
        for h in c.itertuples():
            if all(min(h.seq_to, t.seq_to) - max(h.seq_from, t.seq_from) + 1 <= overlap for t in tiles):
                tiles.append(h)
        tiles.sort(key=lambda t: t.seq_from)
        cen = np.array([[t.cx, t.cy, t.cz] for t in tiles], dtype=float)

        def own_extent(t, axis):
            e1 = np.array([t.long_x, t.long_y, t.long_z], dtype=float)
            e3 = np.array([t.short_x, t.short_y, t.short_z], dtype=float)
            e2 = np.cross(e1, e3)
            return 2 * np.sqrt(sum((d / 2 * (axis @ e)) ** 2
                                   for d, e in ((t.max_dim, e1), (t.mid_dim, e2), (t.min_dim, e3))))

        for i, t in enumerate(tiles):
            nb = []
            if i > 0 and t.seq_from - tiles[i - 1].seq_to <= gap:
                nb.append(i - 1)
            if i + 1 < len(tiles) and tiles[i + 1].seq_from - t.seq_to <= gap:
                nb.append(i + 1)
            if not nb:
                continue
            axis = cen[nb[-1]] - cen[nb[0]] if len(nb) == 2 else cen[nb[0]] - cen[i]
            n = np.linalg.norm(axis)
            if n == 0:
                continue
            axis /= n
            halves = [abs((cen[j] - cen[i]) @ axis) / 2 for j in nb]
            if len(nb) == 1:
                halves.append(own_extent(t, axis) / 2)
            rows.append((t.PFAM_ID, pdb, chain, len(nb), sum(halves)))
    t = pd.DataFrame(rows, columns=["Pfam", "PDB", "CHAIN", "n_neighbours", "pitch"])
    if t.empty:
        return pd.DataFrame(columns=["Pfam", "neighbour_pitch", "pitch_n", "pitch_structures", "pitch_interior", "inner_domains"])
    return t.groupby("Pfam").agg(neighbour_pitch=("pitch", "median"), pitch_n=("pitch", "size"),
                                 pitch_structures=("PDB", "nunique"),
                                 pitch_interior=("n_neighbours", lambda s: int((s == 2).sum()))).reset_index()


def composite_families(m: pd.DataFrame, min_coverage: float = 0.8, min_inner: int = 2,
                       inside: float = 0.8, max_inner_share: float = 0.6) -> pd.DataFrame:
    """Families whose model spans several domains of other families.

    Pfam has models covering a whole array of repeats (``EGF_TEN`` over six to eight EGF
    domains). Counted as one domain with one height they undercount badly, and because they
    score higher than the repeats they cover, the pipeline's overlap resolution would prefer
    them. A family is composite when, over its instances, the median number of whole-model
    hits of other families lying ``inside`` it, side by side (not overlapping each
    other, each at most ``max_inner_share`` of its span, so alternative models of the same
    region do not count), is at least ``min_inner``; it then carries no height, so the domains
    it spans are the ones counted.
    """
    whole = (m["n_res"] >= min_coverage * m["model_len"]).values
    counts = []
    for _, c in m.assign(whole=whole).groupby(["PDB", "CHAIN"]):
        c = c[["PFAM_ID", "seq_from", "seq_to", "whole"]].values
        for fam, a, b, _ in c:          # the outer hit at any coverage: a partly resolved array is still an array
            inner = []
            for fam2, a2, b2, w2 in c:
                if fam2 == fam or not w2:
                    continue
                ov = min(b, b2) - max(a, a2) + 1
                if ov >= inside * (b2 - a2 + 1) and (b2 - a2 + 1) <= max_inner_share * (b - a + 1):
                    inner.append((a2, b2))
            kept = []
            for a2, b2 in sorted(inner, key=lambda x: x[1] - x[0]):
                if all(min(b2, k2) - max(a2, k1) + 1 <= 4 for k1, k2 in kept):
                    kept.append((a2, b2))
            counts.append((fam, len(kept)))
    t = pd.DataFrame(counts, columns=["Pfam", "inner"]).groupby("Pfam").inner.median()
    comp = t[t >= min_inner]
    return pd.DataFrame({"Pfam": comp.index, "inner_domains": comp.values})


def clan_majority_calls(calls: pd.DataFrame, meta: pd.DataFrame, families: pd.Series,
                        min_agreement: float = 0.6) -> pd.DataFrame:
    """Orientation for uncalled families from the calls of their clan-mates (rule 3)."""
    voters = calls.merge(meta[["Pfam", "Clan_x"]], on="Pfam", how="left").dropna(subset=["Clan_x"])
    out = []
    todo = meta[meta["Pfam"].isin(families) & ~meta["Pfam"].isin(calls["Pfam"])].dropna(subset=["Clan_x"])
    for fam, clan in zip(todo["Pfam"], todo["Clan_x"]):
        votes = voters.loc[voters["Clan_x"] == clan, "orientation"].value_counts()
        if votes.empty:
            continue
        top, n = votes.idxmax(), int(votes.sum())
        if votes[top] / n >= min_agreement:
            out.append((fam, top, f"clan majority: {int(votes[top])} of {n} called {clan} members are {top}"))
    return pd.DataFrame(out, columns=["Pfam", "orientation", "basis"])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--measurements", type=Path, required=True)
    ap.add_argument("--orientations", type=Path, required=True,
                    help="hand calls: Pfam, orientation (max|min), basis; read, and rewritten with stacking calls added")
    ap.add_argument("--hmm-dat", type=Path, required=True, help="Pfam-A.hmm.dat, for names and clans")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--min-instances", type=int, default=3,
                    help="instances a family needs before its box dimension is used as its height")
    ap.add_argument("--max-extent-per-residue", type=float, default=4.0)
    ap.add_argument("--height-from-stack-pitch", action="store_true",
                    help="use the stack pitch as the height where one exists (off: box dimension for every family)")
    ap.add_argument("--height-from-neighbour-pitch", action="store_true",
                    help="use the neighbour pitch as the height where one exists, the box dimension along the "
                         "called axis only for families no structure shows beside a neighbour")
    ap.add_argument("--min-pitch-instances", type=int, default=3,
                    help="instances a family needs before its neighbour pitch is used as its height")
    ap.add_argument("--min-coverage", type=float, default=0.8,
                    help="share of the Pfam model an instance must cover to count for any height")
    a = ap.parse_args()

    m = pd.read_csv(a.measurements, dtype={"PDB": str}, low_memory=False)
    bad = m["max_dim"] > a.max_extent_per_residue * m["n_res"]
    print(f"{len(m):,} instances; {int(bad.sum())} dropped as longer than "
          f"{a.max_extent_per_residue} A per residue")
    m = m[~bad]
    partial = m["n_res"] < a.min_coverage * m["model_len"]
    print(f"{int(partial.sum()):,} partial instances (under {a.min_coverage:.0%} of the model) "
          f"kept for the composite test but not for any height")
    m_all = m
    m = m[~partial]

    g = m.groupby("PFAM_ID")
    fam = pd.DataFrame({
        "max_dim_mean": g["max_dim"].mean(), "max_dim_std": g["max_dim"].std(),
        "max_dim_sem": g["max_dim"].sem(), "max_dim_len": g["max_dim"].size(),
        "min_dim_mean": g["min_dim"].mean(), "min_dim_std": g["min_dim"].std(),
        "min_dim_sem": g["min_dim"].sem(), "min_dim_len": g["min_dim"].size(),
        "n_structures": g["PDB"].nunique(), "median_n_res": g["n_res"].median(),
    }).reset_index().rename(columns={"PFAM_ID": "Pfam"})
    few = fam["max_dim_len"] < a.min_instances

    hand = pd.read_csv(a.orientations) if a.orientations.is_file() else \
        pd.DataFrame(columns=["Pfam", "orientation", "basis"])
    hand = hand.dropna(subset=["orientation"])
    # stacking and clan-majority calls from an earlier run are recomputed, not carried over
    derived = hand["basis"].astype(str).str.startswith(("stacking", "clan majority"))
    hand = hand[~derived]
    stack = stacking_calls(m)
    stack["basis"] = stack.apply(lambda r: f"stacking: {r.n_pairs} pairs, "
                                           f"{r.frac_along_long:.0%} along the longest axis", axis=1)
    calls = pd.concat([hand[["Pfam", "orientation", "basis"]],
                       stack[~stack["Pfam"].isin(hand["Pfam"])][["Pfam", "orientation", "basis"]]],
                      ignore_index=True)
    meta = pfam_metadata(a.hmm_dat)
    clan = clan_majority_calls(calls, meta, fam["Pfam"])
    calls = pd.concat([calls, clan], ignore_index=True)
    calls.to_csv(a.orientations, index=False)

    fam = fam.merge(calls, on="Pfam", how="left").rename(columns={"basis": "orientation_basis"})
    fam["height"] = np.where(fam["orientation"] == "max", fam["max_dim_mean"],
                             np.where(fam["orientation"] == "min", fam["min_dim_mean"], np.nan))
    fam["height_basis"] = np.where(fam["orientation"].isna(), None, "box " + fam["orientation"].astype(str))
    fam.loc[few.values, "height"] = np.nan
    fam.loc[few.values, "height_basis"] = "too few whole instances"
    fam = fam.merge(stack_pitch(m), on="Pfam", how="left")
    fam = fam.merge(neighbour_pitch(m), on="Pfam", how="left")
    if a.height_from_neighbour_pitch:
        use = fam["neighbour_pitch"].notna() & (fam["pitch_n"] >= a.min_pitch_instances)
        fam.loc[use, "height"] = fam.loc[use, "neighbour_pitch"]
        fam.loc[use, "height_basis"] = "neighbour pitch"
        fam.loc[~use & fam["height"].notna(), "height_basis"] += ", no neighbour"
    if a.height_from_stack_pitch:
        use = fam["stack_pitch"].notna()
        fam.loc[use, "height"] = fam.loc[use, "stack_pitch"]
        fam.loc[use, "height_basis"] = "stack pitch"
    comp = composite_families(m_all, min_coverage=a.min_coverage)
    is_comp = fam["Pfam"].isin(comp["Pfam"])
    fam = fam.merge(comp, on="Pfam", how="left")
    fam.loc[is_comp, "height"] = np.nan
    fam.loc[is_comp, "height_basis"] = "composite: spans " + fam.loc[is_comp, "inner_domains"].astype(int).astype(str) + " domains"
    fam = fam.merge(meta, on="Pfam", how="left")
    cols = ["Pfam", "max_dim_mean", "max_dim_std", "max_dim_sem", "max_dim_len",
            "min_dim_mean", "min_dim_std", "min_dim_sem", "min_dim_len", "query_name_x",
            "orientation", "Clan_x", "height", "orientation_basis", "n_structures", "median_n_res",
            "height_basis", "stack_pitch", "stack_runs", "neighbour_pitch", "pitch_n", "pitch_structures", "pitch_interior", "inner_domains"]
    fam = fam[cols].sort_values("Pfam")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    fam.to_csv(a.out, index=False)

    print(f"{int(is_comp.sum())} composite families carry no height: "
          + ", ".join(fam.loc[is_comp, "query_name_x"].astype(str).head(12)) + (" ..." if is_comp.sum() > 12 else ""))
    need = fam[fam["orientation"].isna() & ~is_comp & fam["neighbour_pitch"].isna()]
    need.to_csv(Path(str(a.out).replace(".csv", "")).with_suffix(".needs_call.csv"), index=False)
    by = fam["orientation_basis"].fillna("none").str.split(":").str[0].value_counts()
    print(f"{len(fam)} families measured; orientation by basis: {by.to_dict()}")
    print(f"stack pitch recorded for {int(fam['stack_pitch'].notna().sum())} families"
          + (" and used as their height" if a.height_from_stack_pitch else ""))
    print(f"neighbour pitch recorded for {int(fam['neighbour_pitch'].notna().sum())} families"
          + (f"; height by basis: {fam['height_basis'].value_counts(dropna=False).to_dict()}"
             if a.height_from_neighbour_pitch else ""))
    print(f"{len(need)} families need a call by hand -> {a.out.with_suffix('.needs_call.csv')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
