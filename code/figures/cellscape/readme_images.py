#!/usr/bin/env python
"""Draw the README's CellScape illustrations from a run's output.

    conda env create -f code/figures/cellscape/environment.yml
    conda activate cellscape
    python code/figures/cellscape/readme_images.py RUN_DIR --alphafold DIR

Two images, each in a light and a dark version, under ``docs/img/``:

* ``cd45_{light,dark}.svg``: CD45's ectodomain standing on the membrane, with each piece of
  its height marked;
* ``monocyte_surface_{light,dark}.svg``: twenty surface proteins sampled from a classical
  monocyte's proteome by abundance, each drawn to scale; a protein with no AlphaFold model (LRP1's
  ectodomain is too long to predict) is a capsule of its whole height.
* ``contact_{light,dark}.svg``: an NK cell against a HER2-positive breast cancer cell with an
  anti-HER2 antibody (the Figure 6 panel ``nk_aHER2``): for each gap band of the bullseye, the
  pair of single proteins holding most of that band among pairs whose height is the sum of their
  partners' (less the pairs in ``LEFT_OUT_OF_CONTACT``), the NK cell's partner hanging from its
  membrane;
* ``bridge_{light,dark}.svg``: that contact's antibody bridge, drawn the way
  ``surfaceomeTopography.interface.antibody_bridge_heights`` models it: HER2 on the cancer cell,
  trastuzumab lying across HER2's membrane-proximal epitope (an intact human IgG1, PDB 1HZH, its
  thinnest axis vertical, Fab towards HER2, Fc away), and the Fc receptor holding most of the
  contact's bridges hanging from the NK cell onto the Fc; the numbers printed are the table's;
* ``bullseye_{light,dark}.svg``: that contact in section and in plan: the plan is the pipeline's
  own bullseye (``surfaceomeTopography.interface.plot_contact_bullseye``, from the panel's pairs
  table, in the README's blue ramp); the section cuts through it, each ring at the share-weighted
  mean gap of the pairs holding it.

Each protein is the span of its AlphaFold v6 model the height uses, outlined with CellScape
(Silvestre-Ryan et al., https://github.com/jordisr/cellscape) and oriented as CellScape orients a
membrane protein, by its topology: the N-to-C vector vertical and the membrane-side terminus at
the membrane. That side comes from the neighbouring transmembrane segment in the UniProt GFF (the
C-terminal end when the segment follows the ectodomain, as for a type I protein or a GPI anchor).
Whatever height a protein takes from disorder, Pfam domains or sequence is a plain grey capsule
of that length on the side of the structure where those residues lie. Everything is drawn to one
scale, but a structure's drawn extent is its extent along the N-to-C axis, not the dimension the
height table measures, so drawn heights are illustrative; the numbers printed are the table's.

``RUN_DIR`` is a run's output (``python code/run_notebooks.py RUN_DIR``); ``--alphafold`` holds
the ``AF-<accession>-F1-model_v6.pdb`` files (``code/database/inputs/download_alphafold.py``).
"""

from __future__ import annotations

import argparse
import gzip
import importlib.util
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import cellscape
from cellscape.cartoon import plot_polygon, shade_from_color
from cellscape.scene import Membrane
from shapely.geometry import LineString
from shapely.affinity import scale, translate

REPO = Path(__file__).resolve().parents[3]
DATA = REPO / "data"
OUT = REPO / "docs" / "img"

#: A single-hue sequential ramp (blue, light to dark) from the dataviz reference palette:
#: a protein's colour steps with its height. Capsules are neutral grey, as in CellScape.
RAMP = ["#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281"]
CAPSULE = "#e4e3de"
#: Which cell a protein is on, in the contact figure: the first two categorical slots.
IMMUNE_COLOR, CANCER_COLOR = "#2a78d6", "#eb6834"
ANTIBODY_COLOR = "#1baf7a"     # the third categorical slot
IGG_PDB = "https://files.rcsb.org/download/1HZH.pdb"     # intact human IgG1 (Saphire et al. 2001)
#: The Figure 6 contact the README shows, and what to call its two cells.
PANEL = "nk_aHER2"
#: Pairs the pipeline keeps that the contact illustration leaves out, with why. CD44-CD44 stays in
#: every contact as curated, but its binding across cells has been shown only between tumour cells
#: (Liu et al. 2019; Kawaguchi et al. 2020), so it is not the example a reader should take away.
LEFT_OUT_OF_CONTACT = {("P16070", "P16070"): "CD44-CD44: trans binding shown only between tumour cells"}
IMMUNE_CELL, CANCER_CELL = "NK cell", "HER2+ breast\ncancer cell"
THEMES = {"light": {"text": "#0b0b0b", "text2": "#52514e", "rule": "#8a8984", "surface": "#ffffff"},
          "dark": {"text": "#ffffff", "text2": "#c3c2b7", "rule": "#8f8e87", "surface": "#0d1117"}}
#: The bullseye's rings, 50+ nm outermost to 10-20 nm, as an ordinal step of the same blue ramp:
#: darker is taller. Light mode starts no lighter than step 250 and dark mode goes no darker than
#: step 600, so every ring clears 2:1 against its page; the centre (under 10 nm) is the page.
BULLSEYE_RINGS = {"light": ["#104281", "#1c5cab", "#2a78d6", "#5598e7", "#86b6ef"],
                  "dark": ["#184f95", "#256abf", "#3987e5", "#6da7ec", "#9ec5f4"]}
NM = 10.0                      # angstroms per nm; CellScape works in angstroms
CAPSULE_RADIUS = 12.0          # half-width of a capsule, angstroms
SEQ_RATE = 0.04                # nm per residue nothing else covers (METHODS.md §2)

plt.rcParams["svg.fonttype"] = "none"   # keep text as text in the SVG
plt.rcParams["svg.hashsalt"] = "surfaceome"  # stable element ids, so a redraw diffs only if the picture moved
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "DejaVu Sans"]


# ---- what the run says about each protein ------------------------------------------------

def load_run(run: Path):
    heights = pd.read_csv(run / "database" / "height_estimates.csv")
    parts = pd.read_csv(run / "database" / "domain_disorder_assignments_alphafold.csv", low_memory=False)
    return heights, parts


def ecd_range(ecd_id: str) -> tuple[int, int]:
    s, e = re.search(r"\[(\d+):(\d+)\]", ecd_id).groups()
    return int(s), int(e)


def membrane_end(acc: str, start: int, end: int, gff: pd.DataFrame) -> str:
    """'C' when a transmembrane segment or GPI anchor follows the ectodomain, else 'N' if one
    precedes it, else 'C'."""
    feats = gff[gff["acc"] == acc]
    tm = feats[feats["type"] == "Transmembrane"]
    if ((tm["start"] - end).between(1, 3)).any():
        return "C"
    gpi = feats[(feats["type"] == "Lipidation") & feats["note"].str.contains("GPI", na=False)]
    if ((gpi["start"] - end).abs() <= 3).any():
        return "C"
    if ((start - tm["end"]).between(1, 3)).any():
        return "N"
    return "C"


def load_gff() -> pd.DataFrame:
    rows = []
    with gzip.open(DATA / "inputs/proteome/UP000005640.gff.gz", "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > 8 and f[2] in ("Transmembrane", "Lipidation"):
                rows.append((f[0], f[2], int(f[3]), int(f[4]), f[8]))
    return pd.DataFrame(rows, columns=["acc", "type", "start", "end", "note"])


def layout(row: pd.Series, parts: pd.DataFrame, gff: pd.DataFrame) -> dict:
    """The structured span and the non-structure height on each side of it, in nm."""
    acc, ecd_id, total = row["ID link_first"], row["ID"], float(row["total_height"])
    start, end = ecd_range(ecd_id)
    mine = parts[parts["ID"] == ecd_id]
    af = mine[mine["methods"] == "alphafold"]
    side = membrane_end(acc, start, end, gff)
    out = {"acc": acc, "total": total, "side": side, "af": None, "before": 0.0, "after": total}
    if af.empty:
        return out
    a = af.iloc[0]
    s, e, h = int(a["alphafold_start"]), int(a["alphafold_end"]), float(a["height"])
    out["af"] = (s, e, h)
    # non-structure pieces and uncovered residues, by which side of the structure they lie on
    n_side = c_side = 0.0
    covered = np.zeros(end - start + 1, bool)
    covered[max(s, start) - start:min(e, end) - start + 1] = True
    for _, p in mine[mine["methods"] != "alphafold"].iterrows():
        ps = int(p["disorder_start"] if p["methods"] == "disorder" else p["domain_start"])
        pe = int(p["disorder_end"] if p["methods"] == "disorder" else p["domain_end"])
        covered[max(ps, start) - start:min(pe, end) - start + 1] = True
        if pe <= s:
            n_side += float(p["height"])
        else:
            c_side += float(p["height"])
    idx = np.arange(start, end + 1)
    n_side += SEQ_RATE * ((~covered) & (idx < s)).sum()
    c_side += SEQ_RATE * ((~covered) & (idx > e)).sum()
    rest = total - h
    if n_side + c_side > 0:                      # scale to the table's height exactly
        n_side, c_side = rest * n_side / (n_side + c_side), rest * c_side / (n_side + c_side)
    else:
        c_side = rest if side == "N" else 0.0
        n_side = rest if side == "C" else 0.0
    # 'before' lies between the membrane and the structure, 'after' beyond it
    out["before"], out["after"] = (c_side, n_side) if side == "C" else (n_side, c_side)
    return out


# ---- drawing -------------------------------------------------------------------------------

def protein(acc: str, heights: pd.DataFrame, parts: pd.DataFrame, gff: pd.DataFrame, alphafold: Path) -> dict:
    """A protein's layout; with no AlphaFold file at hand its whole height is one capsule."""
    item = layout(heights[heights["ID link_first"] == acc].iloc[0], parts, gff)
    if item["af"] is not None and not (alphafold / f"AF-{acc}-F1-model_v6.pdb").is_file():
        item.update(af=None, before=0.0, after=item["total"])
    return item


def capsule(x: float, y0: float, length: float):
    """A rounded bar from y0 to y0 + length (angstroms), centred on x."""
    r = min(CAPSULE_RADIUS, length / 2)
    return LineString([(x, y0 + r), (x, y0 + length - r)]).buffer(r)


def build(item: dict, alphafold: Path, color: str) -> dict:
    """Polygons for one protein, bottom at y = 0 and left edge at x = 0 (angstroms)."""
    polys = []
    before, after = item["before"] * NM, item["after"] * NM
    if item["af"] is None:
        polys.append({"polygon": capsule(CAPSULE_RADIUS, 0, item["total"] * NM), "facecolor": CAPSULE,
                      "edgecolor": "black", "linewidth": 0.6, "zorder": 0})
        return {"polygons": polys, "width": 2 * CAPSULE_RADIUS, "height": item["total"] * NM,
                "bottom_x": CAPSULE_RADIUS, "top_x": CAPSULE_RADIUS, "segments": []}
    s, e, h = item["af"]
    mol = cellscape.Structure(str(alphafold / f"AF-{item['acc']}-F1-model_v6.pdb"), name=item["acc"],
                              view=False, res_start=s, res_end=e)
    # CellScape's orientation from topology: the N-to-C vector vertical, the membrane-side
    # terminus down (flip=True puts the N terminus up, for a C-terminal membrane anchor)
    mol.auto_view(flip=item["side"] == "C")
    cart = mol.outline("all", depth="contours", depth_contour_interval=10, back_outline=True)
    cart.plot(do_show=False, colors=[color], depth_shading=True, line_width=0.4)
    plt.close("all")
    d = cart.dimensions
    xb, xt, hs = float(d["bottom"][0]), float(d["top"][0]), float(d["height"])
    for p in cart._styled_polygons:
        polys.append({"polygon": translate(p["polygon"].simplify(0.4), 0, before), "facecolor": p["facecolor"],
                      "edgecolor": p["edgecolor"], "linewidth": p["linewidth"], "zorder": 1})
    if before > 0:
        polys.append({"polygon": capsule(xb, 0, before + CAPSULE_RADIUS), "facecolor": CAPSULE,
                      "edgecolor": "black", "linewidth": 0.6, "zorder": 0})
    if after > 0:
        polys.append({"polygon": capsule(xt, before + hs - CAPSULE_RADIUS, after + CAPSULE_RADIUS),
                      "facecolor": CAPSULE, "edgecolor": "black", "linewidth": 0.6, "zorder": 0})
    minx = min(p["polygon"].bounds[0] for p in polys)
    maxx = max(p["polygon"].bounds[2] for p in polys)
    for p in polys:
        p["polygon"] = translate(p["polygon"], -minx, 0)
    return {"polygons": polys, "width": maxx - minx, "height": before + hs + after,
            "bottom_x": xb - minx, "top_x": xt - minx,
            "segments": [("before", 0, before), ("structure", before, before + hs),
                         ("after", before + hs, before + hs + after)]}


def draw(ax, obj: dict, x0: float):
    for p in sorted(obj["polygons"], key=lambda p: p["zorder"]):
        plot_polygon(p["polygon"], facecolor=p["facecolor"], edgecolor=p["edgecolor"],
                     linewidth=p["linewidth"], axes=ax, translate_pre=np.array([x0, 0]),
                     zorder_mod=p["zorder"])


def membrane(ax, width: float):
    # CellScape centres a flat membrane on base_y; ectodomains start at the outer surface, y = 0
    m = Membrane(width=width, thickness=40, axes=ax, base_y=-20)
    m.flat()
    m.draw(lipids=True)


def color_for(height_nm: float, lo: float, hi: float) -> str:
    f = (np.log(height_nm) - np.log(lo)) / (np.log(hi) - np.log(lo)) if hi > lo else 0.5
    return RAMP[int(round(np.clip(f, 0, 1) * (len(RAMP) - 1)))]


def cd45_figure(item: dict, obj: dict, row: pd.Series, parts: pd.DataFrame, out: Path):
    dis = parts[(parts["ID"] == row["ID"]) & (parts["methods"] == "disorder")]
    glycans = int(dis["glycan_count"].sum())
    ds, de = int(dis["disorder_start"].min()), int(dis["disorder_end"].max()) - 1
    s, e, h = item["af"]
    for mode, t in THEMES.items():
        fig, ax = plt.subplots(figsize=(6.4, 4.6))
        x0 = 70.0
        membrane(ax, x0 * 2 + obj["width"])
        draw(ax, obj, x0)
        bx = x0 * 2 + obj["width"] + 12
        labels = {"structure": (f"Folded region · {h:.1f} nm", f"AlphaFold model, residues {s}–{e}"),
                  "after": (f"Disordered region · {item['after']:.1f} nm",
                            f"residues {ds}–{de}, {glycans} glycosylation sites;\na worm-like chain, stiffened by its glycans")}
        for name, y0, y1 in obj["segments"]:
            if y1 - y0 <= 0 or name not in labels:
                continue
            ax.plot([bx, bx + 8, bx + 8, bx], [y0 + 3, y0 + 3, y1 - 3, y1 - 3], color=t["rule"], lw=1)
            head, sub = labels[name]
            ym = (y0 + y1) / 2
            ax.text(bx + 22, ym + 14, head, color=t["text"], fontsize=11, fontweight="bold", va="bottom")
            ax.text(bx + 22, ym + 8, sub, color=t["text2"], fontsize=9.5, va="top", linespacing=1.4)
        top = obj["height"]
        ax.plot([x0 + obj["top_x"] + 20, bx + 250], [top, top], color=t["rule"], lw=0.8, ls=(0, (2, 3)))
        ax.text(bx + 256, top, f"{item['total']:.1f} nm", color=t["text"], fontsize=13, fontweight="bold", va="center")
        ax.text(bx + 256, top - 22, "height of CD45", color=t["text2"], fontsize=9.5, va="center")
        ax.set_xlim(0, bx + 360)
        ax.set_ylim(-60, top + 30)
        ax.set_aspect("equal")
        ax.axis("off")
        path = out / f"cd45_{mode}.svg"
        fig.savefig(path, transparent=True, bbox_inches="tight", pad_inches=0.05, metadata={"Date": None})
        plt.close(fig)
        print(f"  wrote {path.relative_to(REPO)}")


def broken(obj: dict, limit: float) -> dict:
    """A protein taller than ``limit`` cut off there, its drawn part ending in a break."""
    from shapely.geometry import box
    keep = box(-1e4, -1e4, 1e4, limit)
    out = dict(obj, height=limit, true_height=obj["height"])
    out["polygons"] = [dict(p, polygon=p["polygon"].intersection(keep)) for p in obj["polygons"]
                       if not p["polygon"].intersection(keep).is_empty]
    return out


def scene_figure(objs: list[tuple[str, dict]], out: Path, padding: float = 14.0):
    # a protein more than twice as tall as the next tallest is broken just above it, so the rest
    # stay legible; its true height is printed over the break
    tallest = sorted({o["height"] for _, o in objs}, reverse=True)
    if len(tallest) > 1 and tallest[0] > 2 * tallest[1]:
        limit = tallest[1] + 140
        objs = [(g, broken(o, limit) if o["height"] > 2 * tallest[1] else o) for g, o in objs]
    width = sum(o["width"] for _, o in objs) + padding * (len(objs) + 1)
    top = max(o["height"] for _, o in objs) + 40
    fontsize = 11                                # pt; about 14 px at the README's width
    pt_per_a = 12 * 72 / width                   # the figure is 12 in wide
    label_a = max(len(g) for g, _ in objs) * 0.62 * fontsize / pt_per_a
    bottom = -58 - label_a - 10                  # room under the membrane for the longest label
    for mode, t in THEMES.items():
        fig, ax = plt.subplots(figsize=(12, 12 * (top - bottom + 10) / width))
        membrane(ax, width)
        x = padding
        for gene, o in objs:
            draw(ax, o, x)
            if "true_height" in o:                           # the break, and the true height
                cx, yb = x + o["width"] / 2, o["height"] - 70
                for dy in (0, 22):
                    ax.plot([cx - 22, cx + 22], [yb + dy - 8, yb + dy + 8], color="black", lw=0.8, zorder=5)
                ax.fill([cx - 22, cx + 22, cx + 22, cx - 22], [yb - 8, yb + 8, yb + 30, yb + 14],
                        color=t["surface"], zorder=4, lw=0)
                ax.text(cx, o["height"] + 12, f"{o['true_height'] / NM:.0f} nm", color=t["text"],
                        fontsize=10, fontweight="bold", ha="center", va="bottom")
            ax.text(x + o["width"] / 2, -58, gene, color=t["text2"], fontsize=fontsize, rotation=90,
                    ha="center", va="top")
            x += o["width"] + padding
        ax.set_xlim(0, width)
        ax.set_ylim(bottom, top + 10)
        ax.set_aspect("equal")
        ax.axis("off")
        path = out / f"monocyte_surface_{mode}.svg"
        fig.savefig(path, transparent=True, bbox_inches="tight", pad_inches=0.05, metadata={"Date": None})
        plt.close(fig)
        print(f"  wrote {path.relative_to(REPO)}")


def bilayer(ax, x0: float, x1: float, top: float, thickness: float = 40.0, r: float = 4.0):
    """A lipid bilayer from x0 to x1 whose outer surface is at y = top (CellScape's colours)."""
    ax.fill_between([x0, x1], top - r, top - thickness + r, color="#C4E7EF", zorder=1.6, lw=0)
    for x in np.arange(x0 + r, x1 - r + 1e-6, 2 * r):
        for y in (top - r, top - thickness + r):
            ax.add_patch(plt.Circle((x, y), r, facecolor="#D6D1EF", ec="k", lw=0.3, zorder=2))


def flipped(obj: dict, gap: float) -> dict:
    """The object hanging from a membrane at y = gap instead of standing on one at y = 0."""
    out = dict(obj)
    out["polygons"] = [dict(p, polygon=translate(scale(p["polygon"], yfact=-1, origin=(0, 0)), 0, gap))
                       for p in obj["polygons"]]
    return out


def contact_pairs(run: Path, heights: pd.DataFrame, parts: pd.DataFrame, alphafold: Path) -> pd.DataFrame:
    """Per gap band, the single-protein pair holding most of it, among pairs whose height is the
    sum of their partners'."""
    pairs = pd.read_csv(run / "tables" / f"F6_{PANEL}" / f"{PANEL}_pairs.csv")
    pairs["share"] = pairs["interface_shared"] / pairs["interface_shared"].sum()
    single = pairs["complex_name_immune"].isna() & pairs["complex_name_cancer"].isna()
    additive = (pairs["interaction_dim"] - pairs["total_height_prot1"] - pairs["total_height_prot2"]).abs() < 0.01
    left_out = pd.Series([(i, c) in LEFT_OUT_OF_CONTACT for i, c in zip(pairs["Entry_immune"], pairs["Entry_cancer"])],
                         index=pairs.index)
    ok = pairs[single & additive & ~left_out].copy()
    edges = [0, 10, 20, 30, 40, 50, np.inf]
    ok["band"] = pd.cut(ok["interaction_dim"], edges, right=False)
    best = ok.sort_values("share", ascending=False).groupby("band", observed=True).head(1)
    return best.sort_values("interaction_dim")


def contact_figure(chosen: pd.DataFrame, heights: pd.DataFrame, parts: pd.DataFrame, gff: pd.DataFrame,
                   alphafold: Path, out: Path, padding: float = 40.0, margin: float = 230.0):
    built = []
    for _, r in chosen.iterrows():
        im = build(protein(r["Entry_immune"], heights, parts, gff, alphafold), alphafold, IMMUNE_COLOR)
        ca = build(protein(r["Entry_cancer"], heights, parts, gff, alphafold), alphafold, CANCER_COLOR)
        drawn = im["height"] + ca["height"]              # the partners meet; the label gives the table's height
        built.append((r, drawn, flipped(im, drawn), ca))
    # a pair more than twice as tall as the next is shortened through its middle, with a break
    gaps = sorted((g for _, g, _, _ in built), reverse=True)
    breaks = {}
    if len(gaps) > 1 and gaps[0] > 2 * gaps[1]:
        limit = gaps[1] + 250
        from shapely.geometry import box
        for k, (r, g, im, ca) in enumerate(built):
            if g > 2 * gaps[1]:
                cut = limit / 2
                low, high = box(-1e4, -1e4, 1e4, cut), box(-1e4, cut, 1e4, 1e4)
                im = dict(im, polygons=[dict(q, polygon=translate(q["polygon"], 0, limit - g).intersection(high))
                                        for q in im["polygons"]])
                ca = dict(ca, polygons=[dict(q, polygon=q["polygon"].intersection(low)) for q in ca["polygons"]])
                im["polygons"] = [q for q in im["polygons"] if not q["polygon"].is_empty]
                ca["polygons"] = [q for q in ca["polygons"] if not q["polygon"].is_empty]
                built[k] = (r, limit, im, ca)
                breaks[k] = cut
    top = max(g for _, g, _, _ in built)
    width = margin + sum(max(i["width"], c["width"]) + padding for _, _, i, c in built)
    for mode, t in THEMES.items():
        fig, ax = plt.subplots(figsize=(7, 7 * (top + 260) / width))
        x = margin
        for k, (r, gap, im, ca) in enumerate(built):
            w = max(im["width"], ca["width"])
            bilayer(ax, x - padding / 2, x + w + padding / 2, gap + 40)          # the immune cell's
            draw(ax, im, x + (w - im["width"]) / 2)
            draw(ax, ca, x + (w - ca["width"]) / 2)
            if k in breaks:
                cx, yb = x + w / 2, breaks[k] - 11
                ax.fill([cx - 22, cx + 22, cx + 22, cx - 22], [yb - 8, yb + 8, yb + 30, yb + 14],
                        color=t["surface"], zorder=4, lw=0)
                for dy in (0, 22):
                    ax.plot([cx - 22, cx + 22], [yb + dy - 8, yb + dy + 8], color="black", lw=0.8, zorder=5)
            ax.text(x + w / 2, gap + 58, r["Gene Name_immune"], color=t["text"], fontsize=9,
                    ha="center", va="bottom")
            ax.text(x + w / 2, -58, r["Gene Name_cancer"], color=t["text"], fontsize=9, ha="center", va="top")
            ax.text(x + w / 2, -84, f"{r['interaction_dim']:.0f} nm", color=t["text2"], fontsize=8.5,
                    ha="center", va="top")
            x += w + padding
        bilayer(ax, margin - padding / 2, x - padding / 2, 0)                    # the cancer cell's
        ax.text(margin - padding, built[0][1] + 20, IMMUNE_CELL, color=t["text"], fontsize=10,
                fontweight="bold", ha="right", va="center")
        ax.text(margin - padding, -20, CANCER_CELL, color=t["text"], fontsize=10,
                fontweight="bold", ha="right", va="center")
        ax.set_xlim(0, x)
        ax.set_ylim(-130, top + 110)
        ax.set_aspect("equal")
        ax.axis("off")
        path = out / f"contact_{mode}.svg"
        fig.savefig(path, transparent=True, bbox_inches="tight", pad_inches=0.05, metadata={"Date": None})
        plt.close(fig)
        print(f"  wrote {path.relative_to(REPO)}")


def antibody(pdb: Path, color: str) -> dict:
    """An intact IgG lying flat: Fab-to-Fc along x, its thinnest axis vertical, Fab end at x = 0."""
    mol = cellscape.Structure(str(pdb), name="IgG", view=False)
    X = mol.coord - mol.coord.mean(0)
    heavy_fc = np.zeros(len(X), bool)                    # Fc: the heavy chains from residue 240
    for chain in "HK":
        for res_id, res in mol.residues.get(chain, {}).items():
            if res_id >= 240:
                heavy_fc[slice(*res["coord"])] = True
    u = X[heavy_fc].mean(0) - X[~heavy_fc].mean(0)
    u /= np.linalg.norm(u)                               # Fab -> Fc
    rest = X - np.outer(X @ u, u)
    _, vec = np.linalg.eigh(np.cov(rest.T))
    v = vec[:, 0] - (vec[:, 0] @ u) * u                  # thinnest direction across the Fab-Fc axis
    v /= np.linalg.norm(v)
    mol.set_view_matrix(np.column_stack([u, v, np.cross(u, v)]))
    cart = mol.outline("all", depth="contours", depth_contour_interval=10, back_outline=True)
    cart.plot(do_show=False, colors=[color], depth_shading=True, line_width=0.4)
    plt.close("all")
    polys = [{"polygon": p["polygon"].simplify(0.4), "facecolor": p["facecolor"], "edgecolor": p["edgecolor"],
              "linewidth": p["linewidth"], "zorder": 1} for p in cart._styled_polygons]
    minx = min(p["polygon"].bounds[0] for p in polys)
    miny = min(p["polygon"].bounds[1] for p in polys)
    for p in polys:
        p["polygon"] = translate(p["polygon"], -minx, -miny)
    return {"polygons": polys, "width": max(p["polygon"].bounds[2] for p in polys),
            "height": max(p["polygon"].bounds[3] for p in polys)}


def bridge_figure(run: Path, heights: pd.DataFrame, parts: pd.DataFrame, gff: pd.DataFrame,
                  alphafold: Path, igg: Path, out: Path, margin: float = 230.0):
    pairs = pd.read_csv(run / "tables" / f"F6_{PANEL}" / f"{PANEL}_pairs.csv")
    bridges = pairs[pairs["_merge"] == "FcR"].sort_values("interface_shared", ascending=False)
    br = bridges.iloc[0]                                  # the bridge holding most of the contact
    receptor_acc, antigen_acc = br["Entry_immune"], br["Entry_cancer"]
    receptor_h, antigen_h = float(br["total_height_prot2"]), float(br["total_height_prot1"])
    gap = float(br["interaction_dim"])
    antigen = build(protein(antigen_acc, heights, parts, gff, alphafold), alphafold, CANCER_COLOR)
    receptor = build(protein(receptor_acc, heights, parts, gff, alphafold), alphafold, IMMUNE_COLOR)
    ab = antibody(igg, ANTIBODY_COLOR)
    x_ab = margin + antigen["width"] - 6                  # Fab against the antigen's base
    x_fc = x_ab + ab["width"] - receptor["width"] / 2 - 25    # receptor over the Fc
    top = max(ab["height"] + receptor["height"], antigen["height"] + 5)   # receptor meets the Fc
    right = max(x_ab + ab["width"], x_fc + receptor["width"]) + 40
    names = {"antigen": f"{br['Gene Name_cancer']}", "receptor": f"{br['Gene Name_immune']}"}
    for mode, t in THEMES.items():
        fig, ax = plt.subplots(figsize=(7, 7 * (top + 260) / (right + 330)))
        bilayer(ax, margin - 30, right, 0)
        bilayer(ax, margin - 30, right, top + 40)
        draw(ax, antigen, margin)
        draw(ax, ab, x_ab)
        draw(ax, flipped(receptor, top), x_fc)
        ax.text(margin + antigen["width"] / 2, -58, names["antigen"], color=t["text"], fontsize=9,
                ha="center", va="top")
        ax.text(x_fc + receptor["width"] / 2, top + 58, names["receptor"], color=t["text"], fontsize=9,
                ha="center", va="bottom")
        ax.text(x_ab + ab["width"] / 2, -58, "anti-HER2 IgG", color=t["text2"], fontsize=9, ha="center", va="top")
        bx = right + 20
        ax.plot([bx, bx + 8, bx + 8, bx], [3, 3, top - 3, top - 3], color=t["rule"], lw=1)
        ax.text(bx + 22, top / 2 + 12, f"{gap:.1f} nm", color=t["text"], fontsize=11, fontweight="bold", va="bottom")
        ax.text(bx + 22, top / 2 + 4,
                f"{names['receptor']} {receptor_h:.1f} nm\n+ antibody {gap - receptor_h:.2f} nm,\n"
                f"epitope at the membrane\n({names['antigen']} {antigen_h:.1f} nm fits beneath)",
                color=t["text2"], fontsize=8.5, va="top", linespacing=1.4)
        ax.text(margin - 40, top + 20, IMMUNE_CELL, color=t["text"], fontsize=10, fontweight="bold",
                ha="right", va="center")
        ax.text(margin - 40, -20, CANCER_CELL, color=t["text"], fontsize=10,
                fontweight="bold", ha="right", va="center")
        ax.set_xlim(0, bx + 330)
        ax.set_ylim(-110, top + 110)
        ax.set_aspect("equal")
        ax.axis("off")
        path = out / f"bridge_{mode}.svg"
        fig.savefig(path, transparent=True, bbox_inches="tight", pad_inches=0.05, metadata={"Date": None})
        plt.close(fig)
        print(f"  wrote {path.relative_to(REPO)}")


def bullseye_figure(run: Path, out: Path):
    """The contact in section and in plan: a cut through it above, the pipeline's bullseye below.

    The plan is the pipeline's own ``plot_contact_bullseye``: each circle's area is the share of the
    contact held at gaps below a band's upper edge. The section cuts through its centre, so each
    ring keeps its radius, and the monocyte's membrane stands above each ring at the mean gap of
    the pairs holding it, weighted by their shares. Radius is share of the contact, not distance.
    """
    spec = importlib.util.spec_from_file_location("sT_interface", REPO / "code/surfaceomeTopography/interface.py")
    interface = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(interface)
    pairs = pd.read_csv(run / "tables" / f"F6_{PANEL}" / f"{PANEL}_pairs.csv")
    gap, weight = pairs["interaction_dim"], pairs["interface_shared"]
    edges = [0, *sorted(interface.BULLSEYE_BANDS), np.inf]
    for mode, t in THEMES.items():
        colors = [*BULLSEYE_RINGS[mode], t["surface"]]
        plan = interface.plot_contact_bullseye(pairs, colors=colors)
        circles = [c for c in plan.axes[0].get_children() if isinstance(c, matplotlib.patches.Circle)]
        legend = [x.get_text() for x in plan.axes[0].get_legend().get_texts()]
        plt.close(plan)
        radius = np.array([c.get_radius() for c in circles]) / circles[0].get_radius()   # outside in
        # the rings from the centre out: radius, colour, and the mean gap of the pairs holding them
        uppers = [np.inf, *sorted(interface.BULLSEYE_BANDS, reverse=True)]   # circle i holds gaps < uppers[i]
        rings = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            i = uppers.index(hi)
            band = (gap >= lo) & (gap < hi)
            if weight[band].sum() > 0:
                rings.append((radius[i], colors[i], float(np.average(gap[band], weights=weight[band])), lo, hi))
        top = max(g for *_, g, _, _ in rings)

        fig = plt.figure(figsize=(7.2, 7.6))
        side = fig.add_axes([0.20, 0.56, 0.44, 0.38])
        plan_ax = fig.add_axes([0.20, 0.06, 0.44, 0.44 * 7.2 / 7.6])
        # section: the gap above each ring, between the two membranes
        r_in = 0.0
        for r_out, color, g, lo, hi in rings:
            for sgn in (-1, 1):
                x0, x1 = sorted((sgn * r_in, sgn * r_out))
                side.fill_between([x0, x1], 0, g, color=color, lw=0)
                side.fill_between([x0, x1], g, g + 4, color="#C4E7EF", lw=0)
                side.plot([x0, x1], [g, g], color="#9C97C9", lw=1.2)
                side.plot([x0, x1], [g + 4, g + 4], color="#9C97C9", lw=1.2)
            if r_out - r_in > 0.07:
                side.text((r_in + r_out) / 2, g + 7, f"{g:.0f} nm", color=t["text2"], fontsize=8.5,
                          ha="center", va="bottom")
            r_in = r_out
        side.fill_between([-1, 1], -4, 0, color="#C4E7EF", lw=0)
        side.plot([-1, 1], [0, 0], color="#9C97C9", lw=1.2)
        side.plot([-1, 1], [-4, -4], color="#9C97C9", lw=1.2)
        side.text(-1.04, top / 2, IMMUNE_CELL, color=t["text"], fontsize=10, fontweight="bold", ha="right", va="center")
        side.text(-1.04, -2, CANCER_CELL, color=t["text"], fontsize=10, fontweight="bold",
                  ha="right", va="center")
        side.set_xlim(-1, 1)
        side.set_ylim(-6, top + 16)
        for sp in ("top", "right", "bottom"):
            side.spines[sp].set_visible(False)
        side.spines["left"].set_color(t["rule"])
        side.spines["left"].set_position(("axes", 1.03))
        side.yaxis.tick_right()
        side.yaxis.set_label_position("right")
        side.tick_params(axis="y", colors=t["text2"], labelsize=8.5)
        side.set_ylabel("gap between the membranes, nm", color=t["text2"], fontsize=9)
        side.set_xticks([])
        side.patch.set_alpha(0)
        # plan: the pipeline's circles, redrawn at the section's scale
        for c, color in zip(circles, colors):
            plan_ax.add_patch(plt.Circle((0, 0), c.get_radius() / circles[0].get_radius(), facecolor=color,
                                         edgecolor=t["surface"], lw=1.5))
        plan_ax.set_xlim(-1, 1)
        plan_ax.set_ylim(-1, 1)
        plan_ax.set_aspect("equal")
        plan_ax.axis("off")
        handles = [matplotlib.patches.Patch(facecolor=col) for col in colors[:len(legend)]]
        leg = plan_ax.legend(handles, legend, loc="center left", bbox_to_anchor=(1.08, 0.5), frameon=False,
                             fontsize=8.5, handlelength=1.2, labelspacing=0.8)
        for x in leg.get_texts():
            x.set_color(t["text"])
        fig.text(0.42, 0.955, "section through the contact", color=t["text2"], fontsize=9, ha="center")
        fig.text(0.42, 0.515, "the contact seen from above", color=t["text2"], fontsize=9, ha="center")
        path = out / f"bullseye_{mode}.svg"
        fig.savefig(path, transparent=True, bbox_inches="tight", pad_inches=0.05, metadata={"Date": None})
        plt.close(fig)
        print(f"  wrote {path.relative_to(REPO)}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", type=Path, help="a run's output directory")
    ap.add_argument("--alphafold", type=Path, required=True, help="directory of AF-<acc>-F1-model_v6.pdb files")
    ap.add_argument("--n", type=int, default=20, help="proteins in the monocyte scene")
    ap.add_argument("--seed", type=int, default=3, help="random seed for sampling the scene")
    ap.add_argument("--igg", type=Path, default=REPO / "downloads" / "pdb" / "1HZH.pdb",
                    help=f"an intact IgG structure (downloaded from {IGG_PDB} if absent)")
    a = ap.parse_args()

    heights, parts = load_run(a.run_dir)
    gff = load_gff()
    OUT.mkdir(parents=True, exist_ok=True)
    lo, hi = 2.0, 60.0                           # nm: the colour ramp's ends

    # CD45
    row = heights[heights["Entry name_first"] == "PTPRC_HUMAN"].iloc[0]
    item = layout(row, parts, gff)
    cd45_figure(item, build(item, a.alphafold, RAMP[4]), row, parts, OUT)

    # a classical monocyte: proteins sampled by abundance, the weighting the cell-surface panels use
    mono = pd.read_csv(a.run_dir / "tables" / "ravenhill_2020_monocyte_expression_dataset.csv")
    cols = [c for c in mono.columns if c.startswith("Classical monocytes, donation")]
    mono["abundance"] = mono[cols].mean(axis=1)
    mono = mono[mono["abundance"] > 0].merge(heights[["ID link_first", "ID", "total_height", "Entry name_first"]],
                                             left_on="ID link", right_on="ID link_first", suffixes=("_m", ""))
    rng = np.random.default_rng(a.seed)
    picks = rng.choice(len(mono), size=a.n, replace=True, p=(mono["abundance"] / mono["abundance"].sum()).to_numpy())
    objs, cache = [], {}
    for i in picks:
        r = mono.iloc[i]
        acc = r["ID link"]
        if acc not in cache:
            it = protein(acc, heights, parts, gff, a.alphafold)
            cache[acc] = build(it, a.alphafold, color_for(it["total"], lo, hi))
        objs.append((str(r["Gene names  (primary )"]).split(";")[0], cache[acc]))
    scene_figure(objs, OUT)

    # a monocyte against a HER2+ cancer cell with anti-HER2: the pairs, and the bullseye
    contact_figure(contact_pairs(a.run_dir, heights, parts, a.alphafold), heights, parts, gff, a.alphafold, OUT)
    if not a.igg.is_file():
        import urllib.request
        a.igg.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(IGG_PDB, a.igg)
    bridge_figure(a.run_dir, heights, parts, gff, a.alphafold, a.igg, OUT)
    bullseye_figure(a.run_dir, OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
