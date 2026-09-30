#!/usr/bin/env python
"""Draw the README's CellScape illustrations from a run's output.

    conda env create -f code/figures/cellscape/environment.yml
    conda activate cellscape
    python code/figures/cellscape/readme_images.py RUN_DIR --alphafold DIR

Two images, each in a light and a dark version, under ``docs/img/``:

* ``cd45_{light,dark}.svg``: CD45's ectodomain standing on the membrane, with each piece of
  its height marked;
* ``monocyte_surface_{light,dark}.svg``: twenty surface proteins sampled from a classical
  monocyte's proteome by abundance, among those with an AlphaFold model, each drawn to scale.
* ``contact_{light,dark}.svg``: a monocyte against a HER2-positive breast cancer cell with an
  anti-HER2 antibody (the Figure 6 panel ``monocyte_aHER2``): for each gap band of the bullseye,
  the pair of single proteins holding most of that band among pairs whose height is the sum of
  their partners', the monocyte's partner hanging from its membrane, the gap between the
  membranes the pair's height;
* ``bullseye_{light,dark}.svg``: that contact's bullseye, drawn by the pipeline's own
  ``surfaceomeTopography.interface.plot_contact_bullseye`` from the panel's pairs table.

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
PANEL = "monocyte_aHER2"
THEMES = {"light": {"text": "#0b0b0b", "text2": "#52514e", "rule": "#8a8984"},
          "dark": {"text": "#ffffff", "text2": "#c3c2b7", "rule": "#8f8e87"}}
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


def scene_figure(objs: list[tuple[str, dict]], out: Path, padding: float = 14.0):
    width = sum(o["width"] for _, o in objs) + padding * (len(objs) + 1)
    top = max(o["height"] for _, o in objs)
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
    sum of their partners' and whose partners both have AlphaFold models."""
    pairs = pd.read_csv(run / "tables" / f"F6_{PANEL}" / f"{PANEL}_pairs.csv")
    pairs["share"] = pairs["interface_shared"] / pairs["interface_shared"].sum()
    single = pairs["complex_name_immune"].isna() & pairs["complex_name_cancer"].isna()
    additive = (pairs["interaction_dim"] - pairs["total_height_prot1"] - pairs["total_height_prot2"]).abs() < 0.01
    has_model = set(parts.loc[parts["methods"] == "alphafold", "ID"])
    ecd = heights.set_index("ID link_first")["ID"]
    modelled = lambda acc: ecd.get(acc) in has_model and (alphafold / f"AF-{acc}-F1-model_v6.pdb").is_file()
    ok = pairs[single & additive
               & pairs["Entry_immune"].map(modelled) & pairs["Entry_cancer"].map(modelled)].copy()
    edges = [0, 10, 20, 30, 40, 50, np.inf]
    ok["band"] = pd.cut(ok["interaction_dim"], edges, right=False)
    best = ok.sort_values("share", ascending=False).groupby("band", observed=True).head(1)
    return best.sort_values("interaction_dim")


def contact_figure(chosen: pd.DataFrame, heights: pd.DataFrame, parts: pd.DataFrame, gff: pd.DataFrame,
                   alphafold: Path, out: Path, padding: float = 40.0, margin: float = 230.0):
    row = lambda acc: heights[heights["ID link_first"] == acc].iloc[0]
    built = []
    for _, r in chosen.iterrows():
        gap = float(r["interaction_dim"]) * NM
        im = build(layout(row(r["Entry_immune"]), parts, gff), alphafold, IMMUNE_COLOR)
        ca = build(layout(row(r["Entry_cancer"]), parts, gff), alphafold, CANCER_COLOR)
        drawn = im["height"] + ca["height"]              # the partners meet; the label gives the table's height
        built.append((r, drawn, flipped(im, drawn), ca))
    top = max(g for _, g, _, _ in built)
    width = margin + sum(max(i["width"], c["width"]) + padding for _, _, i, c in built)
    for mode, t in THEMES.items():
        fig, ax = plt.subplots(figsize=(7, 7 * (top + 260) / width))
        x = margin
        for r, gap, im, ca in built:
            w = max(im["width"], ca["width"])
            bilayer(ax, x - padding / 2, x + w + padding / 2, gap + 40)          # the monocyte's
            draw(ax, im, x + (w - im["width"]) / 2)
            draw(ax, ca, x + (w - ca["width"]) / 2)
            ax.text(x + w / 2, gap + 58, r["Gene Name_immune"], color=t["text"], fontsize=9,
                    ha="center", va="bottom")
            ax.text(x + w / 2, -58, r["Gene Name_cancer"], color=t["text"], fontsize=9, ha="center", va="top")
            ax.text(x + w / 2, -84, f"{r['interaction_dim']:.0f} nm", color=t["text2"], fontsize=8.5,
                    ha="center", va="top")
            x += w + padding
        bilayer(ax, margin - padding / 2, x - padding / 2, 0)                    # the cancer cell's
        ax.text(margin - padding, built[0][1] + 20, "monocyte", color=t["text"], fontsize=10,
                fontweight="bold", ha="right", va="center")
        ax.text(margin - padding, -20, "HER2+ breast\ncancer cell", color=t["text"], fontsize=10,
                fontweight="bold", ha="right", va="center")
        ax.set_xlim(0, x)
        ax.set_ylim(-130, top + 110)
        ax.set_aspect("equal")
        ax.axis("off")
        path = out / f"contact_{mode}.svg"
        fig.savefig(path, transparent=True, bbox_inches="tight", pad_inches=0.05, metadata={"Date": None})
        plt.close(fig)
        print(f"  wrote {path.relative_to(REPO)}")


def bullseye_figure(run: Path, out: Path):
    """The pipeline's own bullseye for the panel; only its text colour changes between modes."""
    spec = importlib.util.spec_from_file_location("sT_interface", REPO / "code/surfaceomeTopography/interface.py")
    interface = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(interface)
    pairs = pd.read_csv(run / "tables" / f"F6_{PANEL}" / f"{PANEL}_pairs.csv")
    for mode, t in THEMES.items():
        fig = interface.plot_contact_bullseye(pairs)
        for text in fig.axes[0].get_legend().get_texts():
            text.set_color(t["text"])
            text.set_fontsize(10)
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
    # a protein with no AlphaFold model (too long to predict, e.g. LRP1) would be a bare capsule
    has_model = set(parts.loc[parts["methods"] == "alphafold", "ID"])
    mono = mono[mono["ID"].isin(has_model)
                & mono["ID link"].map(lambda acc: (a.alphafold / f"AF-{acc}-F1-model_v6.pdb").is_file())]
    rng = np.random.default_rng(a.seed)
    picks = rng.choice(len(mono), size=a.n, replace=True, p=(mono["abundance"] / mono["abundance"].sum()).to_numpy())
    objs, cache = [], {}
    for i in picks:
        r = mono.iloc[i]
        acc = r["ID link"]
        if acc not in cache:
            hr = heights[heights["ID link_first"] == acc].iloc[0]
            it = layout(hr, parts, gff)
            cache[acc] = build(it, a.alphafold, color_for(it["total"], lo, hi))
        objs.append((str(r["Gene names  (primary )"]).split(";")[0], cache[acc]))
    scene_figure(objs, OUT)

    # a monocyte against a HER2+ cancer cell with anti-HER2: the pairs, and the bullseye
    contact_figure(contact_pairs(a.run_dir, heights, parts, a.alphafold), heights, parts, gff, a.alphafold, OUT)
    bullseye_figure(a.run_dir, OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
