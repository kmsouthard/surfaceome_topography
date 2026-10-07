#!/usr/bin/env python
"""Every protein's CellScape cartoon, as outlines the web explorer draws.

    conda activate cellscape        # code/figures/cellscape/environment.yml
    python code/explorer/build_cartoons.py RUN_DIR --alphafold DIR [--dest explorer/data] [--jobs N]

The explorer draws a protein the way the README's illustrations do
(``code/figures/cellscape/readme_images.py``, whose layout and outline this reuses): the span
of its AlphaFold model the height uses, outlined with CellScape and oriented by its topology,
with a grey capsule for the height it takes from disorder, domains or sequence on each side.
An SVG per protein would run to hundreds of megabytes, so what is written is the outline
itself, and the page draws it in its own colours:

    cartoons/NNN.json   {accession: {w, h, before, structure, after, bx, tx, polys}}, in 256
                        files by a hash of the accession

Lengths are in angstroms, y upward from the membrane. ``polys`` are the model's
depth contours from back to front, each ``[shade, rings]``: ``shade`` 0 to 1 from the darkest
contour to the lightest, a ring a flat list of x, y. A protein with no model has no entry; the
page draws its whole height as a capsule.

As in the README, a structure's drawn extent is its extent along its N-to-C axis, not the
dimension the height table measures, so a drawn height is illustrative.
"""

from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "code" / "figures" / "cellscape"))

SHARDS = 256
BASE = "#808080"        # outlined in grey; the page maps each contour's shade onto its own colour
_STATE = {}


def shard(accession: str) -> int:
    """The cartoon file an accession is in; `app.js` computes the same."""
    n = 7
    for c in accession:
        n = (n * 31 + ord(c)) % 65521
    return n % SHARDS


def _start(run: str, alphafold: str) -> None:
    import readme_images as R
    heights, parts = R.load_run(Path(run))
    _STATE.update(R=R, heights=heights, parts=parts, gff=R.load_gff(), alphafold=Path(alphafold))


def _rings(geometry) -> list:
    polygons = getattr(geometry, "geoms", [geometry])
    out = []
    for polygon in polygons:
        if polygon.is_empty or polygon.geom_type != "Polygon":
            continue
        for ring in (polygon.exterior, *polygon.interiors):
            out.append([int(round(v)) for xy in ring.coords[:-1] for v in xy])
    return out


def cartoon(accession: str):
    from matplotlib.colors import to_rgb
    R = _STATE["R"]
    try:
        item = R.protein(accession, _STATE["heights"], _STATE["parts"], _STATE["gff"], _STATE["alphafold"])
        if item["af"] is None:
            return accession, None, ""
        obj = R.build(item, _STATE["alphafold"], BASE)
    except Exception as problem:  # one model CellScape cannot outline should not stop the rest
        return accession, None, f"{type(problem).__name__}: {problem}"
    structure = [p for p in obj["polygons"] if p["zorder"] == 1]
    light = [sum(to_rgb(p["facecolor"])) / 3 for p in structure]
    lo, hi = min(light), max(light)
    polys = [[round((l - lo) / (hi - lo), 2) if hi > lo else 0.5, _rings(p["polygon"].simplify(1.2))]
             for p, l in zip(structure, light)]
    (_, _, before), (_, _, top), (_, _, total) = obj["segments"]
    whole = lambda v: int(round(v))
    return accession, {"w": whole(obj["width"]), "h": whole(total), "before": whole(before),
                       "structure": whole(top - before), "after": whole(total - top),
                       "bx": whole(obj["bottom_x"]), "tx": whole(obj["top_x"]),
                       "polys": [p for p in polys if p[1]]}, ""


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", type=Path, help="a finished run's output directory")
    ap.add_argument("--alphafold", required=True, type=Path, help="directory of AlphaFold v6 models")
    ap.add_argument("--dest", type=Path, default=REPO / "explorer" / "data")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--only", nargs="+", metavar="ACCESSION", help="draw these, and print instead of writing")
    a = ap.parse_args()

    import pandas as pd
    accessions = a.only or list(pd.read_csv(a.run / "database" / "height_estimates.csv")["ID link_first"])
    shards, failed = [dict() for _ in range(SHARDS)], []
    with ProcessPoolExecutor(a.jobs, initializer=_start, initargs=(str(a.run), str(a.alphafold))) as pool:
        for i, (accession, drawn, problem) in enumerate(pool.map(cartoon, accessions, chunksize=8), 1):
            if drawn is not None:
                shards[shard(accession)][accession] = drawn
            if problem:
                failed.append((accession, problem))
            if i % 250 == 0:
                print(f"  {i:,} of {len(accessions):,}", flush=True)
    if a.only:
        for content in shards:
            for accession, drawn in content.items():
                print(accession, {k: v for k, v in drawn.items() if k != "polys"}, len(drawn["polys"]), "contours,",
                      sum(len(r) for _, rings in drawn["polys"] for r in rings) // 2, "points,",
                      len(json.dumps(drawn, separators=(",", ":"))), "bytes")
    else:
        folder = a.dest / "cartoons"
        folder.mkdir(parents=True, exist_ok=True)
        for i, content in enumerate(shards):
            (folder / f"{i:03d}.json").write_text(json.dumps(content, separators=(",", ":")))
        size = sum(f.stat().st_size for f in folder.glob("*.json"))
        print(f"  {sum(len(s) for s in shards):,} cartoons in {SHARDS} files, {size / 1e6:.1f} MB, under {folder}")
    for accession, problem in failed:
        print(f"  not drawn: {accession}: {problem}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
