#!/usr/bin/env python
"""Assemble a bound pair from AlphaFold ectodomains placed on a structure of the complex.

    python code/figures/validation/build_interaction_models.py --templates CSV --models DIR \\
        --af-confidence CSV --structures DIR --out-dir DIR

The 2018 assembled models (``data/README.md``) put full ectodomain models into a
crystal structure of the bound pair by hand in PyMOL, mostly on mouse complexes. This does
the same step scripted, on the human AlphaFold v6 ectodomains and a human structure of the
complex, so the model can be rebuilt and cited:

1. load the template entry and the two AlphaFold models;
2. cut each model to its confident ectodomain span (``af_confidence.csv``: ``ecd_start`` plus
   ``n_trim`` to ``ecd_start`` plus ``c_trim``, the span the height measurement uses);
3. superpose each ectodomain on its partner's chain in the template with PyMOL ``super``
   (sequence-independent, so a template fragment such as an integrin I domain still places
   the whole ectodomain), and keep the RMSD and atom count of the fit;
4. box the union of the two placed ectodomains with the inertia-axis bounding box used for
   every other measurement (``measure_structures.iabb``): ``dim_1`` is the model's length,
   which ``interaction_height_rules.csv`` records as the pair's height.

``--templates`` has one row per pair: ``pair``, ``uniprot_1``, ``uniprot_2``, ``pdb``,
``chain_1``, ``chain_2``, ``pmid``. The placed model is written as ``<pair>.pdb.gz`` beside a
table of fits and dimensions. A poor fit (RMSD above a few angstrom, or few atoms aligned)
means the template does not hold that partner's ectodomain and the row should be re-examined
rather than used.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "database" / "inputs"))
from measure_structures import _setup, iabb  # noqa: E402


def ecd_span(conf: pd.DataFrame, acc: str) -> tuple[int, int]:
    r = conf[conf["ID"] == acc]
    if r.empty:
        raise KeyError(f"{acc} not in af_confidence")
    r = r.iloc[0]
    return int(r["ecd_start"] + r["n_trim"]), int(r["ecd_start"] + r["c_trim"])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--templates", type=Path, required=True)
    ap.add_argument("--models", type=Path, required=True, help="AlphaFold models, AF-<acc>-F1-model_v6.pdb")
    ap.add_argument("--af-confidence", type=Path, required=True)
    ap.add_argument("--structures", type=Path, required=True, help="template mmCIFs, <pdbid>.cif")
    ap.add_argument("--out-dir", type=Path, required=True)
    a = ap.parse_args()

    conf = pd.read_csv(a.af_confidence, low_memory=False)
    tpl = pd.read_csv(a.templates, dtype=str)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    cmd = _setup()
    rows = []
    for t in tpl.itertuples():
        cmd.delete("all")
        cmd.load(str(a.structures / f"{t.pdb.lower()}.cif"), "tpl")
        cmd.remove("tpl and not polymer")
        fits = []
        parts = []
        for k, (acc, chain) in enumerate(((t.uniprot_1, t.chain_1), (t.uniprot_2, t.chain_2)), 1):
            model = a.models / f"AF-{acc}-F1-model_v6.pdb"
            if not model.is_file():
                print(f"  {t.pair}: no AlphaFold model for {acc}", file=sys.stderr)
                break
            s, e = ecd_span(conf, acc)
            name = f"p{k}"
            cmd.load(str(model), name)
            cmd.remove(f"{name} and not (resi {s}-{e})")
            target = f"tpl and chain {chain} and polymer"
            try:
                rms, n_atoms, *_ = cmd.super(name, target)
            except Exception as exc:  # no alignment
                print(f"  {t.pair}: super failed for {acc} on chain {chain}: {exc}", file=sys.stderr)
                break
            fits.append((acc, chain, s, e, round(rms, 2), int(n_atoms)))
            parts.append(name)
        if len(parts) < 2:
            continue
        cmd.create("assembly", " or ".join(parts))
        n_atoms = cmd.count_atoms("assembly")
        dims = iabb(cmd, "assembly")
        out = a.out_dir / f"{t.pair.replace('/', '_')}.pdb"
        cmd.save(str(out), "assembly")
        with open(out, "rb") as fh, gzip.open(str(out) + ".gz", "wb") as gz:
            gz.write(fh.read())
        out.unlink()
        (acc1, ch1, s1, e1, r1, n1), (acc2, ch2, s2, e2, r2, n2) = fits
        rows.append({"pair": t.pair, "pdb": t.pdb.upper(), "pmid": getattr(t, "pmid", ""),
                     "uniprot_1": acc1, "chain_1": ch1, "ecd_1": f"{s1}-{e1}", "rmsd_1": r1, "atoms_1": n1,
                     "uniprot_2": acc2, "chain_2": ch2, "ecd_2": f"{s2}-{e2}", "rmsd_2": r2, "atoms_2": n2,
                     "n_atoms": n_atoms, "dim_1": round(dims[0], 2), "dim_2": round(dims[1], 2),
                     "dim_3": round(dims[2], 2), "length_nm": round(dims[0] / 10, 2),
                     "model": out.name + ".gz"})
        print(f"  {t.pair:<16} {t.pdb}  fit {r1} A/{n1} atoms, {r2} A/{n2} atoms  ->  {dims[0]/10:.1f} nm", flush=True)
    table = pd.DataFrame(rows)
    table.to_csv(a.out_dir / "interaction_models_2026.csv", index=False)
    print(f"wrote {len(table)} models to {a.out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
