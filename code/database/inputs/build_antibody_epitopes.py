#!/usr/bin/env python
"""Each therapeutic antibody's epitope, read off a structure of it bound to its antigen.

    python code/database/inputs/build_antibody_epitopes.py --pdb DIR [--table CSV] [--check]

``data/curated/antibody_epitopes.csv`` names, per antibody, the PDB entry of the antibody bound
to its antigen and the paper that reports it. This fills in what follows from the entry: the
antigen chain, the two antibody chains on it, and the epitope, the antigen residues with a heavy
atom within ``cutoff_A`` of the antibody, in UniProt numbering. An antibody with no entry keeps
its row: with the residues a paper maps by other means and how (``basis``: ``peptide mapping``,
``domain``), or without residues (``none``). Only ``structure`` rows are derived here.

The numbering is the entry's own alignment of the antigen chain to its UniProt sequence
(``_struct_ref_seq``) where it has one. A peptide antigen without one is aligned to the UniProt
sequence here instead, and ``numbering`` says which was used. Where an entry holds several copies of
the complex, the antigen chain with the most contacts is read.

Entries are fetched from RCSB into ``--pdb`` as ``<id>.cif`` unless already there. ``--check``
rebuilds without writing and reports any row that differs from the committed table.
"""

from __future__ import annotations

import argparse
import gzip
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from scipy.spatial import cKDTree

REPO = Path(__file__).resolve().parents[3]
TABLE = REPO / "data/curated/antibody_epitopes.csv"
PROTEOME = REPO / "data/inputs/proteome/UP000005640.fasta.gz"
DERIVED = ["antigen_chain", "antibody_chains", "numbering", "n_residues", "residues"]
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
         "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
         "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def fetch(pdb: str, folder: Path) -> Path:
    path = folder / f"{pdb.lower()}.cif"
    if not path.exists():
        folder.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb.upper()}.cif", path)
    return path


def atoms(path: Path) -> pd.DataFrame:
    """The polymer heavy atoms of an entry's first model, with each one's UniProt accession and number."""
    cif = MMCIF2Dict(str(path))
    n = len(cif["_atom_site.id"])
    get = lambda name: cif.get(f"_atom_site.{name}", ["?"] * n)
    a = pd.DataFrame({"group": get("group_PDB"), "element": get("type_symbol"), "chain": get("auth_asym_id"),
                      "entity": get("label_entity_id"), "seq": get("label_seq_id"), "resname": get("label_comp_id"),
                      "model": get("pdbx_PDB_model_num"), "x": get("Cartn_x"), "y": get("Cartn_y"), "z": get("Cartn_z")})
    polymer = {e for e, t in zip(cif["_entity.id"], cif["_entity.type"]) if t == "polymer"}
    a = a[(a["group"] == "ATOM") & (a["model"] == a["model"].iloc[0]) & (a["element"] != "H") & a["entity"].isin(polymer)]
    a = a.astype({"x": float, "y": float, "z": float, "seq": int}).reset_index(drop=True)

    a["acc"], a["num"] = "", 0              # the entry's own alignment of each chain to its reference
    ref = lambda name: cif.get(f"_struct_ref_seq.{name}", [])
    for chain, acc, first, last, db_first in zip(ref("pdbx_strand_id"), ref("pdbx_db_accession"), ref("seq_align_beg"),
                                                 ref("seq_align_end"), ref("db_align_beg")):
        here = (a["chain"] == chain) & a["seq"].between(int(first), int(last))
        a.loc[here, "acc"] = acc
        a.loc[here, "num"] = a.loc[here, "seq"] - int(first) + int(db_first)
    return a


def uniprot_sequence(accession: str) -> str:
    keep, out = False, []
    with gzip.open(PROTEOME, "rt") as f:
        for line in f:
            if line.startswith(">"):
                keep = f"|{accession}|" in line
            elif keep:
                out.append(line.strip())
    return "".join(out)


def aligned_numbers(chain: pd.DataFrame, accession: str) -> dict:
    """UniProt number per ``seq`` of a chain with no SIFTS mapping, from a local alignment."""
    residues = chain.drop_duplicates("seq").sort_values("seq")
    aligner = PairwiseAligner()
    aligner.mode, aligner.open_gap_score, aligner.extend_gap_score, aligner.mismatch_score = "local", -5, -1, -1
    best = aligner.align("".join(THREE.get(r, "X") for r in residues["resname"]), uniprot_sequence(accession))[0]
    seqs = residues["seq"].to_numpy()
    return {int(seqs[i]): j + 1 for (a0, a1), (b0, _) in zip(*best.aligned) for i, j in zip(range(a0, a1), range(b0, b0 + a1 - a0))}


def epitope(path: Path, accession: str, cutoff: float, antigen_is_peptide: bool) -> dict:
    a = atoms(path)
    mapped = a[a["acc"] == accession]
    if mapped.empty and not antigen_is_peptide:
        raise SystemExit(f"{path.name}: no residue maps to {accession}")
    if mapped.empty:                        # the shortest polymer is the peptide
        sizes = a.groupby("chain")["seq"].nunique()
        antigen_chains, numbering = [sizes.idxmin()], "aligned"
    else:
        antigen_chains, numbering = list(mapped["chain"].unique()), "entry"
    carries = set(a.loc[a["acc"] == accession, "chain"]) | set(antigen_chains)
    others = a[~a["chain"].isin(carries)]
    tree = cKDTree(others[["x", "y", "z"]].to_numpy())

    best = None
    for chain in antigen_chains:
        own = a[a["chain"] == chain] if numbering == "aligned" else mapped[mapped["chain"] == chain]
        near = tree.query_ball_point(own[["x", "y", "z"]].to_numpy(), cutoff)
        partner = pd.Series([c for hit in near for c in set(others["chain"].to_numpy()[hit])]).value_counts()
        top = list(partner.index[:2])
        touching = [bool(set(others["chain"].to_numpy()[hit]) & set(top)) for hit in near]
        found = own[touching]
        if best is None or found["seq"].nunique() > best[1]["seq"].nunique():
            best = (chain, found, top, own)
    chain, found, top, own = best
    if numbering == "aligned":
        number = aligned_numbers(own, accession)
        numbers = sorted({number[s] for s in found["seq"] if s in number})
    else:
        numbers = sorted(set(found["num"].astype(int)))
    return {"antigen_chain": chain, "antibody_chains": ";".join(sorted(top)), "numbering": numbering,
            "n_residues": str(len(numbers)), "residues": ranges(numbers)}


def ranges(numbers: list) -> str:
    """``[1, 2, 3, 7]`` as ``1-3,7``."""
    out, start = [], None
    for i, n in enumerate(numbers):
        start = n if start is None else start
        if i + 1 == len(numbers) or numbers[i + 1] != n + 1:
            out.append(str(start) if start == n else f"{start}-{n}")
            start = None
    return ",".join(out)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb", required=True, type=Path, help="directory of mmCIF files, filled as needed")
    ap.add_argument("--table", type=Path, default=TABLE)
    ap.add_argument("--check", action="store_true", help="compare with the table instead of writing it")
    a = ap.parse_args()

    table = pd.read_csv(a.table, dtype=str, keep_default_na=False)
    built = table.copy()
    for i, r in table.iterrows():
        if r["basis"] != "structure":           # a mapped or domain-level epitope is recorded, not derived
            continue
        found = dict.fromkeys(DERIVED, "")
        if r["pdb"]:
            found = epitope(fetch(r["pdb"], a.pdb), r["accession"], float(r["cutoff_A"]), r["antigen_form"] == "peptide")
            print(f"  {r['antibody']:<14} {r['gene']:<7} {r['pdb']}  chain {found['antigen_chain']} against "
                  f"{found['antibody_chains']}  {found['n_residues']:>3} residues ({found['numbering']})  {found['residues']}")
        for column in DERIVED:
            built.at[i, column] = found[column]
    if a.check:
        differ = (built != table).any(axis=1)
        print(f"  {differ.sum()} row(s) differ from {a.table}")
        return int(differ.any())
    built.to_csv(a.table, index=False)
    print(f"  wrote {a.table}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
