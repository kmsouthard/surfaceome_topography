#!/usr/bin/env python
"""Measure every Pfam domain instance on the selected solved-structure chains.

    python code/database/inputs/measure_domain_instances.py --instances CSV --structures DIR \\
        --hmm Pfam-A.hmm --work DIR --out CSV [--evalue 1e-3] [--limit N]

Replaces the 2020 run of ``measure_domain_sizes.py`` in PyMOL, which took each domain's
residue range from the RCSB search table. Here the range comes from HMMER: each family's
profile is fetched from the pressed Pfam release with ``hmmfetch`` and searched against the
sequences of the chains selected for it (``select_domain_instances.py``), and every hit with
an independent E-value under ``--evalue`` is one domain instance. That keeps the measurement
to the domain's own residues -- the 2020 Sushi measurement was contaminated by whole chains,
which cannot happen when the range is a profile hit.

Three passes, each cached under ``--work`` so an interrupted run resumes:

1. **sequences** -- load each structure once in PyMOL, read the observed residues of each
   selected chain (CA atoms, first alternate location, state 1) and write them per family as
   FASTA, with the residue identifiers kept beside the sequence so a hit maps back to the
   structure by index rather than by residue number. Insertion codes and gaps are therefore
   harmless.
2. **domains** -- ``hmmsearch --domtblout`` per family (or, with ``--all-families``, one search of the whole library over every chain, so every domain on a chain is found); ``env_from``/``env_to`` are 1-based
   indices into the observed sequence.
3. **measure** -- load each structure once more and box every instance with the same
   inertia-axis-aligned bounding box as every other measurement (``measure_structures.iabb``),
   reloading between instances because the box transforms the object in place. Only the
   first deposited model of a multi-model entry is boxed, and entries whose coordinates are
   not an atomic model (solution scattering) are skipped and listed.

The output has one row per instance: family, PDB, chain, the hit's residue span (structure
numbering and sequence index) and length, its E-value, the three box dimensions in angstrom,
largest first, and the instance's centroid and longest and shortest inertia axes in the
structure's frame, from which ``build_domain_heights.py`` orients each family by how its
consecutive copies stack along a chain.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from measure_structures import _setup, iabb  # noqa: E402


def principal_axes(cmd, selection: str):
    """Mass-weighted centroid and inertia axes of a selection in the structure's own frame.

    Returned as (centroid, longest-extent axis, shortest-extent axis), unit vectors; the
    longest extent lies along the smallest inertia eigenvalue. Used by build_domain_heights.py
    to orient a family from how consecutive domains stack along a chain (METHODS.md §2).
    """
    model = cmd.get_model(selection)
    m = np.array([a.get_mass() for a in model.atom])
    x = np.array([a.coord for a in model.atom])
    c = (x * m[:, None]).sum(0) / m.sum()
    d = x - c
    I = np.zeros((3, 3))
    I[0, 0] = (m * (d[:, 1] ** 2 + d[:, 2] ** 2)).sum()
    I[1, 1] = (m * (d[:, 0] ** 2 + d[:, 2] ** 2)).sum()
    I[2, 2] = (m * (d[:, 0] ** 2 + d[:, 1] ** 2)).sum()
    I[0, 1] = I[1, 0] = -(m * d[:, 0] * d[:, 1]).sum()
    I[0, 2] = I[2, 0] = -(m * d[:, 0] * d[:, 2]).sum()
    I[1, 2] = I[2, 1] = -(m * d[:, 1] * d[:, 2]).sum()
    val, vec = np.linalg.eigh(I)
    return c, vec[:, 0], vec[:, 2]

def _tool(name: str) -> str:
    """HMMER binaries: beside this interpreter (the conda env) or on PATH."""
    beside = Path(sys.executable).parent / name
    return str(beside) if beside.is_file() else (shutil.which(name) or name)


THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M",
    "SEC": "U", "PYL": "O",
}


def chain_residues(cmd, chain: str) -> list[tuple[str, str]]:
    """[(resi, resn)] of the observed polymer residues of a chain, in file order."""
    seen: list[tuple[str, str]] = []
    cmd.iterate_state(1, f"s and polymer and chain {chain} and name CA and alt ''+A",
                      "seen.append((resi, resn))", space={"seen": seen})  # noqa: F841
    return seen


def pass_sequences(cmd, inst: pd.DataFrame, structures: Path, work: Path) -> pd.DataFrame:
    """FASTA per family and a residue table per chain; returns the chains that loaded."""
    fasta_dir = work / "fasta"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    resi_path = work / "chain_residues.tsv"
    done: set[tuple[str, str]] = set()
    if resi_path.is_file():
        prev = pd.read_csv(resi_path, sep="\t", dtype=str)
        done = set(zip(prev.PDB, prev.CHAIN))
    per_family: dict[str, list[str]] = defaultdict(list)
    rows = []
    n_pdb = inst.PDB.nunique()
    for n, (pdb, grp) in enumerate(inst.groupby("PDB", sort=True), 1):
        chains = sorted(set(grp.CHAIN))
        if all((pdb, c) in done for c in chains):
            continue
        path = structures / f"{pdb.lower()}.cif"
        if not path.is_file():
            continue
        cmd.delete("all")
        try:
            cmd.load(str(path), "s")
        except Exception:
            continue
        for chain in chains:
            if (pdb, chain) in done:
                continue
            res = chain_residues(cmd, chain)
            if len(res) < 3:
                continue
            seq = "".join(THREE_TO_ONE.get(r, "X") for _, r in res)
            rows.append((pdb, chain, seq, ",".join(i for i, _ in res)))
            for fam in set(grp[grp.CHAIN == chain].PFAM_ID):
                per_family[fam].append(f">{pdb}_{chain}\n{seq}\n")
        if n % 200 == 0:
            print(f"  sequences: {n:,}/{n_pdb:,} entries", flush=True)
    if rows:
        new = pd.DataFrame(rows, columns=["PDB", "CHAIN", "seq", "resi"])
        new.to_csv(resi_path, sep="\t", index=False, mode="a", header=not resi_path.is_file())
    for fam, recs in per_family.items():
        with open(fasta_dir / f"{fam}.fasta", "a") as fh:
            fh.writelines(recs)
    return pd.read_csv(resi_path, sep="\t", dtype=str) if resi_path.is_file() else pd.DataFrame()


def family_names(hmm_dat: Path) -> dict[str, str]:
    """{PFxxxxx: family NAME} from Pfam-A.hmm.dat; hmmfetch keys on the name, not the accession."""
    names, ident = {}, None
    with open(hmm_dat) as fh:
        for line in fh:
            if line.startswith("#=GF ID"):
                ident = line.split(None, 2)[2].strip()
            elif line.startswith("#=GF AC") and ident:
                names[line.split(None, 2)[2].strip().split(".")[0]] = ident
    return names


def pass_domains(inst: pd.DataFrame, hmm: Path, hmm_dat: Path, work: Path,
                 evalue: float) -> pd.DataFrame:
    """hmmsearch per family over its chains' sequences -> one row per domain hit."""
    fasta_dir, hmm_dir, dom_dir = work / "fasta", work / "hmm", work / "domtbl"
    hmm_dir.mkdir(exist_ok=True)
    dom_dir.mkdir(exist_ok=True)
    names = family_names(hmm_dat)
    rows = []
    fams = sorted(set(inst.PFAM_ID))
    for k, fam in enumerate(fams, 1):
        fa = fasta_dir / f"{fam}.fasta"
        if not fa.is_file():
            continue
        prof = hmm_dir / f"{fam}.hmm"
        if not prof.is_file():
            if fam not in names:
                print(f"  {fam}: not in {hmm_dat.name}", file=sys.stderr)
                continue
            r = subprocess.run([_tool("hmmfetch"), str(hmm), names[fam]], capture_output=True, text=True)
            if r.returncode or not r.stdout:
                print(f"  {fam}: not in {hmm.name}", file=sys.stderr)
                continue
            prof.write_text(r.stdout)
        tbl = dom_dir / f"{fam}.domtbl"
        if not tbl.is_file():
            subprocess.run([_tool("hmmsearch"), "--domtblout", str(tbl), "--noali", "-o", "/dev/null",
                            "--cpu", "4", str(prof), str(fa)], check=True)
        for line in tbl.read_text().splitlines():
            if line.startswith("#"):
                continue
            f = line.split()
            # target, -, tlen, query, acc, qlen, E, score, bias, #, of, cE, iE, score, bias,
            # hmm_from, hmm_to, ali_from, ali_to, env_from, env_to, acc, desc
            ie = float(f[12])
            if ie > evalue:
                continue
            pdb, chain = f[0].rsplit("_", 1)
            rows.append((fam, pdb, chain, int(f[19]), int(f[20]), int(f[5]), ie, float(f[13])))
        if k % 100 == 0:
            print(f"  domains: {k}/{len(fams)} families, {len(rows):,} hits", flush=True)
    return pd.DataFrame(rows, columns=["PFAM_ID", "PDB", "CHAIN", "env_from", "env_to",
                                       "model_len", "i_evalue", "score"])


def pass_domains_all(residues: pd.DataFrame, hmm: Path, work: Path, evalue: float,
                     cpu: int = 8) -> pd.DataFrame:
    """One hmmsearch of the whole Pfam-A library over every chain -> every domain on every chain.

    Used with ``--all-families``: the neighbours of a measured domain, whatever their family,
    are what defines its stacking axis (``build_domain_heights.py --height-from-neighbour-pitch``).
    """
    fa = work / "all_chains.fasta"
    if not fa.is_file():
        with open(fa, "w") as fh:
            for r in residues.itertuples():
                fh.write(f">{r.PDB}_{r.CHAIN}\n{r.seq}\n")
    tbl = work / "all_chains.domtbl"
    if not tbl.is_file():
        print(f"  hmmsearch of {hmm.name} over {len(residues):,} chains ...", flush=True)
        subprocess.run([_tool("hmmsearch"), "--domtblout", str(tbl), "--noali", "-o", "/dev/null",
                        "--cpu", str(cpu), str(hmm), str(fa)], check=True)
    rows = []
    with open(tbl) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            ie = float(f[12])
            if ie > evalue:
                continue
            pdb, chain = f[0].rsplit("_", 1)
            fam = f[4].split(".")[0]
            rows.append((fam, pdb, chain, int(f[19]), int(f[20]), int(f[5]), ie, float(f[13])))
    return pd.DataFrame(rows, columns=["PFAM_ID", "PDB", "CHAIN", "env_from", "env_to",
                                       "model_len", "i_evalue", "score"])


NON_ATOMIC = ("SCATTERING", "THEORETICAL")


def experimental_method(path: Path) -> str:
    """``_exptl.method`` of an mmCIF entry (the first, if several)."""
    with open(path) as fh:
        for line in fh:
            if line.startswith("_exptl.method"):
                return line.split(None, 1)[1].strip().strip("'\"")
            if line.startswith("_atom_site."):
                break
    return ""


def atomic(method: str) -> bool:
    """Whether an entry's coordinates are an experimental atomic model.

    Solution-scattering entries are rigid-body or bead models fitted to SAXS curves, often
    deposited as several conformers; their domain arrangements are not measurements.
    """
    return not any(k in method.upper() for k in NON_ATOMIC)


def load_first_state(cmd, path: Path, name: str = "s") -> None:
    """Load an entry as a single-state object.

    NMR ensembles and multi-conformer entries load as several states; a box taken over all
    of them spans the ensemble's spread rather than one domain, so only the first deposited
    model is kept.
    """
    cmd.load(str(path), name)
    if cmd.count_states(name) > 1:
        cmd.create(name + "_1", name, 1, 1)
        cmd.delete(name)
        cmd.set_name(name + "_1", name)


def pass_measure(cmd, hits: pd.DataFrame, residues: pd.DataFrame, structures: Path,
                 work: Path) -> pd.DataFrame:
    """Box every hit; cached per PDB so a rerun only measures new entries."""
    out_path = work / "measurements.tsv"
    done = set()
    if out_path.is_file():
        done = set(pd.read_csv(out_path, sep="\t", dtype={"PDB": str})["PDB"])
    resi = {(r.PDB, r.CHAIN): r.resi.split(",") for r in residues.itertuples()}
    rows = []
    todo = hits[~hits.PDB.isin(done)]
    n_pdb = todo.PDB.nunique()
    skipped = []
    for n, (pdb, grp) in enumerate(todo.groupby("PDB", sort=True), 1):
        path = structures / f"{pdb.lower()}.cif"
        if not path.is_file():
            continue
        method = experimental_method(path)
        if not atomic(method):
            skipped.append((pdb, method))
            continue
        for h in grp.itertuples():
            ids = resi.get((pdb, h.CHAIN))
            if not ids:
                continue
            span = ids[h.env_from - 1:h.env_to]
            sel = f"s and polymer and chain {h.CHAIN} and resi " + "+".join(span)
            cmd.delete("all")
            try:
                load_first_state(cmd, path)
                n_res = cmd.count_atoms(f"{sel} and name CA and alt ''+A", state=1)
                if n_res >= 3:
                    cen, ax_long, ax_short = principal_axes(cmd, sel)
                dims = iabb(cmd, sel) if n_res >= 3 else None
            except Exception:
                dims = None
                n_res = 0
            if dims is None:
                continue
            rows.append((h.PFAM_ID, pdb, h.CHAIN, span[0], span[-1], h.env_from, h.env_to,
                         n_res, h.model_len, h.i_evalue, h.score, *dims,
                         *np.round(cen, 2), *np.round(ax_long, 4), *np.round(ax_short, 4)))
        if n % 100 == 0:
            print(f"  measure: {n:,}/{n_pdb:,} entries, {len(rows):,} instances", flush=True)
            _flush(rows, out_path)
            rows = []
    _flush(rows, out_path)
    if skipped:
        print(f"  skipped {len(skipped)} entries whose model is not atomic: "
              + ", ".join(f"{p} ({m})" for p, m in skipped), flush=True)
    return pd.read_csv(out_path, sep="\t", dtype={"PDB": str}) if out_path.is_file() else pd.DataFrame()


COLS = ["PFAM_ID", "PDB", "CHAIN", "resi_from", "resi_to", "seq_from", "seq_to", "n_res",
        "model_len", "i_evalue", "score", "max_dim", "mid_dim", "min_dim",
        "cx", "cy", "cz", "long_x", "long_y", "long_z", "short_x", "short_y", "short_z"]


def _flush(rows, out_path: Path) -> None:
    if rows:
        pd.DataFrame(rows, columns=COLS).to_csv(out_path, sep="\t", index=False, mode="a",
                                                 header=not out_path.is_file())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--instances", type=Path, required=True, help="from select_domain_instances.py")
    ap.add_argument("--structures", type=Path, required=True, help="directory of <pdbid>.cif")
    ap.add_argument("--hmm", type=Path, required=True, help="Pfam-A.hmm (pressed, for hmmfetch)")
    ap.add_argument("--hmm-dat", type=Path, help="Pfam-A.hmm.dat (default: beside --hmm), for accession to name")
    ap.add_argument("--work", type=Path, required=True, help="cache directory")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--evalue", type=float, default=1e-3, help="max independent E-value per domain hit")
    ap.add_argument("--limit", type=int, help="first N PDB entries only (a smoke test)")
    ap.add_argument("--all-families", action="store_true",
                    help="search the whole Pfam-A library over every selected chain and measure every hit, "
                         "not only the families each chain was selected for")
    a = ap.parse_args()

    inst = pd.read_csv(a.instances, dtype=str)
    have = {p.stem for p in a.structures.glob("*.cif")}
    inst = inst[inst.PDB.str.lower().isin(have)]
    if a.limit:
        inst = inst[inst.PDB.isin(sorted(set(inst.PDB))[:a.limit])]
    print(f"{len(inst):,} chains over {inst.PDB.nunique():,} entries on disk, "
          f"{inst.PFAM_ID.nunique()} families")
    a.work.mkdir(parents=True, exist_ok=True)

    cmd = _setup()
    residues = pass_sequences(cmd, inst, a.structures, a.work)
    print(f"  sequences for {len(residues):,} chains")
    hmm_dat = a.hmm_dat or a.hmm.with_suffix(".hmm.dat")
    if a.all_families:
        hits = pass_domains_all(residues, a.hmm, a.work, a.evalue)
    else:
        hits = pass_domains(inst, a.hmm, hmm_dat, a.work, a.evalue)
    print(f"  {len(hits):,} domain hits in {hits.PFAM_ID.nunique()} families")
    meas = pass_measure(cmd, hits, residues, a.structures, a.work)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    meas.to_csv(a.out, index=False)
    print(f"wrote {a.out}: {len(meas):,} instances, "
          f"{meas.PFAM_ID.nunique() if len(meas) else 0} families")
    return 0


if __name__ == "__main__":
    sys.exit(main())
