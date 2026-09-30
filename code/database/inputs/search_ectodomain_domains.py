#!/usr/bin/env python
"""Search the Pfam library over the ectodomain sequences: the table domain counting starts from.

    python code/database/inputs/search_ectodomain_domains.py --fasta FASTA --hmm Pfam-A.hmm --out DOMTBL

Runs ``hmmsearch --domtblout`` of every Pfam-A profile against the ectodomain FASTA that
``build_surfaceome.py`` writes, so a Pfam refresh needs no hand step. ``build_domain_assignments.py``
turns the table into the assignments the height notebook reads. HMMER's binaries are looked for
beside the interpreter (the conda environment) and then on PATH.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def _tool(name: str) -> str:
    beside = Path(sys.executable).with_name(name)
    if beside.is_file():
        return str(beside)
    found = shutil.which(name)
    if not found:
        raise SystemExit(f"{name} not found beside {sys.executable} or on PATH; install HMMER")
    return found


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", type=Path, required=True, help="ectodomain sequences (build_surfaceome.py)")
    ap.add_argument("--hmm", type=Path, required=True, help="Pfam-A.hmm")
    ap.add_argument("--out", type=Path, required=True, help="the --domtblout table")
    ap.add_argument("--cpu", type=int, default=8)
    a = ap.parse_args()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    cmd = [_tool("hmmsearch"), "--domtblout", str(a.out), "--noali", "-o", "/dev/null",
           "--cpu", str(a.cpu), str(a.hmm), str(a.fasta)]
    print("  " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)
    n = sum(1 for l in open(a.out) if not l.startswith("#"))
    print(f"wrote {a.out}: {n:,} domain hits")
    return 0


if __name__ == "__main__":
    sys.exit(main())
