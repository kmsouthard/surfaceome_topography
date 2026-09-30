#!/usr/bin/env python
"""Record a completed run in ``data/snapshot/``, with a manifest that goes stale loudly.

    python code/build_snapshot.py RUN_DIR [--dry-run]

``data/snapshot/`` holds the tables the last full run produced (``database/`` and the
``tables/`` that are read back), so one notebook can run on its own and a clone can be checked
against the code that made the numbers.  This copies each ``origin: run`` table out of a
finished run directory and rewrites ``MANIFEST.json``:

    built_utc        when
    code_commit      the commit, and whether the tree was dirty at the time
    code_trees       git hashes of the code that decides what the numbers are: the notebooks,
                     the package, and ``surfaceome_config.py``.  ``check_snapshot.py``
                     compares these to say the snapshot is out of date.  The runner and the
                     checkers are recorded but not compared: editing a checker cannot change
                     an output.
    run              the command, the output directory, and how many inputs the run read
                     from each tier
    databases        the upstream versions behind the run
    files[]          a sha256 per committed file under ``data/``, keyed by its data-relative
                     path; generated tables carry ``produced_by``, the notebook and cell that
                     wrote them, parsed by the same code that keeps ``GENERATED`` honest

Entries whose origin is not ``run`` are inputs: re-hashed and carried over, never replaced.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from file_hashes import content_sha256, sha256
from check_generated_paths import write_targets
from surfaceome_config import DATA, REPO_ROOT, SNAPSHOT, GENERATED

#: Upstream versions behind the committed run.  Update these with the inputs, not after.
DATABASES = {
    "uniprot": "2026_02",
    "alphafold": "v6",
    "pfam": "38.2",
    "mobidb": "2026 (prediction-disorder-mobidb_lite / -priority)",
    "glygen": "2.11.1",
    "glycomine": "2014, frozen",
    "sifts": "2026/09/08 (PDB 36.26, UniProt 2026.04); domain instances and validation structures",
    "string": "v12.0 links between surfaceome proteins at combined score >= 400; pair labels from data/curated/interaction_types.csv (2020 calls) and string_v12_pair_review.csv",
    "cellphonedb": "v5.0.0 (cellphonedb-data); 2022-02-25 snapshot under the baseline",
    "david": "knowledgebase v2025_1, clustering run 2026-09-16",
    "expression_atlas": "E-PROT-1 and E-PROT-27, re-downloaded 2026-09-16 (values unchanged, symbols current)",
    "domain_heights": "measured 2026-09 on current PDB structures: neighbour pitch per family, "
                      "clan averages by build_clan_heights.py; 2020 values only where no structure",
}


def git(*args: str) -> str:
    return subprocess.run(["git", "-C", str(REPO_ROOT), *args],
                          capture_output=True, text=True, check=False).stdout.strip()


#: The code whose content decides the outputs; ``check_snapshot.check_stale`` compares these.
DECIDES_OUTPUT = ("code/database/notebooks", "code/figures/notebooks", "code/surfaceomeTopography",
                  "code/surfaceome_config.py")

#: Recorded for the record, never compared.
RECORDED_ONLY = ("code/run_notebooks.py",)


def code_trees() -> dict[str, str]:
    return {p: git("rev-parse", f"HEAD:{p}") for p in DECIDES_OUTPUT + RECORDED_ONLY}


def provenance(run_dir: Path) -> dict:
    dirty = bool(git("status", "--porcelain", "--", *DECIDES_OUTPUT))
    prov = {
        "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "code_commit": git("rev-parse", "HEAD"),
        "code_commit_subject": git("log", "-1", "--format=%s"),
        "code_trees": code_trees(),
        "code_dirty": dirty,
        "databases": DATABASES,
        "run": {"output_dir": run_dir.name,
                "command": "python code/run_notebooks.py " + run_dir.name},
    }
    trace = run_dir / "input_provenance.tsv"
    if trace.is_file():
        tiers: dict[str, set] = {}
        for line in trace.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) >= 2:
                tiers.setdefault(parts[0], set()).add(parts[1])
        prov["run"]["inputs_by_tier"] = {k: len(v) for k, v in sorted(tiers.items())}
    return prov


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", type=Path, help="output directory of a completed run")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    run_dir = a.run_dir.resolve()
    if not run_dir.is_dir():
        print(f"no such run directory: {run_dir}", file=sys.stderr)
        return 2

    man = json.loads((SNAPSHOT / "MANIFEST.json").read_text())
    producers = write_targets()

    # The snapshot ships every table in GENERATED, not just the ones the manifest already
    # lists: a table added to the pipeline after the manifest was first written used to be
    # silently never shipped, and check_snapshot reported 0 problems while the figure stage
    # could not run on its own.  A run entry the pipeline no longer produces is dropped.
    # A run table lives under snapshot/ (its origin "run", or "run: <how>" for one promoted with
    # a note); a measured input elsewhere may carry a "run: <script>" origin too and is not one.
    is_run = lambda rel, e: rel.startswith("snapshot/") and str(e.get("origin", "")).split(":")[0].strip() == "run"
    for rel in sorted(GENERATED):
        man["files"].setdefault(f"snapshot/{rel}", {"origin": "run"})
    stale = [rel for rel, e in man["files"].items()
             if is_run(rel, e) and rel.split("/", 1)[1] not in GENERATED]
    for rel in stale:
        print(f"    NO LONGER GENERATED  {rel} (dropped from the manifest)")
        del man["files"][rel]

    updated, missing, kept = [], [], []
    for rel, entry in sorted(man["files"].items()):
        if not is_run(rel, entry):
            kept.append(rel)
            continue
        # a snapshot entry "snapshot/database/x.csv" is the run's "database/x.csv"
        src = run_dir / rel.split("/", 1)[1] if rel.startswith("snapshot/") else run_dir / rel
        if not src.is_file():
            missing.append(rel)
            continue
        dest = DATA / rel
        same = dest.is_file() and sha256(src) == sha256(dest)
        if not a.dry_run:
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, dest)
        updated.append((rel, same))

    for rel in [r for r in man["files"] if not (DATA / r).is_file()]:
        print(f"    GONE  {rel} (dropped from the manifest)")
        del man["files"][rel]
    for rel, entry in man["files"].items():
        p = DATA / rel
        entry["sha256"] = sha256(p)
        entry["bytes"] = p.stat().st_size
        if str(rel).endswith(".gz"):
            entry["content_sha256"] = content_sha256(p)
        out_rel = rel.split("/", 1)[1] if rel.startswith("snapshot/") else rel
        if out_rel in producers:
            entry["produced_by"] = producers[out_rel]

    changed = [r for r, same in updated if not same]
    print(f"  {len(updated)} output(s) taken from the run, {len(changed)} of them different")
    for rel in changed:
        print(f"    CHANGED  {rel}")
    for rel in missing:
        print(f"    NOT IN THE RUN  {rel}")
    print(f"  {len(kept)} input(s) carried over unchanged")

    if a.dry_run:
        print("\n  --dry-run: nothing written")
        return 1 if missing else 0

    man["provenance"] = provenance(run_dir)
    (SNAPSHOT / "MANIFEST.json").write_text(json.dumps(man, indent=1) + "\n")
    print(f"\n  manifest written at {man['provenance']['code_commit'][:12]}"
          f"{' (DIRTY TREE)' if man['provenance']['code_dirty'] else ''}")
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
