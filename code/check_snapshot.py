#!/usr/bin/env python
"""Verify the committed data against its manifest, and that the 2026 proteome is what a run reads.

    python code/check_snapshot.py [--bulk-only] [--quick]

``data/snapshot/MANIFEST.json`` records a sha256 for every committed input and every table
the last full run produced, plus the code version that produced them.  Three things are
checked:

* every listed file is present with the recorded hash;
* the 2026 proteome under ``data/inputs/proteome/`` is what ``resolve()`` returns for it
  (unless ``SURFACEOME_BASELINE=1`` asks for the 2022 vintage).  The surfaceome tables are in
  ``GENERATED``, so a run rebuilds them into ``OUT_ROOT`` and every downstream notebook reads
  those: a 2022 proteome read here would quietly turn every height into a mixture of
  vintages that no label distinguishes;
* the snapshot is still what the current code produces (``check_stale``): the manifest
  records the git hash of each part of ``code/`` that decides the numbers, and a change
  since the snapshot was built means a run will not reproduce the committed tables.

``--bulk-only`` answers just the vintage question; ``run_notebooks.py`` calls it before
executing anything.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

from build_snapshot import DECIDES_OUTPUT
from surfaceome_config import (DATA, GENERATED, INPUTS, REFRESH_ROOT, REPO_ROOT, SNAPSHOT,
                               USE_BASELINE, resolve)

MANIFEST = SNAPSHOT / "MANIFEST.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load() -> dict:
    if not MANIFEST.is_file():
        raise SystemExit(f"no manifest at {MANIFEST}")
    return json.loads(MANIFEST.read_text())


#: The proteome files whose vintage decides every downstream table.
PROTEOME = [
    "inputs/proteome/UP000005640.gff",
    "inputs/proteome/UP000005640.fasta",
    "inputs/proteome/UP000005640.tab",
]


def check_bulk(quiet: bool = False) -> list[str]:
    """Problems with the proteome vintage; empty means the 2026 proteome resolves."""
    if USE_BASELINE:
        return []
    problems = []
    for rel in PROTEOME:
        got = resolve(rel)
        if not got.is_file():
            problems.append(f"MISSING   {rel}")
        elif not (str(got).startswith(str(INPUTS))
                  or (REFRESH_ROOT and str(got).startswith(str(REFRESH_ROOT)))):
            problems.append(f"NOT 2026  {rel}\n              resolves to {got}")
    return problems


def check_stale() -> list[str]:
    """Is the snapshot still what current code produces?

    Documentation, the runner and the checkers are outside the compared set on purpose:
    none of them can change an output, and a check that demanded a full rerun to fix a typo
    in a checker is one people learn to ignore.

    A tree without a snapshot (the public release, which ships no result tables) has nothing
    to be stale, so there is nothing to report.
    """
    if not MANIFEST.is_file():
        return []
    man = load()
    prov = man.get("provenance")
    if not prov:
        return ["NO PROVENANCE  the manifest records no code version, so staleness cannot be judged"]

    def git(*args):
        return subprocess.run(["git", "-C", str(REPO_ROOT), *args],
                              capture_output=True, text=True, check=False).stdout.strip()
    was = prov.get("code_trees", {})
    moved = [p for p in DECIDES_OUTPUT
             if git("rev-parse", f"HEAD:{p}") and was.get(p) and git("rev-parse", f"HEAD:{p}") != was[p]]
    problems = []
    if moved:
        problems.append(f"STALE     {', '.join(moved)} changed since the snapshot was built")
        problems.append(f"          built {prov.get('built_utc', '?')} at commit "
                        f"{prov.get('code_commit', '?')[:12]} ({prov.get('code_commit_subject', '')[:60]})")
        problems.append("          rerun the pipeline and: python code/build_snapshot.py RUN_DIR")
    elif git("status", "--porcelain", "--", *DECIDES_OUTPUT):
        problems.append("UNCOMMITTED  the notebooks or package have uncommitted changes; the "
                        "snapshot matches the last commit, not the working tree")
    return problems


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bulk-only", action="store_true", help="only check which proteome vintage resolves")
    ap.add_argument("--quick", action="store_true", help="compare sizes, not hashes")
    a = ap.parse_args()

    man = load()
    problems = []
    if not a.bulk_only:
        for rel, want in man["files"].items():
            p = DATA / rel
            if not p.is_file():
                problems.append(f"MISSING   {rel}")
            elif p.stat().st_size != want["bytes"]:
                problems.append(f"SIZE      {rel}")
            elif not a.quick and sha256(p) != want["sha256"]:
                problems.append(f"CHANGED   {rel}")
        print(f"  {len(man['files'])} committed files checked against the manifest")
        # Every table the pipeline generates must ship, or the figure stage cannot run from
        # the snapshot alone.  The manifest is the list of what ships, so a generated table
        # with no entry is the builder's omission, not a missing file.
        for rel in sorted(GENERATED):
            if f"snapshot/{rel}" not in man["files"]:
                problems.append(f"NOT SHIPPED  snapshot/{rel}: generated but not in the manifest; "
                                f"rerun build_snapshot.py")

    bulk = check_bulk(quiet=a.quick)
    print(f"  {len(PROTEOME)} proteome file(s) checked for vintage")
    stale = [] if a.bulk_only else check_stale()
    if not a.bulk_only:
        prov = man.get("provenance", {})
        print(f"  built {prov.get('built_utc', '?')} at {prov.get('code_commit', '?')[:12]}")
    for msg in problems + bulk + stale:
        print(f"    {msg}")
    total = len(problems) + len(bulk) + len(stale)
    if bulk:
        print("\n  The 2026 proteome is not what resolves; a run would rebuild a 2022 surfaceome")
        print("  that overrides the committed 2026 tables downstream.  Restore data/inputs/proteome/,")
        print("  or set SURFACEOME_BASELINE=1 if a 2022 run is what you want.")
    print(f"\n  {total} problem(s)")
    return 1 if total else 0


if __name__ == "__main__":
    sys.exit(main())
