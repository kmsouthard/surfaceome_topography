#!/usr/bin/env python
"""Refresh the MobiDB disorder annotations the pipeline reads.

    python code/database/inputs/build_mobidb_annotations.py --validate --archive DIR
    python code/database/inputs/build_mobidb_annotations.py --surfaceome CSV --out-dir DIR

The pipeline reads two files per proteome, in a format the retired MobiDB consensus API
produced in 2019:

    <proteome>_disorder_full.txt   "mobidb_consensus.disorder.full.full.regions"
    <proteome>_disorder_lite.txt   "mobidb_consensus.disorder.predictors.mobidb-lite.regions"

Tab separated, every field quoted, the third holding a region list:

    "P25063"\t"mobidb_consensus...mobidb-lite.regions"\t"[[28,58,""D_WC""]]"

`mobidb_annotations.import_disorder` slices that string with `x[2:-2]` and splits on `],[`,
then `parse_disorder` splits each region on commas and keeps those whose third field matches
`D|d`. So the format is load-bearing in a way a schema would not suggest, and this script
writes it rather than inventing a new one — the notebooks are unchanged.

`download_mobidb.to_frame()` deliberately does **not** produce this: its own docstring warns
that its output corresponds to post-filter data and must not be fed to `import_disorder`.
This script bridges that gap.

Track mapping, following `download_mobidb.MOBIDB_KEYS`:

    lite  <- prediction-disorder-mobidb_lite    a direct equivalent
    full  <- prediction-disorder-priority       the closest current key

The legacy `full` track annotated *every* residue with a structural state — `D`/`d`
disordered, `S`/`s` structured, `C`/`c` context-dependent — and `parse_disorder` then kept
only `D|d`. The current key returns the disordered regions directly, so they are written with
the code `D` and the structured states are simply absent. Downstream that is equivalent:
`parse_disorder` would have discarded them.

`--validate` proves the format round-trips before any refreshed data is written: it re-emits
the regions parsed out of the archived files and checks that reading them back reproduces the
archived parse exactly, protein for protein and region for region.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # code/, for surfaceome_config and the package
from surfaceomeTopography.mobidb_annotations import import_disorder, parse_disorder  # noqa: E402

#: Legacy consensus key written into column 2 of each file.
MOBI_KEY = {
    "full": "mobidb_consensus.disorder.full.full.regions",
    "lite": "mobidb_consensus.disorder.predictors.mobidb-lite.regions",
}

#: Current MobiDB API key backing each legacy track.
API_KEY = {
    "full": "prediction-disorder-priority",
    "lite": "prediction-disorder-mobidb_lite",
}


def write_track(path: Path, regions: dict[str, list[tuple[int, int, str]]], method: str) -> int:
    """Write one legacy-format track file. Returns the number of proteins written."""
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", quoting=csv.QUOTE_ALL)
        for acc in sorted(regions):
            rs = regions[acc]
            if not rs:
                continue
            body = ",".join(f'[{a},{b},"{code}"]' for a, b, code in rs)
            w.writerow([acc, MOBI_KEY[method], f"[{body}]"])
            n += 1
    return n


def parsed_regions(directory: str, proteome: str, method: str) -> dict[str, list[tuple]]:
    """{accession: [(start, end, code)]} as the pipeline actually parses them."""
    df = parse_disorder(import_disorder(directory, proteome, method=method))
    out: dict[str, list[tuple]] = {}
    for r in df.itertuples():
        out.setdefault(r._1 if hasattr(r, "_1") else getattr(r, "ID link", None), [])
    out = {}
    for acc, start, end, typ in zip(df["ID link"], df.disorder_start,
                                    df.disorder_end, df.disorder_type):
        out.setdefault(acc, []).append((int(start), int(end), str(typ)))
    return out


def validate(archive_dir: str, proteome: str, tmp: Path) -> int:
    problems = 0
    for method in ("full", "lite"):
        ref = parsed_regions(archive_dir, proteome, method)
        n = write_track(tmp / f"{proteome}_disorder_{method}.txt", ref, method)
        back = parsed_regions(str(tmp) + "/", proteome, method)
        same = ref == back
        print(f"  {method:<5} {len(ref):>5,} proteins re-emitted ({n:,} written), "
              f"round-trip {'IDENTICAL' if same else 'DIFFERS'}")
        if not same:
            problems += 1
            only_ref = set(ref) - set(back)
            only_back = set(back) - set(ref)
            diff = [a for a in set(ref) & set(back) if ref[a] != back[a]]
            print(f"        only before {len(only_ref)}, only after {len(only_back)}, "
                  f"differing {len(diff)}  e.g. {diff[:3]}")
    return problems


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--proteome", default="UP000005640")
    ap.add_argument("--validate", action="store_true",
                    help="round-trip the archived files through the writer and stop")
    ap.add_argument("--archive", help="directory holding the archived *_disorder_*.txt")
    ap.add_argument("--surfaceome", help="surfaceome CSV; its 'ID link' column is fetched")
    ap.add_argument("--out-dir", type=Path)
    ap.add_argument("--cache", type=Path, help="cache the raw MobiDB pull here")
    a = ap.parse_args()

    if a.validate:
        if not a.archive:
            print("--validate needs --archive", file=sys.stderr)
            return 2
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            print("format round-trip against the archived files:")
            return 1 if validate(a.archive, a.proteome, Path(td)) else 0

    if not (a.surfaceome and a.out_dir):
        print("need --surfaceome and --out-dir", file=sys.stderr)
        return 2

    from surfaceomeTopography.download_mobidb import fetch

    accs = sorted(pd.read_csv(a.surfaceome, low_memory=False)["ID link"].dropna().unique())
    print(f"fetching MobiDB for {len(accs):,} surfaceome accessions")

    records = []
    for rec in fetch(accs):
        records.append(rec)
    print(f"  {len(records):,} records returned")
    if a.cache:
        a.cache.parent.mkdir(parents=True, exist_ok=True)
        pd.Series([str(r) for r in records]).to_csv(a.cache, index=False, header=False)

    for method in ("full", "lite"):
        key = API_KEY[method]
        regions: dict[str, list[tuple]] = {}
        for rec in records:
            acc = rec.get("acc")
            block = rec.get(key) or {}
            for r in block.get("regions", []):
                if len(r) >= 2:
                    regions.setdefault(acc, []).append((int(r[0]), int(r[1]), "D"))
        path = a.out_dir / f"{a.proteome}_disorder_{method}.txt"
        n = write_track(path, regions, method)
        total = sum(len(v) for v in regions.values())
        print(f"  {method:<5} {n:,} proteins, {total:,} regions -> {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
