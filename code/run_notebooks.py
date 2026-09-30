#!/usr/bin/env python
"""Execute the manuscript pipeline notebooks headlessly, in dependency order.

    python code/run_notebooks.py [OUTPUT_DIR] [--stage database|figures|all] [--only SUBSTRING ...]

Tables land in OUTPUT_DIR/database/ (what ships) and OUTPUT_DIR/tables/ (intermediates and
per-panel tables), figures in OUTPUT_DIR/figures/, and executed copies of the notebooks (with
fresh outputs) in OUTPUT_DIR/notebooks/<stage>/; the notebooks in ``code/database/notebooks/``
and ``code/figures/notebooks/`` are left untouched.
Exit status is non-zero if any notebook failed.  ``--stage database`` builds the height and
interaction tables only; ``--stage figures`` draws the figures from them.

Why nbclient's async API rather than ``jupyter nbconvert --execute``: on the pinned
stack (jupyter_core 5.8.1 / Python 3.9) nbconvert intermittently aborts with
``RuntimeError: no running event loop`` before running a single cell.  It is a race in
the tooling, not a problem with the notebooks - roughly every other invocation failed.
Driving nbclient inside ``asyncio.run()`` guarantees a running loop and removes it.
"""

from __future__ import annotations

import argparse
import asyncio
import collections
import json
import os
import sys
import traceback
from pathlib import Path

import nbformat
from nbclient import NotebookClient
from nbclient.exceptions import CellExecutionError

import check_snapshot
from surfaceome_config import DATA, DATA_ROOT, GENERATED

REPO = Path(__file__).resolve().parent.parent
NBDIR = REPO / "code"   # notebook names below are relative to this

# The archive checks live under code/scripts/archive/ and are for the working repository:
# they compare the vendored data with the Box archive of record.  A release tree can drop
# that directory, and the pipeline runs without them.
sys.path.insert(0, str(REPO / "code" / "scripts" / "archive"))
try:
    import archive_guard
    import check_archive_subsets
except ImportError:  # release tree without the archive tooling
    archive_guard = None
    check_archive_subsets = None

#: The database stage: the three notebooks that build the height table and the interaction
#: table.  The height notebook produces the tables every other notebook reads, so it runs
#: first.
DATABASE = [
    "database/notebooks/01_height_estimates.ipynb",
    "database/notebooks/02_cellphonedb_interactions.ipynb",
    "database/notebooks/03_interaction_heights.ipynb",
]

#: The figure stage: reads the database stage's tables (from this run, or from the vendored
#: snapshot when run on its own) and writes figures and per-panel tables.
FIGURES = [
    "figures/notebooks/F1_domain_sizes.ipynb",
    "figures/notebooks/F2_height_histograms.ipynb",
    "figures/notebooks/F3_F4_expression_weighted_topography.ipynb",
    "figures/notebooks/S2_methods_vs_structures.ipynb",
    "figures/notebooks/F5_david_clusters.ipynb",
    "figures/notebooks/F4_ravenhill_monocyte_subtypes.ipynb",
    "figures/notebooks/F6_cell_contacts.ipynb",
]

#: Dependency order over both stages.  ``check_generated_paths.py`` enforces that no
#: notebook reads a table before the notebook that writes it has run.
ORDER = DATABASE + FIGURES

STAGES = {"database": DATABASE, "figures": FIGURES, "all": ORDER}


def prepare_output_tree(out_dir: Path) -> None:
    """Create the run's output directories: ``database/``, ``tables/``, ``figures/``.

    Every write goes through ``out_path()``, which creates parents as it goes, so this only
    makes the layout visible before the first notebook runs.
    """
    for sub in ("database", "tables", "figures", "notebooks"):
        (out_dir / sub).mkdir(parents=True, exist_ok=True)


def _first_error(exc: CellExecutionError) -> str:
    text = str(exc)
    for line in reversed(text.splitlines()):
        line = line.strip()
        if line and ("Error" in line or "Exception" in line) and ":" in line:
            return line[:110]
    return (text.splitlines()[-1][:110] if text else "?")


async def run_one(name: str, out_dir: Path, timeout: int) -> tuple[str, str]:
    nb = nbformat.read(NBDIR / name, as_version=4)
    # Run with cwd = the output directory: OUT_ROOT defaults to it, so every out_path()
    # lands there, and every read goes through surfaceome_config.resolve().
    prepare_output_tree(out_dir)
    client = NotebookClient(
        nb,
        timeout=timeout,
        kernel_name="python3",
        allow_errors=False,
        resources={"metadata": {"path": str(out_dir)}},
    )
    try:
        await client.async_execute()
        result = ("PASS", "")
    except CellExecutionError as exc:
        result = ("FAIL", _first_error(exc))
    except Exception as exc:  # kernel death, timeout, ...
        result = ("ERROR", f"{type(exc).__name__}: {exc}"[:110])
    executed = out_dir / "notebooks" / Path(name).parts[0] / Path(name).name   # out/notebooks/<stage>/<name>
    executed.parent.mkdir(parents=True, exist_ok=True)
    nbformat.write(nb, executed)
    return result


async def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("output_dir", nargs="?", default="nbrun")
    ap.add_argument("--stage", choices=sorted(STAGES), default="all",
                    help="database: build the height and interaction tables; figures: draw "
                         "the figures from them; all (default): both, in order")
    ap.add_argument("--only", nargs="+", default=None,
                    help="run only notebooks whose name contains one of these substrings")
    ap.add_argument("--timeout", type=int, default=2400, help="per-cell timeout, seconds")
    args = ap.parse_args()

    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    todo = STAGES[args.stage]
    if args.only:
        todo = [n for n in todo if any(s.lower() in n.lower() for s in args.only)]
        if not todo:
            print("no notebooks matched --only", file=sys.stderr)
            return 2

    # The archive is read-only. Notebook cells that build write paths by string
    # concatenation bypass out_path() and land in it; snapshot before and check after so
    # a new bypass is caught on the run that introduces it.
    # Refuse to run a mixed-vintage pipeline.  The surfaceome tables are in GENERATED, so a
    # run rebuilds them into OUT_ROOT and resolve() prefers OUT_ROOT over every stored copy --
    # so running without the 2026 proteome quietly rebuilds a 2022 surfaceome that then
    # overrides the vendored 2026 tables for every downstream notebook, with nothing to show
    # for it.  See code/check_snapshot.py.
    # Inputs whose licences forbid redistribution are not committed to the public release;
    # they are downloaded and rebuilt by fetch_restricted_inputs.py.  Without them the height
    # notebook would fail several minutes in, so say so now.
    restricted = DATA / "inputs" / "RESTRICTED.json"
    if restricted.is_file():
        absent = [rel for src in json.loads(restricted.read_text())["sources"].values()
                  for rel in src["files"] if not (DATA / rel).is_file()]
        if absent:
            print("=== REFUSING TO RUN: inputs that cannot be redistributed are not fetched ===",
                  file=sys.stderr)
            for rel in absent:
                print(f"  missing data/{rel}", file=sys.stderr)
            print("\n  Fetch and rebuild them (about 45 MB of downloads, checked by hash):",
                  file=sys.stderr)
            print("    python code/database/fetch_restricted_inputs.py", file=sys.stderr)
            return 2

    bulk = check_snapshot.check_bulk(quiet=True)
    if bulk:
        print("=== REFUSING TO RUN: the 2026 bulk inputs are not reachable ===", file=sys.stderr)
        for msg in bulk:
            print(f"  {msg}", file=sys.stderr)
        print("\n  The surfaceome build would produce a 2022 surfaceome that overrides the",
              file=sys.stderr)
        print("  vendored 2026 tables downstream, mixing vintages silently.", file=sys.stderr)
        print("  Point SURFACEOME_REFRESH at the 2026 overlay, or set SURFACEOME_BASELINE=1",
              file=sys.stderr)
        print("  to run the 2022 baseline deliberately.", file=sys.stderr)
        return 2

    # The GlycoMine subsets under data/inputs/glycosylation/glycomine/ stand in for 166 MB
    # of archive files; if they have been edited the glycan densities, and so every
    # disordered height, are quietly not what their manifest says they are.
    subsets = check_archive_subsets.check_integrity() if check_archive_subsets else []
    if subsets:
        print("=== REFUSING TO RUN: the GlycoMine subsets do not match their manifest ===",
              file=sys.stderr)
        for msg in subsets:
            print(f"  {msg}", file=sys.stderr)
        print("\n  Rebuild with: python code/scripts/archive/build_archive_subsets.py", file=sys.stderr)
        return 2

    guard = archive_guard.take() if (archive_guard and DATA_ROOT.is_dir()) else None

    # Record where every input actually came from.  The kernels inherit this, and
    # surfaceome_config.resolve() appends one line per resolution.
    trace = out_dir / "input_provenance.tsv"
    trace.unlink(missing_ok=True)
    os.environ["SURFACEOME_TRACE"] = str(trace)

    print(f"{'NOTEBOOK':<52} RESULT", flush=True)
    results = []
    for name in todo:
        try:
            status, detail = await run_one(name, out_dir, args.timeout)
        except Exception:
            status, detail = "ERROR", traceback.format_exc().splitlines()[-1][:110]
        results.append((name, status, detail))
        print(f"{name[:50]:<52} {status}  {detail}", flush=True)

    print("\n=== SUMMARY ===")
    for state in ("PASS", "FAIL", "ERROR"):
        print(f"  {state:<6} {sum(1 for _, s, _ in results if s == state)}")
    print(f"\nexecuted copies written to {out_dir}")

    if trace.is_file():
        counts = collections.Counter()
        paths = collections.defaultdict(set)
        for line in trace.read_text().splitlines():
            src, rel, _, *_ = line.split("\t") + [""]
            counts[src] += 1
            paths[src].add(rel)
        print("\n=== INPUT PROVENANCE ===")
        for src in ("generated", "refresh", "baseline", "snapshot", "committed", "archive", "missing"):
            if src in counts:
                print(f"  {src:<18} {counts[src]:>5} reads   {len(paths[src]):>4} distinct files")

        # An archive read is an input that is neither vendored nor hash-checked: it is
        # whatever that Box folder holds today, and a run elsewhere cannot reproduce it.
        # Name them rather than leaving a count, so the list stays short deliberately.
        if paths.get("archive"):
            print(f"\n  {len(paths['archive'])} input(s) read from the archive only -- not "
                  "vendored, not hash-checked:")
            for r in sorted(paths["archive"]):
                print(f"    {r}")
        stale = {r for r in paths.get("snapshot", set()) | paths.get("baseline", set()) | paths.get("archive", set())
                 if r in GENERATED}
        if stale:
            print(f"\n  WARNING: {len(stale)} pipeline-generated table(s) read from a stored")
            print("  copy rather than this run's output -- stage order or a failed upstream stage:")
            for r in sorted(stale):
                print(f"    {r}")

    # Coverage can only be judged against the surfaceome this run actually built, which is
    # why it is checked here rather than up front: a refreshed proteome adds accessions that
    # no stored table shows, and an accession outside the GlycoMine subsets loses its
    # glycosites silently.
    subsets_covered = True
    for msg in (check_archive_subsets.check_coverage() if check_archive_subsets else []):
        subsets_covered = False
        print(f"\n  {msg}")

    # Not fatal: this run wrote its own tables into OUT_ROOT and read those.  It matters for
    # anyone reading data/snapshot/ directly, or running one notebook on its own.
    stale = check_snapshot.check_stale()
    if stale:
        print("\n=== data/snapshot/ IS NOT WHAT THIS CODE PRODUCES ===")
        for msg in stale:
            print(f"  {msg}")

    archive_clean = True
    if guard is not None:
        modified, created, deleted = archive_guard.diff(guard, archive_guard.take())
        touched = modified + created + deleted
        if touched:
            archive_clean = False
            print(f"\n=== ARCHIVE WRITTEN TO: {len(touched)} file(s) ===")
            for k in touched[:20]:
                print(f"  {k}")
            if len(touched) > 20:
                print(f"  ... and {len(touched) - 20} more")
            print("  Generated files belong in OUT_ROOT; see code/scripts/archive/archive_guard.py.")

    ok = all(s == "PASS" for _, s, _ in results)
    return 0 if (ok and archive_clean and subsets_covered) else 1


if __name__ == "__main__":
    sys.exit(asyncio.run(main()))
