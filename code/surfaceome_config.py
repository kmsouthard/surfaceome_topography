"""Single source of truth for every data path in the surfaceome topography pipeline.

Historically the notebooks hardcoded absolute paths from whichever laptop they were last run
on.  They have all been replaced by ``resolve()`` calls against a path relative to ``data/``::

    from surfaceome_config import resolve, out_path
    df = pd.read_csv(resolve('inputs/proteome/UP000005640.tab'))
    df.to_csv(out_path('database/height_estimates.csv'))

The committed tree under ``data/``:

    curated/         decision tables and the citation sheet, edited by hand
    inputs/          every external input a run reads, grouped by source (proteome, mapping,
                     disorder, glycosylation, domains, alphafold, structures, string,
                     cellphonedb, expression, david, phyre2); the 2026 vintage
    measurements/    tables measured here on structures: domain heights, PDB structure sizes,
                     the assembled interaction models
    snapshot/        the tables the last full run produced, with ``MANIFEST.json``, so one
                     notebook can run alone and a clone can be checked against the code
    baseline_2022/   the submitted 2022 vintage, at the same relative names, read under
                     ``SURFACEOME_BASELINE=1``; ``superseded/`` holds files nothing reads
    sources/         complete upstream files the derived inputs are rebuilt from; never read

``ARCHIVE_PATHS.csv`` beside them maps each committed file to its path in the Box archive of
record, for the archive checks.

Resolution order for a read:

0. ``OUT_ROOT/<path>`` for a table the pipeline produces (``GENERATED``): what this run just
   wrote beats every stored copy, which is what makes the stages chain.  Set
   ``SURFACEOME_STORED_INPUTS=1`` to read stored copies throughout instead.
1. ``REFRESH_ROOT/<path>`` when ``SURFACEOME_REFRESH`` names an overlay directory: files there
   win, everything else falls through, so a difference between two runs is attributable to
   the overlaid data alone.
2. ``data/baseline_2022/<path>`` under ``SURFACEOME_BASELINE=1``.
3. ``data/snapshot/<path>`` for a ``GENERATED`` table, else ``data/<path>``.
4. ``DATA_ROOT/<archive path>``, the Box archive, mapped through ``ARCHIVE_PATHS.csv``; never
   needed to run anything.  A missing file resolves there without creating a directory.

Writes go through ``out_path()`` (or ``resolve(..., for_write=True)``) and always land under
``OUT_ROOT``: ``database/`` for the tables that ship, ``tables/`` for intermediates and
per-panel tables, ``figures/`` for figure files.
"""

from __future__ import annotations

import csv
import os
from pathlib import Path

__all__ = [
    "REPO_ROOT", "DATA", "DATA_ROOT", "INPUTS", "CURATED", "MEASUREMENTS", "SNAPSHOT",
    "BASELINE", "SOURCES", "REFRESH_ROOT", "OUT_ROOT", "USE_BASELINE", "USE_STORED_INPUTS",
    "GENERATED", "resolve", "out_path", "archive_path", "committed_files",
]

REPO_ROOT = Path(__file__).resolve().parent.parent   # <repo>/code/ -> <repo>

#: The committed data tree and its tiers.
DATA = REPO_ROOT / "data"
INPUTS = DATA / "inputs"
CURATED = DATA / "curated"
MEASUREMENTS = DATA / "measurements"
SNAPSHOT = DATA / "snapshot"
BASELINE = DATA / "baseline_2022"
SOURCES = DATA / "sources"

#: Root of the full Box archive of record, for the archive tools under code/scripts/archive/.
#: Set SURFACEOME_DATA to point at it.  Not needed to run anything: every input is committed,
#: and the default is a path that does not exist.
DATA_ROOT = Path(os.environ.get("SURFACEOME_DATA", REPO_ROOT / "archive")).expanduser()

#: Read the submitted 2022 vintage instead of the 2026 one.
USE_BASELINE = os.environ.get("SURFACEOME_BASELINE", "") not in ("", "0")
if USE_BASELINE and not BASELINE.is_dir():
    # Without this, every baseline lookup would miss and fall through to the 2026 files,
    # and a run labelled 2022 would quietly be a 2026 run.
    raise RuntimeError(f"SURFACEOME_BASELINE is set but {BASELINE} does not exist; the 2022 "
                       "vintage is kept in the working repository, not in this release")

#: Read stored copies of the generated tables throughout instead of this run's own output.
USE_STORED_INPUTS = os.environ.get("SURFACEOME_STORED_INPUTS", "") not in ("", "0")

#: Optional overlay consulted before everything committed, for trying a refreshed input.
_refresh = os.environ.get("SURFACEOME_REFRESH")
REFRESH_ROOT = Path(_refresh).expanduser() if _refresh else None

#: Where generated files go.  Never the archive and never ``data/``.  Override with
#: SURFACEOME_OUT; ``run_notebooks.py`` sets it to the run's output directory.
OUT_ROOT = Path(os.environ.get("SURFACEOME_OUT", Path.cwd())).expanduser()

#: The tables the pipeline itself produces and reads back, by their ``OUT_ROOT`` path.  Kept
#: explicit so an unrelated file in ``OUT_ROOT`` can never shadow an input.
#: ``check_generated_paths.py`` fails if this list and the notebooks' writes disagree.
GENERATED = frozenset({
    # 01_height_estimates: the surfaceome and its ectodomains
    "database/UP000005640_ecds_total.csv",
    "database/UP000005640_surfaceome_ecds_classified.csv",
    "database/UP000005640_surfaceome_largest_ecds_classified.csv",
    "database/UP000005640_surfaceome_disorder_domain_counting_disorder.csv",
    "database/UP000005640_surfaceome_largest_ecds_domain_counting_disorder.csv",
    "database/UP000005640_disorder_within_largest_ecd_domain_counting_only.csv",
    "database/glycosites.csv",
    # 01_height_estimates: the height tables, read by every figure notebook
    "database/height_estimates.csv",
    "database/height_estimates_with_domain_counts.csv",
    "database/height_estimates_alphafold.csv",
    "database/height_estimates_homology.csv",
    "database/height_estimates_domain_counting.csv",
    "database/domain_disorder_assignments_alphafold.csv",
    "database/domain_disorder_assignments_homology.csv",
    "database/complete_models_alphafold.csv",
    "database/complete_models_homology.csv",
    "database/incomplete_models_alphafold.csv",
    "database/incomplete_models_homology.csv",
    "database/unmodeled_ecds_longest_seq.csv",
    # 02_cellphonedb_interactions
    "database/protein_heights_for_interactions.csv",
    "database/cellphonedb_interaction_sizes.csv",
    "database/cellphonedb_tm_complexes_sizes.csv",
    "database/cellphonedb_tm_interaction_sizes.csv",
    # 03_interaction_heights, read by the Figure 5 and Figure 6 notebooks
    "database/interaction_heights.csv",
    "database/interaction_heights_trans_all.csv",
    # figure-stage tables the database ships
    "database/david_interaction_clustering_long.csv",
    "database/surfaceome_pdbs.csv",
    "database/pdb_structures_measured.csv",
})

#: When ``SURFACEOME_TRACE`` names a file, every read is appended to it as
#: ``<source>\t<path>\t<resolved path>``; ``run_notebooks.py`` summarises it so a run shows
#: which tier each input came from.  Sources: generated, refresh, baseline, snapshot,
#: committed, archive, missing.
_TRACE_PATH = os.environ.get("SURFACEOME_TRACE")


def _trace(source: str, rel, resolved) -> None:
    if not _TRACE_PATH:
        return
    try:
        with open(_TRACE_PATH, "a") as fh:
            fh.write(f"{source}\t{rel}\t{resolved}\n")
    except OSError:
        pass


def _stored(root, rel: Path) -> Path | None:
    """``root/rel``, or its ``.gz`` sibling, if either is a file.

    Text too large to commit raw is committed gzipped; a reader gets the ``.gz`` path back,
    pandas decompresses by extension and ``database_parsing`` opens FASTA through ``gzip``.
    """
    if root is None:
        return None
    plain = root / rel
    if plain.is_file():
        return plain
    packed = root / (str(rel) + ".gz")
    return packed if packed.is_file() else None


_ARCHIVE_PATHS: dict[str, str] | None = None


def archive_path(relpath: str | os.PathLike) -> str:
    """The Box-archive path of a committed file, from ``data/ARCHIVE_PATHS.csv``.

    Files without an entry (everything made in 2026) have no archive copy; their own path
    is returned so the archive tier can still be asked.
    """
    global _ARCHIVE_PATHS
    if _ARCHIVE_PATHS is None:
        _ARCHIVE_PATHS = {}
        table = DATA / "ARCHIVE_PATHS.csv"
        if table.is_file():
            with open(table, newline="") as fh:
                for row in csv.DictReader(fh):
                    _ARCHIVE_PATHS[row["path"]] = row["archive_path"]
    rel = str(relpath)
    return _ARCHIVE_PATHS.get(rel) or _ARCHIVE_PATHS.get(rel + ".gz") or rel


def resolve(relpath: str | os.PathLike, *, for_write: bool = False) -> Path:
    """Resolve a ``data/``-relative path to a concrete file (see the module docstring)."""
    rel = Path(relpath)
    if rel.is_absolute():
        raise ValueError(f"resolve() takes a data-relative path, got absolute {rel!r}.")
    if for_write:
        return out_path(rel)

    generated = str(rel) in GENERATED
    if generated and not USE_STORED_INPUTS:
        produced = OUT_ROOT / rel
        if produced.is_file():
            _trace("generated", rel, produced)
            return produced

    overlaid = _stored(REFRESH_ROOT, rel)
    if overlaid is not None:
        _trace("refresh", rel, overlaid)
        return overlaid

    if USE_BASELINE:
        baseline = _stored(BASELINE, rel)
        if baseline is not None:
            _trace("baseline", rel, baseline)
            return baseline

    if generated:
        snap = _stored(SNAPSHOT, rel)
        if snap is not None:
            _trace("snapshot", rel, snap)
            return snap

    committed = _stored(DATA, rel)
    if committed is not None:
        _trace("committed", rel, committed)
        return committed

    # Nothing committed matches.  Return the archive path without creating anything: a read
    # of a missing file fails where it is read, and a write must go through out_path().
    target = DATA_ROOT / archive_path(rel)
    _trace("archive" if target.exists() else "missing", rel, target)
    return target


def out_path(relpath: str | os.PathLike) -> Path:
    """Path for a generated file, under ``OUT_ROOT``, with parent directories created."""
    target = OUT_ROOT / Path(relpath)
    target.parent.mkdir(parents=True, exist_ok=True)
    return target


def committed_files(root: Path = DATA) -> list[Path]:
    """Every committed data file under ``root``, relative to it."""
    if not root.is_dir():
        return []
    return sorted(p.relative_to(root) for p in root.rglob("*")
                  if p.is_file() and not p.name.startswith("."))


def describe() -> str:
    """One-line summary, handy to print at the top of a notebook."""
    ok = "present" if DATA_ROOT.is_dir() else "not mounted"
    return (
        f"REPO_ROOT={REPO_ROOT}\n"
        f"DATA={DATA}  [{len(committed_files())} files; "
        f"{'baseline 2022' if USE_BASELINE else '2026'} vintage]\n"
        f"DATA_ROOT={DATA_ROOT}  [{ok}]\n"
        f"OUT_ROOT={OUT_ROOT}"
    )


if __name__ == "__main__":
    print(describe())
