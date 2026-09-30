#!/usr/bin/env python
"""Download and rebuild the inputs whose licences do not allow them to be redistributed.

    python code/database/fetch_restricted_inputs.py [--cache DIR] [--check]

Three sources the pipeline reads cannot be committed to a public repository
(``data/inputs/RESTRICTED.json`` says why for each):

    glycomine     GlycoMine glycosylation predictions, reduced to the surfaceome
    surfy         the SURFY surface classification, per surfaceome protein
    cellphonedb   the CellphoneDB v5.0.0 tables

For each, this downloads the public release, rebuilds the file the pipeline reads with the
same code that built the committed copy, and checks the result against the sha256 recorded
in ``RESTRICTED.json`` before writing it under ``data/inputs/``. A file that is already
present with the recorded hash is left alone, so a second run downloads nothing.

Downloads are cached under ``--cache`` (default ``downloads/`` at the repository root, which
git ignores). ``--check`` only reports which files are present and correct.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import subprocess
import sys
import tempfile
import urllib.request
import zipfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "data"
INPUTS_CODE = REPO / "code" / "database" / "inputs"
RESTRICTED = DATA / "inputs" / "RESTRICTED.json"


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def present(rel: str, want: dict) -> bool:
    p = DATA / rel
    return p.is_file() and p.stat().st_size == want["bytes"] and sha256(p.read_bytes()) == want["sha256"]


def download(url: str, dest: Path, want_sha: str | None = None) -> Path:
    if not dest.is_file():
        dest.parent.mkdir(parents=True, exist_ok=True)
        print(f"    downloading {url}")
        with urllib.request.urlopen(url, timeout=300) as r:
            data = r.read()
        if want_sha and sha256(data) != want_sha:
            raise SystemExit(f"{url}: the download does not match its recorded sha256; the source has changed")
        dest.write_bytes(data)
    return dest


def build_glycomine(src: dict, cache: Path, tmp: Path) -> None:
    """The header and every row whose accession is in the recorded list, in file order."""
    z = zipfile.ZipFile(download(src["url"], cache / "glycomine" / "proteome_research.zip",
                                 src.get("download_sha256")))
    accs = set((DATA / "inputs/glycosylation/glycomine/accessions.txt").read_text().split())
    for rel in src["files"]:
        name = Path(rel).name
        fh = io.TextIOWrapper(io.BytesIO(z.read(name)))  # text mode: CRLF in the release becomes LF
        rows = [fh.readline()] + [ln for ln in fh if (ln.split(None, 1) or [""])[0] in accs]
        (tmp / rel).parent.mkdir(parents=True, exist_ok=True)
        (tmp / rel).write_text("".join(rows))


def build_surfy(src: dict, cache: Path, tmp: Path) -> None:
    table = download(src["url"], cache / "surfy" / "table_S3_surfaceome.xlsx")
    (rel,) = src["files"]
    (tmp / rel).parent.mkdir(parents=True, exist_ok=True)
    subprocess.run([sys.executable, str(INPUTS_CODE / "build_surfaceome_classification.py"),
                    "--surfaceome", str(DATA / "inputs/annotation/surfaceome_accessions.csv"),
                    "--proteome", str(DATA / "inputs/proteome/UP000005640.tab.gz"),
                    "--surfy", str(table), "--out", str(tmp / rel)],
                   check=True, capture_output=True)


def build_cellphonedb(src: dict, cache: Path, tmp: Path) -> None:
    bundle = download(src["url"], cache / "cellphonedb" / "cellphonedb.zip")
    out_dir = tmp / Path(next(iter(src["files"]))).parent
    subprocess.run([sys.executable, str(INPUTS_CODE / "build_cellphonedb_tables.py"),
                    "--bundle", str(bundle), "--version", "v5.0.0", "--out-dir", str(out_dir)],
                   check=True, capture_output=True)


BUILDERS = {"glycomine": build_glycomine, "surfy": build_surfy, "cellphonedb": build_cellphonedb}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cache", type=Path, default=REPO / "downloads", help="where downloads are kept")
    ap.add_argument("--check", action="store_true", help="report what is present, fetch nothing")
    a = ap.parse_args()

    sources = json.loads(RESTRICTED.read_text())["sources"]
    failed = 0
    for name, src in sources.items():
        missing = [rel for rel, want in src["files"].items() if not present(rel, want)]
        if not missing:
            print(f"  {name:12s} present")
            continue
        if a.check:
            print(f"  {name:12s} MISSING  {', '.join(missing)}")
            failed += 1
            continue
        print(f"  {name:12s} rebuilding from {src['url']}")
        with tempfile.TemporaryDirectory() as t:
            tmp = Path(t)
            BUILDERS[name](src, a.cache, tmp)
            for rel, want in src["files"].items():
                data = (tmp / rel).read_bytes()
                if sha256(data) != want["sha256"]:
                    print(f"    MISMATCH {rel}: the rebuilt file is not the one the pipeline was run on")
                    failed += 1
                    continue
                (DATA / rel).parent.mkdir(parents=True, exist_ok=True)
                (DATA / rel).write_bytes(data)
                print(f"    wrote {rel}")
    print(f"\n  {failed} problem(s)")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
