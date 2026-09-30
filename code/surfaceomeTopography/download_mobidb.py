#!/usr/bin/env python3
"""Download MobiDB disorder annotations.

Rewritten 2026 for the current MobiDB API.

Replaces two near-identical scripts that both called
``http://mobidb.bio.unipd.it/ws/{acc}/consensus``, an endpoint MobiDB has retired:

* ``parse_mobidb_annotations.py``      - kept ``mobidb_consensus.disorder.predictors.mobidb-lite.regions``
* ``parse_mobidb_annotations_full.py`` - kept ``mobidb_consensus.disorder.full.full.regions``

Both ran their work at import time against hardcoded absolute paths on one laptop.

The current API is ``https://mobidb.org/api/download``, returning newline-delimited JSON.
Each value is ``{"regions": [[start, end], ...]}``, 1-based inclusive.

Mapping the old keys across is not a straight rename:

``mobidb_consensus.disorder.predictors.mobidb-lite.regions``
    Disordered regions from MobiDB-lite.  Direct equivalent:
    ``prediction-disorder-mobidb_lite``.

``mobidb_consensus.disorder.full.full.regions``
    **Not** simply a region list.  The legacy "full" track annotated *every* residue with
    a structural state - ``D``/``d`` disordered, ``S``/``s`` structured, ``C``/``c``
    context-dependent, lowercase meaning lower confidence - e.g.
    ``[[1,17,"s"],[18,398,"S"],[399,407,"D"]]``.  ``mobidb_annotations.parse_disorder()``
    then filtered it to ``D|d``.  The closest current key is
    ``prediction-disorder-priority``, which returns those disordered regions directly.

Consequence: ``to_frame()`` output corresponds to the **post-filter** legacy data, i.e. what
``mobidb_annotations.parse_disorder()`` returns - not to the raw archived
``*_disorder_full.txt``.  Do not feed this module's CSV to ``import_disorder()``; that
function parses the legacy raw format.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Iterable, Iterator

import pandas as pd
import requests

__all__ = ["MOBIDB_KEYS", "fetch", "fetch_proteome", "to_frame", "download"]

API = "https://mobidb.org/api/download"

#: Legacy consensus name -> current API key.
MOBIDB_KEYS = {
    "mobidb_lite": "prediction-disorder-mobidb_lite",
    "priority": "prediction-disorder-priority",
    "curated": "curated-disorder-priority",
    "alphafold": "prediction-disorder-alphafold",
}

#: MobiDB accepts a comma-separated accession list; keep batches modest.
BATCH = 100


def _get(params: dict, *, retries: int = 4, timeout: int = 300) -> requests.Response:
    last = None
    for attempt in range(retries):
        try:
            r = requests.get(API, params=params, timeout=timeout, stream=True)
            if r.status_code == 200:
                return r
            if r.status_code in (429, 500, 502, 503, 504):
                last = f"HTTP {r.status_code}"
            else:
                r.raise_for_status()
        except requests.RequestException as exc:
            last = str(exc)
        time.sleep(2 ** attempt)
    raise RuntimeError(f"MobiDB request failed after {retries} attempts: {last} ({params})")


def _iter_ndjson(resp: requests.Response) -> Iterator[dict]:
    for line in resp.iter_lines():
        if line:
            yield json.loads(line)


def fetch(accessions: Iterable[str], *, batch_size: int = BATCH) -> Iterator[dict]:
    """Yield one MobiDB record per accession."""
    accs = [a for a in dict.fromkeys(accessions) if a]
    for i in range(0, len(accs), batch_size):
        chunk = accs[i:i + batch_size]
        print(f"  MobiDB {i + 1}-{i + len(chunk)} of {len(accs)}")
        yield from _iter_ndjson(_get({"acc": ",".join(chunk), "format": "json"}))


def fetch_proteome(proteome_id: str) -> Iterator[dict]:
    """Yield every MobiDB record for a whole proteome, e.g. 'UP000005640'."""
    yield from _iter_ndjson(_get({"proteome": proteome_id, "format": "json"}))


def to_frame(records: Iterable[dict], keys: Iterable[str] = ("mobidb_lite", "priority")) -> pd.DataFrame:
    """Flatten MobiDB records into one row per disordered region.

    Columns: ``acc``, ``source``, ``disorder_start``, ``disorder_end``, ``disorder_len``.
    Coordinates are 1-based inclusive, matching the legacy consensus output.
    """
    wanted = {k: MOBIDB_KEYS[k] for k in keys}
    rows = []
    for rec in records:
        acc = rec.get("acc")
        for label, api_key in wanted.items():
            for start, end in rec.get(api_key, {}).get("regions", []):
                rows.append({
                    "acc": acc,
                    "source": label,
                    "disorder_start": start,
                    "disorder_end": end,
                    "disorder_len": end - start + 1,
                })
    return pd.DataFrame(rows, columns=["acc", "source", "disorder_start",
                                       "disorder_end", "disorder_len"])


def download(
    accessions: Iterable[str] | None = None,
    *,
    proteome_id: str | None = None,
    out_path: str | Path,
    keys: Iterable[str] = ("mobidb_lite", "priority"),
) -> Path:
    """Fetch disorder annotations and write a tidy CSV.

    Give either ``accessions`` or ``proteome_id``.
    """
    if (accessions is None) == (proteome_id is None):
        raise ValueError("pass exactly one of accessions= or proteome_id=")
    records = fetch_proteome(proteome_id) if proteome_id else fetch(accessions)
    df = to_frame(records, keys=keys)
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"wrote {len(df):,} regions for {df['acc'].nunique():,} proteins -> {out}")
    return out


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--proteome", help="e.g. UP000005640")
    g.add_argument("--surfaceome-csv", help="CSV with an 'ID link' or 'Entry' column")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--keys", nargs="+", default=["mobidb_lite", "priority"],
                    choices=sorted(MOBIDB_KEYS))
    a = ap.parse_args()

    if a.proteome:
        download(proteome_id=a.proteome, out_path=a.out, keys=a.keys)
    else:
        df = pd.read_csv(a.surfaceome_csv)
        col = "ID link" if "ID link" in df.columns else "Entry"
        download(df[col].dropna().unique().tolist(), out_path=a.out, keys=a.keys)
