#!/usr/bin/env python3
"""Download whole proteomes from UniProt.

Rewritten 2026 for the current UniProt REST API.

The original version (Jamie Heather, 2018; adapted 2019) called
``https://www.uniprot.org/uniprot/?query=proteome:...&columns=...``, which UniProt
retired in mid-2022.  Requests to it now fail.  This module targets
``https://rest.uniprot.org/uniprotkb/`` instead.

One column was silently renamed in the migration::

    legacy 'Gene names'  ->  current 'Gene Names'

Eleven notebooks in this project index that column by its legacy name, so by default
``download_proteome`` renames it back.  Pass ``legacy_headers=False`` to keep UniProt's
current spelling.
"""

from __future__ import annotations

import re
import time
from pathlib import Path
from typing import Iterator

import requests

__all__ = ["TAB_FIELDS", "LEGACY_HEADER_MAP", "download_proteome"]

BASE = "https://rest.uniprot.org/uniprotkb"

#: Columns the pipeline expects in the ``tab``/``tsv`` proteome download.
#: Legacy names, in order, were: id, protein_names, genes, length,
#: comment(SUBCELLULAR LOCATION), feature(TOPOLOGICAL DOMAIN), feature(TRANSMEMBRANE).
TAB_FIELDS = (
    "accession",
    "protein_name",
    "gene_names",
    "length",
    "cc_subcellular_location",
    "ft_topo_dom",
    "ft_transmem",
    "reviewed",
)

#: Header text that changed between the legacy and current APIs.
LEGACY_HEADER_MAP = {"Gene Names": "Gene names"}

#: 'tab' is the legacy name for what the current API calls 'tsv'.
_FORMAT_ALIASES = {"tab": "tsv"}


def _request(url: str, params: dict | None, *, retries: int = 5, timeout: int = 180) -> requests.Response:
    """GET with retry on transient network and server errors."""
    last = None
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, timeout=timeout)
            if r.status_code == 200:
                return r
            if r.status_code in (429, 500, 502, 503, 504):
                last = f"HTTP {r.status_code}"
            else:
                r.raise_for_status()
        except requests.RequestException as exc:  # network hiccup or SSL EOF
            last = f"{type(exc).__name__}: {exc}"
        time.sleep(2 ** attempt)
    raise RuntimeError(
        f"UniProt request failed after {retries} attempts: {last}\n{url} {params}"
    )


_NEXT = re.compile(r'<(?P<url>[^>]+)>;\s*rel="next"')


def _paged(params: dict, *, page_size: int = 500) -> Iterator[str]:
    """Yield response bodies page by page, following UniProt's cursor links.

    The /stream endpoint drops the connection on large proteomes, so the paginated
    /search endpoint is used instead - this is UniProt's documented approach for
    result sets of this size.
    """
    url = f"{BASE}/search"
    p = dict(params, size=page_size)
    total = None
    fetched = 0
    while url:
        r = _request(url, p)
        p = None  # subsequent URLs already carry their query string
        if total is None:
            total = r.headers.get("x-total-results")
            if total:
                print(f"  {int(total):,} entries to fetch")
        body = r.text
        fetched += max(0, body.count("\n") - (1 if fetched == 0 else 0))
        if total:
            print(f"  ... {fetched:,}/{int(total):,}", end="\r", flush=True)
        yield body
        m = _NEXT.search(r.headers.get("Link", ""))
        url = m.group("url") if m else None
    print()


def _apply_legacy_headers(text: str) -> str:
    """Rename the header line back to the legacy column spellings."""
    nl = text.find("\n")
    if nl == -1:
        return text
    header = text[:nl]
    cols = header.split("\t")
    renamed = [LEGACY_HEADER_MAP.get(c, c) for c in cols]
    return "\t".join(renamed) + text[nl:]


def download_proteome(
    proteome_id: str,
    output_file_type: str = "tsv",
    output_path: str | Path = ".",
    *,
    fields: tuple[str, ...] = TAB_FIELDS,
    legacy_headers: bool = True,
    include_isoforms: bool = False,
    filename: str | None = None,
) -> Path:
    """Download one proteome and write it to ``output_path``.

    Parameters
    ----------
    proteome_id
        e.g. ``'UP000005640'`` (human).
    output_file_type
        ``'tsv'`` (or its legacy alias ``'tab'``), ``'fasta'``, or ``'gff'``.
    fields
        Only used for tsv output.
    legacy_headers
        Rename current UniProt column headers back to their pre-2022 spellings so that
        existing notebook code keeps working.  See ``LEGACY_HEADER_MAP``.

    Returns the path written.
    """
    fmt = _FORMAT_ALIASES.get(output_file_type, output_file_type)
    params: dict[str, object] = {
        "query": f"proteome:{proteome_id}",
        "format": fmt,
        "includeIsoform": str(bool(include_isoforms)).lower(),
    }
    if fmt == "tsv":
        params["fields"] = ",".join(fields)

    print(f"retrieving {proteome_id} as {fmt} ...")
    chunks = []
    for i, body in enumerate(_paged(params)):
        if i and fmt == "tsv":
            body = body.split("\n", 1)[1] if "\n" in body else ""  # drop repeated header
        chunks.append(body)
    text = "".join(chunks)
    if fmt == "tsv" and legacy_headers:
        text = _apply_legacy_headers(text)

    out_dir = Path(output_path)
    out_dir.mkdir(parents=True, exist_ok=True)
    name = filename or f"proteomes{proteome_id}.{output_file_type}"
    out = out_dir / name
    out.write_text(text, encoding="utf-8")
    print(f"  wrote {out} ({out.stat().st_size:,} bytes)")
    return out


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("proteome_id", help="e.g. UP000005640")
    ap.add_argument("-f", "--format", default="tsv", choices=["tsv", "tab", "fasta", "gff"])
    ap.add_argument("-o", "--out", default=".")
    ap.add_argument("--current-headers", action="store_true",
                    help="keep UniProt's current column names instead of the legacy ones")
    a = ap.parse_args()
    download_proteome(a.proteome_id, a.format, a.out, legacy_headers=not a.current_headers)
