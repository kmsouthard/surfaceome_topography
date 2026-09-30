#!/usr/bin/env python
"""Keep ``surfaceome_config.GENERATED`` honest against what the notebooks actually write.

    python code/check_generated_paths.py [--list]

``GENERATED`` is the set of tables the pipeline produces and reads back, and ``resolve()``
consults ``OUT_ROOT`` first for exactly those.  That is what makes the stages chain, so the
set has to match reality: a table missing from it is silently read from the committed
snapshot instead of from this run, and nothing says so.

Every write in the notebooks goes through ``out_path('<database|tables|figures>/...')``.
This parses those calls in the notebooks ``run_notebooks.ORDER`` lists and checks:

* every ``database/`` table is in ``GENERATED``, and so is every ``tables/`` file some
  notebook reads back; nothing else is;
* **one writer per table**: two notebooks writing the same file means whichever ran last
  silently won;
* **no table read before the notebook that writes it has run**: otherwise the read falls
  back to the snapshot copy from an earlier run.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

from surfaceome_config import GENERATED

REPO = Path(__file__).resolve().parent.parent
NBDIR = REPO / "code"   # notebook paths in run_notebooks.ORDER are relative to this

OUT_PATH = re.compile(r"out_path\(\s*(.+?)\s*\)")
RESOLVE = re.compile(r"resolve\(\s*(f?'[^']+'|f?\"[^\"]+\")")

#: Constants the notebooks concatenate into paths.
CONSTS = {"proteome": "UP000005640"}

FIGURE_EXT = (".pdf", ".svg", ".png", ".html", ".fasta")


def _literal(expr: str) -> str | None:
    """The path an ``out_path``/``resolve`` argument names, with the proteome id filled in."""
    e = expr.strip()
    e = (e.replace("'+proteome+'", CONSTS["proteome"]).replace('"+proteome+"', CONSTS["proteome"])
         .replace("{proteome}", CONSTS["proteome"]))
    if e[:1] in "fF":
        e = e[1:]
    if len(e) >= 2 and e[0] == e[-1] and e[0] in "'\"" and e[0] not in e[1:-1]:
        return e[1:-1]
    return None


def _code(cell: dict) -> str:
    return "\n".join(l for l in "".join(cell["source"]).splitlines() if not l.strip().startswith("#"))


def write_targets() -> dict[str, str]:
    """{out path written: "notebook.ipynb:cell"} over every notebook, figures included.

    Shared with ``build_snapshot.py``, which records the producer of each snapshot table.
    """
    from run_notebooks import ORDER  # here, not at import: run_notebooks imports this module's users
    found: dict[str, str] = {}
    for name in ORDER:
        nb = json.loads((NBDIR / name).read_text())
        for i, c in enumerate(nb["cells"]):
            if c["cell_type"] != "code":
                continue
            for m in OUT_PATH.finditer(_code(c)):
                rel = _literal(m.group(1))
                if rel:
                    found.setdefault(rel, f"{name}:{i}")
    return found


def pipeline_graph() -> tuple[list[str], dict[str, set], dict[str, set]]:
    """(run order, {notebook: files written}, {notebook: files read}), from the notebooks."""
    from run_notebooks import ORDER
    writes: dict[str, set] = {}
    reads: dict[str, set] = {}
    for name in ORDER:
        nb = json.loads((NBDIR / name).read_text())
        w, r = set(), set()
        for c in nb["cells"]:
            if c["cell_type"] != "code":
                continue
            src = _code(c)
            for m in OUT_PATH.finditer(src):
                rel = _literal(m.group(1))
                if rel:
                    w.add(rel)
            for m in RESOLVE.finditer(src):
                rel = _literal(m.group(1))
                if rel:
                    r.add(rel)
        writes[name], reads[name] = w, r
    return ORDER, writes, reads


def check_order() -> list[str]:
    """One writer per table, and no table read before the notebook that writes it has run."""
    order, writes, reads = pipeline_graph()
    pos = {nb: i for i, nb in enumerate(order)}
    owners: dict[str, list[str]] = {}
    for nb in order:
        for rel in writes[nb]:
            if not rel.endswith(FIGURE_EXT):
                owners.setdefault(rel, []).append(nb)
    problems = []
    for rel, who in sorted(owners.items()):
        if len(who) > 1:
            problems.append(f"WRITTEN BY {len(who)} NOTEBOOKS  {rel}  -- " + " | ".join(who))
    for nb in order:
        for rel in sorted(reads[nb]):
            if rel not in GENERATED:
                continue
            who = owners.get(rel, [])
            if not who:
                problems.append(f"READ, NEVER WRITTEN  {nb} reads {rel}, which no notebook produces")
            elif min(pos[w] for w in who) >= pos[nb]:
                problems.append(f"READ BEFORE WRITTEN  {nb} reads {rel}, written by " + ", ".join(who))
    return problems


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true", help="print every write target found")
    args = ap.parse_args()

    found = write_targets()
    _, _, reads = pipeline_graph()
    read_back = set().union(*reads.values())

    if args.list:
        for k in sorted(found):
            print(f"  {k:<80} {found[k]}")
        print(f"\n  {len(found)} write targets")

    tables = {k for k in found if not k.endswith(FIGURE_EXT)}
    should = {k for k in tables if k.startswith("database/") or k in read_back}
    unlisted = should - set(GENERATED)
    extra = {k for k in GENERATED if k in tables and k not in should}
    stale = set(GENERATED) - tables
    misplaced = {k for k in found if not k.startswith(("database/", "tables/", "figures/"))}

    for k in sorted(unlisted):
        print(f"  NOT IN GENERATED   {k}   ({found[k]})")
    for k in sorted(extra):
        print(f"  IN GENERATED BUT NEITHER database/ NOR READ BACK   {k}")
    for k in sorted(stale):
        print(f"  DECLARED, NOT WRITTEN BY ANY NOTEBOOK   {k}")
    for k in sorted(misplaced):
        print(f"  OUTSIDE database/, tables/, figures/   {k}   ({found[k]})")

    order_problems = check_order()
    for msg in order_problems:
        print(f"  {msg}")

    problems = len(unlisted) + len(extra) + len(stale) + len(misplaced) + len(order_problems)
    print(f"\n  {len(found)} write targets, {len(GENERATED)} declared, "
          f"{len(order_problems)} ordering problem(s), {problems} problem(s)")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
