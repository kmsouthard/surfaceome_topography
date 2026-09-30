#!/usr/bin/env python3
"""How much of an ectodomain its annotations account for.

`annotation_coverage` adds three columns, in percent of the ectodomain's sequence length, to a
per-ectodomain frame aggregated from an assignment table -- the ``<column>_sum`` and
``seq_len_first`` columns the height notebook's ``groupby('ID').agg(...)`` produces:

    dom_coverage       every assigned residue, whichever method assigned it (``dom_len_sum``)
    dom_coverage_each  domains + disorder + the structural model, summed per method
    best_coverage      the larger of the two; an ectodomain is "well modelled" at >= 70

``dom_len`` already carries the disorder lengths (``dom_len == disorder_len`` on every disorder
row), so adding ``disorder_len_sum`` to ``dom_coverage`` would count disordered residues twice.
"""

from __future__ import annotations

import pandas as pd

__all__ = ["annotation_coverage"]


def annotation_coverage(heights: pd.DataFrame, model_len_col: str) -> pd.DataFrame:
    """Return ``heights`` with ``dom_coverage``, ``dom_coverage_each`` and ``best_coverage``.

    ``model_len_col`` is the summed length of the structural model the assignment used:
    ``alphafold_len_sum`` or ``homology_len_sum``.
    """
    seq_len = heights["seq_len_first"]
    heights = heights.assign(
        dom_coverage=heights["dom_len_sum"] / seq_len * 100,
        dom_coverage_each=(heights["disorder_len_sum"] + heights["domain_len_sum"]
                           + heights[model_len_col]) / seq_len * 100)
    return heights.assign(best_coverage=heights[["dom_coverage", "dom_coverage_each"]].max(axis=1))
