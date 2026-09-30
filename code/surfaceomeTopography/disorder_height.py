#!/usr/bin/env python3
"""Height of a disordered segment from its length and glycan grafting density.

The pipeline's original rule switched between two persistence lengths at a single
glycan-density threshold::

    glycan_density > 0.03  ->  l_p = 8 nm   (bottle-brush)
    glycan_density <= 0.03 ->  l_p = 1 nm   (bare chain)

That threshold does not sit at a break in the data. Of 2835 annotated disordered
segments, 1973 (70%) carry *no* predicted or annotated glycosite at all, and the
remaining 862 spread smoothly over two orders of magnitude with no trough near 0.03;
a third of them fall within a factor of two of the cutoff. The real structure in the
data is "glycosylated or not", after which density varies continuously -- so an 8x
step in stiffness at one density is the wrong shape for the model.

This module keeps the same two physical limits, which come from Kuo et al.,
Nat. Phys. 14:658-669 (2018), and interpolates between them instead:

    sigma -> 0     l_p -> l_p_bare    (bare disordered polypeptide)
    sigma -> large l_p -> l_p_brush   (persistence length ~ glycan side-chain length)

The crossover is set by when side chains begin to crowd along the backbone. Glycans
grafted at density `sigma` per residue sit `b / sigma` apart, with b = 0.4 nm the
per-residue rise; they start to overlap once that spacing falls to the glycan's own
size R. So the natural crossover density is `sigma_50 = b / R`.

R is not a free parameter. The brush limit `l_p_brush` is itself Kuo's statement that a
densely grafted bottle-brush has a persistence length of about its side-chain length, so
self-consistency requires the *same* R in both places: R = l_p_brush, hence

    sigma_50 = b / l_p_brush        = 0.4 / 8 = 0.050
    l_p(sigma) = l_p_bare + (l_p_brush - l_p_bare) * sigma / (sigma + sigma_50)

The model therefore has no fitted parameter -- only l_p_bare and l_p_brush, both taken
from the literature with stated ranges.

Two independent checks, neither used to set anything:

* The sigma_50 that leaves the pipeline's *total* disorder height budget unchanged across
  all 2835 segments is 0.0485, within 3% of the 0.050 the crowding argument gives. So the
  change redistributes stiffness between segments rather than rescaling the population.
* At CD45's density (sigma = 0.0354) the model returns l_p = 3.9 nm. The EM-derived
  ectodomain height for CD45RABC (~40 nm total, ~15 nm folded) implies l_p = 4.0 nm. CD45
  was not used to fit anything.

Every parameter is a literature quantity with a stated range rather than a fitted
value, because the available anchors cannot pin them down: CD45's own two isoforms,
which share a folded region and glycan chemistry, imply l_p = 1.8 nm (CD45RO, 22 nm)
and 4.1 nm (CD45RABC, 40 nm) from the same measurements. Report the interval.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .polymer_models import AA_RISE_NM, wlc_end_to_end

__all__ = ["L_P_BARE_NM", "L_P_BRUSH_NM",
           "persistence_length", "disorder_height", "disorder_height_interval",
           "merge_runs", "combine", "heterogeneous_wlc", "plddt_disorder_regions",
           "assign_disorder_heights"]

#: Bare disordered polypeptide, nm. Pipeline used 1.0; the code's own note records
#: "reported values vary from 0.5-3 nm".
L_P_BARE_NM = (0.5, 1.0, 3.0)

#: Densely grafted bottle-brush, nm. Kuo: mucin glycan side chains 5-10 nm, and at high
#: grafting density l_p is about the side-chain length. Pipeline used 8.0.
L_P_BRUSH_NM = (5.0, 8.0, 10.0)

#: The crowding crossover density is rise / l_p_brush; it is not an independent parameter.


def persistence_length(glycan_density, *, l_p_bare=1.0, l_p_brush=8.0,
                       rise_nm=AA_RISE_NM):
    """Effective persistence length, nm, for a segment at this grafting density.

    The crowding crossover is fixed at ``rise_nm / l_p_brush`` -- see the module
    docstring; there is no free shape parameter.
    """
    sigma = np.clip(np.asarray(glycan_density, dtype=float), 0.0, None)
    sigma_50 = float(rise_nm) / float(l_p_brush)
    return l_p_bare + (l_p_brush - l_p_bare) * sigma / (sigma + sigma_50)


def disorder_height(n_residues, glycan_density, **kw):
    """Kratky-Porod end-to-end distance, nm, with a density-dependent l_p."""
    lp = persistence_length(glycan_density, **kw)
    return wlc_end_to_end(n_residues, persistence_nm=lp,
                          rise_nm=kw.get("rise_nm", AA_RISE_NM))


def disorder_height_interval(n_residues, glycan_density):
    """(low, central, high) height in nm, propagating the literature ranges above.

    The bounds come from the extremes of L_P_BARE_NM, L_P_BRUSH_NM and GLYCAN_SIZE_NM,
    not from a fit, so they describe how far the answer moves across values the
    literature actually supports.
    """
    n = np.asarray(n_residues, dtype=float)
    s = np.asarray(glycan_density, dtype=float)
    lo_b, mid_b, hi_b = L_P_BARE_NM
    lo_r, mid_r, hi_r = L_P_BRUSH_NM
    corners = [disorder_height(n, s, l_p_bare=b, l_p_brush=r)
               for b in (lo_b, hi_b) for r in (lo_r, hi_r)]
    stack = np.stack(corners)
    central = disorder_height(n, s, l_p_bare=mid_b, l_p_brush=mid_r)
    return stack.min(axis=0), central, stack.max(axis=0)


# ---------------------------------------------------------------------------
# Combining pieces along an ectodomain
# ---------------------------------------------------------------------------
#
# An ECD is a linear chain of alternating rigid and flexible pieces. Each piece
# contributes its *typical extent*: a bounding-box dimension for a folded domain
# (rigid, so there is no ensemble spread) or a Kratky-Porod end-to-end distance for a
# disordered stretch. Those are the same kind of quantity, so they can be combined.
#
# How they combine depends on what separates them:
#
# * Within ONE continuous flexible stretch there is nothing to align -- a chain does not
#   align with itself. It has a single end-to-end distance. Splitting it into annotated
#   sub-segments and summing their heights double-counts, and is simply wrong.
# * Across pieces separated by a folded domain there is a real hinge. Summing linearly is
#   the *aligned* limit (an upper bound); adding in quadrature is the *orientationally
#   independent* limit (a lower bound). The manuscript's quantity -- fully extended height
#   as the maximum membrane separation -- is the aligned limit, so linear is the intended
#   choice, but it should be labelled as such and bracketed by the quadrature value.
#
# Merging is deliberately conservative: only across gaps that are BOTH short and free of
# any assigned structural element. A long unassigned tract is ignorance, not evidence of
# continuity, and merging across one would silently invent a single long chain. That is a
# coverage problem and belongs in the coverage accounting, not here.

MERGE_GAP_MAX_AA = 5


def merge_runs(segments, occupied=(), gap_max=MERGE_GAP_MAX_AA):
    """Merge annotated disordered segments into continuous flexible runs.

    `segments` is an iterable of (start, end, glycan_count). `occupied` is an iterable of
    (start, end) spans assigned to any structural element. Two consecutive segments merge
    only when the gap between them is at most `gap_max` residues and no occupied span
    intersects it.

    Returns a list of (start, end, glycan_count) runs. Lengths follow the pipeline's
    convention of ``end - start``.
    """
    segs = sorted(segments, key=lambda s: s[0])
    occ = list(occupied)
    out = []
    for start, end, gly in segs:
        if out:
            p_start, p_end, p_gly = out[-1]
            lo, hi = p_end + 1, start - 1
            gap = start - p_end - 1
            clear = not any(not (e < lo or s > hi) for s, e in occ) if hi >= lo else True
            if 0 <= gap <= gap_max and clear:
                out[-1] = (p_start, max(p_end, end), p_gly + gly)
                continue
        out.append((start, end, gly))
    return out


def combine(heights, mode="linear"):
    """Combine piece extents along an ectodomain.

    ``linear``     aligned limit -- the sum. What the manuscript reports.
    ``quadrature`` orientationally independent limit -- sqrt of the sum of squares.
    """
    h = np.asarray(list(heights), dtype=float)
    if h.size == 0:
        return 0.0
    if mode == "linear":
        return float(h.sum())
    if mode == "quadrature":
        return float(np.sqrt((h ** 2).sum()))
    raise ValueError(f"unknown mode {mode!r}")


def heterogeneous_wlc(lengths, persistences, rise_nm=AA_RISE_NM):
    """RMS end-to-end distance, nm, of a chain whose persistence length varies along it.

    A disordered run is rarely uniform: glycan density changes along it, so stiffness does
    too. Pooling to a mean density and applying plain Kratky-Porod is an approximation.
    This is the exact result for piecewise-constant persistence length.

    For a worm-like chain the tangent correlation decays through whatever it crosses::

        <t(s).t(s')> = exp( -integral ds'' / l_p(s'') )

    Carrying that through <R^2> = int int <t(s).t(s')> ds ds' gives, for pieces i of contour
    length L_i and persistence l_i::

        <R^2> = sum_i [ 2 l_i L_i - 2 l_i^2 (1 - e^-L_i/l_i) ]
                + 2 sum_{i<j} mu_i mu_j prod_{i<k<j} gamma_k

        mu_i    = l_i (1 - e^-L_i/l_i)     mean projection of piece i on the entry tangent
        gamma_k = e^-L_k/l_k               tangent survival across piece k

    The first term is the usual per-piece Kratky-Porod; the second carries the orientational
    correlation between pieces, which is exactly what quadrature throws away and what pooling
    approximates. With a single persistence length throughout it reduces to plain
    Kratky-Porod over the total length, to machine precision.
    """
    L = rise_nm * np.asarray(lengths, dtype=float)
    lp = np.asarray(persistences, dtype=float)
    if L.size == 0:
        return 0.0
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        r2 = float(np.sum(2.0 * lp * L - 2.0 * lp**2 * (1.0 - np.exp(-L / lp))))
        mu = lp * (1.0 - np.exp(-L / lp))
        gamma = np.exp(-L / lp)
    for i in range(L.size):
        surviving = 1.0
        for j in range(i + 1, L.size):
            r2 += 2.0 * mu[i] * mu[j] * surviving
            surviving *= gamma[j]
    return float(np.sqrt(max(r2, 0.0)))


# ---------------------------------------------------------------------------
# Calling disorder from AlphaFold confidence
# ---------------------------------------------------------------------------

#: pLDDT at or above which a residue is treated as ordered. Piovesan, Monzon & Tosatto,
#: Protein Science 2022;31(11):e4466 report an F1-optimal disorder threshold of pLDDT
#: < 68.8% on the CAID DisProt benchmark; 70 is that value at the resolution of the
#: pipeline's per-residue confidence string, which stores pLDDT // 10.
PLDDT_DISORDER_MAX = 70

#: Shortest run of low-confidence residues accepted as a disordered region. Below roughly
#: this length a segment cannot support its own glycan-density estimate -- one predicted
#: glycosite changes the density, and hence the persistence length, by more than the whole
#: literature range -- so short runs contribute noise rather than signal.
MIN_DISORDER_RUN_AA = 20

#: Minimum for a run in the confidence-trimmed tails, where nothing else covers the residues.
#: Discarding a short internal run costs nothing -- it lies inside the AlphaFold bounding box
#: and the overlap resolution drops it anyway -- but discarding a short terminal run leaves
#: those residues contributing no height at all. See `plddt_disorder_regions`.
MIN_TERMINAL_RUN_AA = 5


def plddt_disorder_regions(confidence_string, ecd_start=0,
                           plddt_max=PLDDT_DISORDER_MAX,
                           min_length=MIN_DISORDER_RUN_AA,
                           measured_span=None,
                           terminal_min_length=MIN_TERMINAL_RUN_AA):
    """Disordered regions called from an AlphaFold per-residue confidence string.

    `confidence_string` holds one digit per residue, pLDDT // 10, indexed from `ecd_start`.
    Returns a list of (start, end) in the same coordinates as `ecd_start`.

    Two things differ from the original inline implementation::

        re.finditer('[0-' + str(cutoff) + ']{1,}', confidence_string)   # cutoff = 5

    * `[0-5]` is inclusive of the digit 5, i.e. pLDDT 50-59, so that call thresholded at
      pLDDT < 60 while its comment said 50. Here the threshold is stated in pLDDT units and
      converted once, so the code and its description cannot drift apart.
    * `{1,}` accepted single-residue regions. The median region it produced was 4 residues
      and 85% were shorter than 20, which is what fragmented the disorder annotation and
      left most segments too short to estimate a glycan density from.

    Pass `measured_span` as the (lo, hi) residue range the structural bounding box covers --
    for AlphaFold that is ``(ecd_start + n_trim, ecd_start + c_trim)``. Runs falling wholly
    outside it are in the confidence-trimmed tails, where no structural measurement reaches,
    so they are held to `terminal_min_length` instead of `min_length`. Dropping a short
    *internal* run costs nothing, because it lies inside the bounding box and the overlap
    resolution discards it regardless; dropping a short *terminal* run would leave those
    residues contributing no height at all. With the asymmetric defaults the annotation ends
    up covering more of the trimmed tails than the original implementation did (4464 orphaned
    residues against 8441) while still raising the median region from 7 to 24 residues.
    """
    s = str(confidence_string)
    digit_max = int(plddt_max) // 10          # pLDDT < 70  ->  digits 0-6

    def keep(start, end):
        if measured_span is None:
            return end - start >= min_length
        lo, hi = measured_span
        terminal = (end + ecd_start <= lo) or (start + ecd_start >= hi)
        return end - start >= (terminal_min_length if terminal else min_length)

    out, run_start = [], None
    for i, ch in enumerate(s):
        low = ch.isdigit() and int(ch) < digit_max
        if low and run_start is None:
            run_start = i
        elif not low and run_start is not None:
            if keep(run_start, i):
                out.append((run_start + ecd_start, i + ecd_start))
            run_start = None
    if run_start is not None and keep(run_start, len(s)):
        out.append((run_start + ecd_start, len(s) + ecd_start))
    return out


def assign_disorder_heights(heights, id_col="target_name"):
    """Merge continuous disordered runs within each ECD and give each one a height.

    Drop-in replacement for the inline block in
    `database/notebooks/01_height_estimates.ipynb` cell 55, which applied a
    binary persistence length row by row::

        PL_gly = 4.0; PL_idp = 0.5; glycan_cutoff = 0.03
        heights['height'] = heights.apply(lambda row: sqrt(4*PL*L*(...)), axis=1)

    Three changes, each argued in README_disorder_height.md:

    1. `l_p` varies continuously with glycan density instead of stepping 1 -> 8 nm at 0.03.
    2. Disordered rows that are contiguous, with no structural row between them, are one
       polymer and get one calculation rather than being summed as separate pieces.
    3. Within a merged run the persistence length is allowed to vary piece by piece, using
       the exact heterogeneous worm-like chain rather than a pooled mean density.

    Structural rows keep their measured bounding-box height, converted from Angstrom to nm.
    That conversion lived in the final `else` of the original lambda --
    ``... else row['height']/10`` -- reached only when `glycan_density` was NaN, i.e. on
    structural rows. It is easy to miss inside a disorder-height expression, and dropping it
    silently leaves every domain height a factor of ten too large: EGF (PF00008) is stored as
    32.674972, which is Angstrom, not nm.

    Returns a new frame; merged runs carry `n_merged` and the rows they replaced are dropped.
    """
    df = heights.copy()
    is_dis = df["disorder_len"].notna()
    struct = df[~is_dis].copy()
    struct["height"] = struct["height"] / 10.0        # Angstrom -> nm, as the original lambda did
    dis = df[is_dis].sort_values([id_col, "disorder_start"])

    occupied = {}
    for r in struct.itertuples():
        s = getattr(r, "dom_start", np.nan)
        e = getattr(r, "dom_end", np.nan)
        if pd.notna(s) and pd.notna(e):
            occupied.setdefault(getattr(r, id_col), []).append((s, e))

    out_rows = []
    for key, grp in dis.groupby(id_col, sort=False):
        segs = [(r.disorder_start, r.disorder_end, (r.glycan_count if pd.notna(r.glycan_count) else 0.0))
                for r in grp.itertuples()]
        runs = merge_runs(segs, occupied.get(key, []))
        for run_start, run_end, run_gly in runs:
            members = grp[(grp.disorder_start >= run_start) & (grp.disorder_end <= run_end)]
            row = members.iloc[0].copy()
            run_len = max(run_end - run_start, 1.0)
            if len(members) == 1:
                row["height"] = float(disorder_height(run_len, run_gly / run_len))
            else:
                lengths, lps, cursor = [], [], run_start
                for m in members.itertuples():
                    if m.disorder_start > cursor:                     # unassigned gap residues
                        lengths.append(m.disorder_start - cursor)
                        lps.append(float(persistence_length(run_gly / run_len)))
                    lengths.append(max(m.disorder_end - m.disorder_start, 1.0))
                    lps.append(float(persistence_length(
                        m.glycan_density if pd.notna(m.glycan_density) else 0.0)))
                    cursor = m.disorder_end
                if run_end > cursor:
                    lengths.append(run_end - cursor)
                    lps.append(float(persistence_length(run_gly / run_len)))
                row["height"] = heterogeneous_wlc(lengths, lps)
            row["disorder_start"] = run_start
            row["disorder_end"] = run_end
            row["disorder_len"] = run_len
            row["dom_len"] = run_len
            row["glycan_count"] = run_gly
            row["glycan_density"] = run_gly / run_len
            row["n_merged"] = len(members)
            out_rows.append(row)

    merged = pd.DataFrame(out_rows) if out_rows else dis.iloc[0:0].copy()
    return pd.concat([struct, merged], sort=False).reset_index(drop=True)
