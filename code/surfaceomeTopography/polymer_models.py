#!/usr/bin/env python3
"""Polymer models for the extension of disordered segments.

The manuscript estimates the height contributed by a disordered stretch from its sequence
length rather than from any predicted structure, using the Kratky-Porod (worm-like chain)
model for the mean end-to-end distance::

    <R^2> = 2 L l_p - 2 l_p^2 (1 - e^(-L / l_p))

with `L` the contour length and `l_p` the persistence length, and `L` = 0.4 nm per residue.

The pipeline's own implementation lives inline in a notebook cell and is written in a
half-persistence-length parameterisation::

    <R^2> = 4 P L (1 - (2P/L)(1 - e^(-L/2P)))        with PL_gly = 4.0, PL_idp = 0.5

That is algebraically identical to the Kratky-Porod form above with ``l_p = 2P`` -- they agree
to 5e-6 nm over 1-2000 residues. So the effective persistence lengths actually used are

    glycosylated (bottle-brush) : l_p = 8.0 nm   (PL_gly = 4.0)
    unglycosylated IDP          : l_p = 1.0 nm   (PL_idp = 0.5)

l_p = 8.0 nm follows Kuo, Gandhi, Zia & Paszek, *Physical biology of the cancer cell
glycocalyx*, Nat. Phys. 14:658-669 (2018), which reports mucin glycan side chains of 5-10 nm
and notes that a densely grafted bottle-brush has a persistence length of about the same size
as its side chains. 8 nm is a mid-range choice within that estimate.

Which of the two values applies is decided per segment by glycan density, with a cutoff of
0.03: denser segments are treated as bottle-brushes (l_p = 8 nm), sparser ones as bare chains
(l_p = 1 nm). Most annotated disordered segments -- 933 of 1463 -- fall on the bare-chain side,
so applying the bottle-brush model everywhere overestimates total disorder height by ~56%.
The caller is responsible for choosing the branch; this module only supplies the constants.
"""

from __future__ import annotations

import numpy as np

__all__ = ["AA_RISE_NM", "PERSISTENCE_LENGTH_NM", "PERSISTENCE_LENGTH_IDP_NM",
           "contour_length", "wlc_end_to_end"]

#: Contour length contributed per amino acid, nm (Ainavarapu et al.).
AA_RISE_NM = 0.4

#: Persistence length, nm, for a glycosylated (bottle-brush) segment. The pipeline uses
#: PL_gly = 4.0 as a *half* persistence length, so l_p = 8.0 nm.
PERSISTENCE_LENGTH_NM = 8.0

#: Persistence length, nm, for an unglycosylated intrinsically disordered segment
#: (pipeline PL_idp = 0.5, again a half persistence length).
PERSISTENCE_LENGTH_IDP_NM = 1.0


def contour_length(n_residues, rise_nm: float = AA_RISE_NM):
    """Contour length L in nm for a segment of `n_residues`."""
    return np.asarray(n_residues, dtype=float) * rise_nm


def wlc_end_to_end(n_residues, persistence_nm: float = PERSISTENCE_LENGTH_NM,
                   rise_nm: float = AA_RISE_NM):
    """Root-mean-square end-to-end distance, nm, of a disordered segment.

    >>> round(float(wlc_end_to_end(100)), 2)   # 100 aa, L = 40 nm, l_p = 8 nm
    22.65
    >>> round(float(wlc_end_to_end(0)), 2)
    0.0

    For L >> l_p this approaches the random-coil limit sqrt(2 L l_p); for L << l_p it
    approaches L, a rigid rod.
    """
    L = contour_length(n_residues, rise_nm)
    lp = np.asarray(persistence_nm, dtype=float)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        r2 = 2.0 * L * lp - 2.0 * lp**2 * (1.0 - np.exp(-L / lp))
    return np.sqrt(np.clip(np.nan_to_num(r2, nan=0.0), 0.0, None))
