"""One cell-cell contact, from two surfaces and the interaction table to its tables and figures.

The Figure 6 panels were five copies of one notebook, each with its immune cell, its cancer cell
and whether an antibody bridges them typed into the code. This module is that notebook's
computation once; ``data/curated/contact_panels.csv`` names the panels, and the F6 notebook
loops over them. Every quantity is what the copies computed: the pairs that can form across the
contact (`interface.pair_interactions`), the antibody bridges' gaps, the shared-abundance weight,
the contact's exclusion on each side and the cis pairs on each side, and the clearance matrix
behind the interaction-by-protein heatmaps.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.colors import SymLogNorm

from .interface import (antibody_bridge_heights, cis_pairs_at_contact, clearance, contact_exclusion,
                        interface_weights, pair_interactions, select_interactions, surface_units)

__all__ = ["contact_pairs", "contact_layout", "plot_clearance", "plot_contact_matrix"]

#: The columns of a pair the summary views show, in this order.
SUMMARY = ["interface_shared", "interface_prob", "combined_expression", "interaction_dim",
           "percent_expression_cancer", "percent_expression_immune", "limiting_cell",
           "Gene Name_immune", "Gene Name_cancer", "complex_name_immune", "complex_name_cancer"]


def contact_pairs(immune: pd.DataFrame, cancer: pd.DataFrame, interactions: pd.DataFrame,
                  antibody_target: str | None = None) -> pd.DataFrame:
    """Every trans pair that can form between two surfaces, with its gap and its share of the contact.

    ``interactions`` is the interaction table; ``antibody_target`` keeps the antibody bridges to
    that antigen (an entry-name substring, ``'ERBB2'`` for trastuzumab) and ``None`` drops every
    bridge. The bridges' gaps are set from the epitope (`antibody_bridge_heights`), which leaves
    a table without bridges unchanged, and each pair's share of the contact follows from the two
    cells' abundances by mass action (`interface_weights`).
    """
    pairs = pair_interactions(immune, cancer, select_interactions(interactions, antibody_target))
    return interface_weights(antibody_bridge_heights(pairs))


def contact_layout(pairs: pd.DataFrame, immune: pd.DataFrame, min_share: float = 0.1,
                   min_expression: float = 1.0):
    """The clearance matrix: the pairs holding at least ``min_share`` % of the contact, by gap,
    against the immune cell's units above ``min_expression`` %, by height.

    Returns ``(td, units, vs)``: the selected pairs indexed by ``immune,cancer`` gene names and
    sorted by gap; the selected units indexed by gene (complex) name and sorted by height; and
    `interface.clearance` of the two, gap minus height, positive where the protein fits beneath.
    """
    td = pairs[pairs["interface_shared"] >= min_share].sort_values("interaction_dim").reset_index()
    td["names"] = td["Gene Name_immune"] + "," + td["Gene Name_cancer"]
    td = td.set_index("names")
    units = immune[immune["percent_expression"] > min_expression].copy()
    units.loc[:, "Gene Name"] = np.where(units["complex_name"].isnull(), units["Gene Name"], units["complex_name"])
    units = units.sort_values("total_height").set_index("Gene Name")
    return td, units, clearance(td["interaction_dim"], units["total_height"])


def plot_clearance(vs: pd.DataFrame):
    """The clearance matrix alone, interactions by proteins, on a symmetric log scale."""
    fig = plt.figure(figsize=(20, 15))
    ax = sns.heatmap(vs, cmap="RdGy", norm=SymLogNorm(linthresh=3, linscale=1),
                     cbar_kws={"ticks": [-60, -40, -20, -10, -8, -6, -4, -2, -1, 0, 1, 2, 4, 6, 8, 10, 20, 40, 60]})
    ax.invert_yaxis()
    return fig


def plot_contact_matrix(vs: pd.DataFrame, units: pd.DataFrame, td: pd.DataFrame):
    """The clearance matrix with its margins: each protein's height and share of its surface
    along the bottom, each pair's gap, share of the contact and the two partners' shares down
    the side. ``units`` and ``td`` are `contact_layout`'s."""
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    palette_r = sns.diverging_palette(321, 179, s=80, l=50, n=99)
    greys = sns.light_palette((206, 19, 43), input="husl")
    light = sns.light_palette((4, 5, 50), input="husl")

    fig = plt.figure(figsize=(15, 10))
    ax1 = plt.subplot2grid((20, 20), (0, 0), colspan=15, rowspan=17)
    ax2 = plt.subplot2grid((20, 20), (18, 0), colspan=12, rowspan=1)
    ax3 = plt.subplot2grid((20, 20), (0, 16), colspan=1, rowspan=17)
    ax4 = plt.subplot2grid((20, 20), (19, 0), colspan=12, rowspan=1)
    ax5 = plt.subplot2grid((20, 20), (0, 17), colspan=1, rowspan=17)
    ax6 = plt.subplot2grid((20, 20), (0, 18), colspan=1, rowspan=17)
    ax7 = plt.subplot2grid((20, 20), (0, 19), colspan=1, rowspan=17)

    sns.heatmap(vs, ax=ax1, cmap=palette_r, norm=SymLogNorm(linthresh=3, linscale=1),
                cbar_kws={"ticks": [-40, -20, -10, -8, -6, -4, -2, -1, 0, 1, 2, 4, 6, 8, 10, 20, 40]})
    margin = dict(annot=True, cbar=False, xticklabels=False, yticklabels=False)
    sns.heatmap(pd.DataFrame(units["total_height"].round(0)).transpose(), ax=ax2, cmap=greys,
                **{**margin, "xticklabels": True})
    sns.heatmap(pd.DataFrame(td["interaction_dim"].round(0)), ax=ax3, cmap=greys, **margin)
    sns.heatmap(pd.DataFrame(units["percent_expression"].round(0)).transpose(), ax=ax4, cmap=greys, **margin)
    sns.heatmap(pd.DataFrame(td["interface_shared"].round(1)), ax=ax5, cmap=light, **margin)
    sns.heatmap(pd.DataFrame(td["percent_expression_immune"].round(2)), ax=ax6, cmap=light, **margin)
    sns.heatmap(pd.DataFrame(td["percent_expression_cancer"].round(2)), ax=ax7, cmap=light, **margin)
    for ax in (ax1, ax3, ax5, ax6, ax7):
        ax.invert_yaxis()
    ax2.xaxis.tick_bottom()
    ax2.set_xticklabels(labels=units.index, rotation=40)
    return fig
