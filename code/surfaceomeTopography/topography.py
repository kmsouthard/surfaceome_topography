"""A cell's surface topography: the heights of the proteins it expresses, weighted by expression.

A cell type's surface is the surfaceome proteins it expresses. Each protein's share of that surface
is its expression over the cell's total, in percent (``percent_expression``), and the cell's
topography is the height distribution weighted by those shares.

`cell_surface` builds one cell type's surface from a table with one row per protein and a column per
cell type (or per replicate). `percent_expression` does the same share calculation for a long table
holding many cell types at once, as the Cell-Cell panels use. Both give each protein's share of its
own cell's total.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = ["HEIGHT_BINS", "OVERLAY_STYLE", "cell_surface", "percent_expression", "mean_height",
           "plot_surfaces"]

#: Log-spaced height bins, 1 nm to 1.4 um, shared by every topography histogram.
HEIGHT_BINS = np.logspace(np.log10(1), np.log10(1400), 30)

#: Overlaid surfaces, in order: histogram colour, label colour, mean-height box colour.
OVERLAY_STYLE = (("#59656F", "#59656F", (0.8, 0.85, 0.85)),
                 ("#CD489A", "#CD489A", (0.9, 0.5, 0.7)),
                 ("#00A8A8", "#0B8481", (0.0, 0.7, 0.67)))


def cell_surface(table: pd.DataFrame, expression, columns=None) -> pd.DataFrame:
    """One cell type's surface: the proteins it expresses, each with its share of the total.

    ``table`` has one row per protein. ``expression`` names the cell type's expression column, or a
    list of replicate columns, which are combined by their median (``expression_median``, with
    ``expression_std`` beside it). A protein is on the surface only if every named column has a
    value. ``columns`` selects and orders the columns kept; by default all of them.
    """
    expression = [expression] if isinstance(expression, str) else list(expression)
    surface = (table if columns is None else table[list(columns)]).dropna(subset=expression)
    if len(expression) > 1:
        replicates = surface[expression]
        surface = surface.assign(expression_median=replicates.median(axis="columns"),
                                 expression_std=replicates.std(axis="columns"))
        level = surface["expression_median"]
    else:
        level = surface[expression[0]]
    return surface.assign(percent_expression=level / level.sum() * 100)


def percent_expression(long: pd.DataFrame, sample: str = "sample", expression: str = "expression",
                       protein: str = "ID link") -> pd.Series:
    """Each row's share of its sample's total surface expression, in percent.

    For a long table with one row per protein and sample. A protein counts once toward its sample's
    total however many rows it has.
    """
    totals = long.drop_duplicates([sample, protein]).groupby(sample)[expression].sum()
    return long[expression] / long[sample].map(totals) * 100


def mean_height(surface: pd.DataFrame, height: str = "total_height",
                weight: str = "percent_expression") -> float:
    """Expression-weighted mean height of a surface, nm."""
    return np.average(surface[height], weights=surface[weight])


def _histogram(ax, surface, color, height="total_height", weight="percent_expression"):
    ax.hist(surface[height], bins=HEIGHT_BINS, histtype="step", color=color,
            weights=surface[weight])


def _rug(ax, surface, color, height="total_height"):
    x = surface[height]
    ax.plot(x, np.zeros(x.shape), "|", ms=10, color=color)


def _log_height_axis(ax):
    # after everything is drawn: switching to log first changes how the data limits are taken
    ax.set_xscale("log")
    ax.set_xlabel("log height (nm)", fontsize=14)
    ax.set_ylabel("count")


def _tag(ax, text, xy, fc, color=None):
    ax.annotate(text, xy=xy, xycoords="axes fraction", xytext=xy, textcoords="offset points",
                size=14, color=color, bbox=dict(boxstyle="round", fc=fc, ec="none"))


def plot_surfaces(surfaces: dict, overlay: bool = False, height: str = "total_height",
                  weight: str = "percent_expression"):
    """Topography figure for one or more cell surfaces; ``surfaces`` maps a label to a surface.

    Stacked (the default) gives each surface its own panel, all sharing both axes. ``overlay``
    draws up to three in one panel, coloured by `OVERLAY_STYLE`, with the first surface's rug.
    Each surface is labelled with its expression-weighted mean height. ``height`` and ``weight``
    name the columns plotted, so the same figure shows a contact's gaps (``interaction_dim``
    weighted by an interface weight). Returns the figure.
    """
    import matplotlib.pyplot as plt

    if overlay:
        if len(surfaces) > len(OVERLAY_STYLE):
            raise ValueError(f"overlay draws at most {len(OVERLAY_STYLE)} surfaces")
        fig = plt.figure(figsize=(7, 4), constrained_layout=True)
        ax = fig.add_subplot(fig.add_gridspec(1, 1, figure=fig)[0:1, 0])
        for surface, (color, _, _) in zip(surfaces.values(), OVERLAY_STYLE):
            _histogram(ax, surface, color, height, weight)
        _rug(ax, next(iter(surfaces.values())), OVERLAY_STYLE[0][0], height)
        _log_height_axis(ax)
        for i, ((label, surface), (_, label_color, mean_fc)) in enumerate(
                zip(surfaces.items(), OVERLAY_STYLE)):
            y = 0.9 - 0.1 * i
            _tag(ax, mean_height(surface, height, weight).round(1), (0.86, y), mean_fc)
            _tag(ax, label, (0.05, y), (0.95, 0.95, 0.95), color=label_color)
        return fig

    n = len(surfaces)
    fig = plt.figure(figsize=(5, 2.5 * n), constrained_layout=True)
    gs = fig.add_gridspec(n, 1, figure=fig)
    first = None
    for i, (label, surface) in enumerate(surfaces.items()):
        ax = fig.add_subplot(gs[i:i + 1, 0], sharex=first, sharey=first)
        first = first or ax
        _histogram(ax, surface, "#59656F", height, weight)
        _rug(ax, surface, "#59656F", height)
        _log_height_axis(ax)
        _tag(ax, mean_height(surface, height, weight).round(1), (0.86, 0.9), (0.9, 0.5, 0.7))
        _tag(ax, label, (0.1, 0.8), (0.9, 0.9, 0.9))
    return fig
