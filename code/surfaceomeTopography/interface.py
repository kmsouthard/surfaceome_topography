"""The interface between two cells: which interactions can form across it, how much of it each
can occupy, and the gap each one sets between the membranes.

Each cell's surface is taken in *interacting units* -- the proteins it expresses, with a
CellphoneDB complex counted as one unit -- and paired through the trans interaction table
(`human_interaction_heights_string_cellphoneDB_well_20220228_fcr.csv`). An interaction occupies
the interface only as far as both partners are present, so it is weighted by their shares of
their own cell's surface (``percent_expression``). Three weights are reported:

    interface_limiting   the scarcer partner's share (``limiting_expression``), normalised so the
                         interface sums to 100 with each limiting protein counted once per side
    interface_prob       the product of the two shares (``interface_mult``), normalised over
                         unique interactions: how often an encounter between the two surfaces
                         is this pair
    interface_shared     the pair's share of bound complexes when each protein's abundance is
                         shared among all of its partners, so none is counted in full by each:
                         the competitive binding equilibrium, one binding strength for every
                         pair, at strong binding

An interaction row names each partner by an accession. A protein matches its own; a complex
matches any of its subunits', so an interaction recorded against either chain of an integrin
reaches the integrin. A CellphoneDB row that names a complex (``partner_a``/``partner_b``) matches
that complex only, never another that shares the subunit; if that complex is not on the cell, the
row falls back to the protein it names, alone. A complex can be missing only because a small
subunit went undetected -- DAP12, SIRPB1's signalling partner, has no row in the E-PROT-1
proteomics -- and the partner that binds is still there. Where several rows reach the same pair
of units the tallest is kept, so each pair is counted once and a complex stands at its tallest
subunit -- the rule the CellphoneDB mapping uses for a complex's height.

`clearance` sets each interaction's height against each protein's: a trans interaction holds the
membranes about its own height apart, and a protein taller than that gap is pushed out of it.
`contact_exclusion` turns that into where on the contact each protein fits.

Cis pairs, two proteins bound on the same membrane, are not sized: `cis_pairs_at_contact` reports
only which of them have a partner bound across the contact, and which the contact splits, one
partner fitting beneath the gap and the other excluded.

Two assumptions stand behind every weight here: all pairs bind equally well (affinities are
unknown for most pairs, and at a membrane they shift with crowding), and heights are the extended
estimates, upper bounds for flexible proteins, so gaps and exclusion are both overstated.

Column names follow the Cell-Cell notebooks: the two cells are named by ``names`` and their
columns carry it as a suffix (``percent_expression_immune``, ``Entry_cancer``).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = ["expressed_complexes", "surface_units", "select_interactions", "pair_interactions",
           "antibody_bridge_heights", "interface_weights", "clearance", "contact_exclusion",
           "cis_pairs_at_contact", "plot_contact_bullseye"]

_SUBUNIT = ["ID link", "expression", "percent_expression", "sample"]
_SUBUNITS = ["uniprot_1", "uniprot_2", "uniprot_3", "uniprot_4"]
_UNIT = ["Entry", "sample", "expression", "percent_expression", "Gene Name", "complex_name"]


def expressed_complexes(long: pd.DataFrame, complexes: pd.DataFrame) -> pd.DataFrame:
    """CellphoneDB complexes expressed in each sample of a long expression table.

    A complex is expressed when its first two subunits are, and its third if it has one; its
    ``expression`` and ``percent_expression`` are the lower of its first two subunits'. A fourth
    subunit is not required. Each row carries subunit 1's own columns (``Entry``, height, names)
    and every subunit's accession (``uniprot_1`` to ``uniprot_4``); ``complex_height`` is its
    tallest subunit's height.
    """
    units = complexes[["complex_name", "uniprot_1", "uniprot_2", "uniprot_3", "uniprot_4"]]
    c = pd.merge(units, long, right_on="ID link", left_on="uniprot_1")
    c = pd.merge(c, long[_SUBUNIT], right_on=["ID link", "sample"], left_on=["uniprot_2", "sample"],
                 suffixes=("_1", "_2"))
    c = pd.merge(c, long[_SUBUNIT], right_on=["ID link", "sample"], left_on=["uniprot_3", "sample"],
                 how="left")
    c = pd.merge(c, long[_SUBUNIT], right_on=["ID link", "sample"], left_on=["uniprot_4", "sample"],
                 suffixes=("_3", "_4"), how="left")
    c = c.dropna(subset=["uniprot_1", "uniprot_2", "expression_1", "expression_2"])
    c = c[~(~c.uniprot_3.isna() & c.expression_3.isna())]
    c = c.assign(expression=c[["expression_1", "expression_2"]].min(axis=1))
    c = c.assign(percent_expression=c[["percent_expression_1", "percent_expression_2"]].min(axis=1))
    height = long.drop_duplicates("ID link").set_index("ID link")["total_height"]
    return c.assign(complex_height=c[_SUBUNITS].apply(lambda subunit: subunit.map(height)).max(axis=1))


def surface_units(long: pd.DataFrame, complexes_expressed: pd.DataFrame, sample: str,
                  include=()) -> pd.DataFrame:
    """One cell's surface in interacting units.

    Its expressed complexes, and the proteins it expresses that are not the first or second
    subunit of one. ``include`` adds proteins back by a substring of their entry name even when
    they are in a complex -- the HER2 panels add ERBB2 so that an antibody against it can bind.
    """
    cplx = complexes_expressed[complexes_expressed["sample"] == sample]
    expressed = long[long["sample"] == sample].sort_values("percent_expression", na_position="first")
    in_complex = expressed["ID link"].isin(cplx["uniprot_1"]) | expressed["ID link"].isin(cplx["uniprot_2"])
    single = expressed[~in_complex].dropna(subset=["percent_expression"])
    extra = [expressed[expressed["Entry name"].str.contains(name)] for name in include]
    return pd.concat([cplx, single, *extra], sort=False)


def select_interactions(interactions: pd.DataFrame, antibody_target: str | None = None) -> pd.DataFrame:
    """The interactions that can form across the interface.

    Rows with ``_merge == 'FcR'`` are antibody bridges: an Fc-gamma receptor on the immune cell
    bound through an antibody to an antigen on the other cell. Without ``antibody_target`` there
    is no antibody and every bridge is dropped; with it, only bridges to that antigen (a
    substring of its entry name, such as ``'ERBB2'``) are kept.
    """
    bridge = interactions["_merge"] == "FcR"
    if antibody_target is None:
        return interactions[~bridge]
    target = (interactions["Entry_prot1"].str.contains(antibody_target)
              | interactions["Entry_prot2"].str.contains(antibody_target))
    return interactions[~(bridge & ~target)]


def _match_keys(units: pd.DataFrame) -> pd.DataFrame:
    """One row per accession a unit answers to: a protein its own, a complex each subunit's."""
    subunits = units.reindex(columns=_SUBUNITS).to_numpy()
    keys = [[acc for acc in row if isinstance(acc, str)] if isinstance(name, str) else [entry]
            for row, name, entry in zip(subunits, units["complex_name"], units["Entry"])]
    return units.assign(matched=keys).explode("matched")


def _named_complex_only(pairs: pd.DataFrame, side: str, complex_col: str,
                        complexes_on_cell: set) -> pd.DataFrame:
    """Keep a CellphoneDB row's match on ``side`` only if it is the complex the row names.

    When that complex is not on the cell, the single protein the row names matches instead.
    """
    partner = pairs[f"partner_{'a' if side == 'prot1' else 'b'}"]
    names_complex = partner.notna() & (partner != pairs[f"Human ID link_{side}"])
    named = pairs[complex_col] == partner
    fallback = ~partner.isin(complexes_on_cell) & pairs[complex_col].isna()
    return pairs[~names_complex | named | fallback]


def pair_interactions(cell_a: pd.DataFrame, cell_b: pd.DataFrame, interactions: pd.DataFrame,
                      names=("immune", "cancer")) -> pd.DataFrame:
    """Every interaction with one partner on each cell, once per pair of units.

    The table is taken in both orientations. See the module docstring for how a row matches a
    complex, and which row is kept when several reach the same pair. ``matched_<name>`` records
    the accession each side matched through.
    """
    a, b = names
    unit_a = _match_keys(cell_a)[_UNIT + ["matched"]].dropna(subset=["expression"])
    unit_b = _match_keys(cell_b)[_UNIT + ["matched"]].dropna(subset=["expression"])
    on_a, on_b = set(unit_a["complex_name"].dropna()), set(unit_b["complex_name"].dropna())
    sides = []
    for side_a, side_b in (("prot1", "prot2"), ("prot2", "prot1")):
        with_a = pd.merge(unit_a, interactions, left_on="matched", right_on=f"Human ID link_{side_a}")
        with_a = _named_complex_only(with_a, side_a, "complex_name", on_a)
        pairs = pd.merge(unit_b, with_a, left_on="matched", right_on=f"Human ID link_{side_b}",
                         suffixes=(f"_{b}", f"_{a}"))
        sides.append(_named_complex_only(pairs, side_b, f"complex_name_{b}", on_b))
    pairs = pd.concat(sides, ignore_index=True)
    units = pd.DataFrame({side: pairs[f"complex_name_{side}"].fillna(pairs[f"Entry_{side}"])
                          for side in (a, b)})
    tallest_first = pairs["interaction_dim"].sort_values(ascending=False, kind="mergesort").index
    keep = units.loc[tallest_first].drop_duplicates().index
    return pairs.loc[sorted(keep)]


def antibody_bridge_heights(pairs: pd.DataFrame, epitope_height=0.0) -> pd.DataFrame:
    """The gap an antibody bridge holds: the Fc receptor, the antibody, and where the epitope is.

    A bridge's height in the interaction table is antigen + Fc receptor + 3.24 nm, the antibody's
    extra length (PDB 5DK3 and 4X4M; 03_interaction_heights). The antibody holds the antigen at
    its epitope, not at its top, so the gap is the epitope's height above the membrane plus the
    receptor and the antibody -- and never less than the antigen itself, which has to fit beneath:

        gap = max(epitope_height + receptor + 3.24 nm, antigen height)

    ``epitope_height`` is one height for every bridge, or a mapping from the antigen's accession
    (``Human ID link_prot1``) to its antibody's, as ``antibody_epitope_heights.csv`` records them
    (`epitope.curated_epitope_heights`); an antigen the mapping lacks is taken at 0, a
    membrane-proximal epitope, which was the 2020 rule for every bridge. Trastuzumab binds HER2's
    domain IV, 1.6 nm up, and its gap is 11.6-15.3 nm across the Fc receptors, near the measured
    bridge (trastuzumab on HER2 7.9 nm plus FcgR 6.2 nm, Son et al. 2020) and Bakalar et al.
    2018's epitope height + 11.5 nm.
    """
    if "total_height_prot1" not in pairs:
        raise KeyError("antibody_bridge_heights needs the antigen's height, total_height_prot1, "
                       "from the interaction table; keep that column when loading it")
    bridge = pairs["_merge"] == "FcR"
    antigen = pairs["total_height_prot1"]        # on bridge rows prot1 is the antigen, prot2 the receptor
    if not np.isscalar(epitope_height):
        epitope_height = pairs["Human ID link_prot1"].map(epitope_height).fillna(0.0)
    held = pairs["interaction_dim"] - antigen + epitope_height
    return pairs.assign(interaction_dim=np.where(bridge, np.maximum(held, antigen), pairs["interaction_dim"]))


def _shared_abundance(pairs: pd.DataFrame, a: str, b: str, strength: float = 1e6,
                      tol: float = 1e-12, max_iter: int = 100000) -> pd.Series:
    """Each pair's share of bound complexes, percent, with proteins shared among their partners.

    Mass action with conservation: a complex forms as ``strength`` x free a x free b, and each
    unit's free amount plus its complexes equals its share of its own cell's surface. One strength
    for every pair; at 1e6, with shares as fractions, the result no longer depends on it to four
    significant figures -- the strong-binding limit.
    """
    if pairs.empty:
        return pd.Series(dtype=float, index=pairs.index)
    ia, _ = pd.factorize(pairs[f"complex_name_{a}"].fillna(pairs[f"Entry_{a}"]))
    ib, _ = pd.factorize(pairs[f"complex_name_{b}"].fillna(pairs[f"Entry_{b}"]))
    total_a = pd.Series(pairs[f"percent_expression_{a}"].to_numpy() / 100).groupby(ia).first().to_numpy()
    total_b = pd.Series(pairs[f"percent_expression_{b}"].to_numpy() / 100).groupby(ib).first().to_numpy()
    free_a, free_b = total_a.copy(), total_b.copy()
    for _ in range(max_iter):
        next_a = total_a / (1 + strength * np.bincount(ia, weights=free_b[ib], minlength=len(total_a)))
        next_b = total_b / (1 + strength * np.bincount(ib, weights=next_a[ia], minlength=len(total_b)))
        next_a, next_b = np.sqrt(free_a * next_a), np.sqrt(free_b * next_b)
        done = (np.allclose(next_a, free_a, rtol=tol, atol=0)
                and np.allclose(next_b, free_b, rtol=tol, atol=0))
        free_a, free_b = next_a, next_b
        if done:
            break
    else:
        raise RuntimeError("binding equilibrium did not converge")
    bound = strength * free_a[ia] * free_b[ib]
    return pd.Series(bound / bound.sum() * 100, index=pairs.index)


def interface_weights(pairs: pd.DataFrame, names=("immune", "cancer")) -> pd.DataFrame:
    """Weight each paired interaction by how much of the interface it can occupy.

    Adds ``limiting_expression`` (the scarcer partner's share), ``limiting_cell``,
    ``combined_expression``, ``interface_limiting``, ``interface_mult``, ``interface_prob`` and
    ``interface_shared``; see the module docstring. A partner that is a complex is named by its complex name.
    """
    a, b = names
    shares = [f"percent_expression_{b}", f"percent_expression_{a}"]
    for side in (b, a):
        gene, cplx = f"Gene Name_{side}", f"complex_name_{side}"
        pairs = pairs.assign(**{gene: np.where(pairs[cplx].isnull(), pairs[gene], pairs[cplx])})
    pairs = pairs.assign(limiting_expression=pairs[shares].min(axis=1))
    pairs = pairs.assign(limiting_cell=pairs[shares].idxmin(axis=1))
    pairs = pairs.assign(combined_expression=pairs[shares].sum(axis=1))

    limiting = [pairs[pairs["limiting_cell"] == f"percent_expression_{side}"]
                .drop_duplicates(subset=[f"Gene Name_{side}"]).limiting_expression.sum()
                for side in (b, a)]
    pairs = pairs.assign(interface_limiting=pairs.limiting_expression / (limiting[0] + limiting[1]) * 100)

    pairs = pairs.assign(interface_mult=pairs[shares[0]] * pairs[shares[1]])
    total = (pairs.sort_values("interface_mult", ascending=False)
             .drop_duplicates(subset=["Human interaction_id", f"complex_name_{a}", f"complex_name_{b}",
                                      f"Entry_{a}", f"Entry_{b}"], keep="first")
             .interface_mult.sum())
    pairs = pairs.assign(interface_prob=pairs.interface_mult / total * 100)
    return pairs.assign(interface_shared=_shared_abundance(pairs, a, b))


def clearance(interaction_heights: pd.Series, protein_heights: pd.Series) -> pd.DataFrame:
    """Interaction height minus protein height, for every interaction (rows) and protein (columns).

    Positive: the protein fits beneath the gap the interaction sets. Negative: it is taller, and
    is pushed out of that contact. Proteins sharing a label share a column, holding the last.
    """
    gaps = pd.DataFrame()
    for label, height in protein_heights.items():
        gaps[label] = interaction_heights - height
    return gaps


def contact_exclusion(surface: pd.DataFrame, pairs: pd.DataFrame, weight: str = "interface_shared",
                      margin: float = 5.0) -> pd.DataFrame:
    """Where on the contact each of one cell's proteins fits, by height.

    Each bound pair holds the membranes at its own height (``interaction_dim``) over its share of
    the contact (``weight``). A protein fits where it is no taller than that gap, is partly
    excluded where it is up to ``margin`` taller, and is excluded beyond: about 5 nm drives
    exclusion at reconstituted membrane interfaces (Schmid et al. 2016). Heights are extended
    estimates, so exclusion is overstated for flexible proteins.

    ``surface`` is one cell's units (`surface_units`); a complex stands at its tallest subunit.
    Returns one row per unit, with the percent of the contact where it fits, is partly excluded
    and is excluded.
    """
    units = _unit_heights(surface)
    gaps, share = _contact_shares(pairs, weight)
    h = units["unit_height"].to_numpy()[:, None]
    fits = ((h <= gaps) * share).sum(axis=1)
    excluded = ((h >= gaps + margin) * share).sum(axis=1)
    return pd.DataFrame({"name": units["unit_name"].to_numpy(),
                         "Entry": units["Entry"].to_numpy(), "height": units["unit_height"].to_numpy(),
                         "percent_expression": units["percent_expression"].to_numpy(),
                         "fits": fits, "partly_excluded": 100 - fits - excluded,
                         "excluded": excluded}).sort_values("height", ignore_index=True)


def _unit_heights(surface: pd.DataFrame) -> pd.DataFrame:
    """One cell's units, once each, with ``unit`` (its complex name or accession), ``unit_name``
    (its complex or gene name) and ``unit_height``: a complex stands at its tallest subunit."""
    units = surface.dropna(subset=["percent_expression"])
    units = units.assign(unit=units["complex_name"].fillna(units["Entry"])).drop_duplicates("unit")
    height = units["total_height"]
    if "complex_height" in units:
        height = units["complex_height"].where(units["complex_name"].notna() & units["complex_height"].notna(), height)
    return units.assign(unit_name=units["complex_name"].fillna(units["Gene Name"]), unit_height=height)


def _contact_shares(pairs: pd.DataFrame, weight: str):
    """Each bound pair's gap, nm, and its share of the contact, percent."""
    return pairs["interaction_dim"].to_numpy(), pairs[weight].to_numpy() / pairs[weight].sum() * 100


def cis_pairs_at_contact(surface: pd.DataFrame, cis_pairs: pd.DataFrame, pairs: pd.DataFrame, side: str,
                         weight: str = "interface_shared", margin: float = 5.0) -> pd.DataFrame:
    """The cis pairs on one cell at a contact: which are bound across it, and which it splits.

    ``cis_pairs`` lists two proteins that bind on the same membrane (``Human ID link_prot1``,
    ``Human ID link_prot2``, ``subtype``). A pair is on the cell when both partners are, each
    matching a unit as in `pair_interactions`: a protein by its own accession, a complex by any
    subunit's. Two subunits of one complex on the cell are already that complex and are not listed;
    a pair reaching the same two units through several rows is listed once.

    ``pairs`` are the trans pairs at the contact (`interface_weights`), and ``side`` names this
    cell among their columns (``'immune'`` or ``'cancer'``). Partner 1 is the shorter. For each pair:

        bound_1, bound_2   the share of the contact, percent, held by trans pairs that partner is
                           in. Above 0, it is bound across the contact, and holds its cis partner
                           there with it.
        split              the share of the contact, percent, where partner 1 fits beneath the gap
                           and partner 2 is excluded from it, by the rule of `contact_exclusion`: a
                           cis pair held together there would have to bend, tilt or come apart.

    Nothing else is modelled: the cis bond does not change either partner's height, its abundance
    or any trans pair's weight.
    """
    units = _unit_heights(surface)
    keys = _match_keys(units)[["matched", "unit", "unit_name", "unit_height"]]
    on_cell = (cis_pairs.merge(keys.add_suffix("_1"), left_on="Human ID link_prot1", right_on="matched_1")
                        .merge(keys.add_suffix("_2"), left_on="Human ID link_prot2", right_on="matched_2"))
    pair = on_cell["Human ID link_prot1"] + "," + on_cell["Human ID link_prot2"]
    on_cell = on_cell[~pair.isin(pair[on_cell["unit_1"] == on_cell["unit_2"]])]
    swap = on_cell["unit_height_1"] > on_cell["unit_height_2"]
    for col in ("unit", "unit_name", "unit_height"):
        first, second = on_cell[f"{col}_1"], on_cell[f"{col}_2"]
        on_cell = on_cell.assign(**{f"{col}_1": first.where(~swap, second), f"{col}_2": second.where(~swap, first)})
    on_cell = on_cell.drop_duplicates(["unit_1", "unit_2"])

    gaps, share = _contact_shares(pairs, weight)
    held = pd.Series(share, index=pairs.index).groupby(pairs[f"complex_name_{side}"].fillna(pairs[f"Entry_{side}"])).sum()
    short, tall = on_cell["unit_height_1"].to_numpy()[:, None], on_cell["unit_height_2"].to_numpy()[:, None]
    split = (((short <= gaps) & (tall >= gaps + margin)) * share).sum(axis=1)
    return pd.DataFrame({"name_1": on_cell["unit_name_1"].to_numpy(), "name_2": on_cell["unit_name_2"].to_numpy(),
                         "subtype": on_cell["subtype"].to_numpy(),
                         "height_1": on_cell["unit_height_1"].to_numpy(), "height_2": on_cell["unit_height_2"].to_numpy(),
                         "bound_1": on_cell["unit_1"].map(held).fillna(0).to_numpy(),
                         "bound_2": on_cell["unit_2"].map(held).fillna(0).to_numpy(),
                         "split": split}).sort_values(["split", "name_1", "name_2"], ascending=[False, True, True],
                                                      ignore_index=True)


#: Gap bands of the contact bullseye, nm, and its ring colours from the outside in.
BULLSEYE_BANDS = (10, 20, 30, 40, 50)
BULLSEYE_COLORS = ("#8C7976", "#D8CECD", "#8CA592", "#CDD9D3", "#698C79", "#FFFFFF")


def plot_contact_bullseye(pairs: pd.DataFrame, weight: str = "interface_shared",
                          bands=BULLSEYE_BANDS, names=("immune", "cancer"), colors=BULLSEYE_COLORS):
    """The contact as nested circles, one ring per band of gap heights.

    Each circle's area is the share of the contact held at gaps below a band's upper edge, so the
    outer circle is the whole contact and each ring is one band: 50 nm and over, 40-50, ..., 10-20,
    with the contact under 10 nm left white at the centre. The legend names each ring by the pair
    holding most of it, with the ring's share of the contact and that pair's share of the ring.
    ``colors`` gives the rings from the outside in, then the centre. Returns the figure.
    """
    import matplotlib.pyplot as plt

    a, b = names
    share = pairs[weight] / pairs[weight].sum()
    gap = pairs["interaction_dim"]
    label = pairs[f"Gene Name_{a}"].astype(str) + "–" + pairs[f"Gene Name_{b}"].astype(str)
    edges = [np.inf, *sorted(bands, reverse=True)]
    fig, ax = plt.subplots(figsize=(5, 5))
    handles, legend = [], []
    for i, (upper, color) in enumerate(zip(edges, colors)):
        circle = plt.Circle((0.5, 0.5), np.sqrt(share[gap < upper].sum() / np.pi), color=color)
        ax.add_artist(circle)
        if i + 1 < len(edges):
            lower = edges[i + 1]
            band = (gap >= lower) & (gap < upper)
            span = f"{lower:g}+ nm" if np.isinf(upper) else f"{lower:g}–{upper:g} nm"
            top = share[band].groupby(label[band]).sum()
            handles.append(circle)
            legend.append(f"{span}, {share[band].sum():.0%} of contact: {top.idxmax()} "
                          f"({top.max() / top.sum():.0%})" if len(top) else f"{span}, none")
    ax.set_xlim((-0.1, 1.1))
    ax.set_ylim((-0.1, 1.1))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.legend(handles, legend, loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False)
    return fig
