"""Resolve overlapping domain annotations into one non-overlapping assignment per ECD.

HMMER reports every profile that matches, so a single stretch of sequence often carries several
competing Pfam assignments -- and Pfam 38.2 makes this worse than Pfam 34 did, because it adds
long composite families that span the individual repeats they contain. Choosing among them is
a maximum-weight independent set problem on the overlap graph, solved here by depth-first search
over connected assignments, keeping the highest bit score.

Two consequences worth knowing:

* The selection is by **score**, not by length or by whether the family has a measured height.
  A long heightless family can therefore displace several short ones that do have heights --
  which is why `build_domain_assignments.py` drops heightless families before this runs.
* `allowed_overlap` (default 4 residues) tolerates small alignment-edge overlaps rather than
  treating them as genuine competition.
"""

## Domain assigment for Human Extracellular domains

#Using the domain level output from the hmmscan of the Pfam-A hmm profiles, assign each extracellular domain, non overlapping domains.

#Steps:
#1. Find best domain assignment for each ECD
#2. pick that domain
#3. find next best, check for overlap
#4. if no overlap keep, if overlap detected discard

#Finding the best domains with a little graph therory:
#https://stackoverflow.com/questions/54969074/python-3-remove-overlaps-in-table
import numpy as np
import pandas as pd


def connections(graph, id):
    def dict_to_df(d):
        df = pd.DataFrame(data=[d.keys(), d.values()], index=['ID', 'Subgraph']).T
        df['id'] = id
        return df[['id', 'Subgraph', 'ID']]

    def dfs(node, num):
        visited[node] = num
        for _node in graph.loc[node].iloc[0]:
            if _node not in visited:
                dfs(_node, num)

    visited = {}
    graph = graph.loc[id]
    for (num, node) in enumerate(graph.index):
        if node not in visited:
            dfs(node, num)

    return dict_to_df(visited)

def find_overlaps(data, allowed_overlap = 4, columns = ('align_from', 'align_to'), verbose = True):
    """Shrink each annotation by ``allowed_overlap`` at both ends, so that annotations
    overlapping by no more than that are not treated as competing.

    The tolerance is applied **per row**. It used to be applied to the whole frame inside a
    bare ``try/except``::

        try:                 x[start] + allowed_overlap, x[end] - allowed_overlap
        except: try:         +2 / -2
                except: try: +1 / -1
                        except: +0 / -0

    ``df.apply`` is all-or-nothing, so one row too short to survive the shrink -- it needs
    ``end - start > 2 * allowed_overlap`` -- sent **every** row in the frame down a level, and
    the cascade bottoms out at zero. A two-residue multipass loop is enough to trigger it, and
    the surfaceome has 84 ectodomains under 0.5 nm made of exactly those. The stated
    ``allowed_overlap`` of 4 (cells 19, 31, 46, 62) and 1 (cell 70) was therefore **in effect
    nowhere in the pipeline** -- every overlap resolution ran at zero tolerance.

    That is not a cosmetic difference. At zero tolerance a single shared residue makes two
    annotations compete, and the loser is discarded outright: it is why a disordered region
    reaching one residue into an AlphaFold model lost its entire polymer contribution
    (see ``trim_to_unmodelled``), and it applies equally to domain-versus-domain
    resolution everywhere else.

    Here each row takes the largest tolerance it can support, ``min(allowed_overlap,
    (span - 1) // 2)``, so a short row degrades alone instead of dragging the frame with it.
    Rows that can carry the intended tolerance always get it. Equal endpoints are legal --
    ``pd.Interval(5, 5, closed='neither')`` is empty and overlaps nothing, which is the right
    reading of an annotation too short to have an interior.
    """
    #make copy of data so not mutating it
    df = data.copy()
    #set an id for each row
    df.loc[:,'ID'] = range(df.shape[0])

    lo = pd.to_numeric(df[columns[0]], errors='coerce')
    hi = pd.to_numeric(df[columns[1]], errors='coerce')
    span = hi - lo

    #largest tolerance each row can carry, capped at the one asked for
    tol = np.minimum(allowed_overlap, np.maximum(0, (span - 1) // 2))
    tol = tol.fillna(0)

    if verbose:
        short = int((tol < allowed_overlap).sum())
        if short:
            print(f"  find_overlaps: {short:,} of {len(df):,} rows are too short for "
                  f"allowed_overlap={allowed_overlap} and use a smaller one "
                  f"(min {tol.min():.0f}); the rest keep {allowed_overlap}")

    a = (lo + tol).to_numpy()
    b = (hi - tol).to_numpy()
    df.loc[:,'Interval'] = [
        pd.Interval(x, y, closed='neither') if (pd.notna(x) and pd.notna(y) and y >= x)
        else pd.Interval(0.0, 0.0, closed='neither')
        for x, y in zip(a, b)
    ]

    return df

def overlap_graph(overlap_df, grouping_var = 'target_name'):

    columns = [grouping_var, 'Interval', 'ID']
    connected = overlap_df[columns].merge(overlap_df[columns], on=grouping_var)
    connected['Overlap'] = connected.apply(lambda x: x['Interval_x'].overlaps(x['Interval_y']), axis=1)
    connected = connected.loc[connected['Overlap'] == True, [grouping_var, 'ID_x', 'ID_y']]

    graph = connected.groupby([grouping_var, 'ID_x']).agg(list)

    return graph

def find_connections(graph):
    dfs = []
    for id in graph.index.get_level_values(0).unique():
        dfs.append(connections(graph, id))

    conns = pd.concat(dfs)

    return conns

def annotate_subgraphs(connect, overlaps):
    #merge the connections and overlaps dataframes so each subgraph is annotated
    data = overlaps.merge(connect[['Subgraph', 'ID']], on=['ID'])

    return data


def select_max(x, col):
    m = x[col].min()
    if len(x) > 1 and (x[col] == m).all():
        return -1
    else:
        return x[col].idxmax()


#selected = data.groupby(['target_name', 'Subgraph'])['score_dom', 'ID'].apply(select_max)
#selected = selected[selected >= 0]


def trim_to_unmodelled(df, *, id_col="ID", start_col="dom_start", end_col="dom_end",
                       len_col="dom_len", model_mask_col="mean_confidence_trimmed",
                       disorder_mask_col="disorder_start"):
    """Clip disordered regions off the structurally modelled span rather than losing them.

    A disordered region that runs a few residues into an AlphaFold model's confident span is
    **discarded outright** by the overlap resolution, because the model row carries the
    maximum ``score_dom`` and wins its whole subgraph. The region contributes nothing at all
    and its residues fall through to the sequence-length estimate of 0.04 nm/aa.

    The size of that cliff is not the tolerance it appears to be. ``find_overlaps`` wraps its
    interval construction in a bare ``try/except`` and, if *any* row in the frame has
    ``end - start <= 2 * allowed_overlap`` -- a two-residue multipass loop will do it --
    silently degrades the tolerance for the **entire frame**, ultimately to zero. So the
    stated ``allowed_overlap`` of 1 or 4 is not in effect, and a **single residue** of overlap
    is enough to discard a region of any length.

    NPFF2 is the worked example. Its ectodomain is 1-147 and its model is confident from 133:

        2019 MobiDB   disorder (1, 133)   132 aa   no overlap    polymer, 17.78 nm
        2026 MobiDB   disorder (1, 135)   134 aa   overlaps by 2 discarded, 0 nm

    Two residues of annotation cost 132 residues of polymer and 12.4 nm of height. Refreshing
    MobiDB puts 100% of one affected group across that boundary, by a median of 4 residues,
    which is how the cliff came to light -- but it sits under the 2022 published numbers too.

    This trims instead: a region overlapping the modelled span by k residues loses k residues,
    not all of them. Regions that would be split in two by an interior model keep their longer
    piece. Regions swallowed whole are dropped, which is correct -- those residues *are*
    modelled.

    Returns a copy with ``start_col``/``end_col``/``len_col`` updated on disorder rows.
    """
    import numpy as np
    import pandas as pd

    out = df.copy()
    is_model = out[model_mask_col].notna() if model_mask_col in out.columns else pd.Series(False, index=out.index)
    is_disorder = out[disorder_mask_col].notna() if disorder_mask_col in out.columns else pd.Series(False, index=out.index)

    spans = {}
    for key, lo, hi in zip(out.loc[is_model, id_col],
                           out.loc[is_model, start_col],
                           out.loc[is_model, end_col]):
        if pd.notna(lo) and pd.notna(hi):
            spans.setdefault(key, []).append((float(lo), float(hi)))

    trimmed = dropped = 0
    for i in out.index[is_disorder]:
        key = out.at[i, id_col]
        a, b = out.at[i, start_col], out.at[i, end_col]
        if key not in spans or pd.isna(a) or pd.isna(b):
            continue
        a, b = float(a), float(b)
        for lo, hi in spans[key]:
            if b <= lo or a >= hi:
                continue                      # no overlap with this model span
            left, right = (a, min(b, lo)), (max(a, hi), b)
            pieces = [(x, y) for x, y in (left, right) if y - x >= 1]
            if not pieces:
                a = b = np.nan                # the model covers it entirely
                break
            a, b = max(pieces, key=lambda p: p[1] - p[0])
        if pd.isna(a):
            out.at[i, len_col] = 0.0
            dropped += 1
        elif (a, b) != (float(out.at[i, start_col]), float(out.at[i, end_col])):
            out.at[i, start_col] = a
            out.at[i, end_col] = b
            out.at[i, len_col] = b - a
            trimmed += 1

    print(f"  trim_to_unmodelled: {trimmed:,} disordered regions clipped off a modelled span, "
          f"{dropped:,} fully covered and dropped")
    return out
