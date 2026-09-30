#!/usr/bin/env python
"""Extend the Pfam clan-height table to the families a Pfam refresh adds.

    python code/database/inputs/build_clan_heights.py --validate --oriented CSV --archived CSV
    python code/database/inputs/build_clan_heights.py --oriented CSV --archived CSV --pfam-dat FILE \
        --assignments CSV --out CSV

Domain heights come from `human_domain_clans_averaged_20200129.csv`, measured in 2020. Pfam
38.2 assigns **476 families it does not cover**, and `build_domain_assignments.py` drops them
— correctly, since keeping a heightless family lets it displace the repeats it spans and
collapse a protein's height. But dropping them is not free: they are **30% of all Pfam 38.2
assignment rows**, touching 815 of 1,554 proteins, and 123 proteins lose every domain they
have.

### The derivation, recovered

Not documented anywhere, so it was recovered from the two archived tables and reproduces them
exactly:

1. `height` per family = `max_dim_mean` if `orientation == 'max'` else `min_dim_mean`
   — **514 of 514 exact**. This is the one manual step: 231 families called `max`, 195 `min`,
   88 never called.
2. **PF00084 (Sushi) is excluded**, and nothing else is. Its measurement is contaminated —
   98.2 Å over 234 instances for a domain that is really ~4 nm, almost certainly whole chains
   measured rather than isolated domains.
3. `height_clan` = mean of family heights within a clan, **weighted by instance count**
   (`max_dim_len`). With Sushi excluded this reproduces **95 of 95 clans exactly**; leaving it
   in breaks CL0001 by 34 Å.
4. A family with no measurement of its own still receives its clan's height — `PF12946`
   carries `height_pfam = NaN` and `height_clan = 32.67` in the archived table.

### The extension

Step 4 is the whole point. Most of what Pfam 38.2 added is not new folds but **finer
subdivisions of families that already have heights** — `Ig_CLSTN`, `Ig_VEGFR-1-like_5th`,
`Cadherin_FAT4_N`, `EGF_Teneurin`, `FN3_DSCAM-DSCAML_C`. Of the 476:

| | families | assignment rows |
|---|---:|---:|
| in a clan that already has a measured height | **305** | **5,113 (89%)** |
| no clan, or a clan with no height | 171 | 633 (11%) |

So 89% of the loss is recovered by applying rule 4 to families the 2020 table predates. That
is not a new modelling decision — it is the rule the table already runs on.

PF00084 falls out neatly: still excluded from the *average* so its bad measurement cannot
pollute CL0001, but it now *receives* the CL0001 height of 32.7 Å ≈ 3.3 nm, which is about
right for a Sushi domain and replaces the contaminated 98.2 Å.

### The guard that makes inheritance safe

A clan height is the height of **one domain**, so only a family that really is one domain may
inherit it. Pfam 38.2 also adds *composite* families spanning several repeats, and giving one
of those a single-domain height would undercount badly — `PF25024 EGF_TEN` has a model length
of 234 against CL0001's measured median of 36, so it covers about six EGF repeats. It is also
the exact family that collapsed NOTCH1 from 88.4 nm to 19.8 nm when heightless families were
kept, and handing it 3.3 nm would collapse it again by the opposite route.

Model length separates the two cleanly:

| family | model length | clan median | ratio | |
|---|---:|---:|---:|---|
| `EGF_Teneurin` | 30 | 36 | 0.8 | a real subdivision |
| `Cadherin_FAT4_N` | 85 | 96 | 0.9 | a real subdivision |
| `Ig_CLSTN` | 123 | 96 | 1.3 | a real subdivision |
| **`EGF_TEN`** | **234** | **36** | **6.5** | **a composite** |

So a family inherits only if its model length is within `--max-length-ratio` (default 2.0) of
the median model length of its clan's *measured* members. Composites stay dropped, which is
exactly the status quo for them and therefore never worse.

Scaling a composite by its ratio instead — six EGF repeats at six times the EGF height — would
be defensible, since the domain-counting method sums per-domain heights anyway. That is a
modelling decision rather than a recovery of the existing rule, so it is not done here.

The families left over need structures measured, and only 45 of the 476 appear in the 2018
Pfam→PDB mapping at all, so the rest would need the RCSB search API. They stay dropped, and
`--out` writes a `.uncovered.txt` beside it so the limitation is explicit rather than silent.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SUSHI = "PF00084"


def pfam_clans(path: Path) -> dict[str, dict]:
    """{accession: {'clan':…, 'id':…, 'ml':…}} from Pfam-A.hmm.dat[.gz]."""
    opener = gzip.open if path.suffix == ".gz" else open
    recs, cur = {}, {}
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("//"):
                if cur.get("AC"):
                    recs[cur["AC"].split(".")[0]] = {"clan": cur.get("CL"),
                                                     "id": cur.get("ID"),
                                                     "ml": cur.get("ML")}
                cur = {}
                continue
            m = re.match(r"^#=GF\s+(\S+)\s+(.*)$", line.rstrip("\n"))
            if m:
                cur[m.group(1)] = m.group(2).strip()
    return recs


def family_heights(oriented: pd.DataFrame) -> pd.DataFrame:
    """Per-family height: the table's own ``height`` where it has one, else from the orientation call."""
    out = oriented.copy()
    out["height_calc"] = np.where(out.orientation == "max", out.max_dim_mean,
                          np.where(out.orientation == "min", out.min_dim_mean, np.nan))
    if "height" in out.columns:
        own = out["height"].notna()
        out.loc[own, "height_calc"] = out.loc[own, "height"]
    return out


def clan_heights(oriented: pd.DataFrame, exclude=(SUSHI,)) -> pd.Series:
    """Instance-count-weighted mean family height per clan."""
    g = oriented.dropna(subset=["height_calc"])
    g = g[~g.Pfam.isin(exclude)]
    return g.groupby("Clan_y").apply(
        lambda x: np.average(x.height_calc, weights=x.max_dim_len))


def rebuild(o: pd.DataFrame, arch: pd.DataFrame, a) -> int:
    """The clan-height table from a fresh family table, by the same four rules.

    1. A family's height is the mean of its instances' longest or shortest dimension, by its
       orientation call.  2. Clan height = the instance-count-weighted mean of its measured
       families (no family excluded: the Sushi exclusion answered a contamination the 2020
       whole-chain measurement had, and a profile-hit measurement cannot have).  3. A family
       with no height of its own, measured or not, inherits its clan's height if its model
       length is within --max-length-ratio of the clan's measured members.  Clan names are
       taken from --archived where a clan is already named, else the clan accession.
    """
    if not (a.pfam_dat and a.out):
        print("--rebuild needs --pfam-dat and --out", file=sys.stderr)
        return 2
    o = o.copy()
    if "Clan_y" not in o.columns:
        o["Clan_y"] = o["Clan_x"]
    recs = pfam_clans(a.pfam_dat)
    names = arch.dropna(subset=["Clan_y"]).drop_duplicates("Clan_y").set_index("Clan_y")["clan_name"].to_dict()
    ch = clan_heights(o, exclude=()).to_dict()
    o["ml"] = o.Pfam.map(lambda f: float((recs.get(f, {}) or {}).get("ml") or "nan"))
    clan_ml = o.dropna(subset=["height_calc"]).groupby("Clan_y").ml.median().to_dict()
    fam_name = o.set_index("Pfam")["query_name_x"].to_dict() if "query_name_x" in o else {}

    rows, composites, uncalled = [], [], []
    # Families the height builder found to span several domains of other families carry no
    # height and may not inherit one either: the domains they cover are the ones counted.
    flagged = set(o.loc[o.get("height_basis", pd.Series(dtype=str)).astype(str).str.startswith("composite"), "Pfam"]) \
        if "height_basis" in o.columns else set()
    for fam in sorted(flagged):
        rec = recs.get(fam, {})
        composites.append((fam, rec.get("id"), rec.get("clan"), float(rec.get("ml") or "nan"), np.nan, np.nan))
    o = o[~o.Pfam.isin(flagged)]
    measured = o.dropna(subset=["height_calc"])
    for r in measured.itertuples():
        clan = r.Clan_y if isinstance(r.Clan_y, str) else None
        rows.append({"Clan_y": clan, "clan_name": names.get(clan, clan) if clan else None,
                     "height_clan": ch.get(clan, np.nan) if clan else np.nan,
                     "height_pfam": r.height_calc, "Pfam": r.Pfam,
                     "query_name_y": fam_name.get(r.Pfam, recs.get(r.Pfam, {}).get("id")),
                     "height": r.height_calc})
    done = set(measured.Pfam)
    for fam in sorted(set(recs) | set(o.Pfam)):
        if fam in done or fam in flagged:
            continue
        rec = recs.get(fam, {})
        clan = rec.get("clan") or (o.loc[o.Pfam == fam, "Clan_y"].iloc[0] if fam in set(o.Pfam) else None)
        if not isinstance(clan, str) or clan not in ch:
            if fam in set(o.Pfam):
                uncalled.append(fam)
            continue
        ml = float(rec.get("ml") or "nan")
        cm = clan_ml.get(clan, np.nan)
        if ml == ml and cm == cm and cm > 0 and ml / cm > a.max_length_ratio:
            composites.append((fam, rec.get("id"), clan, ml, cm, ml / cm))
            continue
        rows.append({"Clan_y": clan, "clan_name": names.get(clan, clan), "height_clan": ch[clan],
                     "height_pfam": np.nan, "Pfam": fam, "query_name_y": rec.get("id"),
                     "height": ch[clan]})
    for r in rows:
        r["source"] = "2026 measurement" if r["height_pfam"] == r["height_pfam"] else "clan average"
    # A family the archive measured that no current human or mouse structure covers keeps
    # its archived measurement, marked as such, rather than dropping out of the table.
    kept = arch[arch.Pfam.isin(set(arch.Pfam) - {r["Pfam"] for r in rows} - {c[0] for c in composites})]
    kept = kept.dropna(subset=["height"])
    for r in kept.to_dict("records"):
        r["source"] = ("2020 measurement, no 2026 structure" if r["height_pfam"] == r["height_pfam"]
                       else "2020 clan average, no 2026 structure")
        rows.append(r)
    out = pd.DataFrame(rows)
    print(f"archived families carried over unmeasured: {len(kept)}")
    print(f"rebuilt from {len(measured)} measured families in {len(ch)} clans: "
          f"{len(out):,} families with a height ({len(out) - len(measured):,} inheriting a clan height)")
    print(f"held back as composites (model length > {a.max_length_ratio}x their clan's measured "
          f"members): {len(composites):,}")
    print(f"measured but uncalled and in no clan with a height: {len(uncalled)}")
    was = arch.dropna(subset=["height"]).set_index("Pfam")["height"]
    both = out.set_index("Pfam")["height"].to_frame("now").join(was.rename("was"), how="inner")
    ratio = both.now / both.was
    print(f"against the archived table: {len(both):,} families in both; height within 10%: "
          f"{(abs(ratio - 1) < 0.1).sum():,}, within 25%: {(abs(ratio - 1) < 0.25).sum():,}; "
          f"archived families without a height now: {len(set(was.index) - set(out.Pfam)):,}; "
          f"new families with a height: {len(set(out.Pfam) - set(was.index)):,}")
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, index=False)
    print(f"wrote {a.out}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--oriented", required=True,
                    help="human_domain_sizes_oriented_clans_20200123.csv")
    ap.add_argument("--archived", required=True,
                    help="human_domain_clans_averaged_20200129.csv")
    ap.add_argument("--validate", action="store_true",
                    help="reproduce the archived table and stop")
    ap.add_argument("--pfam-dat", type=Path, help="Pfam-A.hmm.dat.gz for clan membership")
    ap.add_argument("--assignments", help="Pfam 38.2 filtered assignments, to report coverage")
    ap.add_argument("--max-length-ratio", type=float, default=2.0,
                    help="a family may inherit its clan height only if its model length is "
                         "within this factor of the clan's measured members (default 2.0)")
    ap.add_argument("--out")
    ap.add_argument("--rebuild", action="store_true",
                    help="derive the whole table from --oriented (a fresh measurement, e.g. "
                         "build_domain_heights.py's 2026 table) instead of extending --archived; "
                         "--archived then supplies only the clan names")
    a = ap.parse_args()

    o = family_heights(pd.read_csv(a.oriented))
    arch = pd.read_csv(a.archived)
    if a.rebuild:
        return rebuild(o, arch, a)

    exact = np.isclose(o.height_calc, o.height, equal_nan=True).sum()
    print(f"family height from the orientation call: {exact}/{len(o)} exact")

    ch = clan_heights(o)
    t = arch[["Clan_y", "height_clan"]].drop_duplicates().merge(
        ch.rename("calc"), left_on="Clan_y", right_index=True)
    ok = np.isclose(t.calc, t.height_clan).sum()
    print(f"clan height, Sushi excluded and instance-weighted: {ok}/{len(t)} exact")
    if a.validate:
        return 0 if (exact == len(o) and ok == len(t)) else 1

    if not (a.pfam_dat and a.out):
        print("extension needs --pfam-dat and --out", file=sys.stderr)
        return 2

    recs = pfam_clans(a.pfam_dat)
    print(f"Pfam families in the .dat: {len(recs):,}")

    ch_map = ch.to_dict()

    # Median model length of each clan's *measured* members, for the composite guard.
    o["ml"] = o.Pfam.map(lambda f: float((recs.get(f, {}) or {}).get("ml") or "nan"))
    clan_ml = o.dropna(subset=["height_calc"]).groupby("Clan_y").ml.median().to_dict()

    # Every archived row is kept verbatim; this only ever ADDS.  The first version rebuilt
    # the table from clan membership instead, and silently lost the 95 archived families that
    # have no clan at all -- Notch, NOD, NODP, LRRNT, LRRCT, Kunitz_BPTI among them.  They
    # carry a family-level height and no clan, so a clan-driven rebuild cannot see them.
    # NOTCH1 lost its three Notch domains that way, which is how it was caught.
    known = set(arch.Pfam)
    rows = [r._asdict() if hasattr(r, "_asdict") else dict(r)
            for _, r in arch.iterrows()]
    rows = arch.to_dict("records")

    composites = []
    for fam in sorted(set(recs) - known):
        r = recs.get(fam, {})
        clan = r.get("clan")
        if clan not in ch_map:
            continue
        ml = float(r.get("ml") or "nan")
        cm = clan_ml.get(clan, np.nan)
        if ml == ml and cm == cm and cm > 0 and ml / cm > a.max_length_ratio:
            composites.append((fam, r.get("id"), clan, ml, cm, ml / cm))
            continue
        named = arch[arch.Clan_y == clan]
        rows.append({
            "Clan_y": clan,
            "clan_name": named.clan_name.iloc[0] if len(named) else clan,
            "height_clan": ch_map[clan],
            "height_pfam": np.nan,
            "Pfam": fam,
            "query_name_y": r.get("id"),
            "height": ch_map[clan],
        })
    out = pd.DataFrame(rows)

    missing = known - set(out.Pfam)
    if missing:
        print(f"  ERROR: {len(missing)} archived families dropped: {sorted(missing)[:8]}",
              file=sys.stderr)
        return 1
    print(f"  every one of the {len(known):,} archived families preserved")
    print(f"\nclan-height table: {len(arch):,} families -> {len(out):,}")
    print(f"held back as composites (model length > {a.max_length_ratio}x their clan's "
          f"measured members): {len(composites):,}")
    for fam, fid, clan, ml, cm, ratio in sorted(composites, key=lambda x: -x[5])[:6]:
        print(f"    {fam} {str(fid)[:24]:<26} clan {clan}  ML {ml:>5.0f} vs {cm:>5.0f}  x{ratio:.1f}")

    if a.assignments:
        d = pd.read_csv(a.assignments, low_memory=False)
        before = d.Pfam.isin(set(arch.Pfam))
        after = d.Pfam.isin(set(out.Pfam))
        print(f"Pfam 38.2 assignment rows covered: {before.sum():,} ({100*before.mean():.0f}%) "
              f"-> {after.sum():,} ({100*after.mean():.0f}%)")
        still = sorted(set(d.Pfam[~after]))
        print(f"families still without a height: {len(still):,} "
              f"({(~after).sum():,} rows)")
        Path(a.out).with_suffix(".uncovered.txt").write_text(
            "\n".join(f"{f}\t{recs.get(f,{}).get('id','')}\t{recs.get(f,{}).get('clan') or ''}"
                      for f in still) + "\n")

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, index=False)
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
