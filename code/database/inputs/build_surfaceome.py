#!/usr/bin/env python
"""Derive the surfaceome (extracellular domains + GPI-anchored proteins) from a UniProt pull.

    python code/database/inputs/build_surfaceome.py --input-dir DIR --out-dir DIR [--proteome UP000005640]

This is the first stage of the pipeline, ported from `code/raw_code/Homo sapiens ECDs.ipynb`,
which does not run as written: it imports `Bio.Alphabet.IUPAC` (removed in Biopython 1.78),
calls `subset_ecds(gff)` against a signature that now takes two arguments, and hardcodes
absolute paths from a retired machine.

Inputs are the three UniProt downloads for a proteome — `.gff`, `.tab`, `.fasta` — named
`proteomes<ID>.{gff,tab,fasta}` in `--input-dir`. Fetch them with
`surfaceomeTopography.download_proteome`.

Outputs, in `--out-dir`:

    <ID>.csv                 the surfaceome: one row per ECD/GPI/short/long segment
    <ID>_no_long_split.csv   same, with long ECDs left unsplit
    proteomes<ID>_ecds.fasta      ECD sequences within the Phyre2 length window
    proteomes<ID>_gpi.fasta       GPI-anchored chain sequences, same window
    proteomes<ID>_all_ecds.fasta  every ECD and GPI sequence over 30 aa
    proteomes<ID>_long.fasta      long ECDs, split for modelling (if any)
    proteomes<ID>_ecds_<a>_<b>.fasta   the ECD set in 500-sequence batches

Behaviour is preserved exactly, including two quirks worth knowing about: `seq_len` is computed
as `end - start` rather than `end - start + 1`, and `fasta_subset` slices sequences with
0-based Python indices against 1-based GFF coordinates. Both are left as they were so that
output stays comparable with the published dataset; changing them would silently shift every
sequence-scaled height.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # code/, for surfaceome_config and the package

from surfaceomeTopography.database_parsing import (
    import_gff, import_tab, import_fasta, fasta_subset,
)
from surfaceomeTopography.identify_domains import (
    subset_ecds, subset_gpi, split_long, phyre_filter, surfaceome_df,
)

BATCH = 500


def build(input_dir: Path, out_dir: Path, proteome: str) -> dict:
    stem = input_dir / f"proteomes{proteome}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_stem = out_dir / f"proteomes{proteome}"

    print(f"reading {stem}.gff / .tab / .fasta")
    def _in(ext: str) -> str:
        """`<stem>.<ext>`, or its `.gz` sibling -- the vendored proteome is gzipped."""
        plain = Path(f"{stem}.{ext}")
        return str(plain if plain.is_file() else Path(f"{stem}.{ext}.gz"))

    gff = import_gff(_in("gff"))
    tab = import_tab(_in("tab"))
    fasta = import_fasta(_in("fasta"))
    print(f"  gff rows {len(gff):,} | tab rows {len(tab):,} | sequences {len(fasta):,}")

    gpi = subset_gpi(tab, gff)
    ecds = subset_ecds(gff, tab)
    print(f"  GPI segments {len(gpi):,} over {gpi['ID link'].nunique():,} proteins")
    print(f"  ECD segments {len(ecds):,} over {ecds['ID link'].nunique():,} proteins")

    # long ECDs exceed the Phyre2 submission limit and are split for modelling
    long = split_long(ecds) if (ecds["seq_len"] >= 3000).any() else pd.DataFrame()
    if not long.empty:
        SeqIO.write(fasta_subset(long.copy(), fasta), f"{out_stem}_long.fasta", "fasta")
        print(f"  long ECDs split: {len(long):,} segments")

    phyre_ecd = phyre_filter(ecds)
    phyre_gpi = phyre_filter(gpi)
    short = ecds[ecds["seq_len"] < 30]
    print(f"  within Phyre2 window: {len(phyre_ecd):,} ECD, {len(phyre_gpi):,} GPI; "
          f"{len(short):,} short (<30 aa)")

    gpi_seq = fasta_subset(phyre_gpi.copy(), fasta)
    ecd_seq = fasta_subset(phyre_ecd.copy(), fasta)
    SeqIO.write(gpi_seq, f"{out_stem}_gpi.fasta", "fasta")
    SeqIO.write(ecd_seq, f"{out_stem}_ecds.fasta", "fasta")

    for start in range(0, len(ecd_seq), BATCH):
        SeqIO.write(ecd_seq[start:start + BATCH],
                    f"{out_stem}_ecds_{start}_{start + BATCH}.fasta", "fasta")

    total = surfaceome_df(phyre_gpi, phyre_ecd, short, long=long)
    total.to_csv(out_dir / f"{proteome}.csv", index=False)

    # the same table with long ECDs left unsplit
    big = ecds[ecds["seq_len"] >= 3000].assign(g=0)
    surfaceome_df(phyre_gpi, phyre_ecd, short, long=big).to_csv(
        out_dir / f"{proteome}_no_long_split.csv", index=False)

    all_ecds = pd.concat([gpi, ecds])
    not_short = all_ecds[all_ecds["seq_len"] > 30]
    SeqIO.write(fasta_subset(not_short.copy(), fasta), f"{out_stem}_all_ecds.fasta", "fasta")

    summary = {
        "proteins_in_tab": len(tab),
        "ecd_segments": len(ecds),
        "ecd_proteins": int(ecds["ID link"].nunique()),
        "gpi_segments": len(gpi),
        "gpi_proteins": int(gpi["ID link"].nunique()),
        "surfaceome_rows": len(total),
        "surfaceome_proteins": int(total["ID link"].nunique()),
    }
    print("\n  " + "  ".join(f"{k}={v:,}" for k, v in summary.items()))
    return summary


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input-dir", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--proteome", default="UP000005640")
    a = ap.parse_args()
    build(a.input_dir, a.out_dir, a.proteome)
    return 0


if __name__ == "__main__":
    sys.exit(main())
