#!/usr/bin/env python
"""Classify every surfaceome protein by location, topology, function and surface evidence.

    python code/database/inputs/build_surfaceome_classification.py --surfaceome CSV \\
        --proteome UniProt.tab --surfy table_S3_surfaceome.xlsx --out CSV

The surfaceome is defined inclusively: every reviewed human protein UniProt gives an
extracellular topological domain or a GPI anchor. That admits proteins whose annotated
location is an intracellular membrane (ER, Golgi, lysosome, endosome, vesicle), which reach the
plasma membrane only in some states, if at all. Rather than filter them, the height table
carries a classification, so a reader (or a figure) can choose a stricter surfaceome and say so:

    uniprot_location   plasma membrane | GPI-anchored | secreted | membrane, unspecified |
                       intracellular membrane only, from UniProt's subcellular location
    tm_count           transmembrane segments in UniProt (0 = GPI-anchored or peripheral)
    surfy_label        surface | nonsurface | unscored: SURFY's random-forest call
                       (Bausch-Fluck et al. 2018, PNAS, PMID 30373828)
    surfy_score        its score, 0 to 1
    surfy_fpr_class    1, 2 or 3: the false-positive-rate class (1%, 5%, 15%) the surface
                       call was made at; "(nonTM)" marks proteins scored without a TM segment
    almen_main         Receptors | Transporters | Enzymes | Miscellaneous | Unclassified,
                       the Almen membranome classification carried by the SURFY table
                       (Almen et al. 2009, BMC Biology, PMID 19678920: function and
                       evolutionary origin of the human membrane proteome)
    almen_sub          its subclass (GPCR;Rhodopsin, SLC;..., Channels;..., ...)
    cspa_category      Cell Surface Protein Atlas evidence of the protein at a cell surface
                       (Bausch-Fluck et al. 2015, PLoS One, PMID 25894527): high confidence,
                       putative, or none
    cd_number          CD nomenclature, where one exists

``--surfy`` is the SURFY master table, ``table_S3_surfaceome.xlsx`` from
https://wollscheidlab.org/SURFY (20,193 human proteins; its terms are not stated, so only this
derived table for the surfaceome is committed).
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

PLASMA = re.compile(r"cell membrane|cell surface|plasma membrane|apical|basolateral|cell junction|"
                    r"synap|cilium|flagell|sarcolemma|microvill|cell projection|membrane raft|"
                    r"podosome|invadopodium|lateral cell membrane", re.I)
INTRA = re.compile(r"endoplasmic reticulum|golgi|lysosom|endosom|vesicle|mitochondri|nucle|"
                   r"peroxisom|melanosom|vacuole|sarcoplasmic", re.I)


def location_class(text: str) -> str:
    s = str(text) if isinstance(text, str) else ""
    if "GPI-anchor" in s:
        return "GPI-anchored"
    if PLASMA.search(s):
        return "plasma membrane"
    if re.search(r"\bSecreted\b", s):
        return "secreted"
    if INTRA.search(s):
        return "intracellular membrane only"
    if "membrane" in s.lower():
        return "membrane, unspecified"
    return "no location"


def read_surfy(path: Path) -> pd.DataFrame:
    raw = pd.read_excel(path, sheet_name="SurfaceomeMasterTable", header=None)
    hdr = next(i for i in range(5) if raw.iloc[i].astype(str).str.contains("UniProt", case=False).any())
    t = raw.iloc[hdr + 1:].copy()
    t.columns = [str(c).strip() for c in raw.iloc[hdr]]
    t["accession"] = t["UniProt accession"].astype(str).str.split("-").str[0]
    return t.drop_duplicates("accession").set_index("accession")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--surfaceome", type=Path, required=True,
                    help="a surfaceome table with an 'ID link' accession column (stage 0's <proteome>.csv)")
    ap.add_argument("--proteome", type=Path, required=True, help="the UniProt proteome .tab (gz accepted)")
    ap.add_argument("--surfy", type=Path, required=True, help="table_S3_surfaceome.xlsx")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()

    accs = sorted(set(pd.read_csv(a.surfaceome, usecols=["ID link"])["ID link"].astype(str)))
    prot = pd.read_csv(a.proteome, sep="\t", low_memory=False).set_index("Entry")
    surfy = read_surfy(a.surfy)

    out = pd.DataFrame(index=pd.Index(accs, name="accession"))
    loc = prot["Subcellular location [CC]"].reindex(accs)
    out["uniprot_location"] = loc.map(location_class)
    out["tm_count"] = prot["Transmembrane"].reindex(accs).fillna("").str.count("TRANSMEM").astype(int)
    label = surfy["Surfaceome Label"].reindex(accs)
    out["surfy_label"] = label.where(label.notna(), "unscored")
    out["surfy_score"] = pd.to_numeric(surfy["MachineLearning score"].reindex(accs), errors="coerce").round(4)
    out["surfy_fpr_class"] = surfy["MachineLearning FPR class (1=1%, 2=5%, 3=15%)"].reindex(accs).astype(object)
    out["almen_main"] = surfy["Membranome Almen main-class"].reindex(accs)
    out["almen_sub"] = surfy["Membranome Almen sub-class"].reindex(accs)
    out["cspa_category"] = surfy["CSPA category"].reindex(accs)
    out["cd_number"] = surfy["CD number"].reindex(accs)
    out = out.reset_index()

    a.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, index=False)
    print(f"{len(out)} surfaceome proteins classified -> {a.out}")
    print("  location:", out["uniprot_location"].value_counts().to_dict())
    print("  SURFY:", out["surfy_label"].value_counts().to_dict())
    print("  Almen:", out["almen_main"].fillna("not in table").value_counts().to_dict())
    print("  CSPA evidence:", int(out["cspa_category"].notna().sum()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
