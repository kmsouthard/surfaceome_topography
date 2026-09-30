#!/usr/bin/env python
"""DAVID Functional Annotation Clustering of the trans-interaction proteins, through the web service.

    python code/database/inputs/build_david_clustering.py --interactions database/interaction_heights.csv \\
        --email-file DIR/david/registered_email.txt --cache DIR/david --out inputs/david/interaction_clustering.csv
        [--names data/curated/david_cluster_names.csv] [--from-cache]

Figure 5 groups the trans pairs by function: a pair belongs to a cluster when both partners are
in it. The clusters come from DAVID's Functional Annotation Clustering (kappa-linked groups of
enriched terms; Huang et al. 2009, Sherman et al. 2022 and 2026), run on every protein of the
trans interaction table. This was the pipeline's last hand step, a browser session; it is now a
call to the DAVID web service (a SOAP interface, spoken here directly so no client library is
needed) with the settings recorded in ``method_citations.csv``:

    categories   GOTERM_BP_DIRECT, GOTERM_CC_DIRECT, GOTERM_MF_DIRECT, REACTOME_PATHWAY, and the
                 UniProt keyword categories (UP_KW_*)
    stringency   Medium: kappa 0.50, term overlap 3, initial and final group size 3, linkage 0.50

The service needs an email address registered with DAVID (davidbioinformatics.nih.gov/webservice/
register.htm), read from ``--email-file`` so no address is written into the repository. The raw
response is kept under ``--cache`` as ``termClusterReport.xml`` with the knowledgebase version
and the date, and ``--from-cache`` rebuilds the table from it without a call, so the committed
table can be regenerated offline and the call repeated only for a new protein set or release.

Output: one row per term per cluster, in the columns Figure 5 reads (Category, Term, Count, %,
PValue, Genes, List Total, Pop Hits, Pop Total, Fold Enrichment, Bonferroni, Benjamini, FDR,
Annotation Cluster, Enrichment Score, Cluster Name). Clusters are numbered by enrichment score;
a cluster's name is its number and its top term unless ``--names`` maps that top term to a
hand-given name, which survives renumbering because it is keyed on the term, not the number.
The service allows 200 jobs a day and lists of at most 3,000 genes.
"""

from __future__ import annotations

import argparse
import datetime as dt
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import requests

ENDPOINT = ("https://davidbioinformatics.nih.gov/webservice/services/"
            "DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/")
NS_SVC = "http://service.session.sample"
NS_SOAP = "http://schemas.xmlsoap.org/soap/envelope/"

CATEGORIES = ("GOTERM_BP_DIRECT,GOTERM_CC_DIRECT,GOTERM_MF_DIRECT,REACTOME_PATHWAY,"
              "UP_KW_BIOLOGICAL_PROCESS,UP_KW_CELLULAR_COMPONENT,UP_KW_CODING_SEQUENCE_DIVERSITY,"
              "UP_KW_DEVELOPMENTAL_STAGE,UP_KW_DISEASE,UP_KW_DOMAIN,UP_KW_LIGAND,"
              "UP_KW_MOLECULAR_FUNCTION,UP_KW_PTM,UP_KW_TECHNICAL_TERM")
#: Medium stringency, DAVID's default: overlap, initial group size, final group size, linkage, kappa (%)
MEDIUM = (3, 3, 3, 0.50, 50)
ID_TYPE = "UNIPROT_ACCESSION"

COLUMNS = ["Category", "Term", "Count", "%", "PValue", "Genes", "List Total", "Pop Hits", "Pop Total",
           "Fold Enrichment", "Bonferroni", "Benjamini", "FDR", "Annotation Cluster",
           "Enrichment Score", "Cluster Name"]


class David:
    """The DAVID web service as plain SOAP 1.1 over requests; one session cookie carries the list."""

    def __init__(self, timeout: int = 900):
        self.s = requests.Session()
        self.timeout = timeout

    def call(self, op: str, *args) -> ET.Element:
        body = "".join(f"<ns:args{i}>{_xml(a)}</ns:args{i}>" for i, a in enumerate(args))
        envelope = (f'<soapenv:Envelope xmlns:soapenv="{NS_SOAP}" xmlns:ns="{NS_SVC}">'
                    f"<soapenv:Body><ns:{op}>{body}</ns:{op}></soapenv:Body></soapenv:Envelope>")
        r = self.s.post(ENDPOINT, data=envelope.encode(), timeout=self.timeout,
                        headers={"Content-Type": "text/xml; charset=utf-8", "SOAPAction": f'"urn:{op}"'})
        root = ET.fromstring(r.content)
        fault = root.find(f".//{{{NS_SOAP}}}Fault")
        if fault is not None or r.status_code != 200:
            text = "".join(fault.itertext()) if fault is not None else r.text[:300]
            raise RuntimeError(f"DAVID {op}: HTTP {r.status_code}: {text.strip()}")
        return root

    def value(self, op: str, *args) -> str:
        ret = self.call(op, *args).find(f".//{{{NS_SVC}}}return")
        return "" if ret is None else (ret.text or "")


def _xml(v) -> str:
    return str(v).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def parse_report(root: ET.Element) -> pd.DataFrame:
    """The term-cluster report as rows; clusters in the order the service returns them (by score)."""
    rows = []
    clusters = [c for c in root.iter() if _local(c.tag) == "return"]
    for n, c in enumerate(clusters, 1):
        fields = {_local(e.tag): (e.text or "") for e in c if _local(e.tag) != "simpleChartRecords"}
        score = float(fields.get("score", "nan"))
        for rec in (e for e in c if _local(e.tag) == "simpleChartRecords"):
            f = {_local(e.tag): (e.text or "") for e in rec}
            rows.append({
                "Category": f.get("categoryName", ""), "Term": f.get("termName", ""),
                "Count": int(f.get("listHits", 0)), "%": float(f.get("percent", "nan")),
                "PValue": float(f.get("ease", "nan")), "Genes": f.get("geneIds", ""),
                "List Total": int(f.get("listTotals", 0)), "Pop Hits": int(f.get("popHits", 0)),
                "Pop Total": int(f.get("popTotals", 0)), "Fold Enrichment": float(f.get("foldEnrichment", "nan")),
                "Bonferroni": float(f.get("bonferroni", "nan")), "Benjamini": float(f.get("benjamini", "nan")),
                "FDR": float(f.get("afdr", "nan")), "Annotation Cluster": f"cluster_{n:02d}",
                "Enrichment Score": score,
            })
    return pd.DataFrame(rows, columns=COLUMNS[:-1])


def name_clusters(t: pd.DataFrame, names: Path | None) -> pd.DataFrame:
    """``Cluster Name``: the cluster's number and its top term, or the hand-given name for that term."""
    given = {}
    if names is not None and names.is_file():
        n = pd.read_csv(names, dtype=str, keep_default_na=False)
        given = dict(zip(n["top_term"], n["name"]))
    out = t.copy()
    out["Cluster Name"] = ""
    for cl, g in t.groupby("Annotation Cluster", sort=False):
        top = g.sort_values("PValue").iloc[0]["Term"]
        top = top.split("~", 1)[-1] if "~" in top else top
        num = cl.split("_")[-1]
        out.loc[g.index, "Cluster Name"] = f"{num} {given.get(top, top)}"
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--interactions", type=Path, required=True,
                    help="the trans interaction table (database/interaction_heights.csv)")
    ap.add_argument("--email-file", type=Path, help="a file holding the DAVID-registered email address")
    ap.add_argument("--cache", type=Path, required=True, help="directory for the raw service response")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--names", type=Path, default=None,
                    help="optional CSV with columns top_term,name giving clusters hand names")
    ap.add_argument("--from-cache", action="store_true", help="rebuild from the cached response, no call")
    ap.add_argument("--list-name", default="trans_interaction_proteins")
    a = ap.parse_args()

    i = pd.read_csv(a.interactions, low_memory=False)
    proteins = sorted(set(i["Human ID link_prot1"].dropna()) | set(i["Human ID link_prot2"].dropna()))
    print(f"{len(proteins)} proteins in the trans interaction table")
    if len(proteins) > 3000:
        print("the DAVID web service takes at most 3,000 genes", file=sys.stderr)
        return 2

    a.cache.mkdir(parents=True, exist_ok=True)
    raw = a.cache / "termClusterReport.xml"
    if a.from_cache:
        if not raw.is_file():
            print(f"no cached response at {raw}", file=sys.stderr)
            return 2
        root = ET.fromstring(raw.read_bytes())
        print(f"rebuilding from {raw}")
    else:
        if not (a.email_file and a.email_file.is_file()):
            print("--email-file is required for a call: a file holding an address registered at "
                  "https://davidbioinformatics.nih.gov/webservice/register.htm", file=sys.stderr)
            return 2
        email = a.email_file.read_text().strip()
        d = David()
        if d.value("authenticate", email) != "true":
            print(f"DAVID does not know this address; register it first", file=sys.stderr)
            return 1
        mapped = float(d.value("addList", ",".join(proteins), ID_TYPE, a.list_name, 0))
        print(f"  list uploaded: {mapped:.1%} of ids mapped")
        cats = d.value("setCategories", CATEGORIES)
        print(f"  categories accepted: {cats}")
        root = d.call("getTermClusterReport", *MEDIUM)
        raw.write_bytes(ET.tostring(root))
        (a.cache / "RUN.txt").write_text(
            f"date {dt.date.today().isoformat()}\nproteins {len(proteins)}\nid_type {ID_TYPE}\n"
            f"categories {cats}\nstringency medium overlap={MEDIUM[0]} initial={MEDIUM[1]} "
            f"final={MEDIUM[2]} linkage={MEDIUM[3]} kappa={MEDIUM[4]}\nmapped {mapped}\n")
        print(f"  response cached at {raw}")

    t = name_clusters(parse_report(root), a.names)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    t.to_csv(a.out, index=False)
    print(f"{len(t)} term rows in {t['Annotation Cluster'].nunique()} clusters "
          f"({(t.groupby('Annotation Cluster')['Enrichment Score'].first() >= 2).sum()} with score >= 2) -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
