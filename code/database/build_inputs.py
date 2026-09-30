#!/usr/bin/env python
"""Stage 0 of the database pipeline: rebuild the derived inputs from the public databases.

    python code/database/build_inputs.py --list
    python code/database/build_inputs.py --downloads DIR --out DIR [--only STEP ...] [--dry-run]

The three database notebooks read derived inputs that are committed under ``data/inputs/``.
This driver runs the scripts in ``code/database/inputs/`` that make them, in dependency
order, from a directory of raw downloads, and compares each result with the committed copy. Nothing here is needed to *run* the pipeline: the committed inputs are what a
run reads. It is needed to *refresh* the database against a newer UniProt, AlphaFold, Pfam or
MobiDB release.

``--downloads DIR`` must hold:

    proteomes/proteomes<ID>.{fasta,gff,tab}   surfaceomeTopography.download_proteome
    alphafold/                                 download_alphafold.py (about 5 GB for the surfaceome)
    pfam/Pfam-A.hmm, pfam/Pfam-A.hmm.dat       the Pfam release, pressed (hmmpress) for hmmfetch
    sifts/pdb_chain_pfam.tsv.gz, sifts/pdb_chain_taxonomy.tsv.gz   SIFTS flat files
    pdb/                                       mmCIF files of the selected domain instances
    domain_orientations.csv                    the orientation calls (a copy of data/curated/)
    cellphonedb/cellphonedb.zip                the cellphonedb-data release bundle (v5.0.0)
    string/<taxon>.protein.{links.detailed,aliases,info}.v<v>.txt.gz   download_string.py (v12.0)
    surfy/table_S3_surfaceome.xlsx             the SURFY master table (wollscheidlab.org/SURFY)
    david/registered_email.txt                 an address registered with the DAVID web service

The MobiDB regions and the UniProt names are fetched by their steps and cached under
``--downloads``. The AlphaFold measurements need PyMOL importable in this environment.

A step whose inputs are missing is skipped with the reason, and the run continues with the
steps that can proceed. Two inputs are never rebuilt: the GlycoMine predictions (the 2014
server is offline) are committed as primary data, and GlyGen sites come from
``surfaceomeTopography.download_glygen`` directly. The STRING step joins the release's links
to the interaction-type calls in ``data/curated/``, so a newer STRING is a re-join, not a
re-annotation; pairs new to the release are listed for a call, never labelled by default. The domain heights are measured afresh on
current structures (four steps); a family no current structure covers keeps its 2020 value.

Every step's output lands under ``--out`` at the data-relative path the notebooks read, so
``--out`` can be a scratch directory for comparison, or ``data`` to replace the committed
inputs, after which a full run and ``build_snapshot.py`` re-record the manifest.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

REPO = Path(__file__).resolve().parents[2]
CODE = REPO / "code"
INPUTS = CODE / "database" / "inputs"
DATA = REPO / "data"
BASELINE = DATA / "baseline_2022"

sys.path.insert(0, str(CODE))
from file_hashes import content_sha256  # noqa: E402

PROTEOME = "UP000005640"


def _interaction_table() -> Path:
    """The trans interaction table DAVID clusters: the snapshot's copy, else a run's in ``out/``.

    The public release ships no snapshot, so there the table comes from running the database
    stage first (``python code/run_notebooks.py out --stage database``).
    """
    snap = DATA / "snapshot" / "database" / "interaction_heights.csv"
    return snap if snap.is_file() else REPO / "out" / "database" / "interaction_heights.csv"

#: The surfaceome the input scripts take (`--surfaceome`, an "ID link" column): an
#: intermediate of this stage, not a vendored input.  The height notebook derives the same
#: surfaceome itself; this copy exists so the AlphaFold and name steps can run before it.
SURFACEOME = f"stage0/{PROTEOME}.csv"


@dataclass
class Step:
    name: str
    script: str
    makes: list[str]                                   # outputs, relative to --out
    needs: list[str] = field(default_factory=list)     # relative to --downloads, or to --out
    tool: str = ""                                     # external requirement, for the listing
    note: str = ""
    args: Callable[[Path, Path], list[str]] = field(default=lambda dl, out: [], repr=False)


def _out(out: Path, rel: str) -> Path:
    p = out / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def _committed(rel: str) -> Path | None:
    """The committed copy of a data-relative path, if there is one."""
    for root in (DATA, BASELINE):
        for cand in (root / rel, root / (rel + ".gz")):
            if cand.is_file():
                return cand
    return None


def steps() -> list[Step]:
    """The stage in dependency order."""
    AF_CONF = "inputs/alphafold/af_confidence.csv"
    CLANS = "inputs/domains/domain_clan_heights.csv"
    return [
        Step("cellphonedb", "build_cellphonedb_tables.py",
             makes=[f"inputs/cellphonedb/{t}_curated.csv"
                    for t in ("interaction", "protein", "complex", "hla")],
             needs=["cellphonedb/cellphonedb.zip"],
             note="the CellphoneDB bundle flattened to the four tables the interaction "
                  "notebooks read",
             args=lambda dl, out: ["--bundle", str(dl / "cellphonedb" / "cellphonedb.zip"),
                                   "--version", "v5.0.0",
                                   "--out-dir", str(_out(out, "inputs/cellphonedb/x").parent)]),
        Step("string", "build_string_interactions.py",
             makes=["inputs/string/string_surfaceome_pairs.csv"],
             needs=["string/9606.protein.links.detailed.v12.0.txt.gz",
                    "string/9606.protein.aliases.v12.0.txt.gz", "string/9606.protein.info.v12.0.txt.gz",
                    SURFACEOME],
             tool="download_string.py",
             note="every STRING v12 link between surfaceome proteins, labelled cis/trans/secreted/"
                  "pathway from data/curated/interaction_types.csv and string_v12_pair_review.csv",
             args=lambda dl, out: ["--string-dir", str(dl / "string"),
                                   "--surfaceome", str(out / SURFACEOME),
                                   "--types", str(DATA / "curated" / "interaction_types.csv"),
                                   "--review", str(DATA / "curated" / "string_v12_pair_review.csv"),
                                   "--out", str(_out(out, "inputs/string/string_surfaceome_pairs.csv"))]),
        Step("david", "build_david_clustering.py",
             makes=["inputs/david/interaction_clustering.csv"],
             needs=["david/registered_email.txt"], tool="network (DAVID web service)",
             note="DAVID Functional Annotation Clustering of every protein in the shipped trans "
                  "interaction table (a figure-stage input built from a database-stage output); "
                  "the raw response is cached under downloads/david/",
             args=lambda dl, out: ["--interactions", str(_interaction_table()),
                                   "--email-file", str(dl / "david" / "registered_email.txt"),
                                   "--cache", str(dl / "david"),
                                   "--names", str(DATA / "curated" / "david_cluster_names.csv"),
                                   "--out", str(_out(out, "inputs/david/interaction_clustering.csv"))]),
        Step("surfaceome", "build_surfaceome.py", makes=[SURFACEOME],
             needs=[f"proteomes/proteomes{PROTEOME}.{ext}" for ext in ("fasta", "gff", "tab")],
             note="the surfaceome and its extracellular regions, from the UniProt pull, plus the "
                  "ectodomain FASTA that hmmsearch runs over",
             args=lambda dl, out: ["--input-dir", str(dl / "proteomes"),
                                   "--out-dir", str(_out(out, SURFACEOME).parent),
                                   "--proteome", PROTEOME]),
        Step("gene_names", "build_gene_name_mapping.py",
             makes=["inputs/mapping/uniprot_surfaceome_gene_name.tab"],
             needs=[SURFACEOME], tool="network (UniProt)",
             note="accession to Entry name, Protein names and Gene names",
             args=lambda dl, out: ["--surfaceome", str(out / SURFACEOME),
                                   "--cache", str(dl / "proteomes" / "uniprot_names_cache.tsv"),
                                   "--out", str(_out(out, "inputs/mapping/uniprot_surfaceome_gene_name.tab"))]),
        Step("disorder", "build_mobidb_annotations.py",
             makes=[f"inputs/disorder/{PROTEOME}_disorder_full.txt",
                    f"inputs/disorder/{PROTEOME}_disorder_lite.txt"],
             needs=[SURFACEOME], tool="network (MobiDB)",
             note="MobiDB-lite disorder regions, in the legacy format the parser reads",
             args=lambda dl, out: ["--surfaceome", str(out / SURFACEOME),
                                   "--out-dir", str(_out(out, f"inputs/disorder/{PROTEOME}_disorder_full.txt").parent),
                                   "--cache", str(dl / "mobidb"), "--proteome", PROTEOME]),
        Step("domain_instances", "select_domain_instances.py",
             makes=["stage0/domain_instances_selected.csv"],
             needs=["sifts/pdb_chain_pfam.tsv.gz", "sifts/pdb_chain_taxonomy.tsv.gz", "pfam/ecds.domtbl"],
             note="the PDB chains to measure each Pfam family on: human, then mouse, up to 20 "
                  "entries per family (SIFTS)",
             args=lambda dl, out: ["--sifts-pfam", str(dl / "sifts" / "pdb_chain_pfam.tsv.gz"),
                                   "--sifts-taxonomy", str(dl / "sifts" / "pdb_chain_taxonomy.tsv.gz"),
                                   "--assignments", str(dl / "pfam" / "ecds.domtbl"),
                                   "--out", str(_out(out, "stage0/domain_instances_selected.csv"))]),
        Step("domain_measure", "measure_domain_instances.py",
             makes=["stage0/domain_measurements.csv"],
             needs=["stage0/domain_instances_selected.csv", "pdb", "pfam/Pfam-A.hmm"], tool="PyMOL, HMMER",
             note="every Pfam domain on every selected chain boxed on its structure (one hmmsearch of "
                  "Pfam-A over the chains gives the residue spans); download the entries with "
                  "figures/validation/download_pdb_structures.py first",
             args=lambda dl, out: ["--all-families",
                                   "--instances", str(out / "stage0/domain_instances_selected.csv"),
                                   "--structures", str(dl / "pdb"), "--hmm", str(dl / "pfam" / "Pfam-A.hmm"),
                                   "--work", str(_out(out, "stage0/domain_work/x").parent),
                                   "--out", str(_out(out, "stage0/domain_measurements.csv"))]),
        Step("domain_heights", "build_domain_heights.py",
             makes=["stage0/domain_heights.csv"],
             needs=["stage0/domain_measurements.csv", "domain_orientations.csv"],
             note="per-family heights: the neighbour pitch (how far a domain advances its chain, "
                  "from its neighbours' centroids) where a family is seen beside a neighbour, else the "
                  "box dimension along the called axis (hand table, stacking, clan majority)",
             args=lambda dl, out: ["--height-from-neighbour-pitch",
                                   "--measurements", str(out / "stage0/domain_measurements.csv"),
                                   "--orientations", str(dl / "domain_orientations.csv"),
                                   "--hmm-dat", str(dl / "pfam" / "Pfam-A.hmm.dat"),
                                   "--out", str(_out(out, "stage0/domain_heights.csv"))]),
        Step("clan_heights", "build_clan_heights.py", makes=[CLANS],
             needs=["stage0/domain_heights.csv", "pfam/Pfam-A.hmm.dat"],
             note="the clan-height table rebuilt from the measured families: clan averages, and "
                  "every Pfam family inheriting its clan's height within the composite guard; "
                  "families with no structure keep their 2020 measurement",
             args=lambda dl, out: [
                 "--rebuild", "--oriented", str(out / "stage0/domain_heights.csv"),
                 "--archived", str(BASELINE / "inputs/domains/domain_clan_heights.csv"),
                 "--pfam-dat", str(dl / "pfam" / "Pfam-A.hmm.dat"),
                 "--out", str(_out(out, CLANS))]),
        Step("ecd_search", "search_ectodomain_domains.py", makes=["stage0/ecds.domtbl"],
             needs=[f"stage0/proteomes{PROTEOME}_all_ecds.fasta", "pfam/Pfam-A.hmm"], tool="HMMER",
             note="hmmsearch of the whole Pfam library over the ectodomain sequences the "
                  "surfaceome step wrote",
             args=lambda dl, out: ["--fasta", str(out / f"stage0/proteomes{PROTEOME}_all_ecds.fasta"),
                                   "--hmm", str(dl / "pfam" / "Pfam-A.hmm"),
                                   "--out", str(_out(out, "stage0/ecds.domtbl"))]),
        Step("domains", "build_domain_assignments.py",
             makes=["inputs/domains/domain_assignments.csv"],
             needs=["stage0/ecds.domtbl", CLANS],
             note="Pfam domain assignments for every ectodomain, from the hmmsearch table",
             args=lambda dl, out: ["--domtbl", str(out / "stage0/ecds.domtbl"),
                                   "--clans", str(out / CLANS),
                                   "--out", str(_out(out, "inputs/domains/domain_assignments.csv"))]),
        Step("classification", "build_surfaceome_classification.py",
             makes=["inputs/annotation/surfaceome_classification.csv"],
             needs=[SURFACEOME, f"proteomes/proteomes{PROTEOME}.tab", "surfy/table_S3_surfaceome.xlsx"],
             note="location, topology, SURFY call, Almen class and CSPA evidence per surfaceome protein",
             args=lambda dl, out: ["--surfaceome", str(out / SURFACEOME),
                                   "--proteome", str(dl / "proteomes" / f"proteomes{PROTEOME}.tab"),
                                   "--surfy", str(dl / "surfy" / "table_S3_surfaceome.xlsx"),
                                   "--out", str(_out(out, "inputs/annotation/surfaceome_classification.csv"))]),
        Step("af_confidence", "build_af_confidence.py", makes=[AF_CONF],
             needs=["alphafold", SURFACEOME],
             note="per-residue AlphaFold confidence per ectodomain, and where to trim it",
             args=lambda dl, out: ["--models", str(dl / "alphafold"),
                                   "--surfaceome", str(out / SURFACEOME),
                                   "--out", str(_out(out, AF_CONF))]),
        Step("af_mapping", "build_af_mapping.py",
             makes=["inputs/proteome/surfaceome_current.tsv"],
             needs=[SURFACEOME, AF_CONF],
             note="which AlphaFold model each surfaceome accession uses",
             args=lambda dl, out: ["--surfaceome", str(out / SURFACEOME),
                                   "--af-confidence", str(out / AF_CONF),
                                   "--out", str(_out(out, "inputs/proteome/surfaceome_current.tsv"))]),
        Step("af_dimensions", "build_af_dimensions.py",
             makes=["inputs/alphafold/alphafold_dimensions.txt"],
             needs=["alphafold", AF_CONF], tool="PyMOL",
             note="the AlphaFold height measurements: the inertia-axis bounding box over the "
                  "confident span of each ectodomain (measure_structures.py)",
             args=lambda dl, out: ["--models", str(dl / "alphafold"),
                                   "--af-confidence", str(out / AF_CONF),
                                   "--out", str(_out(out, "inputs/alphafold/alphafold_dimensions.txt"))]),
    ]


def _present(dl: Path, out: Path, need: str) -> bool:
    return (dl / need).exists() or (out / need).exists()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true", help="print the steps and stop")
    ap.add_argument("--downloads", type=Path, help="directory of raw downloads (see above)")
    ap.add_argument("--out", type=Path, help="where rebuilt inputs go, by data-relative path")
    ap.add_argument("--only", nargs="+", metavar="STEP", help="run only these steps, by name")
    ap.add_argument("--dry-run", action="store_true", help="print the commands without running")
    a = ap.parse_args()

    S = steps()
    if a.list or not (a.downloads and a.out):
        print(f"{'step':<14} {'script':<30} {'needs':<28} makes")
        for s in S:
            print(f"{s.name:<14} {s.script:<30} {s.tool:<28} {s.makes[0]}")
            for m in s.makes[1:]:
                print(f"{'':<74} {m}")
        if not a.list:
            print("\n--downloads DIR and --out DIR run the stage; --list only prints it.")
        return 0

    dl, out = a.downloads.resolve(), a.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    todo = [s for s in S if not a.only or s.name in a.only]
    ran = skipped = failed = 0
    for s in todo:
        missing = [n for n in s.needs if not _present(dl, out, n)]
        if missing:
            print(f"SKIP  {s.name:<14} missing: {', '.join(missing)}")
            skipped += 1
            continue
        cmd = [sys.executable, str(INPUTS / s.script), *s.args(dl, out)]
        print(f"RUN   {s.name:<14} {s.script} {' '.join(cmd[2:])}")
        if a.dry_run:
            continue
        if subprocess.run(cmd).returncode:
            print(f"FAIL  {s.name}")
            failed += 1
            continue
        ran += 1
        for rel in s.makes:
            made, was = out / rel, _committed(rel)
            if not made.is_file():
                print(f"      {rel}: not written")
            elif was is None:
                print(f"      {rel}: written (an intermediate; nothing committed to compare)")
            elif content_sha256(made) == content_sha256(was):
                print(f"      {rel}: identical to the committed copy")
            else:
                print(f"      {rel}: DIFFERS from {was.relative_to(REPO)}")
    print(f"\n{ran} ran, {skipped} skipped, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
