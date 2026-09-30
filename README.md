# Surfaceome topography: the heights of human cell-surface proteins

How far does each protein on a human cell stand above the membrane? This repository estimates
that height for every cell-surface protein in the human proteome, 3,146 proteins, and uses the
heights to describe the surface topography of immune and cancer cells and the gaps that bound
protein pairs hold between two cells in contact.

> **Status.** Analysis code for a manuscript in preparation. The result tables and figures are
> not posted here until the manuscript is published; the pipeline regenerates all of them from
> its inputs in about two minutes (below). Figure and table numbers refer to the
> manuscript draft and may change.

## Quick start

```bash
conda env create -f environment.yml
conda activate surfaceome
python code/database/fetch_restricted_inputs.py
python code/run_notebooks.py out
```

Three inputs cannot be redistributed under their licences (GlycoMine, SURFY and CellphoneDB;
see [Data sources](#data-sources)). `fetch_restricted_inputs.py` downloads them, about 45 MB,
rebuilds the files the pipeline reads, and checks each against the hash of the file the
analysis was run on. Everything else is committed.

The run takes about two minutes on a laptop. It writes the height and interaction tables to
`out/database/`, intermediate and per-panel tables to `out/tables/`, the figures to
`out/figures/` and the executed notebooks to `out/notebooks/`. The run fails if any notebook
raises.

## What the height is

A protein's height is how far its extracellular region, its *ectodomain*, can extend from the
membrane: its folded and disordered parts laid end to end along the membrane normal. This is an
upper bound by construction; a tilted or folded-back ectodomain stands lower.

Each ectodomain is covered by annotations, resolved so that no residue is counted twice, and
each piece contributes a height. In order of preference:

| piece | its height |
|---|---|
| AlphaFold v6 model, trimmed to its confident span | longest dimension of its inertia-axis bounding box |
| disordered segment (MobiDB, plus low-confidence AlphaFold runs) | root-mean-square end-to-end distance of a worm-like chain, stiffened by its glycan density |
| Pfam domain, where no model covers the region | its family's height, measured on solved structures as the distance a domain advances along its chain |
| residues nothing else covers | 0.04 nm per residue |

[`docs/METHODS.md`](docs/METHODS.md) gives every rule, its parameters and what it rests on.

## How the pipeline runs

`run_notebooks.py` runs ten notebooks in two stages, each reading what the ones before it wrote:
`--stage database` builds the tables, `--stage figures` draws the figures from them, and the
default runs both.

```
UniProt proteome, AlphaFold models, Pfam domains, MobiDB disorder, glycosylation
        │
        ▼
code/database/notebooks/                     the database stage
   01_height_estimates                       the height table, one row per gene
   02_cellphonedb_interactions               CellphoneDB pairs and complexes
   03_interaction_heights                    the height of each bound pair
        │
        ▼
code/figures/notebooks/                      the figure stage
   F1_domain_sizes                           Figure 1: domain heights by Pfam family
   F2_height_histograms                      Figure 2: heights across the surfaceome
   F3_F4_expression_weighted_topography      Figures 3–4: expression-weighted cell surfaces
   F4_ravenhill_monocyte_subtypes            Figure 4: monocyte subtypes
   F5_david_clusters                         Figure 5: functional clusters of bound pairs
   F6_cell_contacts                          Figure 6: four immune–cancer cell contacts
   S2_methods_vs_structures                  Figure S2: predicted heights against solved structures
```

[`code/figures/figure_map.csv`](code/figures/figure_map.csv) lists which output file is which
manuscript panel.

**1. The height table.** From the UniProt human reference proteome, take every reviewed protein
with an extracellular topological domain or a GPI anchor, keep one entry per gene (its largest
ectodomain), and estimate its height as above. The result, `out/database/height_estimates.csv`,
has one row per protein with its height in nanometres and the sources that went into it
(`methods`); every later step reads it.

**2. Cell surfaces.** A cell type's surface is the surfaceome proteins detected in its proteome,
each weighted by its share of the summed abundance. Figures 3 and 4 compare these weighted
height distributions across immune cell types and monocyte subtypes.

**3. Interactions and contacts.** Pairs of surface proteins that bind across two cells come from
STRING and CellphoneDB. STRING says two proteins associate, not how, so each STRING pair carries
a call, recorded in `data/curated/`, of whether it binds across cells (*trans*), on one cell
(*cis*), through a secreted partner, or only in a pathway. A trans pair's height is the sum of
its partners' heights unless a recorded rule says otherwise (an assembled model of the complex,
or partners that bind side by side). For an immune cell against a cancer cell, the pairs both
cells can form give the distribution of gaps at which the contact holds the two membranes, and
the proteins too tall for those gaps.

## Using the heights for your own cells

After a run, the height table and the library give the surface topography of any cell type
with protein abundance data:

```python
import pandas as pd
from surfaceomeTopography import cell_surface, mean_height, plot_surfaces

heights = pd.read_csv("out/database/height_estimates.csv")
expression = pd.read_csv("my_expression.csv")  # a "uniprot" column, then one column per cell type
table = expression.merge(heights, left_on="uniprot", right_on="ID link_first")

surface = cell_surface(table, "my cell type")  # the proteins it expresses, each with its share
mean_height(surface)                           # expression-weighted mean height, nm
plot_surfaces({"my cell type": surface}).savefig("my_cell_type.pdf")
```

## Every judgement is a table, and every table is cited

Where the pipeline makes a call rather than a measurement (which STRING pairs bind in trans,
which pairs bind side by side, how a domain family orients, what an antibody bridge adds), the
call is a row in a decision table under [`data/curated/`](data/curated/), with the papers it
rests on, and code applies it; nothing is edited by hand downstream. The height model's own
rules are in `method_citations.csv`.
[`data/curated/CITATIONS.md`](data/curated/CITATIONS.md) summarises every table's citations,
each PubMed-verified, and marks the rules that are this work's own assumptions.

## Repository layout

```
code/
  run_notebooks.py        runs the pipeline, by stage or whole
  surfaceome_config.py    resolves every input under data/ and every output under the run directory
  surfaceomeTopography/   the Python library the notebooks use
  database/
    notebooks/            the database stage (01–03)
    inputs/               rebuild the committed inputs from their public sources
    build_inputs.py       runs those rebuild steps in dependency order
    build_citation_sheet.py   gathers the decision tables' citations into data/curated/
  figures/
    notebooks/            the figure stage, one notebook per figure
    figure_map.csv        output file to manuscript panel
    validation/           solved structures and assembled complexes, for checking the heights
  check_generated_paths.py    each table has one writer, and nothing is read before it is written
  check_snapshot.py, build_snapshot.py, file_hashes.py
                          check which proteome vintage a run reads; record a run's tables with
                          the hashes of the code and data that produced them
data/
  curated/                the decision tables and their citations
  inputs/                 every external input a run reads, by source
  measurements/           domain and structure measurements made for this work
docs/METHODS.md           the method as implemented
environment.yml           the pinned software environment
```

[`data/README.md`](data/README.md) describes every input file and its vintage.

## Data sources

Every input a run reads is under `data/inputs/`, at the version below: committed where its
licence allows, and otherwise fetched by `code/database/fetch_restricted_inputs.py`.

| used for | source | version | terms |
|---|---|---|---|
| proteome | [UniProt](https://www.uniprot.org) reference proteome UP000005640 | 2026_02 | CC BY 4.0 |
| structural models | [AlphaFold DB](https://alphafold.ebi.ac.uk) | v6 | CC BY 4.0 |
| domains | [Pfam](https://www.ebi.ac.uk/interpro/entry/pfam/), searched with HMMER | 38.2 | CC0 |
| disorder | [MobiDB](https://mobidb.org) | 2026 | CC BY 4.0 |
| glycosylation, experimental | [GlyGen](https://www.glygen.org) | 2.11.1 | CC BY 4.0 |
| glycosylation, predicted | GlycoMine (Li et al. 2015, *Bioinformatics* 31:1411); server offline, release held by the Internet Archive | 2014 | none stated: fetched |
| surface classification | [SURFY](https://wollscheidlab.org/SURFY) (Bausch-Fluck et al. 2018, *PNAS* 115:E10988), carrying the Cell Surface Protein Atlas and Almén et al. 2009 classes | 2018 | CC BY-NC-ND 4.0: fetched |
| protein abundance | [Expression Atlas](https://www.ebi.ac.uk/gxa/home) E-PROT-1 (Kim et al. 2014) and E-PROT-27 (Tyanova et al. 2016) | 2026 download | CC BY 4.0 |
| monocyte subtype abundance | Ravenhill et al. 2020, *Sci Rep* 10:4560, Supplementary Table S1C | 2020 | CC BY 4.0 |
| interactions | [STRING](https://string-db.org) | v12.0 | CC BY 4.0 |
| interactions and complexes | [CellphoneDB](https://www.cellphonedb.org) data repository | v5.0.0 | none stated: fetched |
| functional clusters | [DAVID](https://davidbioinformatics.nih.gov), over GO, Reactome and UniProt keywords | knowledgebase v2025_1 | results may be published |
| checking the heights | [SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/) and the [PDB](https://www.rcsb.org) | 2026-09-08 | CC BY 4.0; CC0 |

"Fetched" inputs are not in the repository;
`data/inputs/RESTRICTED.json` records each one's source, why it is not committed, and the hash
the rebuilt file must match.

**Rebuilding the inputs.** `python code/database/build_inputs.py --list` shows the steps that
rebuild the committed inputs from these sources, and the header of
[`build_inputs.py`](code/database/build_inputs.py) lists the raw downloads each one needs. The
rebuild also needs HMMER, PyMOL and about 5 GB of AlphaFold models; running the analysis needs
none of them.

## Limitations

- **Heights are the aligned limit.** Folded and disordered parts are laid end to end along the
  membrane normal, so a height is how far a protein *can* extend, not how far it typically does.
  Contact gaps and exclusion use the same heights, so they are upper bounds too.
- **Disordered, heavily glycosylated regions**, such as mucin stalks, take their height from a
  polymer model that few experimental measurements test.
- **Solved structures check the heights and never set them.** They cover a minority of
  ectodomains.
- **The surfaceome is inclusive.** Every reviewed protein with an extracellular segment is in,
  including some whose annotated location is an intracellular membrane. Each entry carries a
  location, topology, function and surface-evidence classification to filter on.
- **Contacts assume every trans pair binds equally well.** Binding strengths are unknown for
  most pairs and change with crowding at the membrane.
- **Antibody bridges assume the antibody binds next to the membrane**, as trastuzumab does on
  HER2; any other antibody needs the height of its epitope.

## Environment

The environment is pinned to Python 3.9 and pandas 1.5 ([`environment.yml`](environment.yml)):
several notebooks use APIs pandas 2 removed, and pandas 2 and 3 change group-by and reduction
defaults in ways that alter results without raising.

## Licence

Code: MIT ([`LICENSE`](LICENSE)). The decision tables and measurements made for this work
(`data/curated/`, `data/measurements/`): CC BY 4.0 ([`data/LICENSE`](data/LICENSE)). Third-party
inputs under `data/inputs/` remain under their providers' terms (the table above); cite the
sources when you use them.

Structural illustrations in the manuscript use [CellScape](https://github.com/jordisr/cellscape).
