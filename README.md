# The cell surface topography

A pipeline to estimate the heights of human cell-surface proteins

**Explore the results online at [kmsouthard.github.io/surfaceome_topography](https://kmsouthard.github.io/surfaceome_topography/)**

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/img/monocyte_surface_dark.svg">
    <img src="docs/img/monocyte_surface_light.svg" width="100%"
         alt="Twenty surface proteins of a classical monocyte drawn to scale on a membrane, from AlphaFold structures, with grey capsules for regions without a structure">
  </picture>
  <br>
  <sub>Twenty proteins sampled by abundance from a classical monocyte's surface (Ravenhill et al.
  2020), drawn to scale with <a href="https://github.com/jordisr/cellscape">CellScape</a> by Jordi
  Silvestre-Ryan and oriented by their topology. Grey capsules are height without a structure:
  disordered regions, domains, or a whole protein with no AlphaFold model, like LRP1, which is
  shortened at the break.</sub>
</p>

> **Status.** Analysis code for a manuscript in preparation. Result tables and figures to follow.

## Explorer

The results can be browsed at
[kmsouthard.github.io/surfaceome_topography](https://kmsouthard.github.io/surfaceome_topography/).

| tab | what it shows |
|---|---|
| Protein | The height of any surface protein, with a drawing of it to scale |
| Epitope | The gap an antibody sets between two cells, from the residues it binds |
| Cell surface | The proteins on a cell type drawn to scale. You can upload your own expression data |
| Contact | The gaps between two cell types, with or without an antibody. You can upload a LIANA result |

## Quick start

```bash
conda env create -f environment.yml
conda activate surfaceome
python code/database/fetch_restricted_inputs.py
python code/run_notebooks.py out
```

`fetch_restricted_inputs.py` downloads three inputs that are not included in the repository
(about 45 MB).

The run writes tables to `out/database/` and
`out/tables/`, figures to `out/figures/`, and the executed notebooks to `out/notebooks/`.

## The method

**1. Heights.** A protein's height is the length of its extracellular region, the
*ectodomain*, with its folded and disordered parts laid end to end along the membrane normal.
It is therefore an upper bound: a tilted or bent ectodomain stands lower. Each residue is
assigned to one source, in this order of preference:

| source | height it contributes |
|---|---|
| AlphaFold v6 model, trimmed to its confident span | longest dimension of the model's inertia-axis bounding box |
| disordered segment (MobiDB, and low-confidence AlphaFold stretches) | end-to-end distance of a worm-like chain whose persistence length is set by its glycan density |
| Pfam domain not covered by a model | the family's height, measured on solved structures as how far one domain advances its chain |
| any residue left | 0.04 nm |

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/img/cd45_dark.svg">
    <img src="docs/img/cd45_light.svg" width="520"
         alt="CD45 on the membrane: a 15.5 nm folded region from its AlphaFold model, and a 25.6 nm glycosylated disordered region above it, 41.1 nm in all">
  </picture>
</p>

The result is one height per gene for every reviewed human protein with an extracellular
domain or a GPI anchor (`out/database/height_estimates.csv`).

**2. Cell surfaces.** A cell type's surface is the surfaceome proteins detected in its
proteome, each weighted by its share of total abundance.

**3. Contacts.** Pairs that bind across two cells come from manual annotations of the STRING database and CellphoneDB. A pair's
height is the sum of its partners' heights, unless a recorded rule gives another (a measured
model of the complex, or partners that bind side by side). An antibody bridge sets the gap from
the height of its epitope.

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/img/contact_dark.svg">
    <img src="docs/img/contact_light.svg" width="520"
         alt="Protein pairs spanning the gap between an NK cell membrane above and a HER2-positive breast cancer cell membrane below, one pair per band of gap heights">
  </picture>
  <br>
  <sub>An NK cell against a HER2-positive breast cancer cell treated with trastuzumab (anti-HER2).
  In each band of gap heights, the pair holding most of it among pairs of two single proteins whose
  height is the sum of the partners' (NK-cell partner blue, cancer-cell partner orange), labelled
  with its height; drawn heights are illustrative.</sub>
</p>

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/img/bridge_dark.svg">
    <img src="docs/img/bridge_light.svg" width="560"
         alt="An antibody bridge: HER2 on the cancer cell, trastuzumab lying across its membrane-proximal epitope, and the Fc receptor CD16 hanging from the NK cell onto the antibody's Fc">
  </picture>
  <br>
  <sub>The same contact's antibody bridge. Trastuzumab binds HER2 1.6 nm above the membrane, so
  the gap is that plus the Fc receptor (CD16, FCGR3A) and the antibody's 3.24 nm. The antibody is
  an intact human IgG1 (PDB 1HZH).</sub>
</p>

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/img/bullseye_dark.svg">
    <img src="docs/img/bullseye_light.svg" width="620"
         alt="The same contact in section and from above: the membranes step apart from about 8 nm at the centre to about 106 nm at the rim, and a bullseye of nested rings for the gap bands, sized by their share of the contact">
  </picture>
  <br>
  <sub>The whole contact (Figure 6). Below, the bullseye: one ring per band of gap heights, each
  ring's area its share of the contact. Above, a section through its centre: over each ring, the
  membranes stand at the mean gap of the pairs holding it. Radius is an estimated share of the contact.</sub>
</p>

[`docs/METHODS.md`](docs/METHODS.md) gives each rule, its parameters and its basis.

## The code

`code/run_notebooks.py` runs the method, then draws the figures from its tables
(`--stage database` or `--stage figures` runs one stage).

The method is three notebooks in `code/database/notebooks/`:

| notebook | builds |
|---|---|
| `01_height_estimates` | the height of every surfaceome protein |
| `02_cellphonedb_interactions` | CellphoneDB pairs and complexes |
| `03_interaction_heights` | the height of each bound pair |

The figures are drawn by one notebook each in `code/figures/notebooks/`:

| notebook | figure |
|---|---|
| `F1_domain_sizes` | 1: domain heights by Pfam family |
| `F2_height_histograms` | 2: heights across the surfaceome |
| `F3_F4_expression_weighted_topography` | 3–4: cell surfaces weighted by abundance |
| `F4_ravenhill_monocyte_subtypes` | 4: monocyte subtypes |
| `F5_david_clusters` | 5: functional clusters of bound pairs |
| `F6_cell_contacts` | 6: four immune–cancer cell contacts |
| `S2_methods_vs_structures` | S2: predicted heights against solved structures |

[`code/figures/figure_map.csv`](code/figures/figure_map.csv) maps each output file to its
manuscript panel. The illustrations in this README are drawn from a run's output with
[CellScape](https://github.com/jordisr/cellscape) by `code/figures/cellscape/readme_images.py`,
in its own environment (`code/figures/cellscape/environment.yml`).

To rebuild the explorer's data from a run:

```bash
python code/explorer/build_explorer.py out --models DIR
conda run -n cellscape python code/explorer/build_cartoons.py out --alphafold DIR
python -m http.server 8765 --directory explorer    # view it at http://localhost:8765
```

## Using the heights for your own cells

```python
import pandas as pd
from surfaceomeTopography import cell_surface, mean_height, plot_surfaces

heights = pd.read_csv("out/database/height_estimates.csv")
expression = pd.read_csv("my_expression.csv")  # a "uniprot" column, then one column per cell type
table = expression.merge(heights, left_on="uniprot", right_on="ID link_first")

surface = cell_surface(table, "my cell type")  # proteins it expresses, each with its share
mean_height(surface)                           # abundance-weighted mean height, nm
plot_surfaces({"my cell type": surface}).savefig("my_cell_type.pdf")
```

## An antibody's epitope

```bash
python code/predict_epitope.py --antibody trastuzumab
python code/predict_epitope.py ERBB2 267-333 --models DIR
```

Gives the epitope's height, the gap the antibody sets to each Fc receptor, and whether CD45,
CD45RO and CD148 fit in that gap. Use a curated antibody, or any residues on any protein. `DIR`
is the folder of AlphaFold models.

## Repository layout

```
code/
  run_notebooks.py          runs the pipeline
  predict_epitope.py        an epitope's height, and the gap an antibody on it holds
  explorer/                 builds the explorer's data files from a run
  surfaceome_config.py      locates inputs under data/ and outputs under the run directory
  surfaceomeTopography/     the Python library the notebooks use
  database/notebooks/       the method
  database/inputs/          rebuild the inputs from their public sources (build_inputs.py)
  figures/notebooks/        the manuscript figures
  figures/validation/       solved structures and complex models, for checking the heights
  figures/cellscape/        the README illustrations
data/
  curated/                  decision tables and their citations
  inputs/                   external inputs, by source
  measurements/             domain and structure measurements made for this work
explorer/                   the web explorer: one page, its script and style, and data/
docs/METHODS.md             the method in full
docs/img/                   the README illustrations
```

[`data/README.md`](data/README.md) describes each input file.

## Data sources

| used for | source | version |
|---|---|---|
| proteome | [UniProt](https://www.uniprot.org) reference proteome UP000005640 | 2026_02 |
| structural models | [AlphaFold DB](https://alphafold.ebi.ac.uk) | v6 |
| domains | [Pfam](https://www.ebi.ac.uk/interpro/entry/pfam/), searched with HMMER | 38.2 |
| disorder | [MobiDB](https://mobidb.org) | 2026 |
| glycosylation, experimental | [GlyGen](https://www.glygen.org) | 2.11.1 |
| glycosylation, predicted | GlycoMine (Li et al. 2015, *Bioinformatics* 31:1411) | 2014 |
| surface classification | [SURFY](https://wollscheidlab.org/SURFY) (Bausch-Fluck et al. 2018, *PNAS* 115:E10988) | 2018 |
| protein abundance | [Expression Atlas](https://www.ebi.ac.uk/gxa/home) E-PROT-1 and E-PROT-27 | 2026 download |
| monocyte subtypes | Ravenhill et al. 2020, *Sci Rep* 10:4560, Table S1C | 2020 |
| interactions | [STRING](https://string-db.org) | v12.0 |
| interactions and complexes | [CellphoneDB](https://www.cellphonedb.org) | v5.0.0 |
| functional clusters | [DAVID](https://davidbioinformatics.nih.gov) | knowledgebase v2025_1 |
| checking the heights | [SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/) and the [PDB](https://www.rcsb.org) | 2026-09-08 |

`python code/database/build_inputs.py --list` shows how each input is rebuilt from its source.
That needs HMMER, PyMOL and about 5 GB of AlphaFold models.

## Limitations

- **Heights are upper bounds.** They assume each ectodomain is fully extended, so the gaps
  between cells in contact are upper bounds too.
- **Disordered, heavily glycosylated regions**, such as mucin stalks, rely on a polymer propagation model.
- **Solved structures are used only to check heights**, and cover a minority of ectodomains.
- **The surfaceome is inclusive**: it includes some proteins annotated only to intracellular
  membranes. Each protein has location and surface-evidence classification available for filtering if desired.
- **Contacts assume every pair binds equally well**, since binding strengths are unknown for
  most pairs.
- **Antibody epitopes are known for 12 of the 14 antibodies**, 10 from a structure and 2 from
  mapping studies. The other two are placed at the membrane.

## Environment

Python 3.9 and pandas 1.5, pinned in [`environment.yml`](environment.yml): later pandas versions
change grouping defaults in ways that alter results without an error.

## Licence

Code: MIT ([`LICENSE`](LICENSE)). Data made for this work (`data/curated/`,
`data/measurements/`): CC BY 4.0 ([`data/LICENSE`](data/LICENSE)).

## References

- Proteome, topology and experimental glycosites: UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2023. *Nucleic Acids Res* 2023;51:D523–D531. [doi:10.1093/nar/gkac1052](https://doi.org/10.1093/nar/gkac1052)
- Structure prediction: Jumper J, Evans R, Pritzel A, et al. Highly accurate protein structure prediction with AlphaFold. *Nature* 2021;596:583–589. [doi:10.1038/s41586-021-03819-2](https://doi.org/10.1038/s41586-021-03819-2)
- Structural models: Varadi M, Bertoni D, Magana P, et al. AlphaFold Protein Structure Database in 2024: providing structure coverage for over 214 million protein sequences. *Nucleic Acids Res* 2024;52:D368–D375. [doi:10.1093/nar/gkad1011](https://doi.org/10.1093/nar/gkad1011)
- Domain families: Mistry J, Chuguransky S, Williams L, et al. Pfam: The protein families database in 2021. *Nucleic Acids Res* 2021;49:D412–D419. [doi:10.1093/nar/gkaa913](https://doi.org/10.1093/nar/gkaa913)
- The domain search: Eddy SR. Accelerated Profile HMM Searches. *PLoS Comput Biol* 2011;7:e1002195. [doi:10.1371/journal.pcbi.1002195](https://doi.org/10.1371/journal.pcbi.1002195)
- Solved structures, for the domain heights: Berman HM, Westbrook J, Feng Z, et al. The Protein Data Bank. *Nucleic Acids Res* 2000;28:235–242. [doi:10.1093/nar/28.1.235](https://doi.org/10.1093/nar/28.1.235)
- Residue mapping between UniProt and the PDB: Dana JM, Gutmanas A, Tyagi N, et al. SIFTS: updated Structure Integration with Function, Taxonomy and Sequences resource allows 40-fold increase in coverage of structure-based annotations for proteins. *Nucleic Acids Res* 2019;47:D482–D489. [doi:10.1093/nar/gky1114](https://doi.org/10.1093/nar/gky1114)
- Disordered regions: Piovesan D, Del Conte A, Clementel D, et al. MobiDB: 10 years of intrinsically disordered proteins. *Nucleic Acids Res* 2023;51:D438–D444. [doi:10.1093/nar/gkac1065](https://doi.org/10.1093/nar/gkac1065)
- Experimental glycosites: York WS, Mazumder R, Ranzinger R, et al. GlyGen: Computational and Informatics Resources for Glycoscience. *Glycobiology* 2020;30:72–73. [doi:10.1093/glycob/cwz080](https://doi.org/10.1093/glycob/cwz080)
- Predicted glycosites: Li F, Li C, Wang M, et al. GlycoMine: a machine learning-based approach for predicting N-, C- and O-linked glycosylation in the human proteome. *Bioinformatics* 2015;31:1411–1419. [doi:10.1093/bioinformatics/btu852](https://doi.org/10.1093/bioinformatics/btu852)

## Acknowledgements

This work was done in the lab of Daniel A. Fletcher in the Department of Bioengineering at the
University of California, Berkeley ([fletchlab.berkeley.edu](https://fletchlab.berkeley.edu/)).

The protein illustrations, here and in the manuscript, are drawn with
[CellScape](https://github.com/jordisr/cellscape), Jordi Silvestre-Ryan's tool for turning
protein structures into vector cartoons and composing them into cell-surface scenes:

> Silvestre-Ryan J, Fletcher DA, Holmes I. CellScape: Protein structure visualization with
> vector graphics cartoons. *bioRxiv* 2022.
> [doi:10.1101/2022.06.14.495869](https://doi.org/10.1101/2022.06.14.495869)

## AI assistance

During the preparation of this work the author used Claude (Anthropic) to assist in writing code
and documentation. After using this tool, the author reviewed and edited the content as needed
and takes full responsibility for the content of this work.
