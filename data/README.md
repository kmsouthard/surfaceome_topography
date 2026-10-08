# Data

Everything a run reads is committed here.
Three inputs whose licences do not allow redistribution are the exception: they are
not committed, and `python code/database/fetch_restricted_inputs.py` downloads their sources and
rebuilds them, checking each against the hash recorded in `inputs/RESTRICTED.json`. They are
marked *fetched* below.

Paths in the notebooks are relative to this directory and go through
`surfaceome_config.resolve()`; run `python code/surfaceome_config.py` to see what resolves.

| directory | what it holds | who changes it |
|---|---|---|
| `curated/` | the decision tables (the interaction-type calls per STRING pair and the review list of pairs awaiting a call, interaction rules, cis pairs, the complexes added to CellphoneDB's, antibody bridges, each antibody's epitope and the structure it is read from, the proteins placed in an antibody's gap, height rules, corrections, the domain orientation calls, the assembled-model templates) and the citation sheet behind every decision | hand edits, then `python code/database/build_citation_sheet.py` |
| `inputs/` | every external input a run reads, grouped by source, at the 2026 vintage | `python code/database/build_inputs.py` (stage 0) rebuilds the derived ones from the public databases |
| `measurements/` | tables measured here on structures: the 2020 and 2026 domain measurements, the PDB structure sizes, the assembled interaction models | the scripts under `code/database/inputs/` and `code/figures/validation/` |

## Inputs

| path | source | vintage |
|---|---|---|
| `inputs/proteome/UP000005640.{tab,gff,fasta}.gz` | UniProt reference proteome, `surfaceomeTopography.download_proteome` | 2026_02 |
| `inputs/proteome/surfaceome_current.tsv` | accession to AlphaFold model (stage 0) | 2026 |
| `inputs/mapping/` | accession to gene and protein names (stage 0, UniProt) plus two accessions the export missed | 2026 |
| `inputs/disorder/` | MobiDB disorder regions (stage 0) | 2026 |
| `inputs/glycosylation/glycomine/` | GlycoMine 2014 predictions, reduced to surfaceome accessions (`MANIFEST.json` records the reduction; *fetched*, reduced to the accessions in `accessions.txt`) | 2014, frozen |
| `inputs/glycosylation/glygen_sites_human.csv` | GlyGen experimental glycosylation sites | 2.11.1 |
| `inputs/domains/domain_assignments.csv.gz` | Pfam domains per ectodomain with their heights (stage 0, HMMER over Pfam 38.2) | 2026 |
| `inputs/domains/domain_clan_heights.csv` | the per-family height table every domain assignment takes its height from (stage 0, from `measurements/domains_2026/`) | 2026 |
| `inputs/domains/domain_disorder_annotations_20210419.csv` | the 2021 domain and disorder annotation the height notebook reconciles against | 2021 |
| `inputs/alphafold/` | per-residue AlphaFold confidence, the measured ectodomain dimensions, and where each residue of an antibody's antigen sits along its model's box axes, for the epitope heights (stage 0, AlphaFold v6) | 2026 |
| `inputs/structures/` | SIFTS chain to UniProt mapping, and the 2020 list of surfaceome proteins with a PDB structure | 2026, 2020 |
| `inputs/string/string_surfaceome_pairs.csv.gz` | every STRING v12.0 link between surfaceome proteins at combined score 400 or more, in both orders, labelled cis/trans/secreted/pathway from `curated/interaction_types.csv` (the 2020 calls, made on v11) and `curated/string_v12_pair_review.csv` (stage 0, `build_string_interactions.py`; `download_string.py` fetches the release).  | v12.0 |
| `inputs/cellphonedb/` | CellphoneDB v5.0.0, flattened to the four tables the notebooks read (stage 0); *fetched* | v5.0.0 |
| `inputs/expression/` | Expression Atlas E-PROT-1 and E-PROT-27, re-downloaded 2026 with current gene symbols; Ravenhill 2020 monocyte abundances | 2026, 2020 |
| `inputs/david/interaction_clustering.csv` | DAVID functional annotation clustering of the trans interaction proteins (stage 0, `build_david_clustering.py` through the DAVID web service; needs a registered address in `downloads/david/registered_email.txt`; cluster names from `curated/david_cluster_names.csv` where given). The committed table is the 2026-09-16 run on knowledgebase v2025_1 | 2026 |
| `inputs/phyre2/` | the Phyre2-era ectodomain sizes with disorder (the homology path) | 2020 |
| `inputs/annotation/surfaceome_classification.csv` | per surfaceome protein: UniProt location class and TM count, the SURFY surface call, score and FPR class (Bausch-Fluck 2018), Almen functional class, CSPA evidence, CD number (stage 0, from the SURFY master table at wollscheidlab.org/SURFY); *fetched*, over the accessions in `surfaceome_accessions.csv` | 2026 |

Text too large to commit raw is stored `.gz`; `resolve()` returns the `.gz` path and pandas
decompresses by extension.

## Measurements

`measurements/domains_2026/` is the re-measurement of every Pfam domain on current PDB
structures: the chains selected from SIFTS (`domain_instances_selected.csv`), every domain
instance boxed on them (`domain_measurements.csv`, 38,859 instances of 1,075 families), and
the per-family heights (`domain_heights.csv`) with their basis: the neighbour pitch, the median
advance of a domain along the axis through its sequence neighbours, or the box dimension along
a called axis for families never seen beside a neighbour. `domains_2020.csv` is the 2020
measurement it replaced, kept for comparison; no step reads it. `pdb_structure_sizes.txt` measures the solved
surfaceome structures for the Supplementary Figure 2 comparison. `models_2018/` and
`models_2026/` are the assembled models of bound pairs whose lengths the interaction height
rules take; their templates are recorded in `curated/model_templates.csv`.

## Outputs

A run writes to its output directory: `database/` for the tables that ship (the surfaceome,
the height estimates, the domain and disorder assignments, the interaction heights), `tables/`
for intermediates and per-panel tables, `figures/` for the figure files, `notebooks/` for the
executed notebooks. `surfaceome_config.GENERATED` lists the tables a later notebook reads
back; `resolve()` prefers this run's copy of those.

