# Methods

How the height of every human cell-surface protein is estimated, how the heights are used, and
what each number rests on. This is the methods-level account of the pipeline as committed; the
manuscript's Methods section is written from it. Every rule below is applied by code and, where
it is a judgement rather than a measurement, recorded in a decision table under `data/curated/`
with the papers it rests on (`data/curated/CITATIONS.md` lists them all).

## 1. The surfaceome

**Definition.** Every reviewed human protein (UniProt reference proteome UP000005640, release
2026_02) with an extracellular topological domain or a GPI anchor. Each extracellular segment is
an *ectodomain*; a protein with several transmembrane passes has several. The height table keeps
one entry per gene, the largest ectodomain, so a multi-pass receptor is represented by its
largest loop. This gives 3,146 proteins.

**Classification, not filtering.** The rule is inclusive: it admits 118 proteins whose UniProt
location names only an intracellular membrane (ER, Golgi, lysosome, endosome, vesicle) and 335
that the SURFY classifier calls nonsurface. Rather than narrow the definition, every protein
carries a classification a reader can filter on: the UniProt location class (plasma membrane,
GPI-anchored, secreted, membrane unspecified, intracellular only) and transmembrane count; the
SURFY surface call with its score and false-positive class (Bausch-Fluck et al. 2018); the Almen
functional class, receptor, transporter, enzyme or miscellaneous, with its subclass (Almén et
al. 2009); Cell Surface Protein Atlas evidence of the protein observed at a surface
(Bausch-Fluck et al. 2015); and the CD number. A stricter set, plasma membrane or GPI or secreted
by UniProt, or SURFY surface, or CSPA evidence, keeps 2,920 proteins; it also drops CD99 and platelet
GP9, which is why no filter is applied by default.

## 2. The height of an ectodomain

The height is the distance an ectodomain can extend from the membrane: the length of its
structured and unstructured parts laid end to end along the membrane normal. It is an upper
bound by construction, the aligned limit; a tilted or folded-back ectodomain stands lower.

Each ectodomain is covered by annotations, resolved so that no residue is counted twice, and
each annotated piece contributes a height. The pieces, in the order they are preferred:

**AlphaFold models.** For every ectodomain with an AlphaFold v6 model, the model is cut to its
confident span, trimming low-confidence residues at either end, and the inertia-axis bounding box of that span is measured in PyMOL; the longest box
dimension is the height of the folded part. Low-confidence runs of 20 residues or more inside the
span, and 5 or more at its ends (pLDDT below 70, the threshold Piovesan et al. 2022 found
F1-optimal for disorder), are treated as disordered rather than folded. 2,520 ectodomains take
their folded height this way.

**Disordered segments.** Regions MobiDB annotates as disordered, merged with the low-confidence
AlphaFold runs, contribute the root-mean-square end-to-end distance of a worm-like chain of their
length, with a persistence length set by glycosylation. Glycosites come from UniProt and GlyGen
experimental annotation unioned with GlycoMine predictions, and the grafting density (sites per
residue) sets the persistence length between the two literature limits, 1 nm for a bare chain and
8 nm for a dense glycan brush (Kuo et al. 2018), by a crowding argument with no fitted parameter.
Consecutive disordered segments with no structure between them are pooled into one chain, and a
chain whose density varies along its length is treated exactly with a piecewise persistence
length. The model reproduces the electron-microscopy height of CD45RABC (39.8 predicted, about
40 measured) and overestimates the short isoform CD45RO (26 against 22); it is validated on those
two proteins and little else, which section 6 returns to.

**Pfam domains, by counting.** Where no AlphaFold model covers a region, its Pfam domains
(HMMER over Pfam 38.2) each contribute the height of their family, and the heights are summed.
A family's height is measured here on structures: every Pfam domain on 9,078 PDB chains
selected through SIFTS (38,859 instances of 1,075 families, first deposited model only,
solution-scattering entries excluded), and for each family the *neighbour pitch*, the median
distance a domain advances its chain along the axis through its sequence neighbours of any
family. This is the height a domain contributes to a string of domains, tilt and overlap
included; summed along a chain it gives the chain's end-to-end length. 378 families take a
pitch. A family never seen beside a neighbour keeps the longest or shortest dimension of its
own bounding box according to an orientation call (the 2020 hand table, then the direction
same-family copies stack, then the majority of its Pfam clan); families with fewer than three
whole instances inherit their clan's instance-weighted mean; a Pfam model that spans several
domains of other families carries no height, so the domains it covers are counted instead.

**Sequence.** Residues no annotation reaches contribute 0.04 nm each, a compact-chain rate.

**Homology models.** The 2020 Phyre2 models are kept as a parallel path (`height_estimates_homology.csv`)
for comparison with the submitted analysis; the shipping height table is the AlphaFold path.

An ectodomain is *well modelled* when its annotations cover at least 70% of its length; the rest
are completed with the sequence rate over the uncovered residues. Which sources went into each
height is recorded in the table's `methods` column.

## 3. Checks against measured structures

No height comes from a solved structure; the structures check the predictions. 393 solved
surfaceome structures (SIFTS 2026-09-08, chains measured with the same bounding box) are
compared with the predicted heights of the same ectodomains in Supplementary Figure 2, and 22
assembled models of bound pairs are compared with the summed heights of section 4. Domain
counting alone, the simplest estimate, is reported as its own table (Supplementary Table 5) and
compared with the model-based heights.

## 4. Interactions and their heights

**Pairs.** Candidate pairs between surface proteins come from STRING v12.0 (every link between
two surfaceome proteins at combined score 400 or more) and CellphoneDB v5.0.0 (interactions and
complexes). STRING says two proteins associate, not how: whether a pair binds across two cells
(*trans*), on one cell (*cis*), through a secreted partner, or only in a pathway is a call made
here per pair and kept in `interaction_types.csv`, the 2020 calls made on STRING v11 and carried
over by pair, so a STRING refresh re-joins the pairs to the calls rather than re-annotating them.
A pair new in v12 has no call and is not trans until it is called; `string_v12_pair_review.csv`
lists the ones worth calling, with their evidence. The trans set is 906 pairs with a height. 32 pairs carry a hand correction of their type recorded in
`interaction_type_corrections.csv`; 494 pairs are classed as cis, binding on the same cell, in
`cis_pairs.csv`, with the pairs left out of that class and why.

**Complexes.** CellphoneDB's complexes are joined by those curated here (`curated_complexes.csv`):
the MHC class II heterodimers HLA-DR (DRA with DRB1, DRB3 or DRB5), HLA-DQ (DQA1 with DQB1, DQA2
with DQB2) and HLA-DP (DPA1 with DPB1). DRA with DRB4 is left out until its pairing is cited. A
trans pair recorded with one chain of a curated complex holds for its other chain, because a
proteome often detects only one chain and the complex then never forms as a unit. With LAG3
added as a class II ligand (`curated_trans_interactions.csv`), this gives CD4 and LAG3 a pair
with every class II chain but DRB4: 16 of the 906 pairs.

**The height of a pair.** By default the two partners' heights are summed, the extended limit.
Three kinds of exception are recorded per pair in `interaction_height_rules.csv`: 22 pairs take
the measured length of an assembled model, the full ectodomains of both partners placed on a
solved structure of the complex (`model_templates.csv` names each template and its paper; the
2026 rebuilds superpose human AlphaFold v6 ectodomains on human complex structures); 26 pairs
whose partners are known to bind extended keep the sum; 6 bind side by side and take the taller
partner (LAG3 with class II among them, from the structure of mouse LAG3 on mouse I-Ab). A complex stands at its tallest subunit. Antibody bridges add the antibody's 3.24 nm
between an Fc receptor and its target, and at a contact hold the gap at the epitope:

    gap = max(epitope height + Fc receptor + 3.24 nm, antigen height)

**Epitope heights.** An antibody's epitope is the antigen residues with a heavy atom within 4 Å of
the antibody in a structure of the two bound, in UniProt numbering (`antibody_epitopes.csv` names
the PDB entry and its paper for 10 of the 14 antibodies; `build_antibody_epitopes.py` derives the
residues). Two more have no structure and a coarser epitope, marked by `basis`: catumaxomab's is
the three EpCAM peptides its binding arm recognises in peptide mapping (residues 49-88 and
175-196), and elotuzumab's is the membrane-proximal IgC2 domain of SLAMF7 it is reported to bind
(131-206), the whole stretch standing in for the epitope. A residue's height is the ectodomain's stack of section 1 read part-way up: the
segments between it and the membrane, plus its C-alpha's position along the box axis its model's
height is measured on, or its share of a disordered chain's or a counted domain's residues. The
membrane end is the ectodomain's C-terminal end when it is the chain's first ectodomain and its
N-terminal end when it is the last; a loop between two helices is placed only inside a model,
with the membrane at the end its two ends lie nearer. An epitope stands at the mean of its
residues (`antibody_epitope_heights.csv`). The two antibodies with no epitope recorded,
tafasitamab and olaratumab, are taken at the membrane. `code/predict_epitope.py` applies the same rule to any residues, and
places CD45, CD45RO and CD148 (`segregation_probes.csv`) in the gap by the exclusion rule of
section 5; it reports an epitope as within reach of phagocytosis up to 10 nm above the membrane
(Bakalar et al. 2018).

**Functional clustering.** The proteins of the trans pairs were clustered by DAVID functional
annotation clustering (knowledgebase v2025_1, run 2026-09-16; Gene Ontology direct terms,
Reactome and UniProt keywords, medium stringency) for Figure 5, called through the DAVID web
service by a stage-0 step, so the clustering is rerun by code whenever the trans set changes;
a cluster's name is its top term unless a hand name is recorded for that term.

## 5. Cell surfaces and contacts

**Expression weighting.** A cell type's surface is the set of surfaceome proteins with a value
in a proteome of that cell type, each weighted by its share of the summed abundance. Immune
cells use Expression Atlas E-PROT-1 (Kim et al. 2014 draft proteome, re-downloaded with current
gene symbols), the monocyte subtypes Ravenhill et al. 2020, breast tumours E-PROT-27 (Tyanova et al. 2016). The
weighted mean height and the height distribution of a cell type are what Figures 3 and 4 show.

**Contacts.** For an immune cell against a cancer cell, every trans pair whose partners both
cells express can form across the contact. Each pair binds equally well, so a protein's
abundance is shared among its partners by mass action, and the contact's gap distribution is the
pair heights weighted by those shares. A protein about 5 nm taller than a gap is excluded from
it (Schmid et al. 2016), which gives each contact its excluded set. Four contacts are drawn in
Figure 6, named in `contact_panels.csv`: monocytes against HER2-positive breast cancer cells with
and without an anti-HER2 antibody, NK cells against HER2-positive cells with the antibody, and NK
cells against triple-negative cells without one.

## 6. Limitations

- Summed heights are the aligned limit. Contact
  gaps and exclusion inherit this.
- The disorder model is validated on two measurements of one protein. Its interval, from the
  literature range of the persistence lengths, is a median 38% of the value, and the interval
  is narrowest for short segments, where the model was wrong for CD45RO.
- Domain counting uses one height per family. β-propellers are counted blade by blade, which
  can overstate a flat propeller; 252 families never seen beside a neighbour still take a box
  dimension along a called axis rather than a measured advance.
- The surfaceome is inclusive; the classification columns say what each entry is, and the
  expression panels contain a few percent by weight of proteins whose surface location is not
  established.
- Expression data cover 970 of 3,146 proteins, so the weighted panels rest on 100 to 170
  proteins each, and every interaction is assumed to bind equally well.
- The breast cancer surfaces are tumour proteomes (Tyanova et al. 2016), not purified cancer
  cells.
- An epitope's height assumes its antigen's model stands on the axis its height is measured on.
  Where the model's membrane end lies in the middle half of that axis the orientation is flagged
  as uncertain (EGFR, for cetuximab), as it is for 376 single-anchored ectodomains in all. Two
  epitopes are no finer than a set of peptides or a domain, and two antibodies have none recorded
  and are taken to bind at the membrane.

## 7. Reproducibility

Every input a run reads is committed under `data/` at the vintage it was used; the pipeline is
10 notebooks run in order by `code/run_notebooks.py`, and a run from a clone reads nothing
outside the repository. `code/check_generated_paths.py` checks that every table has one writer
and that no notebook reads a table before the notebook that writes it has run. The derived inputs rebuild from the public databases
with `code/database/build_inputs.py` (stage 0), byte for byte for the offline steps. The environment is pinned to Python 3.9 and pandas 1.5, because newer
pandas changes group-by defaults in ways that alter results without raising.

## Sources

UniProt 2026_02 (UniProt Consortium 2023). AlphaFold DB v6. Pfam 38.2, HMMER 3. MobiDB 2026.
GlyGen 2.11.1; GlycoMine 2014. SIFTS 2026-09-08 and the PDB. STRING v12.0. CellphoneDB
v5.0.0. Expression Atlas E-PROT-1 and E-PROT-27; Ravenhill et al. 2020. DAVID knowledgebase
v2025_1. SURFY (Bausch-Fluck et al. 2018), the Cell Surface Protein Atlas (Bausch-Fluck et al.
2015) and the Almen classification (Almén et al. 2009). The papers behind each modelling rule
are in `data/curated/CITATIONS.md`.
