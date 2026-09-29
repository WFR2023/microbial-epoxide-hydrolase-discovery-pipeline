# Workflow overview

The computational workflow is organized into sequential analytical layers.

## 1. Genome dataset construction

Target species are defined in `config/target_species.tsv`.

Genome assemblies are queried and prioritized using NCBI Assembly metadata, with RefSeq/GCF assemblies preferred when available and GenBank/GCA assemblies used when necessary.

## 2. Genome acquisition and inventory

Selected assemblies are retrieved from NCBI resources and organized through accession-resolved manifests and inventory files.

## 3. Protein dataset preparation

Native protein FASTA files are used when available.

When necessary, coding sequences and translated proteins are recovered from annotated GenBank files.

## 4. Candidate mining

Protein annotations are screened for epoxide-hydrolase and broader alpha/beta-hydrolase-related terms.

Candidate sequences are deduplicated by exact amino-acid sequence before downstream prioritization.

## 5. Reference comparison

Candidate proteins are compared by BLASTp against curated epoxide hydrolase references representing human EPHX1/mEH, human EPHX2/sEH, and Pseudomonas aeruginosa Cif.

## 6. Candidate prioritization

Annotation specificity and BLASTp support are used to define the prioritized candidate set.

Priority labels are computational categories and should not be interpreted as experimental confidence classes.

## 7. Protein characterization and phylogenetic analysis

Prioritized proteins are characterized using sequence-derived descriptors.

Multiple-sequence alignment is performed with MAFFT, followed by maximum-likelihood phylogenetic inference with IQ-TREE.

## 8. Motif analysis

Sequences are screened for G-X-S-X-G-related motif patterns.

Motif observations are treated as sequence-level compatibility evidence and are not interpreted as proof of a catalytic triad or enzymatic activity.

## 9. Pfam domain analysis

Domain annotation is performed with HMMER/hmmscan against a locally supplied Pfam-A database.

The database used in the original study is identified through checksum-based provenance in `environment/pfam_database_provenance.tsv`.

## 10. Within-species recurrence and genomic context

Candidate recurrence is evaluated across the selected assemblies of the corresponding species.

Local genomic neighborhoods are classified into functional-context categories from neighboring protein-product annotations.

Recurrence values describe only the analyzed assembly set and are not prevalence estimates.

## 11. Functional-context scoring

A heuristic score summarizes three components:

- prior candidate-priority class;
- within-species recurrence;
- selected genomic-neighborhood categories.

Pfam assignments, motif observations, phylogenetic placement, structural modeling, pocket characteristics, and docking results remain analytically separate from the numerical score.

See `docs/scoring_system.md` for the exact scoring framework.

## 12. Structural modeling

Selected candidate proteins and human reference proteins are represented using ColabFold-derived structural models.

The experimental Pseudomonas aeruginosa Cif structure PDB 3KD2 is used as the Cif structural reference.

## 13. Binding-pocket prediction

Putative cavities are predicted with fpocket.

For each receptor, pocket atom-coordinate files are parsed and the cavity containing the largest number of pocket atoms is selected for docking.

The atom count is used only as a reproducible within-receptor selection heuristic and is not treated as an experimentally measured cavity volume.

## 14. Molecular docking

Ethylene oxide, propylene oxide, and glycidol are prepared with Open Babel.

AutoDock Vina docking is performed in a 22 x 22 x 22 Angstrom cubic search box centered on the selected pocket centroid.

The top-ranked Vina score is retained for comparative analysis.

Docking is exploratory and does not establish binding affinity, catalytic turnover, substrate specificity, enzyme kinetics, regioselectivity, stereoselectivity, or biological activity.

## 15. Final analytical outputs

Candidate summaries, sequence-similarity outputs, motif and Pfam results, recurrence records, genomic-context classifications, functional-context scores, structural metadata, pocket coordinates, docking definitions, and docking summaries are retained as machine-readable analytical provenance.
