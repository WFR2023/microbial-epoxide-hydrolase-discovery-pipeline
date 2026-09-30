# Microbial Epoxide Hydrolase Discovery Pipeline

Reproducible computational workflow for genome mining, prioritization, contextual analysis, structural modeling, and exploratory molecular docking of microbial epoxide hydrolase-like proteins.

## Overview

This repository contains the computational workflow developed to identify and prioritize microbial proteins with sequence, domain, phylogenetic, genomic-context, structural, and docking characteristics compatible with epoxide hydrolase-like or related alpha/beta-hydrolase profiles.

The workflow was designed as a candidate-prioritization framework. Computational classifications do not constitute experimental confirmation of enzymatic activity, catalytic function, substrate specificity, or biological activity.

## Main analytical layers

- genome assembly selection and acquisition;
- protein dataset preparation;
- annotation-based candidate mining;
- BLASTp comparison against curated references;
- candidate prioritization;
- protein descriptor calculation;
- MAFFT multiple-sequence alignment;
- IQ-TREE maximum-likelihood phylogeny;
- motif screening;
- Pfam/HMMER domain analysis;
- within-species recurrence analysis;
- genomic-neighborhood profiling;
- functional-context scoring;
- ColabFold structural modeling;
- fpocket cavity prediction;
- AutoDock Vina molecular docking;
- programmatic generation of analytical outputs.

## Repository structure

```text
.
├── config/
│   ├── target_species.tsv
│   └── runtime.env.example
├── docs/
├── environment/
│   ├── software_versions.tsv
│   └── pfam_database_provenance.tsv
├── example/
├── provenance/
│   ├── README.md
│   ├── historical_original_script_sha256.tsv
│   ├── public_sanitized_script_sha256.tsv
│   └── executed_scripts/
└── scripts/
    ├── core/
    ├── fallback/
    ├── optional/
    └── structural/
```

## Runtime configuration

The workflow uses environment variables for values that should not be hard-coded.

For NCBI Entrez queries:

```bash
export NCBI_EMAIL="your.email@example.org"
```

An NCBI API key is optional:

```bash
export NCBI_API_KEY="your_api_key"
```

For the Pfam/HMMER stage:

```bash
export PFAM_DB="/path/to/Pfam-A.hmm"
```

Credentials and API keys must not be committed to this repository.

## Functional-context score

The functional-context score is a heuristic prioritization aid and is not a statistical probability or experimentally calibrated confidence estimate.

It combines:

- prior candidate-priority class;
- within-species recurrence;
- selected genomic-neighborhood functional categories.

Pfam assignments, motif observations, phylogenetic placement, structural confidence, pocket characteristics, and docking results are evaluated separately and are not independently added to the numerical score.

See `docs/scoring_system.md` for the complete scoring definition.

## Structural analysis and docking

Structural modeling was performed with ColabFold for selected candidates and human reference proteins. The experimental Pseudomonas aeruginosa Cif structure PDB 3KD2 was used as the Cif structural reference.

Binding cavities were predicted with fpocket. For each receptor, the cavity containing the largest number of pocket atoms was selected before docking.

Exploratory docking was performed with AutoDock Vina using ethylene oxide, propylene oxide, and glycidol.

Docking scores are computational predictions and should not be interpreted as measured binding affinities or evidence of catalytic activity, substrate specificity, kinetics, regioselectivity, stereoselectivity, or biological function.

## Reproducibility and provenance

This repository separates reusable code from historical execution provenance:

- `scripts/`: portable canonical workflow;
- `provenance/`: sanitized historical execution records;
- `environment/`: software and database provenance;
- `config/`: study configuration.

SHA-256 checksums are retained to distinguish the original internal script snapshot from the sanitized public provenance files.

Large public genome datasets and the Pfam database are not redistributed in this software repository.

Study-level accession manifests, analytical metadata, and selected reproducibility artifacts are publicly available in the Open Science Framework (OSF) reproducibility package: DOI 10.17605/OSF.IO/BKZSP.

## Funding

The development of this software resource and the research underlying its analytical workflow received support from the Fundação de Amparo à Pesquisa do Estado de Minas Gerais (FAPEMIG) through the Science, Technology, and Innovation Development Program (grant RED-00181-23; Call 012/2023, “Structuring Networks for Scientific Research or Technological Development”); from the Brazilian National Council for Scientific and Technological Development (CNPq) through Research Productivity funding (grant 314660/2026-7; Call CNPq No. 23/2025); and in part from the National Institute of Environmental Health Sciences (NIEHS) through the RIVER Award R35ES030443.

The funding agencies had no role in the computational classification of candidates or in the interpretation of the analytical outputs reported by this software resource.


## Citation

Citation metadata are provided in `CITATION.cff`.

The versioned v1.0.0 software release is permanently archived in Zenodo: https://doi.org/10.5281/zenodo.23042986.

## License

This software is released under the MIT License.

## Software creators

- Wellington Francisco Rodrigues
- Sophie Kiss

Wellington Francisco Rodrigues and Sophie Kiss contributed equally to the software resource.

## Disclaimer

This repository supports computational research and candidate prioritization. Its outputs do not constitute clinical, diagnostic, therapeutic, or experimentally validated biochemical conclusions.
