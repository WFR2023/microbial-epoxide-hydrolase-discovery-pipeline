# Reproducibility and provenance

## Reproducibility strategy

This repository separates reusable workflow code from historical execution provenance.

### Canonical workflow

`scripts/` contains the portable, sanitized workflow intended for reuse.

### Historical execution provenance

`provenance/executed_scripts/` contains sanitized archival copies of scripts associated with the original study execution.

Machine-specific absolute paths, local usernames, hostnames, and personal contact information were removed or replaced before public release.

### Checksums

`provenance/historical_original_script_sha256.tsv` contains SHA-256 identifiers for the original internal script snapshot prior to public sanitization.

`provenance/public_sanitized_script_sha256.tsv` contains SHA-256 identifiers for the sanitized archival copies distributed in this repository.

These checksum sets intentionally differ because machine-specific and personal information was removed from the public copies.

## Software provenance

Software and executable-version information is documented in:

`environment/software_versions.tsv`

Where package metadata and executable-reported versions differed, both are retained for transparency.

## Pfam provenance

The Pfam-A database used in the original analysis is not redistributed with this repository.

Its provenance is recorded in:

`environment/pfam_database_provenance.tsv`

The record includes the SHA-256 checksum and HMM profile count of the database used in the analysis.

## External genomic data

Large public genome assemblies are not duplicated in this software repository.

Study-level accession manifests and selected analytical provenance files will be deposited in the associated Open Science Framework reproducibility package.

Genome recurrence estimates refer only to the selected assembly set and must not be interpreted as population prevalence.

## Structural-reference provenance

The Pseudomonas aeruginosa Cif structural reference used in the structural workflow corresponds to PDB 3KD2.

The structural reference is distinct from the Cif sequence reference used during BLASTp-based screening.

## Credentials and local configuration

NCBI credentials are not stored in this repository.

Users should configure a valid NCBI contact email locally:

```bash
export NCBI_EMAIL="your.email@example.org"
```

An NCBI API key may optionally be supplied through the environment:

```bash
export NCBI_API_KEY="your_api_key"
```

The Pfam database path must likewise be provided locally:

```bash
export PFAM_DB="/path/to/Pfam-A.hmm"
```

API keys, credentials, and machine-specific configuration must not be committed to the repository.

## Interpretation limits

Candidate classifications, structural models, pocket predictions, and docking scores are computational evidence layers and do not constitute experimental validation of enzymatic function.
