# Execution provenance

This directory contains sanitized archival copies of scripts used during the
development and execution of the study workflow.

These files are retained for computational provenance and should not be
considered the canonical user-facing pipeline. Machine-specific absolute paths,
local usernames, hostnames, and personal contact information were removed or
replaced with generic placeholders before public release.

The file `historical_original_script_sha256.tsv` records SHA-256 checksums of
the original internal script snapshot before public sanitization.

The file `public_sanitized_script_sha256.tsv` records SHA-256 checksums of the
sanitized archival scripts included in this repository.

Development-only utilities, machine-audit scripts, and superseded placeholder
scripts were excluded from the public provenance package.

The portable and documented workflow intended for reuse is provided separately
under `scripts/`.
