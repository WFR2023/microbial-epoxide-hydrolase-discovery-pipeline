#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from Bio import SeqIO

BASE_DIR = Path("<PROJECT_ROOT>")

REF_FASTA = BASE_DIR / "07_blast_reference" / "epoxide_hydrolase_references.faa"
CAND_FASTA = BASE_DIR / "08_priority_candidates" / "priority_candidates_final.faa"
OUT_DIR = BASE_DIR / "08_priority_candidates"
OUT_FASTA = OUT_DIR / "priority_candidates_plus_references.faa"


def main():
    if not REF_FASTA.exists():
        raise FileNotFoundError(f"Missing reference FASTA: {REF_FASTA}")
    if not CAND_FASTA.exists():
        raise FileNotFoundError(f"Missing candidate FASTA: {CAND_FASTA}")

    ref_records = list(SeqIO.parse(REF_FASTA, "fasta"))
    cand_records = list(SeqIO.parse(CAND_FASTA, "fasta"))

    combined = ref_records + cand_records
    SeqIO.write(combined, OUT_FASTA, "fasta")

    print(f"Combined FASTA written to: {OUT_FASTA}")
    print(f"Reference sequences: {len(ref_records)}")
    print(f"Candidate sequences: {len(cand_records)}")
    print(f"Total sequences: {len(combined)}")


if __name__ == "__main__":
    main()
