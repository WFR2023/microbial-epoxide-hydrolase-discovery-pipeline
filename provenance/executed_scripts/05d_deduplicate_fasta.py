#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


BASE_DIR = Path("<PROJECT_ROOT>")
IN_DIR = BASE_DIR / "06_candidate_mining" / "refined"
OUT_DIR = BASE_DIR / "06_candidate_mining" / "deduplicated"
OUT_DIR.mkdir(parents=True, exist_ok=True)

FILES = {
    "strong": IN_DIR / "candidate_hits_refined_strong.faa",
    "moderate": IN_DIR / "candidate_hits_refined_moderate.faa",
}


def deduplicate_fasta(infile, outfile):
    seen = {}
    total = 0

    for record in SeqIO.parse(infile, "fasta"):
        total += 1
        seq = str(record.seq)
        if seq not in seen:
            seen[seq] = record

    unique_records = list(seen.values())
    SeqIO.write(unique_records, outfile, "fasta")

    return total, len(unique_records)


def main():
    for label, infile in FILES.items():
        if not infile.exists():
            print(f"Skipping missing file: {infile}")
            continue

        outfile = OUT_DIR / f"{label}_deduplicated.faa"
        total, unique_n = deduplicate_fasta(infile, outfile)

        print(f"{label}: total={total} unique={unique_n} outfile={outfile}")


if __name__ == "__main__":
    main()
