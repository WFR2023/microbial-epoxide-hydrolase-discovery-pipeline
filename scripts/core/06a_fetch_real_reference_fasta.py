#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import requests

BASE_DIR = Path(".")
REF_DIR = BASE_DIR / "07_blast_reference"
REF_DIR.mkdir(parents=True, exist_ok=True)

OUT_FASTA = REF_DIR / "epoxide_hydrolase_references.faa"
MANIFEST = REF_DIR / "epoxide_hydrolase_reference_manifest.tsv"

# Curated real references
# Human references from UniProt reviewed entries
# Cif from NCBI Protein
REFERENCES = [
    {
        "label": "human_EPHX2_sEH",
        "source_db": "UniProt",
        "accession": "P34913",
        "description": "Bifunctional epoxide hydrolase 2 (EPHX2), Homo sapiens, reviewed",
        "url": "https://rest.uniprot.org/uniprotkb/P34913.fasta",
    },
    {
        "label": "human_EPHX1_mEH",
        "source_db": "UniProt",
        "accession": "P07099",
        "description": "Epoxide hydrolase 1 (EPHX1), Homo sapiens, reviewed",
        "url": "https://rest.uniprot.org/uniprotkb/P07099.fasta",
    },
    {
        "label": "pseudomonas_Cif",
        "source_db": "NCBI_Protein",
        "accession": "AAG06322.1",
        "description": "CFTR inhibitory factor Cif, Pseudomonas aeruginosa PAO1",
        "url": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=AAG06322.1&rettype=fasta&retmode=text",
    },
]

rows = []
fasta_blocks = []

for ref in REFERENCES:
    r = requests.get(ref["url"], timeout=60)
    r.raise_for_status()

    text = r.text.strip()
    if not text.startswith(">"):
        raise RuntimeError(f"Downloaded content does not look like FASTA for {ref['accession']}")

    # Rewrite header to a standardized format
    lines = text.splitlines()
    seq = "".join(lines[1:]).strip()

    header = f">{ref['label']} accession={ref['accession']} source={ref['source_db']} desc={ref['description']}"
    fasta_blocks.append(header)
    for i in range(0, len(seq), 60):
        fasta_blocks.append(seq[i:i+60])

    rows.append({
        "label": ref["label"],
        "source_db": ref["source_db"],
        "accession": ref["accession"],
        "description": ref["description"],
        "url": ref["url"],
        "sequence_length_aa": len(seq),
    })

with open(OUT_FASTA, "w", encoding="utf-8") as f:
    f.write("\n".join(fasta_blocks) + "\n")

import pandas as pd
pd.DataFrame(rows).to_csv(MANIFEST, sep="\t", index=False)

print(f"Reference FASTA written to: {OUT_FASTA}")
print(f"Reference manifest written to: {MANIFEST}")
print(f"Total references fetched: {len(rows)}")
