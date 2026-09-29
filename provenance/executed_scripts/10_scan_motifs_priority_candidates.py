#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import re
import pandas as pd
from Bio import SeqIO

BASE_DIR = Path("<PROJECT_ROOT>")
IN_FASTA = BASE_DIR / "08_priority_candidates" / "priority_candidates_plus_references.faa"
OUT_DIR = BASE_DIR / "12_domain_motif_screen"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_TABLE = OUT_DIR / "priority_candidates_motif_scan.tsv"
OUT_SUMMARY = OUT_DIR / "priority_candidates_motif_summary.tsv"

# Common alpha/beta-hydrolase-like nucleophile motifs
MOTIF_PATTERNS = {
    "GXSXG": re.compile(r"G[A-Z]S[A-Z]G"),
    "GxSxG_relaxed": re.compile(r"G.{1}S.{1}G"),
    "nucleophile_serine_window": re.compile(r"[GSTACVILMFPYWHKRQEND]{0,3}S[A-Z]{0,3}"),
}

def find_all(pattern, seq):
    return [(m.start() + 1, m.group()) for m in pattern.finditer(seq)]

def find_residues(seq, aa):
    return [i + 1 for i, c in enumerate(seq) if c == aa]

def nearest_after(pos_list, anchor):
    vals = [p for p in pos_list if p > anchor]
    return vals[0] if vals else None

def classify_row(name, seq):
    seq = str(seq)
    length = len(seq)

    gxsxg_hits = find_all(MOTIF_PATTERNS["GXSXG"], seq)
    relaxed_hits = find_all(MOTIF_PATTERNS["GxSxG_relaxed"], seq)

    ser_positions = find_residues(seq, "S")
    asp_positions = find_residues(seq, "D")
    glu_positions = find_residues(seq, "E")
    his_positions = find_residues(seq, "H")

    best_anchor = None
    best_motif = ""
    asp_after = None
    glu_after = None
    his_after = None

    if gxsxg_hits:
        best_anchor = gxsxg_hits[0][0] + 2  # serine position inside GXSXG
        best_motif = gxsxg_hits[0][1]
    elif relaxed_hits:
        best_anchor = relaxed_hits[0][0] + 2
        best_motif = relaxed_hits[0][1]

    if best_anchor:
        asp_after = nearest_after(asp_positions, best_anchor)
        glu_after = nearest_after(glu_positions, best_anchor)
        his_after = nearest_after(his_positions, best_anchor)

    # Heuristic classification
    heuristic_class = "uncertain"

    if "pseudomonas_Cif" in name:
        heuristic_class = "reference_cif"
    elif "human_EPHX1" in name:
        heuristic_class = "reference_ephx1"
    elif "human_EPHX2" in name:
        heuristic_class = "reference_ephx2"
    else:
        if best_anchor and his_after:
            if length >= 280 and length <= 340:
                heuristic_class = "cif_like_or_compact_ab_hydrolase"
            elif length >= 350:
                heuristic_class = "ephx_like_or_large_ab_hydrolase"
            else:
                heuristic_class = "ab_hydrolase_like"
        elif best_anchor:
            heuristic_class = "partial_ab_hydrolase_signal"
        else:
            heuristic_class = "no_clear_motif"

    return {
        "sequence_id": name,
        "sequence_length_aa": length,
        "gxsxg_strict_count": len(gxsxg_hits),
        "gxsxg_relaxed_count": len(relaxed_hits),
        "first_gxsxg_motif": gxsxg_hits[0][1] if gxsxg_hits else "",
        "first_gxsxg_start": gxsxg_hits[0][0] if gxsxg_hits else "",
        "best_anchor_serine_pos": best_anchor if best_anchor else "",
        "nearest_Asp_after_anchor": asp_after if asp_after else "",
        "nearest_Glu_after_anchor": glu_after if glu_after else "",
        "nearest_His_after_anchor": his_after if his_after else "",
        "n_total_ser": len(ser_positions),
        "n_total_asp": len(asp_positions),
        "n_total_glu": len(glu_positions),
        "n_total_his": len(his_positions),
        "heuristic_class": heuristic_class,
    }

def main():
    if not IN_FASTA.exists():
        raise FileNotFoundError(f"Missing FASTA: {IN_FASTA}")

    rows = []
    for rec in SeqIO.parse(IN_FASTA, "fasta"):
        rows.append(classify_row(rec.id, rec.seq))

    df = pd.DataFrame(rows)
    df.to_csv(OUT_TABLE, sep="\t", index=False)

    summary = (
        df.groupby("heuristic_class")
        .size()
        .reset_index(name="n_sequences")
        .sort_values("n_sequences", ascending=False)
    )
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    print(f"Motif scan table written to: {OUT_TABLE}")
    print(f"Summary written to: {OUT_SUMMARY}")
    print(f"Sequences scanned: {len(df)}")

if __name__ == "__main__":
    main()
