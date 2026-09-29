#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import csv
import math
from collections import defaultdict

from pymol import cmd

# =========================================================
# KING STRUCTURAL FIGURE SCRIPT
# =========================================================
# Expected project structure:
# <PROJECT_ROOT>/19_structural_pipeline_real/
#   ├── 02_models/pdb/
#   ├── 02_models/manual_models/
#   ├── 05_visualization/
#   └── 06_tables/docking_summary.tsv
#
# Run with:
#   conda activate seh_struct
#   pymol -cq <PROJECT_ROOT>/scripts/19i_make_king_structural_figures.py
#
# Optional environment variables:
#   STRUCT_BASE=/path/to/19_structural_pipeline_real
#   POCKET_DISTANCE=5.0
#   RAY_W=2400
#   RAY_H=1800
# =========================================================

BASE_DIR = os.environ.get(
    "STRUCT_BASE",
    "<PROJECT_ROOT>/19_structural_pipeline_real"
)

DOCKING_SUMMARY = os.path.join(BASE_DIR, "06_tables", "docking_summary.tsv")
OUT_ROOT = os.path.join(BASE_DIR, "05_visualization", "king_figures")
POCKET_DISTANCE = float(os.environ.get("POCKET_DISTANCE", "5.0"))
RAY_W = int(os.environ.get("RAY_W", "2400"))
RAY_H = int(os.environ.get("RAY_H", "1800"))
DPI = int(os.environ.get("DPI", "300"))

DIRS = {
    "best_overview": os.path.join(OUT_ROOT, "01_best_overview"),
    "best_closeup": os.path.join(OUT_ROOT, "02_best_closeup"),
    "gly_overview": os.path.join(OUT_ROOT, "03_glycidol_overview"),
    "gly_closeup": os.path.join(OUT_ROOT, "04_glycidol_closeup"),
    "sessions": os.path.join(OUT_ROOT, "05_sessions"),
    "contact_sheets": os.path.join(OUT_ROOT, "06_contact_sheets"),
    "tables": os.path.join(OUT_ROOT, "07_tables"),
}

for d in DIRS.values():
    os.makedirs(d, exist_ok=True)


def safe_float(x):
    try:
        return float(str(x).strip())
    except Exception:
        return None


def read_tsv(path):
    rows = []
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def get_first_existing(*paths):
    for p in paths:
        if p and os.path.exists(p):
            return p
    return None


def model_class(model_id):
    if model_id.startswith("human_"):
        return "human_reference"
    if model_id.startswith("pseudomonas_"):
        return "bacterial_reference"
    return "candidate"


def protein_color(model_id):
    cls = model_class(model_id)
    if cls == "human_reference":
        return "forest"
    if cls == "bacterial_reference":
        return "tv_orange"
    return "marine"


def ligand_color(ligand_name):
    name = ligand_name.lower()
    if "glycidol" in name:
        return "magenta"
    if "propylene" in name:
        return "hotpink"
    if "ethylene" in name:
        return "cyan"
    return "deeppurple"


def pretty_caption(model_id, ligand, score):
    return f"{model_id}\n{ligand} | score = {score:.3f} kcal/mol"


def find_model_pdb(model_id, row=None):
    candidates = []

    if row:
        for key in ("model_pdb", "protein_pdb", "receptor_pdb", "pdb_path"):
            if key in row and row[key]:
                candidates.append(row[key])

    candidates.extend([
        os.path.join(BASE_DIR, "02_models", "pdb", f"{model_id}.pdb"),
        os.path.join(BASE_DIR, "02_models", "manual_models", f"{model_id}.pdb"),
    ])
    return get_first_existing(*candidates)


def get_pose_pdb(row):
    for key in ("pose_pdb", "pose_pdb_path", "pose_pdb_file"):
        if key in row and row[key]:
            return row[key]
    return None


def setup_publication_style():
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("orthoscopic", 1)
    cmd.set("depth_cue", 0)
    cmd.set("antialias", 2)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_sampling", 10)
    cmd.set("specular", 0.2)
    cmd.set("shininess", 10)
    cmd.set("ambient", 0.25)
    cmd.set("direct", 0.65)
    cmd.set("reflect", 0.1)
    cmd.set("stick_radius", 0.18)
    cmd.set("sphere_scale", 0.23)
    cmd.set("surface_quality", 1)
    cmd.set("transparency_mode", 3)
    cmd.set("label_font_id", 7)
    cmd.set("label_size", -0.35)
    cmd.set("label_color", "black")
    cmd.set("ray_shadows", 0)
    cmd.set("valence", 0)


def render_scene(model_id, protein_pdb, pose_pdb, ligand_name, score, out_png,
                 out_pse=None, closeup=False):
    cmd.reinitialize()
    setup_publication_style()

    prot = "prot"
    lig = "lig"
    pocket = "pocket"
    pocket_surface = "pocket_surface"
    title = "title"

    cmd.load(protein_pdb, prot)
    cmd.load(pose_pdb, lig)

    cmd.remove("solvent")
    cmd.hide("everything", "all")

    cmd.show("cartoon", prot)
    cmd.color(protein_color(model_id), prot)

    cmd.show("sticks", lig)
    cmd.color(ligand_color(ligand_name), lig)

    cmd.select(pocket, f"byres ({prot} within {POCKET_DISTANCE} of {lig}) and polymer.protein")
    cmd.show("sticks", pocket)
    cmd.color("yelloworange", pocket)

    cmd.select(pocket_surface, f"({prot} within {POCKET_DISTANCE + 1.5} of {lig}) and polymer.protein")
    cmd.show("surface", pocket_surface)
    cmd.set("transparency", 0.35, pocket_surface)
    cmd.color("grey80", pocket_surface)

    # Ligand emphasis
    cmd.show("spheres", lig)
    cmd.set("sphere_scale", 0.22, lig)

    if closeup:
        cmd.orient(lig)
        cmd.zoom(f"{lig} or {pocket}", 8)
        cmd.label(f"({pocket}) and name CA", '"%s%s" % (resn,resi)')
    else:
        cmd.orient(prot)
        cmd.zoom(prot, 3)

    # On-figure caption using pseudoatom
    cmd.pseudoatom(title, pos=[0, 0, 0], label=f"{model_id}\\n{ligand_name} | {score:.3f} kcal/mol")
    cmd.hide("everything", title)
    cmd.show("labels", title)
    cmd.set("label_size", 22, title)
    cmd.set("label_color", "black", title)

    # Move label roughly into empty canvas space
    if closeup:
        cmd.translate([0, 0, 0], title)
    else:
        cmd.translate([0, 0, 0], title)

    cmd.png(out_png, width=RAY_W, height=RAY_H, dpi=DPI, ray=1)

    if out_pse:
        cmd.save(out_pse)


def make_contact_sheet(items, out_path, title=""):
    try:
        from PIL import Image, ImageDraw, ImageFont
    except Exception as e:
        print(f"[WARN] Pillow not available. Skipping contact sheet: {out_path}")
        print(f"       Reason: {e}")
        return

    if not items:
        return

    font = ImageFont.load_default()
    title_font = ImageFont.load_default()

    cell_w = 1200
    cell_h = 980
    text_h = 120
    pad = 20
    n = len(items)
    ncols = 2 if n <= 4 else 3
    nrows = math.ceil(n / ncols)
    top_h = 90 if title else 20

    sheet = Image.new("RGB", (ncols * cell_w, nrows * cell_h + top_h), "white")
    draw_sheet = ImageDraw.Draw(sheet)

    if title:
        draw_sheet.text((20, 20), title, fill="black", font=title_font)

    for idx, (img_path, caption) in enumerate(items):
        r = idx // ncols
        c = idx % ncols
        x0 = c * cell_w
        y0 = r * cell_h + top_h

        cell = Image.new("RGB", (cell_w, cell_h), "white")
        draw = ImageDraw.Draw(cell)

        img = Image.open(img_path).convert("RGB")
        img.thumbnail((cell_w - 2 * pad, cell_h - text_h - 2 * pad))

        ix = (cell_w - img.width) // 2
        iy = pad
        cell.paste(img, (ix, iy))

        draw.rectangle((0, 0, cell_w - 1, cell_h - 1), outline=(180, 180, 180), width=2)
        draw.multiline_text((pad, cell_h - text_h + 15), caption, fill="black", font=font, spacing=4)

        sheet.paste(cell, (x0, y0))

    sheet.save(out_path, quality=95)
    print(f"[OK] Contact sheet written: {out_path}")


def write_manifest(path, rows):
    fields = [
        "series",
        "model_id",
        "ligand_name",
        "vina_best_score_kcal_mol",
        "protein_pdb",
        "pose_pdb",
        "overview_png",
        "closeup_png",
        "session_pse"
    ]
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def main():
    if not os.path.exists(DOCKING_SUMMARY):
        raise FileNotFoundError(f"Docking summary not found: {DOCKING_SUMMARY}")

    rows = read_tsv(DOCKING_SUMMARY)

    usable = []
    for row in rows:
        model_id = row.get("model_id", "").strip()
        ligand_name = row.get("ligand_name", "").strip()
        score = safe_float(row.get("vina_best_score_kcal_mol", ""))
        pose_pdb = get_pose_pdb(row)

        if not model_id or not ligand_name or score is None:
            continue
        if not pose_pdb or not os.path.exists(pose_pdb):
            continue

        protein_pdb = find_model_pdb(model_id, row=row)
        if not protein_pdb:
            print(f"[WARN] Protein PDB not found for {model_id}. Skipping.")
            continue

        row["_score"] = score
        row["_protein_pdb"] = protein_pdb
        row["_pose_pdb"] = pose_pdb
        usable.append(row)

    if not usable:
        raise RuntimeError("No usable docking rows found.")

    rows_by_model = defaultdict(list)
    for row in usable:
        rows_by_model[row["model_id"]].append(row)

    best_rows = {}
    gly_rows = {}

    for model_id, model_rows in rows_by_model.items():
        model_rows = sorted(model_rows, key=lambda x: x["_score"])
        best_rows[model_id] = model_rows[0]

        gly = [r for r in model_rows if r["ligand_name"].strip().lower() == "glycidol"]
        if gly:
            gly_rows[model_id] = sorted(gly, key=lambda x: x["_score"])[0]

    best_manifest = []
    gly_manifest = []
    best_overview_items = []
    best_closeup_items = []
    gly_overview_items = []
    gly_closeup_items = []

    # -----------------------------------------------------
    # Render best-per-model figures
    # -----------------------------------------------------
    for model_id in sorted(best_rows):
        row = best_rows[model_id]
        ligand = row["ligand_name"].strip()
        score = row["_score"]
        protein_pdb = row["_protein_pdb"]
        pose_pdb = row["_pose_pdb"]

        overview_png = os.path.join(DIRS["best_overview"], f"{model_id}__best_overview.png")
        closeup_png = os.path.join(DIRS["best_closeup"], f"{model_id}__best_closeup.png")
        session_pse = os.path.join(DIRS["sessions"], f"{model_id}__best_scene.pse")

        print(f"[BEST] Rendering {model_id} with {ligand} ({score:.3f})")
        render_scene(model_id, protein_pdb, pose_pdb, ligand, score, overview_png, closeup=False)
        render_scene(model_id, protein_pdb, pose_pdb, ligand, score, closeup_png, out_pse=session_pse, closeup=True)

        cap = pretty_caption(model_id, ligand, score)
        best_overview_items.append((overview_png, cap))
        best_closeup_items.append((closeup_png, cap))

        best_manifest.append({
            "series": "best_per_model",
            "model_id": model_id,
            "ligand_name": ligand,
            "vina_best_score_kcal_mol": f"{score:.3f}",
            "protein_pdb": protein_pdb,
            "pose_pdb": pose_pdb,
            "overview_png": overview_png,
            "closeup_png": closeup_png,
            "session_pse": session_pse,
        })

    # -----------------------------------------------------
    # Render glycidol-focused figures
    # -----------------------------------------------------
    for model_id in sorted(gly_rows):
        row = gly_rows[model_id]
        ligand = row["ligand_name"].strip()
        score = row["_score"]
        protein_pdb = row["_protein_pdb"]
        pose_pdb = row["_pose_pdb"]

        overview_png = os.path.join(DIRS["gly_overview"], f"{model_id}__glycidol_overview.png")
        closeup_png = os.path.join(DIRS["gly_closeup"], f"{model_id}__glycidol_closeup.png")
        session_pse = os.path.join(DIRS["sessions"], f"{model_id}__glycidol_scene.pse")

        print(f"[GLYCIDOL] Rendering {model_id} with {ligand} ({score:.3f})")
        render_scene(model_id, protein_pdb, pose_pdb, ligand, score, overview_png, closeup=False)
        render_scene(model_id, protein_pdb, pose_pdb, ligand, score, closeup_png, out_pse=session_pse, closeup=True)

        cap = pretty_caption(model_id, ligand, score)
        gly_overview_items.append((overview_png, cap))
        gly_closeup_items.append((closeup_png, cap))

        gly_manifest.append({
            "series": "glycidol_focus",
            "model_id": model_id,
            "ligand_name": ligand,
            "vina_best_score_kcal_mol": f"{score:.3f}",
            "protein_pdb": protein_pdb,
            "pose_pdb": pose_pdb,
            "overview_png": overview_png,
            "closeup_png": closeup_png,
            "session_pse": session_pse,
        })

    # -----------------------------------------------------
    # Contact sheets
    # -----------------------------------------------------
    make_contact_sheet(
        best_overview_items,
        os.path.join(DIRS["contact_sheets"], "best_overview_contact_sheet.png"),
        title="Best docking pose per model – overview"
    )
    make_contact_sheet(
        best_closeup_items,
        os.path.join(DIRS["contact_sheets"], "best_closeup_contact_sheet.png"),
        title="Best docking pose per model – pocket close-up"
    )
    make_contact_sheet(
        gly_overview_items,
        os.path.join(DIRS["contact_sheets"], "glycidol_overview_contact_sheet.png"),
        title="Glycidol series – overview"
    )
    make_contact_sheet(
        gly_closeup_items,
        os.path.join(DIRS["contact_sheets"], "glycidol_closeup_contact_sheet.png"),
        title="Glycidol series – pocket close-up"
    )

    # -----------------------------------------------------
    # Write manifest tables
    # -----------------------------------------------------
    write_manifest(
        os.path.join(DIRS["tables"], "king_best_per_model_manifest.tsv"),
        best_manifest
    )
    write_manifest(
        os.path.join(DIRS["tables"], "king_glycidol_manifest.tsv"),
        gly_manifest
    )

    print("\n[OK] KING structural figure generation completed.")
    print(f"[OUT] Root directory: {OUT_ROOT}")


main()
cmd.quit()
