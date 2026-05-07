#!/usr/bin/env python
"""Render per-case protein-binding-site images with PyMOL (headless).

Each rendered PNG shows:
  - Predicted binding-site residues as cyan sticks
  - Predicted ligand as green sticks
  - Native (crystal) ligand as gray sticks
  - Protein cartoon in transparent cyan

Runs under the PyMOL-PoseBench env (needs `pymol` module).
Produces a grid image per challenge family (4 models × N cases).
"""
from __future__ import annotations
import os, sys, glob
from pathlib import Path

# Add analysis/src to path so we can import config & helpers.
REPO = Path("/mnt/katritch_lab2/aoxu/contrasCF")
sys.path.insert(0, str(REPO / "analysis" / "src"))

# Set env before importing pymol. The PyMOL-PoseBench env ships its own libs.
_ENV_LIB = "/home/aoxu/miniconda3/envs/PyMOL-PoseBench/lib"
if _ENV_LIB not in os.environ.get("LD_LIBRARY_PATH", ""):
    os.environ["LD_LIBRARY_PATH"] = _ENV_LIB + os.pathsep + os.environ.get("LD_LIBRARY_PATH", "")

import pandas as pd  # noqa: E402
import pymol  # noqa: E402  (boots headless)
from pymol import cmd  # noqa: E402

from config import (  # noqa: E402
    CASES, MODELS, NATIVES, NATIVE_DIR, FIGURES_DIR, RESULTS_DIR,
    CDK2_BS_RESIDUES, case_dir, find_first,
)

# Load the target-ligand chain/resname per (model, case) from results.csv so we
# can select exactly the right residue even when cofactors share chain IDs.
_RES = pd.read_csv(RESULTS_DIR / "results.csv")
_LIG_INFO = {
    (r["model"], r["case"]): (r["lig_chain"], r["lig_resname"])
    for _, r in _RES.iterrows()
    if pd.notna(r.get("lig_chain")) and pd.notna(r.get("lig_resname"))
}


# GDH binding site radius (for glucose-family paper-style rendering).
GDH_POCKET_RADIUS = 5.0    # Å

RENDER_DIR = FIGURES_DIR / "per_case"


def reset():
    cmd.reinitialize()
    cmd.bg_color("white")
    cmd.set("ray_shadows", 0)
    cmd.set("ambient", 0.35)
    cmd.set("ray_opaque_background", 1)
    cmd.set("cartoon_transparency", 0.6)
    cmd.set("stick_radius", 0.14)
    cmd.set("orthoscopic", 1)


def _load_and_align(model: str, case: str, native_obj: str = "nat") -> tuple[str, str, bool]:
    """Common setup: load native + predicted, Cα-align, return (pred_lig_sel, native_lig_sel, ok)."""
    d = case_dir(model, case)
    fp = find_first(d, MODELS[model].structure_globs)
    if fp is None:
        return "", "", False

    spec = CASES[case]
    nref = NATIVES[spec.native_key]
    cmd.load(str(NATIVE_DIR / f"{nref.pdb_id}.cif"), native_obj)
    cmd.load(str(fp), "pred")
    cmd.super("pred and polymer.protein", f"{native_obj} and polymer.protein")

    native_lig_sel = f"{native_obj} and resn {nref.ligand_resname}"
    lig_info = _LIG_INFO.get((model, case))
    if lig_info is not None:
        lig_chain, lig_resname = lig_info
        pred_lig_sel = f"pred and chain {lig_chain} and resn {lig_resname}"
    else:
        pred_lig_sel = ("pred and (not polymer.protein) and (not metals) and "
                         "(not resn HOH+NAP+NDP+NAD+MG+ZN+NA+K+ACE+LIG_B)")
    return pred_lig_sel, native_lig_sel, True


def render_paper_style(model: str, case: str, out_png: Path) -> bool:
    """Close-up in the style of Masters et al. 2025 Fig. 1/5.

    - No cartoon, pure sticks
    - Only the 11 CDK2 binding-site residues (or pocket residues for GDH) as sticks
    - Native ligand gray, predicted green, with atomic coloring on heteroatoms
    - Tight zoom on the native ligand (radius ~4 Å)
    """
    reset()
    pred_lig_sel, native_lig_sel, ok = _load_and_align(model, case)
    if not ok:
        print(f"  [skip] {model}/{case} — no structure file")
        return False
    spec = CASES[case]

    cmd.hide("everything")

    # Binding-site sticks (predicted side, to show mutated chemistry).
    if spec.native_key == "cdk2_atp":
        bs_sel = f"pred and polymer.protein and resi {CDK2_BS_RESIDUES}"
    else:
        cmd.select("nat_lig", native_lig_sel)
        cmd.select("bs_tmp", f"byres (pred and polymer.protein) within {GDH_POCKET_RADIUS} of nat_lig")
        bs_sel = "bs_tmp"
    cmd.show("sticks", bs_sel)
    cmd.color("cyan", f"({bs_sel}) and elem C")
    cmd.color("atomic", f"({bs_sel}) and not elem C")

    # Native ligand as gray sticks.
    cmd.show("sticks", native_lig_sel)
    cmd.color("gray60", native_lig_sel)

    # Predicted ligand as green sticks.
    cmd.show("sticks", pred_lig_sel)
    cmd.color("green", f"({pred_lig_sel}) and elem C")
    cmd.color("atomic", f"({pred_lig_sel}) and not elem C")

    # Tight zoom on the native ligand only (paper style).
    cmd.orient(native_lig_sel)
    cmd.zoom(native_lig_sel, 4.0)

    cmd.ray(900, 900)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    cmd.png(str(out_png), dpi=120)
    return True


def render_wide(model: str, case: str, out_png: Path, native_obj: str = "nat") -> bool:
    d = case_dir(model, case)
    fp = find_first(d, MODELS[model].structure_globs)
    if fp is None:
        print(f"  [skip] {model}/{case} — no structure file")
        return False

    reset()
    # Load native (with ligand).
    spec = CASES[case]
    nref = NATIVES[spec.native_key]
    cmd.load(str(NATIVE_DIR / f"{nref.pdb_id}.cif"), native_obj)
    native_lig_selection = f"{native_obj} and resn {nref.ligand_resname}"

    # Load predicted.
    cmd.load(str(fp), "pred")

    # Align predicted protein onto native protein Cα.
    # Use super for robustness against loop differences.
    cmd.super("pred and polymer.protein", f"{native_obj} and polymer.protein")

    # Hide all, then show what we want.
    cmd.hide("everything")
    cmd.show("cartoon", f"{native_obj} and polymer.protein")
    cmd.color("cyan", f"{native_obj} and polymer.protein")

    # Native ligand as gray sticks.
    cmd.show("sticks", native_lig_selection)
    cmd.color("gray60", native_lig_selection)
    cmd.set("stick_radius", 0.15, native_lig_selection)

    # Predicted ligand: use exact chain/resname recorded by the analysis pipeline.
    lig_info = _LIG_INFO.get((model, case))
    if lig_info is not None:
        lig_chain, lig_resname = lig_info
        pred_lig_sel = f"pred and chain {lig_chain} and resn {lig_resname}"
    else:
        # Fall back to broad filter (should not happen after 02_run_analysis.py).
        pred_lig_sel = "pred and (not polymer.protein) and (not metals) and (not resn HOH+NAP+NDP+NAD+MG+ZN+NA+K+ACE+LIG_B)"
    cmd.show("sticks", pred_lig_sel)
    cmd.color("green", pred_lig_sel)
    cmd.color("atomic", pred_lig_sel + " and not elem C")    # recolor non-C by element

    # Binding-site sticks (predicted): residues from paper for CDK2, pocket for GDH.
    if spec.native_key == "cdk2_atp":
        bs_sel = f"pred and polymer.protein and resi {CDK2_BS_RESIDUES}"
    else:
        # all protein residues within 5Å of native ligand
        cmd.select("nat_lig", native_lig_selection)
        cmd.select("bs_tmp", f"byres (pred and polymer.protein) within {GDH_POCKET_RADIUS} of nat_lig")
        bs_sel = "bs_tmp"
    cmd.show("sticks", f"({bs_sel}) and sidechain")
    cmd.color("cyan", f"({bs_sel}) and sidechain and elem C")
    cmd.color("atomic", f"({bs_sel}) and sidechain and not elem C")

    # Zoom to encompass both the native ligand and the predicted ligand
    # (so displaced predictions — e.g. AF3/Chai glucose_3+ — stay in frame).
    cmd.select("view_region", f"({native_lig_selection}) or ({pred_lig_sel})")
    n_lig_atoms = cmd.count_atoms(pred_lig_sel)
    if n_lig_atoms == 0:
        # Safety: if prediction has no atoms selected, just zoom on native.
        cmd.zoom(native_lig_selection, 8.0)
        cmd.orient(native_lig_selection)
    else:
        cmd.orient(native_lig_selection)
        cmd.zoom("view_region", 3.0)

    cmd.ray(900, 700)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    cmd.png(str(out_png), dpi=120)
    return True


# Families where we want the paper-style close-up (crop on native ligand, no cartoon).
# Glucose needs the wide view because heavily methylated variants get ejected from
# the pocket and would otherwise be off-frame.
PAPER_STYLE_FAMILIES = {"bindingsite", "atp_charge"}

MODELS_FOR_RENDER = ["AF3", "Boltz", "Boltz2", "UniDock2", "GNINA", "SurfDock"]


def render_family(family: str) -> None:
    cases = [c for c, sp in CASES.items() if sp.family == family]
    use_paper = family in PAPER_STYLE_FAMILIES
    style = "paper" if use_paper else "wide"
    render_fn = render_paper_style if use_paper else render_wide
    print(f"\n=== rendering family '{family}' in {style} style "
          f"({len(cases)} cases × {len(MODELS_FOR_RENDER)} models) ===")
    for model in MODELS_FOR_RENDER:
        for case in cases:
            out = RENDER_DIR / family / f"{style}_{model}_{case}.png"
            if out.exists():
                print(f"  [cached] {out.name}")
                continue
            ok = render_fn(model, case, out)
            print(f"  [{'OK' if ok else 'skip'}] {out.name}")


def build_grid(family: str) -> None:
    """Stitch per-cell PNGs into a grid image (one row per model, one col per case)."""
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.image import imread

    cases = [c for c, sp in CASES.items() if sp.family == family]
    style = "paper" if family in PAPER_STYLE_FAMILIES else "wide"
    fig, axes = plt.subplots(len(MODELS_FOR_RENDER), len(cases),
                             figsize=(2.6 * len(cases), 2.3 * len(MODELS_FOR_RENDER)))
    if len(MODELS_FOR_RENDER) == 1:
        axes = np.array([axes])
    if len(cases) == 1:
        axes = axes.reshape(-1, 1)
    for i, m in enumerate(MODELS_FOR_RENDER):
        for j, c in enumerate(cases):
            ax = axes[i, j]
            png = RENDER_DIR / family / f"{style}_{m}_{c}.png"
            if png.exists():
                ax.imshow(imread(png))
            ax.set_xticks([]); ax.set_yticks([])
            if i == 0:
                ax.set_title(c.replace(f"{family}_", "").replace("bindingsite_", ""), fontsize=9)
            if j == 0:
                ax.set_ylabel(m, fontsize=10, rotation=0, labelpad=25, va="center")
    title_suffix = "paper-style close-up" if style == "paper" else "wide view"
    fig.suptitle(
        f"Per-case binding-site poses — family: {family} ({title_suffix})\n"
        f"native ligand gray; predicted ligand green; binding-site sidechains cyan",
        fontsize=11,
    )
    fig.tight_layout()
    out = FIGURES_DIR / f"grid_{family}.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"  wrote {out}")


def main():
    pymol.finish_launching(["pymol", "-cq"])  # command-line, quiet
    for fam in ["bindingsite", "atp_charge", "glucose"]:
        render_family(fam)
        build_grid(fam)


if __name__ == "__main__":
    main()
