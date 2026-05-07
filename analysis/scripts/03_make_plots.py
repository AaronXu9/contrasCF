#!/usr/bin/env python
"""Produce summary plots from results.csv.

MODEL_ORDER groups the co-folding models (AF3, Boltz, Boltz2) and the
pose-prediction / docking methods (UniDock2, GNINA, SurfDock). Chai and
RFAA rows remain in the CSV but are not plotted here.
"""
from __future__ import annotations
import os, sys
from pathlib import Path

_ENV_LIB = "/home/aoxu/miniconda3/envs/rdkit_env/lib"
if _ENV_LIB not in os.environ.get("LD_LIBRARY_PATH", ""):
    os.environ["LD_LIBRARY_PATH"] = _ENV_LIB + os.pathsep + os.environ.get("LD_LIBRARY_PATH", "")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from config import RESULTS_DIR, FIGURES_DIR, CASES  # noqa: E402


MODEL_ORDER = ["AF3", "Boltz", "Boltz2", "UniDock2", "GNINA", "SurfDock"]
CASE_ORDER = list(CASES.keys())
MODEL_COLORS = {
    "AF3":      "#4C72B0",
    "Boltz":    "#C44E52",
    "Boltz2":   "#D62728",
    "UniDock2": "#55A868",
    "GNINA":    "#DD8452",
    "SurfDock": "#8C564B",
    "RFAA":     "#8172B2",
    "Chai":     "#937860",
}
FAMILY_ORDER = ["bindingsite", "atp_charge", "glucose"]


def main() -> None:
    df = pd.read_csv(RESULTS_DIR / "results.csv")
    df = df[df["model"].isin(MODEL_ORDER)].copy()
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # --- 1. Heatmap of ligand_rmsd_common --------------------------------
    pivot = df.pivot(index="model", columns="case", values="ligand_rmsd_common")
    pivot = pivot.reindex(index=MODEL_ORDER, columns=CASE_ORDER)

    fig, ax = plt.subplots(figsize=(14, 3.5))
    sns.heatmap(
        pivot, cmap="rocket_r", vmin=0, vmax=8,
        annot=True, fmt=".1f", annot_kws={"fontsize": 8},
        cbar_kws={"label": "ligand RMSD (Å, core)"},
        ax=ax, linewidth=0.3, linecolor="white",
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(
        "Ligand core-RMSD (Å) to native — per model × case\n"
        "(values >8 Å saturate; the ligand has left the pocket)"
    )
    plt.xticks(rotation=35, ha="right")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "rmsd_heatmap.png", dpi=150)
    plt.close()
    print(f"  wrote {FIGURES_DIR / 'rmsd_heatmap.png'}")

    # --- 2. Fig.3-style bar chart: fraction correct per challenge class ----
    df2 = df.copy()
    df2["correct_2A"] = df2["ligand_rmsd_common"] < 2.0
    df2["correct_4A"] = df2["ligand_rmsd_common"] < 4.0
    frac = df2.groupby(["model", "family"])[["correct_2A", "correct_4A"]].mean().reset_index()

    fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
    for ax, metric, label in zip(axes, ["correct_2A", "correct_4A"], ["RMSD < 2 Å", "RMSD < 4 Å"]):
        pv = frac.pivot(index="family", columns="model", values=metric)
        pv = pv.reindex(index=FAMILY_ORDER, columns=MODEL_ORDER)
        pv.plot.bar(ax=ax, color=[MODEL_COLORS[m] for m in pv.columns], width=0.8, edgecolor="none")
        ax.set_ylabel("fraction of cases")
        ax.set_xlabel("")
        ax.set_title(label)
        ax.set_ylim(0, 1)
        ax.grid(axis="y", alpha=0.3)
        ax.legend(title=None, loc="upper right", fontsize=8)
        plt.setp(ax.get_xticklabels(), rotation=0)
    fig.suptitle("Fraction of near-native predictions per challenge family")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig3_bars.png", dpi=150)
    plt.close()
    print(f"  wrote {FIGURES_DIR / 'fig3_bars.png'}")

    # --- 3. Confidence vs RMSD scatter ------------------------------------
    # AF3 / Boltz use ligand_iptm. UniDock2 / GNINA have no iPTM, so we put
    # their docking score on a second subplot.
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))

    # Left: co-folding ligand_iptm
    ax = axes[0]
    for m in ("AF3", "Boltz", "Boltz2"):
        sub = df[df["model"] == m]
        ax.scatter(sub["ligand_rmsd_common"], sub["ligand_iptm"],
                   color=MODEL_COLORS[m], label=m, s=50, edgecolor="white", alpha=0.85)
    ax.axvline(2.0, color="gray", ls="--", alpha=0.5)
    ax.set_xlabel("ligand core RMSD (Å)")
    ax.set_ylabel("ligand-iPTM (co-folding confidence)")
    ax.set_title("Co-folding confidence vs RMSD")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    # Right: docking scores — UniDock2 Vina, GNINA CNNaffinity, SurfDock conf
    ax = axes[1]
    for m, col in [("UniDock2", "dock_vina_score"),
                   ("GNINA", "dock_cnn_affinity"),
                   ("SurfDock", "surfdock_confidence")]:
        sub = df[df["model"] == m]
        ax.scatter(sub["ligand_rmsd_common"], sub[col],
                   color=MODEL_COLORS[m], label=f"{m} ({col})", s=50, edgecolor="white", alpha=0.85)
    ax.axvline(2.0, color="gray", ls="--", alpha=0.5)
    ax.set_xlabel("ligand core RMSD (Å)")
    ax.set_ylabel("docking score  (Vina kcal/mol, GNINA pK, SurfScore a.u.)")
    ax.set_title("Docking / pose-predictor score vs RMSD")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    fig.suptitle("Confidence/score vs. ligand RMSD")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "conf_vs_rmsd.png", dpi=150)
    plt.close()
    print(f"  wrote {FIGURES_DIR / 'conf_vs_rmsd.png'}")

    # --- 4. Score-per-case (physics signal figure) -----------------------
    # Per family (row) × per case (x-axis), plot each model's primary score.
    # The physics question: within a family, does the score degrade as the
    # ligand is perturbed away from the WT/natural chemistry?
    fig, axes = plt.subplots(3, 2, figsize=(13, 9), sharex=False)
    for i, fam in enumerate(FAMILY_ORDER):
        sub = df[df["family"] == fam]
        cases = [c for c in CASE_ORDER if CASES[c].family == fam]
        x = np.arange(len(cases))

        # Left: ligand core RMSD per model.
        ax = axes[i, 0]
        for m in MODEL_ORDER:
            y = [sub[(sub["case"] == c) & (sub["model"] == m)]["ligand_rmsd_common"].mean()
                 for c in cases]
            ax.plot(x, y, marker="o", color=MODEL_COLORS[m], label=m, linewidth=1.8)
        ax.axhline(2.0, color="gray", ls="--", alpha=0.5)
        ax.set_xticks(x); ax.set_xticklabels(cases, rotation=30, ha="right", fontsize=8)
        ax.set_ylabel("ligand core RMSD (Å)")
        ax.set_title(f"{fam}: pose quality")
        ax.grid(alpha=0.3)
        if i == 0:
            ax.legend(fontsize=8, loc="upper left")

        # Right: score axes — co-folding on left y (ligand_iptm), docking on right y.
        ax = axes[i, 1]
        ax2 = ax.twinx()
        for m in ("AF3", "Boltz", "Boltz2"):
            y = [sub[(sub["case"] == c) & (sub["model"] == m)]["ligand_iptm"].mean()
                 for c in cases]
            ax.plot(x, y, marker="o", color=MODEL_COLORS[m], label=f"{m} (iPTM)",
                    linewidth=1.8, linestyle="-")
        for m, col, ls in [("UniDock2", "dock_vina_score", "--"),
                           ("GNINA", "dock_cnn_affinity", "--"),
                           ("SurfDock", "surfdock_confidence", ":")]:
            y = [sub[(sub["case"] == c) & (sub["model"] == m)][col].mean()
                 for c in cases]
            ax2.plot(x, y, marker="s", color=MODEL_COLORS[m],
                     label=f"{m} ({col.replace('dock_', '')})",
                     linewidth=1.8, linestyle=ls)
        ax.set_xticks(x); ax.set_xticklabels(cases, rotation=30, ha="right", fontsize=8)
        ax.set_ylabel("ligand-iPTM", color="#444")
        ax2.set_ylabel("docking score", color="#444")
        ax.set_title(f"{fam}: model scores")
        ax.grid(alpha=0.3)
        if i == 0:
            l1, lb1 = ax.get_legend_handles_labels()
            l2, lb2 = ax2.get_legend_handles_labels()
            ax.legend(l1 + l2, lb1 + lb2, fontsize=7, loc="lower left")

    fig.suptitle("Score-per-case: does confidence/score track perturbation?")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "score_per_case.png", dpi=150)
    plt.close()
    print(f"  wrote {FIGURES_DIR / 'score_per_case.png'}")

    # --- 5. Print summary tables -----------------------------------------
    print("\nPer-family median RMSD (common):")
    med = df.groupby(["family", "model"])["ligand_rmsd_common"].median().unstack(fill_value=np.nan)
    med = med.reindex(index=FAMILY_ORDER, columns=MODEL_ORDER)
    print(med.round(2).to_string())

    print("\nPer-family median docking Vina score (lower = better):")
    med_v = df.groupby(["family", "model"])["dock_vina_score"].median().unstack(fill_value=np.nan)
    med_v = med_v.reindex(index=FAMILY_ORDER, columns=MODEL_ORDER)
    print(med_v.round(2).to_string())

    print("\nPer-family median GNINA CNNaffinity (higher = better, pK):")
    med_c = df.groupby(["family", "model"])["dock_cnn_affinity"].median().unstack(fill_value=np.nan)
    med_c = med_c.reindex(index=FAMILY_ORDER, columns=MODEL_ORDER)
    print(med_c.round(2).to_string())


if __name__ == "__main__":
    main()
