"""Compare the three ligand-RMSD metrics on CASF-285.

Reads results_full.csv (rank-0 only) and reports, per (model, variant):
  - median, P25, P75 of each metric
  - memorization rates at 2 Å under each metric
A scatter plot of common vs bestfit (and common vs fullCa) is emitted to
figures/rmsd_metric_compare_full.{pdf,png}, with the y=x line and a histogram
of (common − bestfit) to make the gap visible.

Run:
    LD_LIBRARY_PATH=/home/aoxu/miniconda3/envs/rdkit_env/lib:$LD_LIBRARY_PATH \
        /home/aoxu/miniconda3/envs/rdkit_env/bin/python \
        analysis/casf_mutagenesis/scripts/15_compare_rmsd_metrics.py
"""
from __future__ import annotations
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
CSV = REPO_ROOT / "analysis/casf_mutagenesis/outputs/results_full.csv"
OUT_DIR = REPO_ROOT / "analysis/casf_mutagenesis/figures"
OUT_TXT = REPO_ROOT / "analysis/casf_mutagenesis/outputs/rmsd_metric_compare_full.txt"

VARIANTS = ["wt", "rem", "pack", "inv"]
MODELS = ["Boltz2", "AF3", "AF3+MSA"]
METRICS = ["ligand_rmsd_a", "ligand_rmsd_fullca_a", "bestfit_rmsd_a"]
LABELS  = {
    "ligand_rmsd_a":        "common (chain-near)",
    "ligand_rmsd_fullca_a": "fullCa (all chains)",
    "bestfit_rmsd_a":       "bestfit (paper-style)",
}


def main() -> int:
    df = pd.read_csv(CSV)
    rank0 = df[(df.pose_idx == 0) & (df.status == "ok")].copy()

    lines: list[str] = []
    def out(s: str) -> None:
        print(s); lines.append(s)

    out(f"CASF-285 rank-0 ok cells: {len(rank0)}")
    out("")
    out("Median ligand RMSD by metric (Å):")
    out(f"  {'model':<8s} {'variant':<6s} {'common':>8s} {'fullCa':>8s} {'bestfit':>8s} {'gap':>8s}  n")
    for model in MODELS:
        for variant in VARIANTS:
            sub = rank0[(rank0.model == model) & (rank0.variant == variant)]
            n = sub[METRICS].dropna().shape[0]
            if n == 0:
                continue
            medians = {m: sub[m].dropna().median() for m in METRICS}
            gap = medians["ligand_rmsd_a"] - medians["bestfit_rmsd_a"]
            out(f"  {model:<8s} {variant:<6s} "
                f"{medians['ligand_rmsd_a']:>8.2f} {medians['ligand_rmsd_fullca_a']:>8.2f} "
                f"{medians['bestfit_rmsd_a']:>8.2f} {gap:>8.2f}  {n}")

    out("")
    out("Memorization rate at <2 Å under each metric:")
    out(f"  {'model':<8s} {'variant':<6s} {'common':>8s} {'fullCa':>8s} {'bestfit':>8s}  n")
    for model in MODELS:
        for variant in VARIANTS:
            sub = rank0[(rank0.model == model) & (rank0.variant == variant)]
            n = sub[METRICS].dropna().shape[0]
            if n == 0:
                continue
            rates = {m: float((sub[m].dropna() < 2.0).mean()) for m in METRICS}
            out(f"  {model:<8s} {variant:<6s} "
                f"{rates['ligand_rmsd_a']:>8.2%} {rates['ligand_rmsd_fullca_a']:>8.2%} "
                f"{rates['bestfit_rmsd_a']:>8.2%}  {n}")

    # -------------------------------------------------------------------------
    # AF3+MSA-only scatter: common vs bestfit (focus on user's question)
    # -------------------------------------------------------------------------
    af3msa = rank0[rank0.model == "AF3+MSA"].dropna(
        subset=["ligand_rmsd_a", "bestfit_rmsd_a", "ligand_rmsd_fullca_a"]
    )

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.5))

    palette = {"wt": "#2E7D32", "rem": "#F9A825", "pack": "#C62828", "inv": "#6A1B9A"}

    # Panel 1: common vs bestfit (log-log so dynamic range fits)
    ax = axes[0]
    for v in VARIANTS:
        s = af3msa[af3msa.variant == v]
        ax.scatter(s.bestfit_rmsd_a, s.ligand_rmsd_a, s=14, alpha=0.5,
                   label=f"{v} (n={len(s)})", color=palette[v])
    lim = max(af3msa.ligand_rmsd_a.max(), af3msa.bestfit_rmsd_a.max())
    ax.plot([0.05, lim], [0.05, lim], "k--", lw=0.8, label="y=x")
    ax.axvline(2.0, color="grey", lw=0.5)
    ax.axhline(2.0, color="grey", lw=0.5)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("bestfit RMSD (Å) — paper-style, pocket-blind")
    ax.set_ylabel("common RMSD (Å) — protein-aligned")
    ax.set_title("AF3+MSA: common vs bestfit\n(points above y=x: pocket misplaced but conformation OK)")
    ax.legend(fontsize=8, loc="lower right")

    # Panel 2: common vs fullCa — for CASF (monomer predictions) these should
    # collapse onto y=x — sanity check that the multi-chain framework picks the
    # same alignment as the chain-near framework on single-chain inputs.
    ax = axes[1]
    for v in VARIANTS:
        s = af3msa[af3msa.variant == v]
        ax.scatter(s.ligand_rmsd_fullca_a, s.ligand_rmsd_a, s=14, alpha=0.5,
                   label=f"{v} (n={len(s)})", color=palette[v])
    lim2 = max(af3msa.ligand_rmsd_a.max(), af3msa.ligand_rmsd_fullca_a.max())
    ax.plot([0.05, lim2], [0.05, lim2], "k--", lw=0.8, label="y=x")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("fullCa RMSD (Å) — all-chain alignment")
    ax.set_ylabel("common RMSD (Å) — chain-near alignment")
    ax.set_title("AF3+MSA: common vs fullCa\n(should collapse — CASF predictions are monomers)")
    ax.legend(fontsize=8, loc="lower right")

    # Panel 3: histogram of (common - bestfit) by variant — the "metric gap"
    ax = axes[2]
    bins = np.linspace(-2, 30, 33)
    for v in VARIANTS:
        s = af3msa[af3msa.variant == v]
        gap = (s.ligand_rmsd_a - s.bestfit_rmsd_a).values
        ax.hist(gap, bins=bins, histtype="step", lw=1.5,
                label=f"{v} median Δ={np.median(gap):.2f} Å",
                color=palette[v])
    ax.set_xlabel("common − bestfit (Å)")
    ax.set_ylabel("# CASF systems")
    ax.set_title("AF3+MSA: metric gap distribution")
    ax.axvline(0, color="k", lw=0.5)
    ax.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "rmsd_metric_compare_full.pdf")
    fig.savefig(OUT_DIR / "rmsd_metric_compare_full.png", dpi=150)
    out("")
    out(f"Scatter written: {OUT_DIR/'rmsd_metric_compare_full.pdf'}")

    OUT_TXT.write_text("\n".join(lines) + "\n")
    out(f"Numeric report: {OUT_TXT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
