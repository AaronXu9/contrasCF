"""Plot Boltz-2 binding-affinity memorization on the ligand_mutagenesis module.

Sibling of `casf_mutagenesis/scripts/10_plot_affinity.py` but specialized
for the ligand-side variant taxonomy. Variants are collapsed into groups for
readability:

  halo    = halo_F_1 + halo_Cl_1 + halo_Br_1
  chrg-   = chrg_neu_methyl + chrg_neu_ethyl + chrg_neu_propyl
  chrg+   = chrg_pos_1 + chrg_pos_2 + chrg_pos_3
  meth    = meth_1..meth_5

Reads `analysis/ligand_mutagenesis/outputs/paired_affinity_ligand.csv` and
produces a 4-panel figure analogous to the casf side:

  (a) WT vs adversarial log[IC50] scatter (y=x = memorized; positive Δ = recognized)
  (b) WT vs adversarial P(binder) scatter (y=x = memorized; negative Δ = recognized)
  (c) Δ log[IC50] histogram per variant group
  (d) Δ P(binder) histogram per variant group

Output: analysis/ligand_mutagenesis/figures/affinity_memorization_ligand.{png,pdf}
"""
from __future__ import annotations
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
sys.path.insert(0, str(REPO_ROOT / "analysis"))

from ligand_mutagenesis.config import OUTPUT_ROOT  # noqa: E402

# Ligand-side variant groups (same convention as 13_plot_overview.py panel b)
LIG_GROUPS = {
    "halo":  ("halo_F_1", "halo_Cl_1", "halo_Br_1"),
    "chrg-": ("chrg_neu_methyl", "chrg_neu_ethyl", "chrg_neu_propyl"),
    "chrg+": ("chrg_pos_1", "chrg_pos_2", "chrg_pos_3"),
    "meth":  ("meth_1", "meth_2", "meth_3", "meth_4", "meth_5"),
}
LIG_GROUP_ORDER = ("halo", "chrg-", "chrg+", "meth")
LIG_COLORS = {
    "halo":  "#4C72B0",
    "chrg-": "#DD8452",
    "chrg+": "#C44E52",
    "meth":  "#55A868",
}
LIG_LABELS = {
    "halo":  "halo (F/Cl/Br)",
    "chrg-": "chrg→neutral",
    "chrg+": "chrg→positive",
    "meth":  "methylation",
}


def variant_to_group(variant: str) -> str | None:
    for g, members in LIG_GROUPS.items():
        if variant in members:
            return g
    return None


def load_paired_affinity() -> dict[str, list[dict]]:
    path = OUTPUT_ROOT / "paired_affinity_ligand.csv"
    by_g: dict[str, list[dict]] = defaultdict(list)
    with path.open() as f:
        for r in csv.DictReader(f):
            if not r["delta_affinity"]:
                continue
            g = variant_to_group(r["variant"])
            if g is None:
                continue
            by_g[g].append({
                "pid": r["pdbid"],
                "wt_aff": float(r["wt_affinity"]),
                "adv_aff": float(r["adv_affinity"]),
                "d_aff": float(r["delta_affinity"]),
                "wt_prob": float(r["wt_probability"]),
                "adv_prob": float(r["adv_probability"]),
                "d_prob": float(r["delta_probability"]),
            })
    return by_g


def main() -> int:
    by_g = load_paired_affinity()
    n_total = sum(len(v) for v in by_g.values())
    print(f"Loaded {n_total} paired-affinity rows across {len(by_g)} variant groups")
    for g in LIG_GROUP_ORDER:
        print(f"  {g}: n={len(by_g[g])}")

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    (ax_aff, ax_prob), (ax_daff, ax_dprob) = axes

    # ---- Panel (a) — WT vs adversarial log[IC50] ----
    all_aff_xy: list[float] = []
    for g in LIG_GROUP_ORDER:
        d = by_g[g]
        if not d:
            continue
        xs = [x["wt_aff"] for x in d]
        ys = [x["adv_aff"] for x in d]
        ax_aff.scatter(xs, ys, s=28, c=LIG_COLORS[g], alpha=0.6,
                       edgecolors="white", linewidth=0.4,
                       label=f"{LIG_LABELS[g]}  (n={len(d)})")
        all_aff_xy.extend(xs); all_aff_xy.extend(ys)
    if all_aff_xy:
        lo, hi = min(all_aff_xy) - 0.3, max(all_aff_xy) + 0.3
        ax_aff.plot([lo, hi], [lo, hi], "k--", lw=0.8, label="y = x (memorization)")
        ax_aff.plot([lo, hi], [lo + 1, hi + 1], "r:", lw=0.8, alpha=0.5,
                    label="y = x + 1 (10× weaker)")
        ax_aff.set_xlim(lo, hi); ax_aff.set_ylim(lo, hi)
    ax_aff.set_xlabel("WT  log[IC50]  (µM)")
    ax_aff.set_ylabel("Adversarial  log[IC50]  (µM)")
    ax_aff.set_title("(a) Predicted affinity: WT vs adversarial (ligand-side)\n"
                     "tighter ←  →  weaker", fontsize=10.5)
    ax_aff.legend(loc="lower right", fontsize=7.5)
    ax_aff.set_aspect("equal", adjustable="box")
    ax_aff.grid(alpha=0.3)

    # ---- Panel (b) — WT vs adversarial P(binder) ----
    for g in LIG_GROUP_ORDER:
        d = by_g[g]
        if not d:
            continue
        xs = [x["wt_prob"] for x in d]
        ys = [x["adv_prob"] for x in d]
        ax_prob.scatter(xs, ys, s=28, c=LIG_COLORS[g], alpha=0.6,
                        edgecolors="white", linewidth=0.4)
    ax_prob.plot([0, 1], [0, 1], "k--", lw=0.8, label="y = x (memorization)")
    ax_prob.axvline(0.5, color="grey", lw=0.5, ls=":", alpha=0.6)
    ax_prob.axhline(0.5, color="grey", lw=0.5, ls=":", alpha=0.6)
    ax_prob.set_xlabel("WT P(binder)")
    ax_prob.set_ylabel("Adversarial P(binder)")
    ax_prob.set_title("(b) Binding probability: WT vs adversarial (ligand-side)\n"
                      "binder ↑  →  non-binder ↓", fontsize=10.5)
    ax_prob.legend(loc="lower right", fontsize=8)
    ax_prob.set_xlim(0, 1); ax_prob.set_ylim(0, 1)
    ax_prob.set_aspect("equal", adjustable="box")
    ax_prob.grid(alpha=0.3)

    # ---- Panel (c) — Δ log[IC50] histogram per group ----
    bins = np.linspace(-3, 6, 28)
    for g in LIG_GROUP_ORDER:
        d = by_g[g]
        if not d:
            continue
        d_affs = [x["d_aff"] for x in d]
        ax_daff.hist(d_affs, bins=bins, alpha=0.55, color=LIG_COLORS[g],
                     label=f"{LIG_LABELS[g]}  (median = {np.median(d_affs):+.2f})",
                     edgecolor="white", linewidth=0.4)
    ax_daff.axvline(0, color="k", lw=0.9, ls="--", label="memorized (Δ = 0)")
    ax_daff.axvline(1, color="r", lw=0.9, ls=":", alpha=0.7, label="10× weaker (Δ = +1)")
    ax_daff.set_xlabel("Δ log[IC50]  (adversarial − WT)")
    ax_daff.set_ylabel("count")
    ax_daff.set_title("(c) Affinity-shift distribution (ligand-side)\n"
                      "← memorized          recognized →", fontsize=10.5)
    ax_daff.legend(loc="upper right", fontsize=7.5)
    ax_daff.grid(alpha=0.3)

    # ---- Panel (d) — Δ P(binder) histogram per group ----
    bins_p = np.linspace(-1, 1, 21)
    for g in LIG_GROUP_ORDER:
        d = by_g[g]
        if not d:
            continue
        d_probs = [x["d_prob"] for x in d]
        ax_dprob.hist(d_probs, bins=bins_p, alpha=0.55, color=LIG_COLORS[g],
                      label=f"{LIG_LABELS[g]}  (median = {np.median(d_probs):+.3f})",
                      edgecolor="white", linewidth=0.4)
    ax_dprob.axvline(0, color="k", lw=0.9, ls="--", label="memorized (Δ = 0)")
    ax_dprob.axvline(-0.3, color="r", lw=0.9, ls=":", alpha=0.7,
                     label="−0.3 (clear drop)")
    ax_dprob.set_xlabel("Δ P(binder)  (adversarial − WT)")
    ax_dprob.set_ylabel("count")
    ax_dprob.set_title("(d) Binding-probability shift distribution (ligand-side)\n"
                       "recognized ←          → memorized", fontsize=10.5)
    ax_dprob.legend(loc="upper left", fontsize=7.5)
    ax_dprob.grid(alpha=0.3)

    fig.suptitle(
        "Boltz-2 binding-affinity memorization on ligand_mutagenesis (full CASF, halo / chrg / meth)",
        fontsize=12, y=0.995,
    )
    fig.tight_layout()

    out_dir = REPO_ROOT / "analysis/ligand_mutagenesis/figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    png_path = out_dir / "affinity_memorization_ligand.png"
    pdf_path = out_dir / "affinity_memorization_ligand.pdf"
    fig.savefig(png_path, dpi=170, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    print(f"\nWrote {png_path}")
    print(f"Wrote {pdf_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
