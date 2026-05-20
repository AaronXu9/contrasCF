"""Plot Boltz-2 binding-affinity memorization on the CASF-mutagenesis subset.

Reads paired_affinity_<scope>.csv (default: subset20) and produces a 4-panel
figure:
  (a) WT vs adversarial log[IC50]   scatter (y=x diagonal; off-diagonal = recognized)
  (b) WT vs adversarial P(binder)   scatter (y=x; off-diagonal = recognized)
  (c) Δ log[IC50] (adv − wt)        histogram per variant
  (d) Δ P(binder) (adv − wt)        histogram per variant

If Boltz-2 is "doing the physics", we expect:
  panel (a): adversarial points HIGHER than WT (weaker IC50 in log µM units)
  panel (b): adversarial points LOWER than WT (no longer classified as binder)
  panel (c): right-shifted positive distribution (≥+1 to +3 log)
  panel (d): left-shifted negative distribution (≤-0.5)

Saves: analysis/casf_mutagenesis/figures/affinity_memorization_<scope>.{png,pdf}

Run:
    LD_LIBRARY_PATH=/home/aoxu/miniconda3/envs/rdkit_env/lib:$LD_LIBRARY_PATH \\
        /home/aoxu/miniconda3/envs/rdkit_env/bin/python \\
        analysis/casf_mutagenesis/scripts/10_plot_affinity.py [--scope subset20]
"""
from __future__ import annotations
import argparse
import csv
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
sys.path.insert(0, str(REPO_ROOT / "analysis"))
from casf_mutagenesis.config import OUTPUT_ROOT  # noqa: E402

VARIANT_COLORS = {"rem": "#4C72B0", "pack": "#DD8452", "inv": "#55A868"}
VARIANT_LABELS = {"rem": "rem (→G)", "pack": "pack (→F)", "inv": "inv (Miyata)"}


def load_paired_affinity(scope: str) -> dict[str, list[dict]]:
    path = OUTPUT_ROOT / f"paired_affinity_{scope}.csv"
    by_v: dict[str, list[dict]] = {"rem": [], "pack": [], "inv": []}
    with path.open() as f:
        for r in csv.DictReader(f):
            if not r["delta_affinity"]:
                continue
            by_v[r["variant"]].append({
                "pid": r["pdbid"],
                "wt_aff": float(r["wt_affinity"]),
                "adv_aff": float(r["adv_affinity"]),
                "d_aff": float(r["delta_affinity"]),
                "wt_prob": float(r["wt_probability"]),
                "adv_prob": float(r["adv_probability"]),
                "d_prob": float(r["delta_probability"]),
            })
    return by_v


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--scope", default="subset20",
                   help="CSV scope suffix (default subset20)")
    args = p.parse_args()

    by_v = load_paired_affinity(args.scope)
    n_total = sum(len(v) for v in by_v.values())
    print(f"Loaded {n_total} paired-affinity rows across {len(by_v)} variants")

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    (ax_aff, ax_prob), (ax_daff, ax_dprob) = axes

    # ---- Panel (a): WT vs adversarial log[IC50] ----
    all_aff_xy: list[float] = []
    for v in ("rem", "pack", "inv"):
        d = by_v[v]
        if not d:
            continue
        xs = [x["wt_aff"] for x in d]
        ys = [x["adv_aff"] for x in d]
        ax_aff.scatter(xs, ys, s=42, c=VARIANT_COLORS[v], alpha=0.75,
                       edgecolors="white", linewidth=0.5,
                       label=f"{VARIANT_LABELS[v]} (n={len(d)})")
        all_aff_xy.extend(xs); all_aff_xy.extend(ys)
    lo, hi = min(all_aff_xy) - 0.3, max(all_aff_xy) + 0.3
    ax_aff.plot([lo, hi], [lo, hi], "k--", lw=0.8, label="y = x (memorization)")
    # +1 log unit (10x weaker) band
    ax_aff.plot([lo, hi], [lo + 1, hi + 1], "r:", lw=0.8, alpha=0.5,
                label="y = x + 1 (10× weaker)")
    ax_aff.set_xlabel("WT  log[IC50]  (µM)")
    ax_aff.set_ylabel("Adversarial  log[IC50]  (µM)")
    ax_aff.set_title("(a) Predicted affinity: WT vs adversarial\n"
                     "tighter ←  →  weaker", fontsize=11)
    ax_aff.legend(loc="lower right", fontsize=8)
    ax_aff.set_xlim(lo, hi); ax_aff.set_ylim(lo, hi)
    ax_aff.set_aspect("equal", adjustable="box")
    ax_aff.grid(alpha=0.3)

    # ---- Panel (b): WT vs adversarial P(binder) ----
    for v in ("rem", "pack", "inv"):
        d = by_v[v]
        if not d:
            continue
        xs = [x["wt_prob"] for x in d]
        ys = [x["adv_prob"] for x in d]
        ax_prob.scatter(xs, ys, s=42, c=VARIANT_COLORS[v], alpha=0.75,
                        edgecolors="white", linewidth=0.5)
    ax_prob.plot([0, 1], [0, 1], "k--", lw=0.8, label="y = x (memorization)")
    ax_prob.axvline(0.5, color="grey", lw=0.5, ls=":", alpha=0.6)
    ax_prob.axhline(0.5, color="grey", lw=0.5, ls=":", alpha=0.6)
    ax_prob.set_xlabel("WT P(binder)")
    ax_prob.set_ylabel("Adversarial P(binder)")
    ax_prob.set_title("(b) Binding probability: WT vs adversarial\n"
                      "binder ↑  →  non-binder ↓", fontsize=11)
    ax_prob.legend(loc="lower right", fontsize=8)
    ax_prob.set_xlim(0, 1); ax_prob.set_ylim(0, 1)
    ax_prob.set_aspect("equal", adjustable="box")
    ax_prob.grid(alpha=0.3)

    # ---- Panel (c): Δ log[IC50] histogram per variant ----
    bins = np.linspace(-3, 6, 19)
    for v in ("rem", "pack", "inv"):
        d = by_v[v]
        if not d:
            continue
        d_affs = [x["d_aff"] for x in d]
        ax_daff.hist(d_affs, bins=bins, alpha=0.6, color=VARIANT_COLORS[v],
                     label=f"{VARIANT_LABELS[v]}  (median = {np.median(d_affs):+.2f})",
                     edgecolor="white", linewidth=0.5)
    ax_daff.axvline(0, color="k", lw=0.8, ls="--", label="memorization (Δ = 0)")
    ax_daff.axvline(1, color="r", lw=0.8, ls=":", alpha=0.7,
                    label="10× weaker (Δ = +1)")
    ax_daff.set_xlabel("Δ log[IC50]  (adversarial − WT)")
    ax_daff.set_ylabel("count")
    ax_daff.set_title("(c) Affinity-shift distribution\n"
                      "← memorized          recognized →", fontsize=11)
    ax_daff.legend(loc="upper right", fontsize=8)
    ax_daff.grid(alpha=0.3)

    # ---- Panel (d): Δ P(binder) histogram per variant ----
    bins_p = np.linspace(-1, 1, 21)
    for v in ("rem", "pack", "inv"):
        d = by_v[v]
        if not d:
            continue
        d_probs = [x["d_prob"] for x in d]
        ax_dprob.hist(d_probs, bins=bins_p, alpha=0.6, color=VARIANT_COLORS[v],
                      label=f"{VARIANT_LABELS[v]}  (median = {np.median(d_probs):+.3f})",
                      edgecolor="white", linewidth=0.5)
    ax_dprob.axvline(0, color="k", lw=0.8, ls="--", label="memorization (Δ = 0)")
    ax_dprob.axvline(-0.3, color="r", lw=0.8, ls=":", alpha=0.7,
                     label="−0.3 (clear drop)")
    ax_dprob.set_xlabel("Δ P(binder)  (adversarial − WT)")
    ax_dprob.set_ylabel("count")
    ax_dprob.set_title("(d) Binding-probability shift distribution\n"
                       "recognized ←          → memorized", fontsize=11)
    ax_dprob.legend(loc="upper left", fontsize=8)
    ax_dprob.grid(alpha=0.3)

    fig.suptitle(
        f"Boltz-2 binding-affinity memorization on CASF-mutagenesis ({args.scope})",
        fontsize=13, y=0.995,
    )
    fig.tight_layout()

    out_dir = REPO_ROOT / "analysis/casf_mutagenesis/figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    png_path = out_dir / f"affinity_memorization_{args.scope}.png"
    pdf_path = out_dir / f"affinity_memorization_{args.scope}.pdf"
    fig.savefig(png_path, dpi=170, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    print(f"\nWrote {png_path}")
    print(f"Wrote {pdf_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
