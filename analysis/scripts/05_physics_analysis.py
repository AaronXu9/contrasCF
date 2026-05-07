#!/usr/bin/env python
"""Physics-aligned analysis of the docking vs co-folding results.

Framing (per Masters et al. 2025):
  - WT cases (`bindingsite_wt`, `glucose_0`): the ligand SHOULD bind in the
    native pocket. Physics-aware model behaviour = low pose RMSD + strong score.
  - ALL other 14 cases are adversarial: the perturbation (pocket mutation,
    altered ATP tail, methylated sugar) is designed to destroy binding. A
    physics-aware model should either (a) place the ligand AWAY from the native
    pocket, or (b) produce a noticeably weaker score than on the WT baseline.
    A model that keeps the ligand at the native position with undiminished
    confidence is MEMORISING, which is the paper's core finding.

Metrics produced:
  1. Memorisation rate: fraction of adversarial cases per model with
     `ligand_rmsd_common` < {2, 4} Å. Lower is more physics-aware.
     (Note: for physics-based docking this is a weaker diagnostic because the
     search box is constrained to the native pocket — docking's physics signal
     lives mainly in the score, next metric.)
  2. Score degradation vs WT baseline: Δscore = score_adv − score_WT for each
     adversarial case, per tool. Positive Δ (worse-than-WT) is physics-aware;
     near-zero Δ is memorisation.
  3. Score vs perturbation magnitude: for families with an ordinal
     perturbation (atp_charge: methyl→propyl→+1→+2→+3; glucose: 0..5 methyls),
     plot score vs ordinal position. Monotonic degradation = physics-aware.

Run:
    LD_LIBRARY_PATH=/home/aoxu/miniconda3/envs/rdkit_env/lib:$LD_LIBRARY_PATH \
        /home/aoxu/miniconda3/envs/rdkit_env/bin/python \
        analysis/scripts/05_physics_analysis.py
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

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from config import RESULTS_DIR, FIGURES_DIR  # noqa: E402


MODEL_ORDER = ["AF3", "Boltz", "Boltz2", "UniDock2", "GNINA", "SurfDock"]
MODEL_COLORS = {
    "AF3":      "#4C72B0",
    "Boltz":    "#C44E52",
    "Boltz2":   "#D62728",
    "UniDock2": "#55A868",
    "GNINA":    "#DD8452",
    "SurfDock": "#8C564B",
}

# WT baseline per family. bindingsite and atp_charge both use CDK2-WT / native
# ATP (only bindingsite_wt is the true WT co-crystal point in our set).
WT_BASELINE = {
    "bindingsite": "bindingsite_wt",
    "atp_charge":  "bindingsite_wt",
    "glucose":     "glucose_0",
}
WT_CASES = {"bindingsite_wt", "glucose_0"}

# Ordinal position within a family, for "score vs perturbation magnitude" plot.
PERTURBATION_ORDER = {
    "atp_charge": [
        "atp_charge_methyl", "atp_charge_ethyl", "atp_charge_propyl",
        "atp_charge_1", "atp_charge_2", "atp_charge_3",
    ],
    "glucose": [f"glucose_{i}" for i in range(6)],
    "bindingsite": ["bindingsite_wt", "bindingsite_rem", "bindingsite_pack", "bindingsite_inv"],
}

# Primary score column per model (higher = more confident binding).
# For Vina scores, we flip sign (−Vina) so all plots follow "bigger = better".
SCORE_COLUMN = {
    "AF3":      "ligand_iptm",         # [0, 1]
    "Boltz":    "ligand_iptm",
    "Boltz2":   "ligand_iptm",
    "UniDock2": "dock_vina_score",     # kcal/mol, more negative = better
    "GNINA":    "dock_cnn_affinity",   # pK, bigger = better
    "SurfDock": "surfdock_confidence", # SurfScore, bigger = better
}
SCORE_SIGN = {   # multiply by this to make "higher = better binding"
    "AF3":      +1,
    "Boltz":    +1,
    "Boltz2":   +1,
    "UniDock2": -1,
    "GNINA":    +1,
    "SurfDock": +1,
}
SCORE_LABEL = {
    "AF3":      "ligand-iPTM",
    "Boltz":    "ligand-iPTM",
    "Boltz2":   "ligand-iPTM",
    "UniDock2": "−Vina score (kcal/mol)",
    "GNINA":    "CNN affinity (pK)",
    "SurfDock": "SurfScore (a.u.)",
}


def _score(df: pd.DataFrame, model: str, case: str) -> float:
    col = SCORE_COLUMN[model]
    sub = df[(df["model"] == model) & (df["case"] == case)]
    if sub.empty:
        return np.nan
    return float(sub.iloc[0][col]) * SCORE_SIGN[model]


def memorisation_rate(df: pd.DataFrame, threshold: float) -> pd.DataFrame:
    """Fraction of ADVERSARIAL cases with ligand_rmsd_common < threshold, per model."""
    adv = df[~df["case"].isin(WT_CASES)].copy()
    adv["near_native"] = adv["ligand_rmsd_common"] < threshold
    return adv.groupby("model")["near_native"].mean().reindex(MODEL_ORDER)


def score_delta_table(df: pd.DataFrame) -> pd.DataFrame:
    """For each adversarial case, Δscore = score − WT_baseline_score, per tool.

    A *positive* Δ means "worse than WT" (for tools where higher = better), i.e.
    physics is registering the perturbation.
    """
    rows = []
    for _, r in df.iterrows():
        if r["case"] in WT_CASES:
            continue
        fam = r["family"]
        baseline_case = WT_BASELINE[fam]
        wt_score = _score(df, r["model"], baseline_case)
        adv_score = _score(df, r["model"], r["case"])
        if np.isnan(wt_score) or np.isnan(adv_score):
            continue
        rows.append({
            "model": r["model"], "family": fam, "case": r["case"],
            "wt_score": wt_score, "adv_score": adv_score,
            "delta": wt_score - adv_score,   # WT is better; positive = perturbation registered
        })
    return pd.DataFrame(rows)


# --- plots -----------------------------------------------------------------


def plot_memorisation_rate(df: pd.DataFrame) -> None:
    """Paper-Fig-3-analog: fraction of adversarial cases with RMSD < threshold."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, thr in zip(axes, (2.0, 4.0)):
        rate = memorisation_rate(df, thr)
        colors = [MODEL_COLORS[m] for m in rate.index]
        ax.bar(rate.index, rate.values, color=colors, edgecolor="white", width=0.6)
        ax.set_ylabel(f"fraction of 14 adversarial cases\nwith RMSD < {thr:g} Å")
        ax.set_ylim(0, 1)
        ax.grid(axis="y", alpha=0.3)
        for i, v in enumerate(rate.values):
            ax.text(i, v + 0.02, f"{v:.2f}", ha="center", fontsize=9)
        ax.set_title(f"Memorisation rate (RMSD < {thr:g} Å)")
    fig.suptitle(
        "Memorisation rate — higher = more pocket-memorising\n"
        "(adversarial cases are designed to break binding; a physics-aware model should NOT retain the native pose)",
        fontsize=10,
    )
    plt.tight_layout()
    out = FIGURES_DIR / "memorisation_rate.png"
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"  wrote {out}")


def plot_score_delta(df: pd.DataFrame) -> None:
    """ΔScore = WT_baseline − adversarial; higher = perturbation registered."""
    d = score_delta_table(df)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), sharey=False)
    for ax, fam in zip(axes, ("bindingsite", "atp_charge", "glucose")):
        sub = d[d["family"] == fam]
        for i, m in enumerate(MODEL_ORDER):
            ys = sub[sub["model"] == m]["delta"].values
            xs = np.full_like(ys, i, dtype=float) + np.random.uniform(-0.12, 0.12, size=len(ys))
            ax.scatter(xs, ys, color=MODEL_COLORS[m], s=55, edgecolor="white", alpha=0.9, label=m if fam == "bindingsite" else None)
            if len(ys):
                ax.plot([i - 0.25, i + 0.25], [ys.mean()] * 2, color="black", linewidth=1.5)
        ax.axhline(0, color="gray", ls="--", alpha=0.6)
        ax.set_xticks(range(len(MODEL_ORDER))); ax.set_xticklabels(MODEL_ORDER, fontsize=9)
        ax.set_title(f"{fam} adversarial cases")
        ax.grid(axis="y", alpha=0.3)
        if fam == "bindingsite":
            ax.set_ylabel("Δ score = WT − adversarial\n(higher = perturbation registered)")
            ax.legend(fontsize=8, loc="upper left")
    fig.suptitle(
        "Physics signal: does the score drop when we move from WT to adversarial?\n"
        "(score column per model: AF3/Boltz ligand-iPTM; UniDock2 −Vina; GNINA CNNaffinity)",
        fontsize=10,
    )
    plt.tight_layout()
    out = FIGURES_DIR / "score_delta_vs_wt.png"
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"  wrote {out}")


def plot_score_vs_perturbation(df: pd.DataFrame) -> None:
    """For ordinal perturbations, score as a function of perturbation index."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for ax, fam in zip(axes, ("atp_charge", "glucose")):
        cases = PERTURBATION_ORDER[fam]
        x = np.arange(len(cases))
        for m in MODEL_ORDER:
            y = [_score(df, m, c) for c in cases]
            ax.plot(x, y, marker="o", color=MODEL_COLORS[m], label=f"{m} — {SCORE_LABEL[m]}", linewidth=1.8)
        ax.set_xticks(x)
        ax.set_xticklabels([c.replace(f"{fam}_", "") for c in cases], rotation=20, ha="right", fontsize=8)
        ax.set_xlabel("perturbation magnitude →")
        ax.set_ylabel("score (higher = better binding)")
        ax.set_title(f"{fam}: score vs perturbation")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc="best")
    fig.suptitle("Does the score monotonically worsen as the ligand is perturbed further? (physics → yes)", fontsize=11)
    plt.tight_layout()
    out = FIGURES_DIR / "score_vs_perturbation.png"
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"  wrote {out}")


def plot_wt_vs_adv_score(df: pd.DataFrame) -> None:
    """Raw WT vs adversarial score per tool — paired view across cases."""
    d = score_delta_table(df)
    n = len(MODEL_ORDER)
    ncols = 3
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, 3.8 * nrows),
                             sharex=False, sharey=False)
    axes_flat = np.array(axes).reshape(-1)
    fam_color = {"bindingsite": "#222222", "atp_charge": "#8888aa", "glucose": "#aa8844"}
    for ax, m in zip(axes_flat, MODEL_ORDER):
        sub = d[d["model"] == m]
        for _, row in sub.iterrows():
            ax.plot([0, 1], [row["wt_score"], row["adv_score"]],
                    marker="o", color=fam_color.get(row["family"], "#888"), alpha=0.55, linewidth=1)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["WT baseline", "adversarial"])
        ax.set_title(f"{m}\n{SCORE_LABEL[m]}")
        ax.grid(axis="y", alpha=0.3)
    # Hide unused axes
    for ax in axes_flat[n:]:
        ax.set_visible(False)
    # Legend on the last visible axis
    handles = [plt.Line2D([0], [0], marker='o', color=c, linewidth=1.5, label=f)
               for f, c in (("bindingsite", "#222222"), ("atp_charge", "#8888aa"), ("glucose", "#aa8844"))]
    axes_flat[n - 1].legend(handles=handles, fontsize=8, loc="best")
    fig.suptitle("Score degradation from WT → adversarial, per tool (each line = one adversarial case)", fontsize=11)
    plt.tight_layout()
    out = FIGURES_DIR / "wt_vs_adv_score.png"
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"  wrote {out}")


def main() -> None:
    df = pd.read_csv(RESULTS_DIR / "results.csv")
    df = df[df["model"].isin(MODEL_ORDER)].copy()
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # --- numeric summaries --------------------------------------------------
    print("Memorisation rate (fraction of 14 adversarial cases with RMSD < 2 Å):")
    print(memorisation_rate(df, 2.0).round(2).to_string(), "\n")
    print("Memorisation rate (RMSD < 4 Å):")
    print(memorisation_rate(df, 4.0).round(2).to_string(), "\n")

    d = score_delta_table(df)
    print("Median Δscore (WT − adv; positive = perturbation registered):")
    piv = d.groupby(["family", "model"])["delta"].median().unstack().reindex(columns=MODEL_ORDER)
    print(piv.round(3).to_string(), "\n")

    # WT baseline scores (sanity)
    print("WT baseline scores (each model on its baseline case):")
    for fam, wt_case in [(f, WT_BASELINE[f]) for f in ("bindingsite", "atp_charge", "glucose")]:
        scores = [f"{m}={_score(df, m, wt_case):.3f}" for m in MODEL_ORDER]
        print(f"  {fam:12s} (baseline={wt_case}): {'  '.join(scores)}")
    print()

    # --- figures ------------------------------------------------------------
    plot_memorisation_rate(df)
    plot_score_delta(df)
    plot_score_vs_perturbation(df)
    plot_wt_vs_adv_score(df)


if __name__ == "__main__":
    main()
