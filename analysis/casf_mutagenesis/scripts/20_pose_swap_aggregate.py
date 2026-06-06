#!/usr/bin/env python3
"""Aggregate the pose-swap panel (output of 19_pose_swap_affinity.py) into the
per-system + panel-level statistics and figure from the design spec, and emit
the decision-table verdict.

Per system (exit-axis ladder d in {0,2,5,10,20} A):
  rho_spear(d, aff), OLS slope (log units/A), native->far gap = aff(exit+20)-aff(native),
  plus the whole+20 sanity residual (must be ~0).
Panel: median/IQR of rho and gap, fraction with gap > +1, paired Wilcoxon
(native vs exit+20 across systems).

Verdict (spec decision table): flat (|median rho|<0.3 AND median gap<+0.5) =>
functionally pose-insensitive; median rho>=+0.7 AND gap>=+2 => pose-sensitive,
mutation-limited; between => weakly pose-sensitive.

Run (protenix env): .../protenix/bin/python 20_pose_swap_aggregate.py \
    --panel /tmp/poseswap_panel --out /tmp/poseswap_panel
"""
from __future__ import annotations
import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

LADDER = [0.0, 2.0, 5.0, 10.0, 20.0]


def per_system(df):
    lad = df[df.axis == "exit"].sort_values("disp_A")
    d = lad.disp_A.values
    a = lad.boltz_aff.values
    nat = float(lad[lad.disp_A == 0].boltz_aff.iloc[0])
    far = lad[lad.disp_A == 20]
    gap = float(far.boltz_aff.iloc[0]) - nat if len(far) else np.nan
    rho = stats.spearmanr(d, a).correlation if len(set(a)) > 1 else 0.0
    slope = np.polyfit(d, a, 1)[0] if len(d) > 1 else np.nan
    prob_far = df[(df.pose_id == "exit+20")].boltz_prob
    prob_nat = df[(df.pose_id == "native")].boltz_prob
    whole = df[df.pose_id == "whole+20"]
    sanity = float(whole.boltz_aff.iloc[0]) - nat if len(whole) else np.nan
    return dict(
        rho=rho, slope=slope, gap=gap,
        prob_gap=(float(prob_far.iloc[0]) - float(prob_nat.iloc[0])) if len(prob_far) and len(prob_nat) else np.nan,
        native_aff=nat, far_aff=nat + gap, sanity_resid=sanity, n_lig=int(df.n_lig.iloc[0]),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", default="/tmp/poseswap_panel")
    ap.add_argument("--out", default="/tmp/poseswap_panel")
    args = ap.parse_args()

    files = sorted(glob.glob(f"{args.panel}/*/poseswap_*.csv"))
    if not files:
        print(f"no panel CSVs under {args.panel}"); return
    summ = []
    curves = {}
    for f in files:
        df = pd.read_csv(f)
        tag = df.tag.iloc[0]
        s = per_system(df); s["pdbid"] = tag
        summ.append(s)
        lad = df[df.axis == "exit"].sort_values("disp_A")
        curves[tag] = (lad.disp_A.values, lad.boltz_aff.values, lad.boltz_prob.values)
    S = pd.DataFrame(summ).set_index("pdbid")

    # sanity gate
    bad = S[S.sanity_resid.abs() > 0.02]
    if len(bad):
        print(f"WARN: {len(bad)} systems failed whole-complex sanity (|resid|>0.02): {list(bad.index)}")

    print(f"\n=== POSE-SWAP PANEL  (n={len(S)} systems; gap = aff(exit+20) - aff(native), + = weaker) ===")
    print(S[["native_aff", "gap", "rho", "slope", "prob_gap", "sanity_resid"]].round(3).to_string())

    med_rho, med_gap = S.rho.median(), S.gap.median()
    iqr = lambda x: (np.percentile(x, 25), np.percentile(x, 75))
    w = stats.wilcoxon(S.far_aff - S.native_aff, alternative="greater") if len(S) >= 6 else None
    frac_big = float((S.gap > 1.0).mean())
    print("\n--- panel summary ---")
    print(f"  median rho(d,aff)      = {med_rho:+.3f}   IQR {tuple(round(v,3) for v in iqr(S.rho.dropna()))}")
    print(f"  median native->far gap = {med_gap:+.3f}   IQR {tuple(round(v,3) for v in iqr(S.gap.dropna()))} log units")
    print(f"  fraction gap > +1      = {frac_big:.0%}")
    print(f"  paired Wilcoxon native vs exit+20 (greater): p = {w.pvalue:.2e}" if w else "  (n<6, Wilcoxon skipped)")

    if abs(med_rho) < 0.3 and med_gap < 0.5:
        verdict = "FUNCTIONALLY POSE-INSENSITIVE (flat under direct 20A ejection)"
    elif med_rho >= 0.7 and med_gap >= 2.0:
        verdict = "POSE-SENSITIVE, MUTATION-LIMITED (head reads geometry; mutations were the weak lever)"
    else:
        verdict = "WEAKLY POSE-SENSITIVE (partial reading)"
    print(f"\n  >>> VERDICT: {verdict}")

    S.round(4).to_csv(Path(args.out) / "pose_swap_summary.csv")

    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for tag, (d, a, p) in curves.items():
        a0 = a[d == 0][0]
        ax[0].plot(d, a - a0, color="#888", alpha=.5, lw=1)
        ax[1].plot(d, p - p[d == 0][0], color="#888", alpha=.5, lw=1)
    # median curve
    dgrid = np.array(LADDER)
    med_a = np.array([np.median([a[d == dd][0] - a[d == 0][0] for d, a, _ in curves.values() if (d == dd).any() and (d == 0).any()]) for dd in dgrid])
    ax[0].plot(dgrid, med_a, color="#c00", lw=3, marker="o", label="median")
    ax[0].axhline(0, color="k", lw=.6)
    ax[0].axhline(3, color="green", ls=":", lw=1, label="physics floor (+3)")
    ax[0].set_xlabel("ligand displacement along exit axis (A)")
    ax[0].set_ylabel("Δ predicted log[IC50] vs native (+ = weaker)")
    ax[0].set_title(f"Boltz-2 affinity vs ligand ejection (n={len(S)})\nmedian gap@20A = {med_gap:+.2f}  |  {verdict.split('(')[0].strip()}")
    ax[0].legend(fontsize=8)
    ax[1].axhline(0, color="k", lw=.6)
    ax[1].set_xlabel("ligand displacement along exit axis (A)")
    ax[1].set_ylabel("Δ P(binder) vs native")
    ax[1].set_title("P(binder) vs ejection")
    fig.suptitle("Pose-swap test — does Boltz-2's affinity head register a ligand pulled out of the pocket?", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, .96])
    out = Path(args.out) / "pose_swap_affinity.png"
    fig.savefig(out, dpi=130)
    print(f"[wrote] {out}\n[wrote] {Path(args.out)/'pose_swap_summary.csv'}")


if __name__ == "__main__":
    main()
