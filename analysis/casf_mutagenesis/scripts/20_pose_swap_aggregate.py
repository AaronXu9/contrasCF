#!/usr/bin/env python3
"""Aggregate the pose-swap panel (output of 19_pose_swap_affinity.py) into
per-system + panel statistics and the figure, and emit the decision verdict.

The ligand is ejected BEYOND the protein bounding sphere (poses native, eject5,
eject15, eject30 = clearance beyond surface; the achieved min ligand-protein
distance is recorded per pose). Decisive pose = eject30 (ligand tens of A from
all protein, every cross-distance past the 22 A distogram saturation).

Per system: rho_spear(min_lig_prot, aff), slope, native->eject30 gap; whole_far
sanity residual (must ~0). Panel: medians/IQR, fraction gap>+1, paired Wilcoxon.
Verdict: |median rho|<0.3 AND median gap<+0.5 => functionally pose-insensitive;
median rho>=+0.7 AND gap>=+2 => pose-sensitive, mutation-limited; else weak.

Run (protenix env): .../protenix/bin/python 20_pose_swap_aggregate.py \
    --panel /tmp/poseswap_panel2 --out /tmp/poseswap_panel2
"""
from __future__ import annotations
import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

CLEARANCE = [0.0, 5.0, 15.0, 30.0]   # nominal clearance beyond protein surface (figure grid)


def per_system(df):
    lad = df[df.axis.isin(["none", "radial"])].sort_values("disp_A")
    sep = lad.min_lig_prot.values                       # actual ligand-protein separation
    a = lad.boltz_aff.values
    nat = float(df[df.pose_id == "native"].boltz_aff.iloc[0])
    far = df[df.pose_id == "eject30"]
    gap = float(far.boltz_aff.iloc[0]) - nat if len(far) else np.nan
    rho = stats.spearmanr(sep, a).correlation if len(set(a)) > 1 else 0.0
    slope = np.polyfit(sep, a, 1)[0] if len(sep) > 1 else np.nan
    pf, pn = df[df.pose_id == "eject30"].boltz_prob, df[df.pose_id == "native"].boltz_prob
    whole = df[df.pose_id == "whole_far"]
    sanity = float(whole.boltz_aff.iloc[0]) - nat if len(whole) else np.nan
    return dict(rho=rho, slope=slope, gap=gap,
                prob_gap=(float(pf.iloc[0]) - float(pn.iloc[0])) if len(pf) and len(pn) else np.nan,
                native_aff=nat, far_aff=nat + gap, far_mindist=float(far.min_lig_prot.iloc[0]) if len(far) else np.nan,
                sanity_resid=sanity, n_lig=int(df.n_lig.iloc[0]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", default="/tmp/poseswap_panel2")
    ap.add_argument("--out", default="/tmp/poseswap_panel2")
    args = ap.parse_args()

    files = sorted(glob.glob(f"{args.panel}/*/poseswap_*.csv"))
    if not files:
        print(f"no panel CSVs under {args.panel}"); return
    summ, curves = [], {}
    for f in files:
        df = pd.read_csv(f)
        tag = df.tag.iloc[0]
        s = per_system(df); s["pdbid"] = tag
        summ.append(s)
        lad = df[df.axis.isin(["none", "radial"])].sort_values("disp_A")
        curves[tag] = (lad.disp_A.values, lad.boltz_aff.values, lad.boltz_prob.values, lad.min_lig_prot.values)
    S = pd.DataFrame(summ).set_index("pdbid")

    bad = S[S.sanity_resid.abs() > 0.02]
    if len(bad):
        print(f"WARN: {len(bad)} systems failed whole-complex sanity (|resid|>0.02): {list(bad.index)}")

    print(f"\n=== POSE-SWAP PANEL  (n={len(S)}; gap = aff(eject30) - aff(native), + = weaker) ===")
    print(S[["native_aff", "far_mindist", "gap", "rho", "slope", "prob_gap", "sanity_resid"]].round(3).to_string())

    med_rho, med_gap = S.rho.median(), S.gap.median()
    iqr = lambda x: tuple(round(float(np.percentile(x, q)), 3) for q in (25, 75))
    w = stats.wilcoxon(S.far_aff - S.native_aff, alternative="greater") if len(S) >= 6 else None
    print("\n--- panel summary ---")
    print(f"  median far min-dist    = {S.far_mindist.median():.0f} A  (ligand fully ejected)")
    print(f"  median rho(sep,aff)    = {med_rho:+.3f}   IQR {iqr(S.rho.dropna())}")
    print(f"  median native->far gap = {med_gap:+.3f}   IQR {iqr(S.gap.dropna())} log units")
    print(f"  fraction gap > +1      = {(S.gap > 1.0).mean():.0%}")
    print(f"  paired Wilcoxon native vs eject30 (greater): p = {w.pvalue:.2e}" if w else "  (n<6, skipped)")

    # GAP MAGNITUDE is decisive: a ligand 30+ A in solvent should weaken affinity
    # by +3..+6 log units. rho only reports the direction of (here negligible) noise.
    if med_gap < 0.5:
        verdict = "FUNCTIONALLY POSE-INSENSITIVE (flat with ligand ejected tens of A into solvent)"
    elif med_gap >= 2.0:
        verdict = "POSE-SENSITIVE, MUTATION-LIMITED (head reads geometry)"
    else:
        verdict = "WEAKLY POSE-SENSITIVE (partial reading)"
    print(f"\n  >>> VERDICT: {verdict}")
    S.round(4).to_csv(Path(args.out) / "pose_swap_summary.csv")

    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    grid = np.array(CLEARANCE)
    for tag, (d, a, p, sep) in curves.items():
        a0 = a[d == 0][0]; p0 = p[d == 0][0]
        ax[0].plot(d, a - a0, color="#888", alpha=.5, lw=1)
        ax[1].plot(d, p - p0, color="#888", alpha=.5, lw=1)
    med = np.array([np.median([a[d == c][0] - a[d == 0][0] for d, a, _, _ in curves.values() if (d == c).any()]) for c in grid])
    ax[0].plot(grid, med, color="#c00", lw=3, marker="o", label="median")
    ax[0].axhline(0, color="k", lw=.6)
    ax[0].axhline(3, color="green", ls=":", lw=1, label="physics floor (+3)")
    ax[0].set_xlabel("ligand ejection — clearance beyond protein surface (Å)\n(native→eject30: ligand ~3 Å → tens of Å from protein)")
    ax[0].set_ylabel("Δ predicted log[IC50] vs native  (+ = weaker)")
    ax[0].set_title(f"Boltz-2 affinity vs ligand ejection (n={len(S)})\nmedian gap = {med_gap:+.2f}  |  {verdict.split('(')[0].strip()}")
    ax[0].legend(fontsize=8)
    ax[1].axhline(0, color="k", lw=.6)
    ax[1].set_xlabel("ligand ejection — clearance beyond protein surface (Å)")
    ax[1].set_ylabel("Δ P(binder) vs native")
    ax[1].set_title("P(binder) vs ejection")
    fig.suptitle("Pose-swap test — does Boltz-2's affinity head notice a ligand ejected into solvent?", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, .95])
    out = Path(args.out) / "pose_swap_affinity.png"
    fig.savefig(out, dpi=130)
    print(f"[wrote] {out}\n[wrote] {Path(args.out)/'pose_swap_summary.csv'}")


if __name__ == "__main__":
    main()
