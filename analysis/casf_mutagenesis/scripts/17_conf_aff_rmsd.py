#!/usr/bin/env python3
"""Joint analysis of (confidence, affinity, RMSD) on the pocket-mutation axis.

Three coupled signals per (system, variant) for Boltz-2:
  RMSD-to-crystal   per-pose  (5 diffusion samples; high = ligand moved = good physics)
  confidence        per-pose  (iptm interface; ptm/complex_plddt global control)
  affinity          per-SYSTEM (one log[IC50] + P(binder) per case, broadcast over poses)

Pose unit = best-of-5 (min ligand RMSD); top-1 (rank-0) + median reported as
robustness checks. Sign conventions:
  dRMSD = adv - wt   (+ = ligand moved away from native; physics-correct)
  dAff  = adv - wt   (+ = predicted weaker binding; physics-correct)  [logIC50, lower=tighter]
  dConf = wt  - adv  (+ = confidence dropped on the adversarial; "aware")

Q1  Why is affinity invariant? memorized-structure vs insensitive-module.
    Decisive evidence is in RESPONDED cells (best-RMSD >= 4 A): if dAff stays ~0
    there too, the affinity module is pose-insensitive; if dAff grows with the
    structural response, affinity just tracks the (memorized) structure.
Q2  Does confidence track the structural response?
    2a between-case dConf ~ dRMSD; 2b within-case pose-level corr(conf, rmsd)
    (WT vs adv: negative within-adv = confident about the WT-like memorized pose);
    2c pooled-pose AUROC near-native(<2) vs displaced(>=4).
Q3  Joint structure, controlling for RMSD: correlation matrix + PARTIAL spearman
    corr(conf, affinity | rmsd) — does confidence carry affinity info beyond pose?

Runs on the protenix env python (pandas/scipy/matplotlib/sklearn).
"""
from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)
ADV = ["rem", "pack", "inv"]
MEM, RESP = 2.0, 4.0          # memorized < 2 A ; responded >= 4 A
AFF_MODELS = ["Boltz2"]       # only Boltz-2 has an affinity head
CONF_MODELS = ["Boltz2", "AF3+MSA"]


# --------------------------------------------------------------------------- #
def boot_ci(x, stat=np.median, n=2000, seed=0, alpha=0.05):
    x = np.asarray([v for v in x if np.isfinite(v)])
    if len(x) < 3:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    bs = [stat(rng.choice(x, len(x), replace=True)) for _ in range(n)]
    return (float(np.percentile(bs, 100 * alpha / 2)),
            float(np.percentile(bs, 100 * (1 - alpha / 2))))


def partial_spearman(a, b, c):
    """corr(a, b | c) on ranks (residualize ranks of a,b on rank of c)."""
    m = np.isfinite(a) & np.isfinite(b) & np.isfinite(c)
    a, b, c = a[m], b[m], c[m]
    if len(a) < 6:
        return np.nan, np.nan, len(a)
    ra, rb, rc = (stats.rankdata(v) for v in (a, b, c))
    X = np.c_[np.ones(len(rc)), rc]
    ea = ra - X @ np.linalg.lstsq(X, ra, rcond=None)[0]
    eb = rb - X @ np.linalg.lstsq(X, rb, rcond=None)[0]
    r, p = stats.pearsonr(ea, eb)
    return float(r), float(p), len(a)


def sp(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 6:
        return np.nan, np.nan, int(m.sum())
    r, p = stats.spearmanr(a[m], b[m])
    return float(r), float(p), int(m.sum())


# --------------------------------------------------------------------------- #
def load():
    df = pd.read_csv(OUT / "results_full.csv")
    df = df[df.status == "ok"].copy()
    for c in ["ligand_rmsd_a", "iptm", "ptm", "complex_plddt", "confidence_score",
              "ligand_iptm", "affinity_pred_value", "affinity_probability_binary", "pose_idx"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def aggregate(df):
    """One row per (model, pdbid, variant): best-of-5, top-1, median summaries."""
    rows = []
    for (model, pid, var), g in df.groupby(["model", "pdbid", "variant"]):
        g = g.dropna(subset=["ligand_rmsd_a"])
        if len(g) == 0:
            continue
        best = g.loc[g.ligand_rmsd_a.idxmin()]
        t = g[g.pose_idx == 0]
        top = t.iloc[0] if len(t) else best
        aff = g.affinity_pred_value.dropna()
        prob = g.affinity_probability_binary.dropna()
        rows.append(dict(
            model=model, pdbid=pid, variant=var, n_pose=len(g),
            rmsd_best=best.ligand_rmsd_a, rmsd_top1=top.ligand_rmsd_a,
            rmsd_med=g.ligand_rmsd_a.median(),
            iptm=best.iptm, ptm=best.ptm, plddt=best.complex_plddt,
            confscore=best.confidence_score,
            iptm_top1=top.iptm, iptm_mean=g.iptm.mean(),
            affinity=aff.iloc[0] if len(aff) else np.nan,
            prob=prob.iloc[0] if len(prob) else np.nan,
        ))
    return pd.DataFrame(rows)


def pair_wt_adv(agg):
    wt = agg[agg.variant == "wt"].set_index(["model", "pdbid"])
    out = []
    for _, r in agg[agg.variant.isin(ADV)].iterrows():
        key = (r.model, r.pdbid)
        if key not in wt.index:
            continue
        w = wt.loc[key]
        if isinstance(w, pd.DataFrame):
            w = w.iloc[0]
        out.append(dict(
            model=r.model, pdbid=r.pdbid, variant=r.variant,
            rmsd_best=r.rmsd_best, rmsd_wt=w.rmsd_best,
            d_rmsd=r.rmsd_best - w.rmsd_best,
            d_iptm=w.iptm - r.iptm, d_ptm=w.ptm - r.ptm, d_plddt=w.plddt - r.plddt,
            iptm_adv=r.iptm, iptm_wt=w.iptm,
            d_aff=r.affinity - w.affinity, aff_adv=r.affinity, aff_wt=w.affinity,
            d_prob=r.prob - w.prob, prob_adv=r.prob, prob_wt=w.prob,
        ))
    p = pd.DataFrame(out)
    p["stratum"] = np.where(p.rmsd_best < MEM, "memorized",
                    np.where(p.rmsd_best >= RESP, "responded", "middle"))
    return p


# --------------------------------------------------------------------------- #
def q1_affinity(pair):
    """Affinity invariance: memorized-structure vs insensitive-module."""
    p = pair[pair.model == "Boltz2"].dropna(subset=["d_aff"])
    print("\n" + "=" * 88)
    print("Q1  AFFINITY INVARIANCE  —  is it memorized structure, or an insensitive module?")
    print("    dAff = adv-wt logIC50 (+ = predicted WEAKER; physics wants +3..+6 on a broken pocket)")
    print("=" * 88)
    rows = []
    for stratum in ["memorized", "middle", "responded"]:
        s = p[p.stratum == stratum]
        if len(s) == 0:
            continue
        lo, hi = boot_ci(s.d_aff.values)
        plo, phi = boot_ci(s.d_prob.values)
        rows.append(dict(stratum=stratum, n=len(s),
                         med_rmsd_adv=round(s.rmsd_best.median(), 2),
                         med_dAff=round(s.d_aff.median(), 3), dAff_ci=f"[{lo:+.2f},{hi:+.2f}]",
                         med_dProb=round(s.d_prob.median(), 3), dProb_ci=f"[{plo:+.2f},{phi:+.2f}]"))
    rep = pd.DataFrame(rows)
    print(rep.to_string(index=False))
    r_ar, p_ar, n = sp(p.rmsd_best.values, p.d_aff.values)
    r_dr, p_dr, _ = sp(p.d_rmsd.values, p.d_aff.values)
    print(f"\n  Spearman( dAff , adv RMSD )  = {r_ar:+.3f}  (p={p_ar:.1e}, n={n})")
    print(f"  Spearman( dAff , dRMSD )     = {r_dr:+.3f}  (p={p_dr:.1e})")
    print("  READ: responded-stratum dAff ~ 0  AND slope ~ 0  =>  affinity module is")
    print("        POSE-INSENSITIVE (memorization doesn't explain it). dAff growing with")
    print("        RMSD => affinity merely tracks the (memorized) structure.")
    rep.to_csv(OUT / "q1_affinity_strata.csv", index=False)
    return p


def q2_confidence(pair, df, agg):
    print("\n" + "=" * 88)
    print("Q2  DOES CONFIDENCE TRACK THE STRUCTURAL RESPONSE?  (mutate -> dRMSD -> dConf)")
    print("=" * 88)
    # 2a between-case
    print("\n[2a] between-case  dConf(iptm) ~ dRMSD   and stratified median drop")
    for model in CONF_MODELS:
        p = pair[pair.model == model].dropna(subset=["d_iptm", "d_rmsd"])
        if len(p) < 6:
            continue
        r, pv, n = sp(p.d_rmsd.values, p.d_iptm.values)
        line = f"  {model:8s} n={n:3d}  Spearman(dRMSD,dConf)={r:+.3f} (p={pv:.1e}) | medΔiptm: "
        line += "  ".join(f"{s}={p[p.stratum==s].d_iptm.median():+.3f}(n{len(p[p.stratum==s])})"
                          for s in ["memorized", "middle", "responded"] if len(p[p.stratum == s]))
        print(line)
    # 2b within-case pose-level corr(conf, rmsd)
    print("\n[2b] within-case pose-level Spearman(iptm, rmsd) across the 5 samples")
    print("     (WT: negative = confidence is a valid quality signal;")
    print("      ADV: negative = MORE confident about the WT-like/memorized pose)")
    wc = []
    for (model, pid, var), g in df.groupby(["model", "pdbid", "variant"]):
        g = g.dropna(subset=["ligand_rmsd_a", "iptm"])
        if len(g) < 3 or g.ligand_rmsd_a.nunique() < 3:
            continue
        rho, _ = stats.spearmanr(g.iptm, g.ligand_rmsd_a)
        wc.append(dict(model=model, variant=var,
                       cls="wt" if var == "wt" else "adv", rho=rho))
    wc = pd.DataFrame(wc)
    for model in CONF_MODELS:
        for cls in ["wt", "adv"]:
            s = wc[(wc.model == model) & (wc.cls == cls)].rho.dropna()
            if len(s) < 6:
                continue
            w = stats.wilcoxon(s, alternative="less") if (s != 0).any() else None
            print(f"  {model:8s} {cls:3s}  n={len(s):3d}  median ρ={s.median():+.3f}  "
                  f"frac ρ<0={np.mean(s < 0):.0%}  (Wilcoxon ρ<0 p={w.pvalue:.1e})")
    wc.to_csv(OUT / "q2_within_case_rho.csv", index=False)
    # 2c pooled-pose AUROC: is confidence higher on near-native poses?
    print("\n[2c] pooled-pose AUROC: confidence separates near-native(<2) vs displaced(>=4) poses")
    for model in CONF_MODELS:
        for scope, vmask in [("WT", df.variant == "wt"), ("ADV", df.variant.isin(ADV))]:
            g = df[(df.model == model) & vmask].dropna(subset=["ligand_rmsd_a", "iptm"])
            near = g[g.ligand_rmsd_a < MEM].iptm
            far = g[g.ligand_rmsd_a >= RESP].iptm
            if len(near) < 10 or len(far) < 10:
                continue
            y = np.r_[np.ones(len(near)), np.zeros(len(far))]
            auc = roc_auc_score(y, np.r_[near, far])
            print(f"  {model:8s} {scope:3s}  near n={len(near):4d}  far n={len(far):4d}  "
                  f"AUROC(conf higher on near-native)={auc:.3f}")
    return wc


def q3_joint(pair, agg):
    print("\n" + "=" * 88)
    print("Q3  JOINT STRUCTURE (Boltz-2, adversarial cells, best-of-5) — controlling for RMSD")
    print("=" * 88)
    a = agg[(agg.model == "Boltz2") & (agg.variant.isin(ADV))].copy()
    cols = {"rmsd_best": "RMSD", "iptm": "iptm", "plddt": "plddt",
            "affinity": "affinity", "prob": "P(bind)"}
    M = a[list(cols)].rename(columns=cols)
    corr = M.corr(method="spearman")
    print("\nSpearman correlation matrix:")
    print(corr.round(3).to_string())
    # partial: conf vs affinity controlling for rmsd
    print("\nPartial correlations (controlling for RMSD):")
    for cname, col in [("iptm", "iptm"), ("plddt", "plddt"), ("confscore", "confscore")]:
        if col not in a:
            continue
        r_raw, p_raw, n = sp(a[col].values, a.affinity.values)
        r_par, p_par, _ = partial_spearman(a[col].values, a.affinity.values, a.rmsd_best.values)
        print(f"  {cname:9s} ↔ affinity : raw ρ={r_raw:+.3f} (p={p_raw:.1e})  |  "
              f"partial ρ|RMSD={r_par:+.3f} (p={p_par:.1e})  n={n}")
    # deltas
    pp = pair[pair.model == "Boltz2"]
    r_par, p_par, n = partial_spearman(pp.d_iptm.values, pp.d_aff.values, pp.d_rmsd.values)
    r_raw, p_raw, _ = sp(pp.d_iptm.values, pp.d_aff.values)
    print(f"\n  Δiptm ↔ Δaff : raw ρ={r_raw:+.3f} (p={p_raw:.1e})  |  "
          f"partial ρ|ΔRMSD={r_par:+.3f} (p={p_par:.1e})  n={n}")
    print("  READ: if conf↔affinity collapses to ~0 after controlling for RMSD, the two heads")
    print("        share no information beyond the pose; if it survives, confidence carries")
    print("        affinity-relevant signal independent of where the ligand landed.")
    corr.to_csv(OUT / "q3_corr_matrix.csv")
    return a, corr


def make_figure(pair, wc, agg, corr):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    b = pair[pair.model == "Boltz2"]
    fig, ax = plt.subplots(2, 2, figsize=(12, 10))

    # (a) Q1: dAff vs adv RMSD
    cmap = {"memorized": "#2b7", "middle": "#fb3", "responded": "#e34"}
    for s, c in cmap.items():
        d = b[b.stratum == s]
        ax[0, 0].scatter(d.rmsd_best, d.d_aff, s=14, c=c, alpha=.55, label=f"{s} (n={len(d)})")
    ax[0, 0].axhline(0, color="k", lw=.7); ax[0, 0].axvline(2, color="gray", ls=":"); ax[0, 0].axvline(4, color="gray", ls=":")
    ax[0, 0].set_xlabel("adversarial best-of-5 ligand RMSD (Å)"); ax[0, 0].set_ylabel("Δaffinity (adv−wt logIC50, + = weaker)")
    ax[0, 0].set_title("Q1  affinity vs structural response"); ax[0, 0].legend(fontsize=8); ax[0, 0].set_xlim(0, 15)

    # (b) Q2a: dConf vs dRMSD
    ax[0, 1].scatter(b.d_rmsd, b.d_iptm, s=14, c="#36c", alpha=.5)
    r, _, n = sp(b.d_rmsd.values, b.d_iptm.values)
    ax[0, 1].axhline(0, color="k", lw=.7); ax[0, 1].axvline(0, color="k", lw=.7)
    ax[0, 1].set_xlabel("ΔRMSD (adv−wt, Å)"); ax[0, 1].set_ylabel("Δconfidence (wt−adv iptm, + = dropped)")
    ax[0, 1].set_title(f"Q2a  confidence tracks response  (ρ={r:+.2f}, n={n})")

    # (c) Q2b: within-case rho distribution WT vs ADV (Boltz2)
    for cls, c in [("wt", "#888"), ("adv", "#e34")]:
        s = wc[(wc.model == "Boltz2") & (wc.cls == cls)].rho.dropna()
        ax[1, 0].hist(s, bins=20, alpha=.6, color=c, label=f"{cls} (med {s.median():+.2f}, n={len(s)})")
    ax[1, 0].axvline(0, color="k", lw=.8)
    ax[1, 0].set_xlabel("within-case Spearman(iptm, RMSD)"); ax[1, 0].set_ylabel("# cases")
    ax[1, 0].set_title("Q2b  pose-level: confident about which pose?"); ax[1, 0].legend(fontsize=8)

    # (d) Q3: correlation heatmap
    im = ax[1, 1].imshow(corr.values, vmin=-1, vmax=1, cmap="RdBu_r")
    ax[1, 1].set_xticks(range(len(corr))); ax[1, 1].set_yticks(range(len(corr)))
    ax[1, 1].set_xticklabels(corr.columns, rotation=45, ha="right"); ax[1, 1].set_yticklabels(corr.index)
    for i in range(len(corr)):
        for j in range(len(corr)):
            ax[1, 1].text(j, i, f"{corr.values[i, j]:.2f}", ha="center", va="center", fontsize=8)
    ax[1, 1].set_title("Q3  joint structure (Boltz-2 adversarial)"); fig.colorbar(im, ax=ax[1, 1], fraction=.046)

    fig.suptitle("Confidence × Affinity × RMSD — CASF pocket mutagenesis (Boltz-2)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, .98])
    out = FIG / "conf_aff_rmsd_pocket.png"
    fig.savefig(out, dpi=130)
    print(f"\n[wrote figure] {out}")


def main():
    df = load()
    df = df[df.model.isin(CONF_MODELS)].copy()
    agg = aggregate(df)
    pair = pair_wt_adv(agg)
    pair.to_csv(OUT / "paired_conf_aff_rmsd.csv", index=False)
    q1_affinity(pair)
    wc = q2_confidence(pair, df, agg)
    a, corr = q3_joint(pair, agg)
    make_figure(pair, wc, agg, corr)
    print(f"\n[wrote] paired_conf_aff_rmsd.csv  ({len(pair)} adversarial cells)")


if __name__ == "__main__":
    main()
