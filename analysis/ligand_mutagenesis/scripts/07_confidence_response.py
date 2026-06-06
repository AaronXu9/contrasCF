#!/usr/bin/env python3
"""Confidence x affinity x RMSD on the LIGAND-mutation axis (Boltz-2 only).

Mirror of the pocket-axis pair:
  casf_mutagenesis/scripts/16_confidence_response.py  (confidence-drop test)
  casf_mutagenesis/scripts/17_conf_aff_rmsd.py        (joint conf x aff x rmsd)
folded into one script, adapted to the ligand axis.

The ligand axis is the cross-axis counterpoint to the pocket axis. There the
mutation is to the POCKET (protein), here it is to the LIGAND (protein constant):
halogen swaps, charge neutralisation (+alkyl), charge addition (+ammonium), and
serial methylation. Variants are DYNAMIC per system and WT is per system. There
is no AF3 on this axis and no paired_full.csv structure-label file, so
memorized/responded is derived from results_ligand.csv ligand_rmsd_a directly.

The documented ligand-side dichotomy (docs/casf_overview.md, casf_affinity.md):
the STRUCTURE head refuses the adversarial ligands (~0-3% < 2 A vs ~69% WT) but
the AFFINITY head says they bind *tighter* (negative median Daff, e.g. -0.46 for
charge-neutralisation). This script adds the third head: does the model's
self-reported CONFIDENCE register the ligand mutation, and does it track the
ligand-RMSD response, the way it did on the pocket axis?

Signals per (system, variant) for Boltz-2:
  ligand RMSD-to-crystal   per-pose  (5 diffusion samples; high = ligand moved)
  confidence (iptm/ligand_iptm interface; ptm/complex_plddt global control) per-pose
  affinity  (log[IC50], P(binder))  per-SYSTEM (one value, broadcast over poses)

Pose unit = best-of-5 (min ligand RMSD); top-1 (rank-0) is a robustness check.
Sign conventions (adversarial vs that system's own WT):
  dRMSD = RMSD(adv) - RMSD(wt)        (+ = ligand moved away from native)
  dConf = iptm(wt)  - iptm(adv)       (+ = confidence dropped on the adversarial)
  dAff  = aff(adv)  - aff(wt)         (+ = predicted WEAKER; logIC50 lower = tighter)

Adversarial variants are grouped both individually and by FAMILY:
  halo (F/Cl/Br swap), chrg_neu (+alkyl), chrg_pos (+ammonium), meth (1-5 methyls).

What it computes (mirrors 16/17, adapted):
  1. CONFIDENCE DROP WT->variant/family: median Diptm, %dropped, one-sided
     Wilcoxon, AUROC(WT vs adv), bootstrap CI, interface(iptm/ligand_iptm) vs
     global(ptm/complex_plddt) control, conditional on structural response.
  2. CONFIDENCE TRACKS RESPONSE: Spearman(dRMSD, dConf) overall + stratified;
     within-case Spearman(iptm, rmsd) across the 5 poses; pooled-pose AUROC.
  3. AFFINITY INVARIANCE / DICHOTOMY: dAff stratified by structural response,
     with the regression-to-mean control (raw vs partial(dAff, displacement |
     WT affinity)); connect to the "affinity says tighter" dichotomy.
  4. JOINT STRUCTURE: Spearman matrix over {RMSD, iptm, affinity, P(bind)} and
     the decoupling test diptm <-> daff (raw + partial | dRMSD).
  5. A 4-panel figure mirroring 17.

CAVEAT (read before over-claiming): on the ligand axis the "memorized" stratum
(best-of-5 < 2 A) is nearly empty (the structure head refuses these ligands), so
the pocket-axis conditional-on-memorization test is n-starved here; the action
is in the RESPONDED stratum. Several adversarial families also have small n
(chrg ~42 systems, meth_4/5 single digits). n is reported everywhere; p-values
are treated cautiously at large n (effect size / AUROC is the content).

Runs on the protenix env python (pandas/scipy/matplotlib/sklearn):
  /mnt/katritch_lab2/aoxu/envs/protenix/bin/python3
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

MEM, RESP = 2.0, 4.0          # memorized < 2 A ; responded >= 4 A
# Confidence metrics: interface should drop; global is the negative control.
INTERFACE = ["iptm", "ligand_iptm"]
GLOBAL = ["ptm", "complex_plddt"]
COMPOSITE = ["confidence_score", "affinity_probability_binary"]
METRICS = INTERFACE + GLOBAL + COMPOSITE

# Adversarial families (DYNAMIC variant names mapped to a family).
FAMILIES = ["halo", "chrg_neu", "chrg_pos", "meth"]
# Individual variant order for the per-variant table.
VARIANT_ORDER = [
    "halo_F_1", "halo_Cl_1", "halo_Br_1",
    "chrg_neu_methyl", "chrg_neu_ethyl", "chrg_neu_propyl",
    "chrg_pos_1", "chrg_pos_2", "chrg_pos_3",
    "meth_1", "meth_2", "meth_3", "meth_4", "meth_5",
]


def family_of(variant: str) -> str:
    if variant == "wt":
        return "wt"
    if variant.startswith("halo"):
        return "halo"
    if variant.startswith("chrg_neu"):
        return "chrg_neu"
    if variant.startswith("chrg_pos"):
        return "chrg_pos"
    if variant.startswith("meth"):
        return "meth"
    return "other"


# --------------------------------------------------------------------------- #
# Small stats helpers (same conventions as scripts 16/17)                     #
# --------------------------------------------------------------------------- #
def boot_ci(x, stat=np.median, n=2000, seed=0, alpha=0.05):
    x = np.asarray([v for v in x if np.isfinite(v)], dtype=float)
    if len(x) < 3:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    bs = [stat(rng.choice(x, len(x), replace=True)) for _ in range(n)]
    return (float(np.percentile(bs, 100 * alpha / 2)),
            float(np.percentile(bs, 100 * (1 - alpha / 2))))


def sp(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 6:
        return np.nan, np.nan, int(m.sum())
    r, p = stats.spearmanr(a[m], b[m])
    return float(r), float(p), int(m.sum())


def partial_spearman(a, b, c):
    """corr(a, b | c) on ranks (residualise ranks of a,b on rank of c)."""
    a, b, c = (np.asarray(v, float) for v in (a, b, c))
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


def wilcoxon_greater(diffs):
    """One-sided Wilcoxon signed-rank, H1: median(diff) > 0 (confidence dropped).

    Returns (p, rank_biserial_effect_size, n_nonzero).
    """
    d = np.asarray([x for x in diffs if np.isfinite(x) and x != 0.0], float)
    n = len(d)
    if n < 6:
        return np.nan, np.nan, n
    try:
        res = stats.wilcoxon(d, alternative="greater")
        p = float(res.pvalue)
    except Exception:
        return np.nan, np.nan, n
    ranks = stats.rankdata(np.abs(d))
    w_pos = ranks[d > 0].sum()
    w_neg = ranks[d < 0].sum()
    total = n * (n + 1) / 2.0
    rb = (w_pos - w_neg) / total
    return p, float(rb), n


def auroc(pos, neg):
    pos = np.asarray([v for v in pos if np.isfinite(v)], float)
    neg = np.asarray([v for v in neg if np.isfinite(v)], float)
    if len(pos) < 1 or len(neg) < 1:
        return np.nan
    y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
    return float(roc_auc_score(y, np.r_[pos, neg]))


# --------------------------------------------------------------------------- #
# Load + aggregate                                                            #
# --------------------------------------------------------------------------- #
def load():
    df = pd.read_csv(OUT / "results_ligand.csv")
    df = df[df.status == "ok"].copy()
    num = ["ligand_rmsd_a", "iptm", "ptm", "ligand_iptm", "complex_plddt",
           "confidence_score", "affinity_pred_value",
           "affinity_probability_binary", "pose_idx"]
    for c in num:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["family"] = df.variant.map(family_of)
    return df


def aggregate(df):
    """One row per (pdbid, variant): best-of-5 and top-1 summaries."""
    rows = []
    for (pid, var), g in df.groupby(["pdbid", "variant"]):
        g = g.dropna(subset=["ligand_rmsd_a"])
        if len(g) == 0:
            continue
        best = g.loc[g.ligand_rmsd_a.idxmin()]
        t = g[g.pose_idx == 0]
        top = t.iloc[0] if len(t) else best
        aff = g.affinity_pred_value.dropna()
        prob = g.affinity_probability_binary.dropna()
        rec = dict(pdbid=pid, variant=var, family=family_of(var), n_pose=len(g),
                   rmsd_best=best.ligand_rmsd_a, rmsd_top1=top.ligand_rmsd_a,
                   rmsd_med=g.ligand_rmsd_a.median(),
                   affinity=aff.iloc[0] if len(aff) else np.nan,
                   prob=prob.iloc[0] if len(prob) else np.nan)
        # best-of-5 confidence (taken on the min-RMSD pose) + top-1 / mean iptm
        for m in METRICS:
            rec[m] = best[m]
        rec["iptm_top1"] = top.iptm
        rec["iptm_mean"] = g.iptm.mean()
        rows.append(rec)
    return pd.DataFrame(rows)


def pair_wt_adv(agg):
    """Pair each adversarial cell to its own system's WT (both must be present)."""
    wt = agg[agg.variant == "wt"].set_index("pdbid")
    out = []
    for _, r in agg[agg.family != "wt"].iterrows():
        if r.pdbid not in wt.index:
            continue
        w = wt.loc[r.pdbid]
        if isinstance(w, pd.DataFrame):
            w = w.iloc[0]
        rec = dict(pdbid=r.pdbid, variant=r.variant, family=r.family,
                   rmsd_best=r.rmsd_best, rmsd_wt=w.rmsd_best,
                   d_rmsd=r.rmsd_best - w.rmsd_best,
                   rmsd_top1=r.rmsd_top1, rmsd_wt_top1=w.rmsd_top1,
                   d_rmsd_top1=r.rmsd_top1 - w.rmsd_top1,
                   d_aff=r.affinity - w.affinity, aff_adv=r.affinity, aff_wt=w.affinity,
                   d_prob=r.prob - w.prob, prob_adv=r.prob, prob_wt=w.prob)
        for m in METRICS:
            rec[f"d_{m}"] = w[m] - r[m]      # + = dropped on adversarial
            rec[f"{m}_adv"] = r[m]
            rec[f"{m}_wt"] = w[m]
        out.append(rec)
    p = pd.DataFrame(out)
    p["stratum"] = np.where(p.rmsd_best < MEM, "memorized",
                    np.where(p.rmsd_best >= RESP, "responded", "middle"))
    return p


# --------------------------------------------------------------------------- #
# Section 1: confidence drop WT -> variant / family                          #
# --------------------------------------------------------------------------- #
def s1_confidence_drop(pair):
    print("\n" + "=" * 94)
    print("S1  CONFIDENCE DROP  WT -> ligand variant   (delta = conf[WT] - conf[adv]; + = dropped)")
    print("    interface metrics (iptm, ligand_iptm) should drop; global (ptm, complex_plddt) = control")
    print("=" * 94)

    def block(label, sub):
        rows = []
        for m in METRICS:
            col = f"d_{m}"
            d = sub[col].dropna()
            wcol, acol = f"{m}_wt", f"{m}_adv"
            if d.empty:
                continue
            p, rb, nnz = wilcoxon_greater(d.values)
            lo, hi = boot_ci(d.values)
            rows.append(dict(
                group=label, metric=m,
                kind=("interface" if m in INTERFACE
                      else "global" if m in GLOBAL else "composite"),
                n=int(len(d)),
                med_wt=round(sub[wcol].median(), 4),
                med_adv=round(sub[acol].median(), 4),
                med_delta=round(d.median(), 4),
                ci_lo=round(lo, 4), ci_hi=round(hi, 4),
                frac_dropped=round(float((d > 0).mean()), 3),
                wilcoxon_p_greater=p, rank_biserial=rb,
                auroc_wt_vs_adv=round(auroc(sub[wcol].values, sub[acol].values), 4),
            ))
        return rows

    all_rows = []
    # by family
    for fam in FAMILIES:
        all_rows += block(fam, pair[pair.family == fam])
    # pooled adversarial
    all_rows += block("ALL_adv", pair)
    # per individual variant
    per_variant = []
    for var in VARIANT_ORDER:
        sub = pair[pair.variant == var]
        if len(sub) < 6:
            # still record n for transparency, but skip stats-heavy block
            if len(sub) > 0:
                d = sub["d_iptm"].dropna()
                per_variant.append(dict(
                    group=var, metric="iptm", kind="interface", n=int(len(sub)),
                    med_wt=round(sub["iptm_wt"].median(), 4),
                    med_adv=round(sub["iptm_adv"].median(), 4),
                    med_delta=round(d.median(), 4) if len(d) else np.nan,
                    ci_lo=np.nan, ci_hi=np.nan,
                    frac_dropped=round(float((d > 0).mean()), 3) if len(d) else np.nan,
                    wilcoxon_p_greater=np.nan, rank_biserial=np.nan,
                    auroc_wt_vs_adv=round(auroc(sub["iptm_wt"].values, sub["iptm_adv"].values), 4),
                ))
            continue
        d = sub["d_iptm"].dropna()
        p, rb, nnz = wilcoxon_greater(d.values)
        lo, hi = boot_ci(d.values)
        per_variant.append(dict(
            group=var, metric="iptm", kind="interface", n=int(len(sub)),
            med_wt=round(sub["iptm_wt"].median(), 4),
            med_adv=round(sub["iptm_adv"].median(), 4),
            med_delta=round(d.median(), 4),
            ci_lo=round(lo, 4), ci_hi=round(hi, 4),
            frac_dropped=round(float((d > 0).mean()), 3),
            wilcoxon_p_greater=p, rank_biserial=rb,
            auroc_wt_vs_adv=round(auroc(sub["iptm_wt"].values, sub["iptm_adv"].values), 4),
        ))

    fam_df = pd.DataFrame(all_rows)
    var_df = pd.DataFrame(per_variant)
    # print family table (interface + a couple global controls)
    print("\n  By FAMILY (and pooled), all metrics:")
    _print_drop_table(fam_df)
    print("\n  Per individual variant (iptm only; small-n variants flagged):")
    _print_drop_table(var_df)

    fam_df.to_csv(OUT / "confidence_stats_ligand.csv", index=False)
    var_df.to_csv(OUT / "confidence_stats_pervariant_ligand.csv", index=False)

    # per-cell iptm dump (mirrors paired_confidence_full.csv)
    dump = pair[["pdbid", "variant", "family", "iptm_wt", "iptm_adv", "d_iptm",
                 "rmsd_best", "d_rmsd", "stratum"]].copy()
    dump = dump.rename(columns={"iptm_wt": "wt_iptm", "iptm_adv": "adv_iptm",
                                "d_iptm": "delta_iptm"})
    dump.to_csv(OUT / "paired_confidence_ligand.csv", index=False)
    return fam_df, var_df


def _print_drop_table(d):
    if d.empty:
        print("    (no rows)")
        return
    hdr = (f"  {'group':16} {'metric':24} {'kind':9} {'n':>4} {'med_wt':>7} "
           f"{'med_adv':>7} {'medD':>7} {'%drop':>5} {'wilcox_p':>9} "
           f"{'rankbis':>7} {'AUROC':>6}")
    print(hdr)
    for _, r in d.iterrows():
        p = r["wilcoxon_p_greater"]
        ps = ("  n/a" if (p != p) else "<1e-4" if p < 1e-4 else f"{p:.4f}")
        rb = r["rank_biserial"]
        rbs = ("  n/a" if (rb != rb) else f"{rb:+.3f}")
        fd = r["frac_dropped"]
        fds = ("  n/a" if (fd != fd) else f"{fd*100:.0f}%")
        md = r["med_delta"]
        mds = ("  n/a" if (md != md) else f"{md:+.3f}")
        au = r["auroc_wt_vs_adv"]
        aus = ("  n/a" if (au != au) else f"{au:.3f}")
        print(f"  {r['group']:16} {r['metric']:24} {r['kind']:9} {int(r['n']):>4} "
              f"{r['med_wt']:>7.3f} {r['med_adv']:>7.3f} {mds:>7} {fds:>5} "
              f"{ps:>9} {rbs:>7} {aus:>6}")


# --------------------------------------------------------------------------- #
# Section 1b: conditional on structural response                             #
# --------------------------------------------------------------------------- #
def s1b_conditional(pair):
    print("\n" + "=" * 94)
    print("S1b CONDITIONAL: does confidence drop within strata of the STRUCTURAL response?")
    print("    memorized = best-of-5 adv RMSD < 2 A   |   responded = best-of-5 adv RMSD >= 4 A")
    print("    NOTE: on the ligand axis the memorized stratum is nearly empty (structure refuses")
    print("          the ligand) -> that test is n-starved here; the responded stratum carries it.")
    print("=" * 94)
    rows = []
    for scope_name, mask in [("ALL_adv", pair.family != "wt")] + \
                            [(f, pair.family == f) for f in FAMILIES]:
        base = pair[mask]
        for stratum in ["memorized", "middle", "responded", "all"]:
            s = base if stratum == "all" else base[base.stratum == stratum]
            for m in ("iptm", "ptm"):       # interface vs global control
                d = s[f"d_{m}"].dropna()
                if len(d) < 6:
                    continue
                p, rb, _ = wilcoxon_greater(d.values)
                rows.append(dict(
                    scope=scope_name, stratum=stratum, metric=m, n=int(len(d)),
                    med_wt=round(s[f"{m}_wt"].median(), 4),
                    med_adv=round(s[f"{m}_adv"].median(), 4),
                    med_delta=round(d.median(), 4),
                    wilcoxon_p_greater=p, rank_biserial=rb,
                    auroc_wt_vs_adv=round(auroc(s[f"{m}_wt"].values, s[f"{m}_adv"].values), 4),
                ))
    cond = pd.DataFrame(rows)
    if not cond.empty:
        hdr = (f"  {'scope':10} {'stratum':10} {'metric':6} {'n':>4} {'med_wt':>7} "
               f"{'med_adv':>7} {'medD':>7} {'wilcox_p':>9} {'rankbis':>7} {'AUROC':>6}")
        print(hdr)
        for _, r in cond.iterrows():
            p = r["wilcoxon_p_greater"]
            ps = "<1e-4" if p < 1e-4 else f"{p:.4f}"
            print(f"  {r['scope']:10} {r['stratum']:10} {r['metric']:6} {int(r['n']):>4} "
                  f"{r['med_wt']:>7.3f} {r['med_adv']:>7.3f} {r['med_delta']:>+7.3f} "
                  f"{ps:>9} {r['rank_biserial']:>+7.3f} {r['auroc_wt_vs_adv']:>6.3f}")
    cond.to_csv(OUT / "confidence_conditional_ligand.csv", index=False)
    return cond


# --------------------------------------------------------------------------- #
# Section 2: does confidence track the structural response?                  #
# --------------------------------------------------------------------------- #
def s2_tracks_response(pair, df):
    print("\n" + "=" * 94)
    print("S2  DOES CONFIDENCE TRACK THE STRUCTURAL RESPONSE?  (mutate ligand -> dRMSD -> dConf)")
    print("=" * 94)

    # 2a between-case dConf(iptm) ~ dRMSD, overall + by family + stratified
    print("\n[2a] between-case  Spearman(dRMSD, dConf[iptm])  + stratified median drop")
    r_all, p_all, n_all = sp(pair.d_rmsd.values, pair.d_iptm.values)
    print(f"  ALL_adv  n={n_all:4d}  Spearman(dRMSD,dConf)={r_all:+.3f} (p={p_all:.1e})")
    rows2a = [dict(group="ALL_adv", n=n_all, spearman_drmsd_dconf=r_all, p=p_all)]
    for fam in FAMILIES:
        s = pair[pair.family == fam]
        r, pv, n = sp(s.d_rmsd.values, s.d_iptm.values)
        line = f"  {fam:9s} n={n:4d}  Spearman(dRMSD,dConf)={r:+.3f} (p={pv:.1e}) | medDiptm: "
        line += "  ".join(
            f"{st}={s[s.stratum==st].d_iptm.median():+.3f}(n{len(s[s.stratum==st])})"
            for st in ["memorized", "middle", "responded"] if len(s[s.stratum == st]))
        print(line)
        rows2a.append(dict(group=fam, n=n, spearman_drmsd_dconf=r, p=pv))
    pd.DataFrame(rows2a).to_csv(OUT / "s2a_dconf_drmsd_ligand.csv", index=False)

    # 2b within-case pose-level Spearman(iptm, rmsd) across the 5 samples
    print("\n[2b] within-case pose-level Spearman(iptm, rmsd) across the 5 diffusion samples")
    print("     (negative = confidence ranks lower-RMSD poses higher = a valid quality signal)")
    wc = []
    for (pid, var), g in df.groupby(["pdbid", "variant"]):
        g = g.dropna(subset=["ligand_rmsd_a", "iptm"])
        if len(g) < 3 or g.ligand_rmsd_a.nunique() < 3:
            continue
        rho, _ = stats.spearmanr(g.iptm, g.ligand_rmsd_a)
        wc.append(dict(pdbid=pid, variant=var, family=family_of(var),
                       cls="wt" if var == "wt" else "adv", rho=rho))
    wc = pd.DataFrame(wc)
    for cls in ["wt", "adv"]:
        s = wc[wc.cls == cls].rho.dropna()
        if len(s) < 6:
            continue
        w = stats.wilcoxon(s, alternative="less") if (s != 0).any() else None
        wp = f"{w.pvalue:.1e}" if w is not None else "n/a"
        print(f"  {cls:3s}  n={len(s):4d}  median rho={s.median():+.3f}  "
              f"frac rho<0={np.mean(s < 0):.0%}  (Wilcoxon rho<0 p={wp})")
    # per-family adv breakdown
    for fam in FAMILIES:
        s = wc[(wc.cls == "adv") & (wc.family == fam)].rho.dropna()
        if len(s) < 6:
            continue
        print(f"       adv:{fam:9s} n={len(s):4d}  median rho={s.median():+.3f}  "
              f"frac rho<0={np.mean(s < 0):.0%}")
    wc.to_csv(OUT / "s2b_within_case_rho_ligand.csv", index=False)

    # 2c pooled-pose AUROC: confidence higher on near-native poses?
    print("\n[2c] pooled-pose AUROC: confidence separates near-native(<2) vs displaced(>=4) poses")
    rows2c = []
    for scope, mask in [("WT", df.variant == "wt"), ("ADV", df.family != "wt")] + \
                       [(f"ADV_{f}", df.family == f) for f in FAMILIES]:
        g = df[mask].dropna(subset=["ligand_rmsd_a", "iptm"])
        near = g[g.ligand_rmsd_a < MEM].iptm
        far = g[g.ligand_rmsd_a >= RESP].iptm
        if len(near) < 10 or len(far) < 10:
            rows2c.append(dict(scope=scope, n_near=len(near), n_far=len(far), auroc=np.nan))
            continue
        au = auroc(near.values, far.values)
        print(f"  {scope:10s}  near n={len(near):4d}  far n={len(far):5d}  "
              f"AUROC(conf higher on near-native)={au:.3f}")
        rows2c.append(dict(scope=scope, n_near=len(near), n_far=len(far), auroc=au))
    pd.DataFrame(rows2c).to_csv(OUT / "s2c_pooled_pose_auroc_ligand.csv", index=False)
    return wc


# --------------------------------------------------------------------------- #
# Section 3: affinity invariance / the dichotomy (with RTM control)          #
# --------------------------------------------------------------------------- #
def s3_affinity(pair):
    print("\n" + "=" * 94)
    print("S3  AFFINITY INVARIANCE / THE DICHOTOMY  (with the regression-to-mean control)")
    print("    dAff = adv-wt logIC50 (+ = predicted WEAKER; physics wants strongly +).")
    print("    The documented dichotomy: structure REFUSES these ligands, but affinity says TIGHTER.")
    print("=" * 94)
    p = pair.dropna(subset=["d_aff"]).copy()

    # stratified by structural response, overall and by family
    print("\n  dAff stratified by structural response (best-of-5 adv RMSD):")
    rows = []
    for scope_name, mask in [("ALL_adv", p.family != "wt")] + \
                            [(f, p.family == f) for f in FAMILIES]:
        base = p[mask]
        for stratum in ["memorized", "middle", "responded", "all"]:
            s = base if stratum == "all" else base[base.stratum == stratum]
            if len(s) == 0:
                continue
            lo, hi = boot_ci(s.d_aff.values)
            plo, phi = boot_ci(s.d_prob.values)
            rows.append(dict(
                scope=scope_name, stratum=stratum, n=len(s),
                med_rmsd_adv=round(s.rmsd_best.median(), 2),
                med_aff_wt=round(s.aff_wt.median(), 3),
                med_dAff=round(s.d_aff.median(), 3), dAff_ci=f"[{lo:+.2f},{hi:+.2f}]",
                med_dProb=round(s.d_prob.median(), 3), dProb_ci=f"[{plo:+.2f},{phi:+.2f}]"))
    rep = pd.DataFrame(rows)
    show = rep[rep.stratum.isin(["all", "responded"])]
    print(show.to_string(index=False))
    rep.to_csv(OUT / "s3_affinity_strata_ligand.csv", index=False)

    # the RTM control: raw vs partial( dAff , displacement | WT affinity )
    print("\n  REGRESSION-TO-MEAN control (the critical one):")
    print("  Memorized stratum is ~empty here, so the pocket-axis 'responds more when memorized'")
    print("  pattern can't form; the RTM question is whether dAff is driven by WT tightness.")
    rtm_rows = []
    for scope_name, mask in [("ALL_adv", p.family != "wt")] + \
                            [(f, p.family == f) for f in FAMILIES]:
        s = p[mask]
        r_wt, p_wt, n = sp(s.d_aff.values, s.aff_wt.values)        # RTM signature
        r_abs, p_abs, _ = sp(s.rmsd_best.values, s.d_aff.values)   # raw |displacement|
        r_dr, p_dr, _ = sp(s.d_rmsd.values, s.d_aff.values)        # raw dRMSD
        # partials controlling for WT affinity
        pa_abs, pp_abs, _ = partial_spearman(s.d_aff.values, s.rmsd_best.values, s.aff_wt.values)
        pa_dr, pp_dr, _ = partial_spearman(s.d_aff.values, s.d_rmsd.values, s.aff_wt.values)
        print(f"\n  [{scope_name}] n={n}")
        print(f"    Spearman(dAff, WT aff)               = {r_wt:+.3f} (p={p_wt:.1e})   <- RTM signature (expect strongly -)")
        print(f"    raw     Spearman(dAff, |displacement|)= {r_abs:+.3f} (p={p_abs:.1e})")
        print(f"    partial (dAff, |displacement| | WTaff)= {pa_abs:+.3f} (p={pp_abs:.1e})")
        print(f"    raw     Spearman(dAff, dRMSD)        = {r_dr:+.3f} (p={p_dr:.1e})")
        print(f"    partial (dAff, dRMSD | WT aff)       = {pa_dr:+.3f} (p={pp_dr:.1e})")
        rtm_rows.append(dict(scope=scope_name, n=n,
                             sp_dAff_wtAff=round(r_wt, 3),
                             raw_dAff_absdisp=round(r_abs, 3),
                             partial_dAff_absdisp=round(pa_abs, 3),
                             raw_dAff_dRMSD=round(r_dr, 3),
                             partial_dAff_dRMSD=round(pa_dr, 3)))
    pd.DataFrame(rtm_rows).to_csv(OUT / "s3_affinity_rtm_ligand.csv", index=False)
    print("\n  READ: a strongly-negative Spearman(dAff, WT aff) = RTM is in play; if the raw")
    print("        dAff<->displacement coupling collapses after partialling out WT affinity,")
    print("        the 'tighter' verdict is RTM on WT tightness, not pose-reading.")
    return rep


# --------------------------------------------------------------------------- #
# Section 4: joint structure + decoupling                                     #
# --------------------------------------------------------------------------- #
def s4_joint(pair, agg):
    print("\n" + "=" * 94)
    print("S4  JOINT STRUCTURE (Boltz-2 adversarial cells, best-of-5) + head decoupling")
    print("=" * 94)
    a = agg[agg.family != "wt"].copy()
    cols = {"rmsd_best": "RMSD", "iptm": "iptm", "complex_plddt": "plddt",
            "affinity": "affinity", "prob": "P(bind)"}
    M = a[list(cols)].rename(columns=cols)
    corr = M.corr(method="spearman")
    print("\nSpearman correlation matrix (adversarial cells, all families pooled):")
    print(corr.round(3).to_string())
    corr.to_csv(OUT / "s4_corr_matrix_ligand.csv")

    print("\nPartial correlations (conf <-> affinity, controlling for RMSD):")
    for cname, col in [("iptm", "iptm"), ("plddt", "complex_plddt"),
                       ("confscore", "confidence_score")]:
        r_raw, p_raw, n = sp(a[col].values, a.affinity.values)
        r_par, p_par, _ = partial_spearman(a[col].values, a.affinity.values, a.rmsd_best.values)
        print(f"  {cname:9s} <-> affinity : raw rho={r_raw:+.3f} (p={p_raw:.1e})  |  "
              f"partial rho|RMSD={r_par:+.3f} (p={p_par:.1e})  n={n}")

    # decoupling of the RESPONSES: diptm <-> daff (raw + partial | dRMSD)
    pp = pair.dropna(subset=["d_iptm", "d_aff"])
    r_raw, p_raw, n = sp(pp.d_iptm.values, pp.d_aff.values)
    r_par, p_par, _ = partial_spearman(pp.d_iptm.values, pp.d_aff.values, pp.d_rmsd.values)
    print(f"\n  Diptm <-> Daff : raw rho={r_raw:+.3f} (p={p_raw:.1e})  |  "
          f"partial rho|dRMSD={r_par:+.3f} (p={p_par:.1e})  n={n}")
    print("  READ: ~0 here = the confidence head and affinity head respond to the ligand")
    print("        mutation INDEPENDENTLY (getting less confident does not coincide with")
    print("        predicting weaker binding) -> same three-head dissociation as the pocket axis.")
    return a, corr


# --------------------------------------------------------------------------- #
# Figure (4-panel, mirrors 17)                                                #
# --------------------------------------------------------------------------- #
def make_figure(pair, wc, corr):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 2, figsize=(12.5, 10))
    fam_colors = {"halo": "#3a7", "chrg_neu": "#e34", "chrg_pos": "#b59", "meth": "#39c"}

    # (a) S3: dAff vs adv RMSD, coloured by family (the dichotomy panel)
    for fam, c in fam_colors.items():
        d = pair[(pair.family == fam)].dropna(subset=["d_aff"])
        ax[0, 0].scatter(d.rmsd_best, d.d_aff, s=12, c=c, alpha=.5,
                         label=f"{fam} (n={len(d)})")
    ax[0, 0].axhline(0, color="k", lw=.7)
    ax[0, 0].axvline(2, color="gray", ls=":"); ax[0, 0].axvline(4, color="gray", ls=":")
    ax[0, 0].set_xlabel("adversarial best-of-5 ligand RMSD (Angstrom)")
    ax[0, 0].set_ylabel("delta affinity (adv-wt logIC50, + = weaker)")
    ax[0, 0].set_title("(a) affinity vs structural response  (the dichotomy:\n"
                       "ligand ejected >>2A yet many dAff<0 = 'tighter')")
    ax[0, 0].legend(fontsize=8, loc="upper right"); ax[0, 0].set_xlim(0, 15)

    # (b) S2a: dConf vs dRMSD, coloured by family
    for fam, c in fam_colors.items():
        d = pair[(pair.family == fam)].dropna(subset=["d_rmsd", "d_iptm"])
        ax[0, 1].scatter(d.d_rmsd, d.d_iptm, s=12, c=c, alpha=.45, label=fam)
    r, _, n = sp(pair.d_rmsd.values, pair.d_iptm.values)
    ax[0, 1].axhline(0, color="k", lw=.7); ax[0, 1].axvline(0, color="k", lw=.7)
    ax[0, 1].set_xlabel("delta RMSD (adv-wt, Angstrom)")
    ax[0, 1].set_ylabel("delta confidence (wt-adv iptm, + = dropped)")
    ax[0, 1].set_title(f"(b) confidence tracks response  (rho={r:+.2f}, n={n})")
    ax[0, 1].legend(fontsize=8, loc="lower right")

    # (c) S2b: within-case rho distribution WT vs ADV
    for cls, c in [("wt", "#888"), ("adv", "#e34")]:
        s = wc[wc.cls == cls].rho.dropna()
        if len(s):
            ax[1, 0].hist(s, bins=20, alpha=.6, color=c,
                          label=f"{cls} (med {s.median():+.2f}, n={len(s)})")
    ax[1, 0].axvline(0, color="k", lw=.8)
    ax[1, 0].set_xlabel("within-case Spearman(iptm, RMSD) across 5 poses")
    ax[1, 0].set_ylabel("# cases")
    ax[1, 0].set_title("(c) pose-level: confident about which pose?")
    ax[1, 0].legend(fontsize=8)

    # (d) S4: correlation heatmap
    im = ax[1, 1].imshow(corr.values, vmin=-1, vmax=1, cmap="RdBu_r")
    ax[1, 1].set_xticks(range(len(corr))); ax[1, 1].set_yticks(range(len(corr)))
    ax[1, 1].set_xticklabels(corr.columns, rotation=45, ha="right")
    ax[1, 1].set_yticklabels(corr.index)
    for i in range(len(corr)):
        for j in range(len(corr)):
            ax[1, 1].text(j, i, f"{corr.values[i, j]:.2f}", ha="center",
                          va="center", fontsize=8)
    ax[1, 1].set_title("(d) joint structure (Boltz-2 adversarial)")
    fig.colorbar(im, ax=ax[1, 1], fraction=.046)

    fig.suptitle("Confidence x Affinity x RMSD - CASF LIGAND mutagenesis (Boltz-2)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, .98])
    out = FIG / "conf_aff_rmsd_ligand.png"
    fig.savefig(out, dpi=130)
    fig.savefig(FIG / "conf_aff_rmsd_ligand.pdf")
    print(f"\n[wrote figure] {out}")


# --------------------------------------------------------------------------- #
def main():
    df = load()
    agg = aggregate(df)
    pair = pair_wt_adv(agg)
    pair.to_csv(OUT / "paired_conf_aff_rmsd_ligand.csv", index=False)

    print("=" * 94)
    print("LIGAND-AXIS confidence x affinity x RMSD  (Boltz-2 only; mirror of pocket scripts 16/17)")
    print(f"  adversarial cells paired to own-system WT: {len(pair)}")
    print("  family n (best-of-5):  " +
          "  ".join(f"{f}={int((pair.family==f).sum())}" for f in FAMILIES))
    print("  structural-response strata (best-of-5 adv RMSD):  " +
          "  ".join(f"{s}={int((pair.stratum==s).sum())}"
                    for s in ["memorized", "middle", "responded"]))
    print("=" * 94)

    s1_confidence_drop(pair)
    s1b_conditional(pair)
    wc = s2_tracks_response(pair, df)
    s3_affinity(pair)
    a, corr = s4_joint(pair, agg)
    make_figure(pair, wc, corr)

    print(f"\n[wrote] paired_conf_aff_rmsd_ligand.csv  ({len(pair)} adversarial cells)")
    print("[done]")


if __name__ == "__main__":
    main()
