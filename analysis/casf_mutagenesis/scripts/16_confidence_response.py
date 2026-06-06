#!/usr/bin/env python3
"""Does model confidence register a broken pocket? (rem / pack / inv)

Reframes the memorization question. The structure head keeps the ligand near
the WT pose and the affinity head keeps ~WT affinity (both memorize). This
asks a different thing: does the model's SELF-REPORTED CONFIDENCE drop when we
mutate the pocket? If interface confidence (iptm / ligand_iptm) drops
significantly — and especially if it drops on the very cells where the
structure memorized — then the model "memorized something but knows it
memorized," and confidence can serve as a practical flag for untrustworthy
predictions.

Four escalating layers:
  1. MARGINAL    paired WT->adv drop per (model, variant, metric):
                 Wilcoxon signed-rank (one-sided), median delta + bootstrap CI,
                 matched-pairs rank-biserial effect size, fraction dropped.
  2. CONTROL     interface confidence (iptm/ligand_iptm) vs global-fold
                 confidence (ptm/complex_plddt). The protein still folds, so
                 global confidence is the built-in negative control; if the
                 interface drops while global doesn't, the drop is localized
                 signal, not run-to-run noise.
  3. UTILITY     AUROC of each confidence metric as a WT-vs-adversarial flag
                 (0.5 = useless, 1.0 = perfect detector).
  4. CONDITIONAL restrict to cells where the structure memorized (WT correct,
                 adv ligand RMSD < 2 A): does confidence STILL drop there?
                 Plus Spearman(delta_confidence, delta_rmsd) to test whether the
                 confidence drop tracks the structural response magnitude.

Pure stdlib (no pandas/scipy) so it runs on any python3.

Inputs  (analysis/casf_mutagenesis/outputs/):
  results_full.csv   per (pdbid, variant, model, pose_idx): confidence + rmsd
  paired_full.csv    per (pdbid, model, variant, pose_idx): structure labels
Outputs (same dir):
  paired_confidence_full.csv        per-cell WT vs adv confidence + deltas
  confidence_stats_full.csv         per (model, variant, metric) test table
  confidence_conditional_full.csv   conditional-on-memorization table
"""
from __future__ import annotations

import csv
import math
import random
from collections import Counter, defaultdict
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "outputs"
RESULTS = OUT / "results_full.csv"
PAIRED = OUT / "paired_full.csv"

ADV = ("rem", "pack", "inv")
# Confidence metrics to test. Interface metrics should drop; global metrics are
# the negative control. affinity_probability_binary only exists for Boltz2.
INTERFACE = ("iptm", "ligand_iptm")
GLOBAL = ("ptm", "complex_plddt")
COMPOSITE = ("confidence_score", "affinity_probability_binary")
METRICS = INTERFACE + GLOBAL + COMPOSITE
NAN = float("nan")


# --------------------------------------------------------------------------- #
# Small stats toolkit (stdlib only)                                           #
# --------------------------------------------------------------------------- #
def _phi(z: float) -> float:
    """Standard normal CDF."""
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def median(xs):
    s = sorted(xs)
    n = len(s)
    if n == 0:
        return NAN
    m = n // 2
    return s[m] if n % 2 else 0.5 * (s[m - 1] + s[m])


def rankdata(vals):
    """1-based ranks with averaged ties (like scipy.stats.rankdata)."""
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    ranks = [0.0] * len(vals)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0  # average of 1-based ranks i+1..j+1
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def wilcoxon_greater(diffs):
    """One-sided Wilcoxon signed-rank, H1: median(diff) > 0.

    diff convention here is wt - adv, so a positive median = confidence dropped
    on the adversarial variant. Returns p for H1, the z-stat, and the
    matched-pairs rank-biserial correlation (effect size, -1..1).
    """
    d = [x for x in diffs if x != 0.0]
    n = len(d)
    if n < 6:
        return dict(n=n, z=NAN, p=NAN, rank_biserial=NAN)
    absd = [abs(x) for x in d]
    ranks = rankdata(absd)
    w_pos = sum(r for r, x in zip(ranks, d) if x > 0)
    w_neg = sum(r for r, x in zip(ranks, d) if x < 0)
    total = n * (n + 1) / 2.0
    mean = n * (n + 1) / 4.0
    var = n * (n + 1) * (2 * n + 1) / 24.0
    for t in Counter(absd).values():       # tie correction
        if t > 1:
            var -= (t ** 3 - t) / 48.0
    if var <= 0:
        return dict(n=n, z=NAN, p=NAN, rank_biserial=NAN)
    z = (w_pos - mean - 0.5) / math.sqrt(var)   # continuity-corrected
    return dict(n=n, z=z, p=1.0 - _phi(z), rank_biserial=(w_pos - w_neg) / total)


def auroc(pos, neg):
    """P(score_pos > score_neg) via Mann-Whitney U (ties count 0.5)."""
    if not pos or not neg:
        return NAN
    allv = [(v, 1) for v in pos] + [(v, 0) for v in neg]
    ranks = rankdata([v for v, _ in allv])
    r_pos = sum(r for r, (_, lab) in zip(ranks, allv) if lab == 1)
    n1, n0 = len(pos), len(neg)
    u = r_pos - n1 * (n1 + 1) / 2.0
    return u / (n1 * n0)


def boot_median_ci(xs, nboot=2000, seed=0, alpha=0.05):
    if len(xs) < 3:
        return (NAN, NAN)
    rng = random.Random(seed)
    n = len(xs)
    meds = sorted(median([xs[rng.randrange(n)] for _ in range(n)]) for _ in range(nboot))
    return (meds[int((alpha / 2) * nboot)], meds[int((1 - alpha / 2) * nboot)])


def spearman(x, y):
    if len(x) < 3:
        return NAN
    rx, ry = rankdata(x), rankdata(y)
    n = len(rx)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else NAN


def _f(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def _b(s):
    return str(s).strip().lower() in ("true", "1", "1.0", "yes")


# --------------------------------------------------------------------------- #
# Load                                                                        #
# --------------------------------------------------------------------------- #
def load_confidence():
    """conf[(model, pdbid, variant)] = {metric: float}, rank-0 pose only."""
    conf = {}
    with RESULTS.open() as fh:
        for row in csv.DictReader(fh):
            if row["pose_idx"] != "0":
                continue
            key = (row["model"], row["pdbid"], row["variant"])
            conf[key] = {m: _f(row.get(m, "")) for m in METRICS}
    return conf


def load_structure_labels():
    """lab[(model, pdbid, variant)] = {wt_correct, memorized, delta_rmsd, adv_rmsd}."""
    lab = {}
    with PAIRED.open() as fh:
        for row in csv.DictReader(fh):
            if row.get("pose_idx", "0") != "0":
                continue
            lab[(row["model"], row["pdbid"], row["variant"])] = dict(
                wt_correct=_b(row.get("wt_correct_2A")),
                memorized=_b(row.get("memorized_given_wt")),
                delta_rmsd=_f(row.get("delta_rmsd_a")),
                adv_rmsd=_f(row.get("adv_rmsd_a")),
                wt_rmsd=_f(row.get("wt_rmsd_a")),
            )
    return lab


# --------------------------------------------------------------------------- #
# Pairing + analysis                                                          #
# --------------------------------------------------------------------------- #
def paired(conf, model, variant, metric):
    """Return (pdbids, wt_vals, adv_vals) for cells with both present."""
    ids, wt, adv = [], [], []
    for (m, pid, v), d in conf.items():
        if m != model or v != variant:
            continue
        a = d.get(metric)
        w = conf.get((model, pid, "wt"), {}).get(metric)
        if a is not None and w is not None and not math.isnan(a) and not math.isnan(w):
            ids.append(pid); wt.append(w); adv.append(a)
    return ids, wt, adv


def models_present(conf):
    seen = []
    for (m, _, _) in conf:
        if m not in seen:
            seen.append(m)
    order = ["Boltz2", "AF3+MSA", "AF3"]
    return [m for m in order if m in seen] + [m for m in seen if m not in order]


def main():
    conf = load_confidence()
    labels = load_structure_labels()
    models = models_present(conf)

    # ---- per-cell paired confidence dump + stats table --------------------
    stat_rows = []
    pair_rows = []
    for model in models:
        for variant in ADV:
            for metric in METRICS:
                ids, wt, adv = paired(conf, model, variant, metric)
                if len(ids) < 6:
                    continue
                diffs = [w - a for w, a in zip(wt, adv)]   # + = dropped
                wl = wilcoxon_greater(diffs)
                lo, hi = boot_median_ci(diffs)
                stat_rows.append(dict(
                    model=model, variant=variant, metric=metric, n=len(ids),
                    median_wt=round(median(wt), 4), median_adv=round(median(adv), 4),
                    median_delta=round(median(diffs), 4),
                    delta_ci_lo=round(lo, 4), delta_ci_hi=round(hi, 4),
                    frac_dropped=round(sum(1 for x in diffs if x > 0) / len(diffs), 3),
                    wilcoxon_p_greater=wl["p"], rank_biserial=wl["rank_biserial"],
                    auroc_wt_vs_adv=round(auroc(wt, adv), 4),
                    kind=("interface" if metric in INTERFACE
                          else "global" if metric in GLOBAL else "composite"),
                ))
                if metric == "iptm":
                    for pid, w, a in zip(ids, wt, adv):
                        pair_rows.append(dict(model=model, pdbid=pid, variant=variant,
                                              wt_iptm=round(w, 4), adv_iptm=round(a, 4),
                                              delta_iptm=round(w - a, 4)))

    # ---- conditional: within structure-memorized cells (Boltz2) -----------
    cond_rows = []
    for model in models:
        for subset_name, keep in (
            ("memorized", lambda L: L["wt_correct"] and L["memorized"]),
            ("responded", lambda L: L["wt_correct"] and (L["adv_rmsd"] or 0) >= 4.0),
            ("all_wtok", lambda L: L["wt_correct"]),
        ):
            for metric in ("iptm", "ptm"):     # interface vs global control
                wt, adv, drmsd = [], [], []
                for (m, pid, v), L in labels.items():
                    if m != model or v not in ADV or not keep(L):
                        continue
                    a = conf.get((model, pid, v), {}).get(metric)
                    w = conf.get((model, pid, "wt"), {}).get(metric)
                    if a is None or w is None or math.isnan(a) or math.isnan(w):
                        continue
                    wt.append(w); adv.append(a)
                    if L["delta_rmsd"] is not None:
                        drmsd.append(L["delta_rmsd"])
                if len(wt) < 6:
                    continue
                diffs = [w - a for w, a in zip(wt, adv)]
                wl = wilcoxon_greater(diffs)
                # Spearman(delta_conf, delta_rmsd) only meaningful on the full set
                rho = NAN
                if len(diffs) == len(drmsd) and subset_name == "all_wtok":
                    rho = spearman(diffs, drmsd)
                cond_rows.append(dict(
                    model=model, subset=subset_name, metric=metric, n=len(diffs),
                    median_wt=round(median(wt), 4), median_adv=round(median(adv), 4),
                    median_delta=round(median(diffs), 4),
                    wilcoxon_p_greater=wl["p"], rank_biserial=wl["rank_biserial"],
                    auroc_wt_vs_adv=round(auroc(wt, adv), 4),
                    spearman_dconf_drmsd=rho,
                ))

    # ---- write ------------------------------------------------------------
    _write_csv(OUT / "paired_confidence_full.csv", pair_rows,
               ["model", "pdbid", "variant", "wt_iptm", "adv_iptm", "delta_iptm"])
    _write_csv(OUT / "confidence_stats_full.csv", stat_rows,
               ["model", "variant", "metric", "kind", "n", "median_wt", "median_adv",
                "median_delta", "delta_ci_lo", "delta_ci_hi", "frac_dropped",
                "wilcoxon_p_greater", "rank_biserial", "auroc_wt_vs_adv"])
    _write_csv(OUT / "confidence_conditional_full.csv", cond_rows,
               ["model", "subset", "metric", "n", "median_wt", "median_adv",
                "median_delta", "wilcoxon_p_greater", "rank_biserial",
                "auroc_wt_vs_adv", "spearman_dconf_drmsd"])

    _report(stat_rows, cond_rows, models)


def _fmt_p(p):
    if p != p:   # nan
        return "  n/a "
    return "<1e-4" if p < 1e-4 else f"{p:.4f}"


def _report(stat_rows, cond_rows, models):
    print("=" * 92)
    print("CONFIDENCE RESPONSE TO POCKET MUTATION  (delta = conf[WT] - conf[adv]; + = dropped)")
    print("=" * 92)
    for model in models:
        rows = [r for r in stat_rows if r["model"] == model]
        if not rows:
            continue
        print(f"\n### {model}    (n per cell shown; rank-0 pose)")
        print(f"{'variant':7} {'metric':26} {'kind':9} {'n':>4} "
              f"{'med_wt':>7} {'med_adv':>7} {'medΔ':>7} {'95%CI':>16} "
              f"{'%drop':>5} {'wilcox_p':>8} {'rankbis':>7} {'AUROC':>6}")
        for r in rows:
            ci = f"[{r['delta_ci_lo']:+.3f},{r['delta_ci_hi']:+.3f}]"
            rb = f"{r['rank_biserial']:+.3f}" if r["rank_biserial"] == r["rank_biserial"] else "  n/a"
            print(f"{r['variant']:7} {r['metric']:26} {r['kind']:9} {r['n']:>4} "
                  f"{r['median_wt']:>7.3f} {r['median_adv']:>7.3f} {r['median_delta']:>+7.3f} "
                  f"{ci:>16} {r['frac_dropped']*100:>4.0f}% {_fmt_p(r['wilcoxon_p_greater']):>8} "
                  f"{rb:>7} {r['auroc_wt_vs_adv']:>6.3f}")

    print("\n" + "=" * 92)
    print("CONDITIONAL: does confidence drop even where the STRUCTURE memorized?")
    print("  memorized = WT correct AND adv ligand RMSD < 2 A (ligand stayed in broken pocket)")
    print("  responded = WT correct AND adv ligand RMSD >= 4 A (model moved the ligand)")
    print("=" * 92)
    for model in models:
        rows = [r for r in cond_rows if r["model"] == model]
        if not rows:
            continue
        print(f"\n### {model}")
        print(f"{'subset':10} {'metric':6} {'n':>4} {'med_wt':>7} {'med_adv':>7} "
              f"{'medΔ':>7} {'wilcox_p':>8} {'rankbis':>7} {'AUROC':>6} {'ρ(Δc,Δrmsd)':>11}")
        for r in rows:
            rho = f"{r['spearman_dconf_drmsd']:+.3f}" if r["spearman_dconf_drmsd"] == r["spearman_dconf_drmsd"] else "    -"
            rb = f"{r['rank_biserial']:+.3f}" if r["rank_biserial"] == r["rank_biserial"] else "  n/a"
            print(f"{r['subset']:10} {r['metric']:6} {r['n']:>4} {r['median_wt']:>7.3f} "
                  f"{r['median_adv']:>7.3f} {r['median_delta']:>+7.3f} "
                  f"{_fmt_p(r['wilcoxon_p_greater']):>8} {rb:>7} {r['auroc_wt_vs_adv']:>6.3f} {rho:>11}")
    print()


def _write_csv(path, rows, cols):
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"[wrote] {path}  ({len(rows)} rows)")


if __name__ == "__main__":
    main()
