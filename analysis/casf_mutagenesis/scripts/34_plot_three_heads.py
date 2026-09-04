"""The three-head dissociation as ONE claim figure.

Replaces `conf_aff_rmsd_pocket.png` for presentation purposes. That figure is
numerically correct (re-verified 2026-09-03: n=189/166/332 and dAff
+0.439/+0.066/+0.068 reproduce exactly from the current results_full.csv) but it
is an *analysis dump* — four panels showing the work, none of which states the
conclusion. A reader has to assemble the claim themselves.

The claim is that Boltz-2's three outputs give contradictory verdicts on the
same perturbation. Carrying it needs TWO questions kept apart, because they have
different answers and the old figure conflated them:

  (a) DIRECTION — does the head move the right way when the pocket is destroyed?
      Common currency: AUROC (WT vs adversarial). All three beat chance. The
      affinity head scores 0.637, which is *why* a direction-only figure would
      understate the problem and a magnitude-only figure would overstate it.

  (b) MAGNITUDE — does it move FAR ENOUGH? Losing the binding pose demands +3 to
      +6 log units. The affinity head delivers +0.229, and under direct
      intervention (ligand ejected 35 A, zero protein contacts) it delivers
      -0.004 — while a pure-physics scorer loses 100% of its binding energy.

So: the affinity head points the right way and moves ~1/20th of the required
distance; under intervention it does not move at all.

Panel (a) is computed live from results_full.csv. Panel (b)'s strata are live;
the two pose-swap constants are published values from docs/casf_confidence.md
(the raw poseswap_*.csv are not on this host) and are labelled as such.

Colors are the Okabe-Ito triple already validated for this repo in
29_plot_paired_rmsd.py (all five dataviz checks PASS). Assignment is semantic:
green = the head that works, vermillion = the head that fails, blue = the head
the rest of the study measures.

Run:
    source env/lab.sh
    $CONTRASCF_PY analysis/casf_mutagenesis/scripts/34_plot_three_heads.py
"""
from __future__ import annotations

import csv
import os
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import roc_auc_score

REPO = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
OUT = REPO / "analysis/casf_mutagenesis/outputs"
FIG = REPO / "analysis/casf_mutagenesis/figures"

VARIANTS = ("rem", "pack", "inv")
# Okabe-Ito, validated. Semantic: works / fails / measured-elsewhere.
C_STRUCT, C_CONF, C_AFF = "#0072B2", "#009E73", "#D55E00"
INK, MUTED, HAIR = "#141B24", "#697687", "#E1E6EC"

# Published in docs/casf_confidence.md "pose-swap test" (n=29 systems, ligand
# ejected to a median 35 A from the protein — zero contacts). Raw per-system
# CSVs live on CARC, not here, so these are quoted, not recomputed.
POSESWAP_AFF_DELTA = -0.004      # median native->eject gap, log units
POSESWAP_N = 29
POSESWAP_VINA = (-8.7, 0.0)      # GNINA Vina kcal/mol, same ejection
# What losing the binding pose demands, in log units (biophysics).
PHYSICS_FLOOR = (3.0, 6.0)


def load(model: str):
    """{(sysid, variant): (best_of_5_rmsd, row)} plus per-system affinity."""
    best: dict[tuple[str, str], tuple[float, dict]] = defaultdict(
        lambda: (float("inf"), None))
    aff: dict[tuple[str, str], float] = {}
    with (OUT / "results_full.csv").open() as f:
        for r in csv.DictReader(f):
            if r["model"] != model or r["status"] != "ok":
                continue
            try:
                rm = float(r["ligand_rmsd_a"])
            except (TypeError, ValueError):
                continue
            k = (r["pdbid"], r["variant"])
            if rm < best[k][0]:
                best[k] = (rm, r)
            if r.get("affinity_pred_value"):
                try:
                    aff[k] = float(r["affinity_pred_value"])
                except ValueError:
                    pass
    return best, aff


def _fnum(row, col):
    try:
        return float(row[col]) if row and row.get(col) else None
    except ValueError:
        return None


def head_aurocs(model: str):
    """AUROC per head for separating WT from adversarial cells, WT-conditioned.

    WT-conditioned because an adversarial cell on a system the model cannot
    solve at WT says nothing (see the skill's conditioning rule). Signs are set
    so that a HIGHER score always means 'noticed the damage'.
    """
    best, aff = load(model)
    sysids = {s for (s, _v) in best}
    wt_ok = {s for s in sysids
             if (s, "wt") in best and best[(s, "wt")][0] < 2.0}
    pairs = {"structure": [], "confidence": [], "affinity": []}
    for s in wt_ok:
        wt_rm, wt_row = best[(s, "wt")]
        for v in VARIANTS:
            if (s, v) not in best:
                continue
            adv_rm, adv_row = best[(s, v)]
            pairs["structure"].append((wt_rm, adv_rm))
            cw, ca = _fnum(wt_row, "iptm"), _fnum(adv_row, "iptm")
            if cw is not None and ca is not None:
                pairs["confidence"].append((cw, ca))
            aw, aa = aff.get((s, "wt")), aff.get((s, v))
            if aw is not None and aa is not None:
                pairs["affinity"].append((aw, aa))

    sign = {"structure": 1, "confidence": -1, "affinity": 1}
    out = {}
    for head, ps in pairs.items():
        if not ps:
            continue
        w = np.array([p[0] for p in ps])
        a = np.array([p[1] for p in ps])
        y = np.r_[np.zeros(len(w)), np.ones(len(a))]
        out[head] = (roc_auc_score(y, sign[head] * np.r_[w, a]), len(ps))
    return out, len(wt_ok)


def responded_daff(model: str = "Boltz2") -> tuple[float, int]:
    """Median dAffinity on WT-conditioned RESPONDED cells (adv RMSD >= 4 A).

    The responded stratum is the decisive one: the ligand has been ejected past
    any plausible binding pose, so it is the clearest possible non-binder.
    """
    best, aff = load(model)
    sysids = {s for (s, _v) in best}
    wt_ok = {s for s in sysids
             if (s, "wt") in best and best[(s, "wt")][0] < 2.0}
    vals = []
    for s in wt_ok:
        if (s, "wt") not in aff:
            continue
        for v in VARIANTS:
            if (s, v) not in best or (s, v) not in aff:
                continue
            if best[(s, v)][0] >= 4.0:
                vals.append(aff[(s, v)] - aff[(s, "wt")])
    return (float(np.median(vals)) if vals else float("nan")), len(vals)


def main() -> int:
    au_b, n_b = head_aurocs("Boltz2")
    au_a, n_a = head_aurocs("AF3+MSA")
    d_resp, n_resp = responded_daff("Boltz2")

    fig = plt.figure(figsize=(12.4, 4.9))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.15],
                          left=0.085, right=0.985, top=0.735, bottom=0.185,
                          wspace=0.30)
    axA, axB = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])

    # ---------- (a) DIRECTION ----------
    heads = [("structure", C_STRUCT), ("confidence", C_CONF), ("affinity", C_AFF)]
    y = np.arange(len(heads))[::-1]
    h = 0.34
    for i, (head, col) in enumerate(heads):
        yy = y[i]
        b = au_b.get(head)
        if b:
            axA.barh(yy + h / 1.85, b[0] - 0.5, left=0.5, height=h, color=col,
                     zorder=3)
            axA.text(b[0] + .006, yy + h / 1.85, f"{b[0]:.3f}", va="center",
                     fontsize=8.8, color=INK, family="monospace",
                     fontweight="bold")
        a = au_a.get(head)
        if a:
            axA.barh(yy - h / 1.85, a[0] - 0.5, left=0.5, height=h, color=col,
                     alpha=0.42, zorder=3)
            axA.text(a[0] + .006, yy - h / 1.85, f"{a[0]:.3f}", va="center",
                     fontsize=8.4, color=MUTED, family="monospace")
        else:
            axA.text(0.515, yy - h / 1.85, "AF3 has no affinity head",
                     va="center", fontsize=7.8, color=MUTED, style="italic")
    axA.axvline(0.5, color=INK, lw=1.1, zorder=4)
    axA.text(0.5, len(heads) - 0.42, "chance", ha="center", va="bottom",
             fontsize=7.8, color=MUTED)
    axA.set_yticks(y)
    axA.set_yticklabels([h[0] for h in heads], fontsize=10.5)
    axA.set_xlim(0.5, 1.0)
    axA.set_ylim(-0.62, len(heads) - 0.2)
    axA.set_xlabel("AUROC — WT vs destroyed pocket", labelpad=6)
    axA.set_title("(a) Direction — does the head move the right way?",
                  fontsize=10.5, loc="left", pad=9, color=INK)
    for sp in ("top", "right", "left"):
        axA.spines[sp].set_visible(False)
    axA.tick_params(axis="y", length=0)
    axA.text(0.0, -0.215, "solid = Boltz-2 · pale = AF3+MSA.   "
             "All three beat chance — direction is NOT the problem.",
             transform=axA.transAxes, fontsize=8.1, color=MUTED)

    # ---------- (b) MAGNITUDE ----------
    axB.axvspan(*PHYSICS_FLOOR, color=C_STRUCT, alpha=0.13, zorder=0)
    axB.text(np.mean(PHYSICS_FLOOR), 2.63,
             "what losing the pose DEMANDS\n+3 to +6 log units",
             ha="center", va="center", fontsize=8.6, color=C_STRUCT,
             fontweight="bold", linespacing=1.35)

    bars = [
        (f"observed\nresponded cells (n={n_resp})", d_resp, C_AFF),
        (f"under intervention\nligand ejected 35 Å (n={POSESWAP_N})",
         POSESWAP_AFF_DELTA, C_AFF),
    ]
    yb = np.array([1.35, 0.35])
    for (lbl, val, col), yy in zip(bars, yb):
        axB.barh(yy, max(val, 0.0), height=0.42, color=col, zorder=3)
        axB.text(max(val, 0.0) + 0.09, yy, f"{val:+.3f}", va="center",
                 fontsize=9.4, color=INK, family="monospace", fontweight="bold")
    axB.set_yticks(yb)
    axB.set_yticklabels([b[0] for b in bars], fontsize=8.8)
    axB.set_xlim(-0.15, 6.6)
    axB.set_ylim(-0.35, 3.05)
    axB.axvline(0, color=MUTED, lw=0.9)
    axB.set_xlabel("Δ predicted log[IC50]   (+ = predicted weaker binding)",
                   labelpad=6)
    axB.set_title("(b) Magnitude — does it move FAR ENOUGH?",
                  fontsize=10.5, loc="left", pad=9, color=INK)
    for sp in ("top", "right", "left"):
        axB.spines[sp].set_visible(False)
    axB.tick_params(axis="y", length=0)
    axB.text(0.0, -0.215,
             f"Pure-physics reference, same 35 Å ejection: GNINA Vina "
             f"{POSESWAP_VINA[0]} → {POSESWAP_VINA[1]:.1f} kcal/mol\n"
             "— 100% of binding energy lost.",
             transform=axB.transAxes, fontsize=8.1, color=MUTED,
             linespacing=1.45, va="top")

    # ---------- framing ----------
    fig.text(0.006, 0.955,
             "The affinity head points the right way — and moves ~1/20th of the "
             "distance physics requires",
             fontsize=13, fontweight="bold", color=INK)
    fig.text(0.006, 0.885,
             "Boltz-2 emits three outputs per prediction. Under the same destroyed pocket they give "
             "contradictory verdicts: confidence registers the damage, the affinity head does not.",
             fontsize=8.9, color=MUTED)
    fig.text(0.006, 0.828,
             f"WT-conditioned throughout (Boltz-2 n={n_b} systems, AF3+MSA n={n_a}); "
             "best-of-5 diffusion samples. Pose-swap values quoted from docs/casf_confidence.md.",
             fontsize=7.9, color=MUTED)

    out = FIG / "three_heads_dissociate.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")
    print(f"  (a) Boltz-2  " + "  ".join(
        f"{k} {v[0]:.3f} (n={v[1]})" for k, v in au_b.items()))
    print(f"  (a) AF3+MSA  " + "  ".join(
        f"{k} {v[0]:.3f} (n={v[1]})" for k, v in au_a.items()))
    print(f"  (b) responded dAff {d_resp:+.3f} (n={n_resp}); "
          f"pose-swap {POSESWAP_AFF_DELTA:+.3f} (n={POSESWAP_N})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
