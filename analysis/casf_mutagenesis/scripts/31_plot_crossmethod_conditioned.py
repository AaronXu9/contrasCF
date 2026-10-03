"""Cross-method WT-conditioned memorization figure, including ICM.

Renews the panel behind doc-2 §4 with ICM added (scored 2026-09-02 by
`30_analyze_icm.py`). Two panels because one number cannot carry the claim:

  (a) WT ceiling vs WT-CONDITIONED adversarial retention, per method.
      Retention = of the systems this method solved at WT, the fraction it
      STILL places within 2 A after the pocket is destroyed. That is the
      memorization measure.

  (b) The same data as the "gap" (WT - mean retention) used in doc 2, shown
      only to make its weakness visible: gap mixes the WT ceiling into a
      memorization score, so a method with a mediocre ceiling is pushed down
      the ranking even when its retention is low. ICM is the clear case —
      retention 0.139 (docking-like, good) but gap +0.546 (below AF3+MSA's
      +0.631, misleading).

Numbers for the five established methods are doc-2's (regenerated
2026-08-28); ICM is computed here from outputs/icm_results.csv.
"""
from __future__ import annotations
import csv
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
OUT = REPO_ROOT / "analysis/casf_mutagenesis/outputs"
FIG = REPO_ROOT / "analysis/casf_mutagenesis/figures"
FIG.mkdir(parents=True, exist_ok=True)

# doc-2 §4 table, regenerated 2026-08-28 (WT-conditioned adversarial columns)
METHODS = [
    # label,        wt,    wt_correct,  rem,   pack,  inv,   kind
    ("SurfDock",   0.876, "218/249",   0.024, 0.043, 0.030, "docking"),
    ("UniDock2",   0.578, "145/251",   0.141, 0.141, 0.104, "docking"),
    ("GNINA",      0.729, "183/251",   0.180, 0.157, 0.151, "docking"),
    ("ICM",        None,  None,        None,  None,  None,  "docking"),   # filled below
    ("AF3+MSA",    0.824, "196/238",   0.439, 0.367, 0.301, "cofolding"),
    ("Boltz-2",    0.594, "136/229",   0.368, 0.360, 0.257, "cofolding"),
]

TEAL, CORAL, INK, MUTED, HAIR = "#0E9C93", "#E4572E", "#141B24", "#697687", "#E1E6EC"


def icm_row():
    rows = [r for r in csv.DictReader(open(OUT / "icm_results.csv"))
            if r["status"] == "ok" and r["rmsd_a"]]
    by = {(r["system"], r["variant"]): float(r["rmsd_a"]) for r in rows}
    wt_all = {s for (s, v) in by if v == "wt"}
    wt_ok = {s for (s, v), x in by.items() if v == "wt" and x < 2.0}
    cond = {}
    for v in ("rem", "pack", "inv"):
        c = [x for (s, vv), x in by.items() if vv == v and s in wt_ok]
        cond[v] = float((np.array(c) < 2).mean())
    return (len(wt_ok) / len(wt_all), f"{len(wt_ok)}/{len(wt_all)}",
            cond["rem"], cond["pack"], cond["inv"])


def main() -> int:
    wt, wtc, rem, pack, inv = icm_row()
    ms = [(m[0], wt, wtc, rem, pack, inv, m[6]) if m[0] == "ICM" else m for m in METHODS]
    ms.sort(key=lambda m: np.mean([m[3], m[4], m[5]]))     # by retention, low = best

    labels = [m[0] for m in ms]
    wts = np.array([m[1] for m in ms])
    ret = np.array([np.mean([m[3], m[4], m[5]]) for m in ms])
    kinds = [m[6] for m in ms]
    gaps = wts - ret
    y = np.arange(len(ms))[::-1]

    fig, ax1 = plt.subplots(figsize=(9.4, 4.4))
    fig.subplots_adjust(left=0.175, right=0.975, top=0.78, bottom=0.255)

    h = 0.34
    ax1.barh(y + h / 1.8, wts, height=h, color=HAIR, edgecolor=MUTED,
             linewidth=0.6, label="WT ceiling (unconditioned)", zorder=3)
    cols = [CORAL if k == "cofolding" else TEAL for k in kinds]
    ax1.barh(y - h / 1.8, ret, height=h, color=cols, zorder=3,
             label="adversarial retention (WT-conditioned)")
    for yy, w, r in zip(y, wts, ret):
        ax1.text(w + .012, yy + h / 1.8, f"{w:.3f}", va="center", fontsize=8.4,
                 color=MUTED, family="monospace")
        ax1.text(r + .012, yy - h / 1.8, f"{r:.3f}", va="center", fontsize=8.8,
                 color=INK, family="monospace", fontweight="bold")
    ax1.set_yticks(y)
    ax1.set_yticklabels([f"{m[0]}\n{m[2]}" for m in ms], fontsize=9)
    ax1.set_xlim(0, 1.02)
    ax1.set_ylim(-0.6, len(ms) + 0.35)
    ax1.set_xlabel("top-1 ligand RMSD < 2 Å rate", labelpad=7)
    for sp in ("top", "right", "left"):
        ax1.spines[sp].set_visible(False)
    ax1.tick_params(axis="y", length=0)
    ax1.legend(loc="upper right", fontsize=8.4, frameon=False,
               bbox_to_anchor=(1.0, 1.02))

    fig.text(0.006, 0.965, "Cross-method memorization, WT-conditioned (full CASF)",
             fontsize=12.5, fontweight="bold", color=INK)
    fig.text(0.006, 0.905,
             "Retention = of the systems a method solved at WT, the fraction still within 2 Å "
             "after the pocket is destroyed.",
             fontsize=8.6, color=MUTED)
    fig.text(0.006, 0.855,
             "Lower retention = less memorization.  teal = docking · coral = co-folding.",
             fontsize=8.6, color=MUTED)
    fig.text(0.006, 0.015,
             "Five methods from the 2026-08-28 regeneration; ICM scored 2026-09-02 with the same "
             "matcher (30_analyze_icm.py). AF3 without MSA omitted: 0/19 WT-correct → conditional undefined.\n"
             "Retention is shown instead of the 'gap' (WT − retention): gap folds the WT ceiling into a "
             "memorization score, and doc-2's gap column is not computed consistently across rows.",
             fontsize=7.2, color=MUTED, linespacing=1.5)
    out = FIG / "crossmethod_conditioned.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")
    print(f"  ICM row: WT {wt:.3f} ({wtc})  rem {rem:.3f} pack {pack:.3f} inv {inv:.3f} "
          f"retention {np.mean([rem,pack,inv]):.3f}  gap +{wt-np.mean([rem,pack,inv]):.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
