#!/usr/bin/env python
"""Per-datapoint WT-vs-mutant ligand RMSD, one panel per method.

Answers "what happened to *this* system when its pocket was destroyed?", which
the aggregate bar chart (13_plot_overview.py) cannot show: a single memorization
rate is compatible with many different per-system distributions.

Encoding
--------
Each point is one system.  x = that method's WT ligand RMSD, y = its RMSD on the
mutated pocket.  The dashed diagonal is "the mutation changed nothing".

    x < 2 Å  and  y < 2 Å   ->  lower-left   MEMORIZED  (still native despite a
                                             destroyed pocket -- the bad outcome)
    x < 2 Å  and  y >= 2 Å  ->  upper-left   RESPONDED  (ligand moved -- desired)
    x >= 2 Å                ->  shaded band  the method failed on WT, so the
                                             mutant cell says nothing.  This is
                                             exactly the region WT-conditioning
                                             removes, so the conditioned rate is
                                             read off the UNSHADED column only.

Colors are the Okabe-Ito CVD-safe triple, validated with the dataviz skill's
validate_palette.js (all five checks PASS; the repo's older green/orange pair
failed protanopia separation at dE 4.5).  Marker shape repeats the variant
identity so the encoding is not colour-alone.

Run:
    source env/lab.sh
    $CONTRASCF_PY analysis/casf_mutagenesis/scripts/29_plot_paired_rmsd.py
"""
from __future__ import annotations

import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
OUT = REPO_ROOT / "analysis" / "casf_mutagenesis" / "outputs"
FIG_DIR = REPO_ROOT / "analysis" / "casf_mutagenesis" / "figures"

VARIANTS = ("rem", "pack", "inv")
VARIANT_LABEL = {"rem": "rem (→Gly)", "pack": "pack (→Phe)", "inv": "inv (Miyata)"}
# Okabe-Ito: blue / vermillion / bluish-green. Validated CVD-safe.
VARIANT_COLOR = {"rem": "#0072B2", "pack": "#D55E00", "inv": "#009E73"}
VARIANT_MARKER = {"rem": "o", "pack": "s", "inv": "^"}   # secondary encoding
THR = 2.0
# x is clipped hard at 4 A: WT RMSD is <2 for most systems (that is the point),
# so an equal-aspect 0-14 x-axis spends ~85% of its width on the uninformative
# shaded band. y keeps the full range because the mutant response is the signal.
CLIP_X = 4.0
CLIP_Y = 14.0


def load_pairs() -> dict[str, dict[str, list[tuple[float, float]]]]:
    """{method: {variant: [(wt_rmsd, adv_rmsd), ...]}} — top-1 pose throughout."""
    out: dict[str, dict[str, list[tuple[float, float]]]] = defaultdict(
        lambda: defaultdict(list))

    # Docking engines: docking_results.csv is already top-1 (rank-1 pose only).
    dock: dict[tuple[str, str], dict[str, float]] = defaultdict(dict)
    with (OUT / "docking_results.csv").open() as f:
        for r in csv.DictReader(f):
            if r["module"] != "casf" or r["status"] != "ok" or not r["rmsd_a"]:
                continue
            try:
                dock[(r["engine"], r["system"])][r["variant"]] = float(r["rmsd_a"])
            except ValueError:
                continue
    for (engine, _sys), d in dock.items():
        if "wt" not in d:
            continue
        for v in VARIANTS:
            if v in d:
                out[engine][v].append((d["wt"], d[v]))

    # Co-folding: results_full.csv, pose_idx 0 = top-1 by model confidence.
    cof: dict[tuple[str, str], dict[str, float]] = defaultdict(dict)
    with (OUT / "results_full.csv").open() as f:
        for r in csv.DictReader(f):
            if r["status"] != "ok" or r["pose_idx"] not in ("0", "0.0"):
                continue
            if not r["ligand_rmsd_a"]:
                continue
            try:
                cof[(r["model"], r["pdbid"])][r["variant"]] = float(r["ligand_rmsd_a"])
            except ValueError:
                continue
    for (model, _pdb), d in cof.items():
        if "wt" not in d:
            continue
        for v in VARIANTS:
            if v in d:
                out[model][v].append((d["wt"], d[v]))
    return out


def panel(ax, method: str, data: dict[str, list[tuple[float, float]]]) -> None:
    # Shade the region conditioning discards: the method failed on WT there.
    ax.axvspan(THR, CLIP_X, color="0.90", zorder=0)
    ax.plot([0, CLIP_Y], [0, CLIP_Y], ls="--", lw=1.0, color="0.55", zorder=1)
    ax.axhline(THR, lw=1.0, color="0.35", zorder=1)
    ax.axvline(THR, lw=1.0, color="0.35", zorder=1)

    n_wt_ok = n_mem = 0
    for v in VARIANTS:
        pts = data.get(v, [])
        if not pts:
            continue
        x = np.clip([p[0] for p in pts], 0, CLIP_X)
        y = np.clip([p[1] for p in pts], 0, CLIP_Y)
        ax.scatter(x, y, s=11, c=VARIANT_COLOR[v], marker=VARIANT_MARKER[v],
                   alpha=0.55, linewidths=0.3, edgecolors="white",
                   label=VARIANT_LABEL[v], zorder=3)
        n_wt_ok += sum(1 for a, _ in pts if a < THR)
        n_mem += sum(1 for a, b in pts if a < THR and b < THR)

    rate = n_mem / n_wt_ok if n_wt_ok else float("nan")
    ax.text(0.03, 0.97,
            f"memorized {n_mem}/{n_wt_ok} = {rate:.3f}",
            transform=ax.transAxes, va="top", ha="left", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="0.7", lw=0.6))
    ax.set_title(method, fontsize=11)
    ax.set_xlim(0, CLIP_X); ax.set_ylim(0, CLIP_Y)
    ax.grid(alpha=0.25, lw=0.5)
    ax.set_axisbelow(True)


def main() -> int:
    data = load_pairs()
    order = [m for m in ("SurfDock", "surfdock", "gnina", "unidock2", "AF3+MSA", "Boltz2")
             if m in data]
    pretty = {"surfdock": "SurfDock", "gnina": "GNINA", "unidock2": "UniDock2",
              "AF3+MSA": "AF3+MSA", "Boltz2": "Boltz-2"}
    if not order:
        print("no data found", file=sys.stderr)
        return 1

    ncol = len(order)
    fig, axes = plt.subplots(1, ncol, figsize=(2.9 * ncol, 4.2), squeeze=False,
                             sharey=True)
    for ax, m in zip(axes[0], order):
        panel(ax, pretty.get(m, m), data[m])
    axes[0][0].set_ylabel("mutant-pocket ligand RMSD (Å)")
    for ax in axes[0]:
        ax.set_xlabel("WT ligand RMSD (Å)")

    h, l = axes[0][0].get_legend_handles_labels()
    fig.suptitle("Per-system ligand RMSD, wild-type vs mutated pocket  (top-1 pose)",
                 fontsize=13)
    fig.tight_layout()
    # Legend then explainer, both BELOW the axes -- keeping the explainer out of
    # the suptitle's band, which collided when it sat at the top.
    fig.legend(h, l, loc="lower center", ncol=3, frameon=False,
               fontsize=9, bbox_to_anchor=(0.5, -0.07))
    fig.text(0.5, -0.155,
             "Each point is one system.  Under the horizontal line = still native on a "
             "destroyed pocket (MEMORIZED, the bad outcome); above it = ligand moved "
             "(desired).\nGrey band = the method failed on WT, so its mutant cell is "
             "uninformative — that is exactly what WT-conditioning removes.  "
             "x clipped at 4 Å, y at 14 Å.",
             ha="center", va="bottom", fontsize=8.5, color="0.3")
    out = FIG_DIR / "paired_rmsd_wt_vs_mutant.png"
    fig.savefig(out, dpi=170, bbox_inches="tight")
    print(f"Wrote {out}")
    for m in order:
        tot = sum(len(v) for v in data[m].values())
        print(f"  {pretty.get(m, m):9s} {tot:4d} paired points")
    return 0


if __name__ == "__main__":
    sys.exit(main())
