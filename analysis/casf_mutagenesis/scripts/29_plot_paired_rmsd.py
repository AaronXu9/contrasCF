#!/usr/bin/env python
"""Per-datapoint WT-vs-mutant ligand RMSD, one panel per method.

Six methods: SurfDock, UniDock2, ICM, GNINA (docking) then AF3+MSA, Boltz-2
(co-folding), so the family boundary is a single vertical cut in the figure.

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
# LINEAR axes with a BROKEN y. Constraints that force this shape:
#   * adversarial RMSD reaches 91.6 A but is 90% under 20 A, so one linear
#     0-92 axis squashes the decision region (0-4 A) into ~4% of the height
#     and hides the 2 A threshold that defines the metric;
#   * clipping at 14 A (the first version of this figure) silently dropped
#     16.2% of points and piled them into a false band on the top edge;
#   * a log axis fixes both but compresses the high end, and readers think
#     about RMSD linearly -- 3 vs 30 A should not look like a small step.
# So: a tall linear 0-20 panel carrying the 2 A line and 90% of the data, and
# a short compressed linear 20-95 panel above it for the remaining 10%. Both
# segments are linear, nothing is clipped in y.
Y_SPLIT = 20.0        # break point
Y_TOP = 95.0          # above the 91.6 A max
X_MAX = 15.0          # covers 97% of WT; the rest sits in the discarded band


def load_pairs() -> dict[str, dict[str, list[tuple[float, float]]]]:
    """{method: {variant: [(wt_rmsd, adv_rmsd), ...]}} — top-1 pose throughout."""
    out: dict[str, dict[str, list[tuple[float, float]]]] = defaultdict(
        lambda: defaultdict(list))

    # Docking engines: already top-1 (rank-1 pose only). ICM lives in its own
    # file (30_analyze_icm.py) but writes the same schema, so it merges here and
    # is scored by the same matcher as GNINA / UniDock2 / SurfDock.
    dock: dict[tuple[str, str], dict[str, float]] = defaultdict(dict)
    for fname in ("docking_results.csv", "icm_results.csv"):
        path = OUT / fname
        if not path.exists():
            continue
        with path.open() as f:
            for r in csv.DictReader(f):
                if r.get("module") != "casf" or r["status"] != "ok" or not r["rmsd_a"]:
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


def panel(ax_hi, ax_lo, method: str, data: dict[str, list[tuple[float, float]]],
          variants: tuple[str, ...] = VARIANTS) -> None:
    """Draw one method across a broken y-axis: ax_lo = 0-20 A, ax_hi = 20-95 A.

    `variants` restricts which mutation cases are drawn.  With a single variant
    the panel carries one series, so the annotated rate is that variant's own
    WT-conditioned memorization rate and matches docking_memorization.csv /
    memorization_full.csv row-for-row.
    """
    for ax in (ax_hi, ax_lo):
        # Shade the region conditioning discards: the method failed on WT there.
        ax.axvspan(THR, X_MAX, color="0.90", zorder=0)
        ax.set_xlim(0, X_MAX)
        ax.grid(alpha=0.25, lw=0.5)
        ax.set_axisbelow(True)
    # Thresholds and the "mutation changed nothing" diagonal live in the lower
    # panel; above 20 A both are off-scale or meaningless.
    ax_lo.plot([0, Y_SPLIT], [0, Y_SPLIT], ls="--", lw=1.0, color="0.55", zorder=1)
    ax_lo.axhline(THR, lw=1.0, color="0.35", zorder=1)
    ax_lo.axvline(THR, lw=1.0, color="0.35", zorder=1)
    ax_hi.axvline(THR, lw=1.0, color="0.35", zorder=1)
    ax_lo.set_ylim(0, Y_SPLIT)
    ax_hi.set_ylim(Y_SPLIT, Y_TOP)
    ax_hi.set_yticks([40, 70])

    n_wt_ok = n_mem = 0
    for v in variants:
        pts = data.get(v, [])
        if not pts:
            continue
        # x beyond X_MAX is pinned to the edge -- those points are all inside the
        # discarded (grey) band, where exact position carries no information.
        x = np.clip([p[0] for p in pts], 0, X_MAX)
        y = np.array([p[1] for p in pts])
        for ax in (ax_hi, ax_lo):
            ax.scatter(x, y, s=11, c=VARIANT_COLOR[v], marker=VARIANT_MARKER[v],
                       alpha=0.55, linewidths=0.3, edgecolors="white",
                       label=VARIANT_LABEL[v], zorder=3)
        n_wt_ok += sum(1 for a, _ in pts if a < THR)
        n_mem += sum(1 for a, b in pts if a < THR and b < THR)

    # Break marks on the facing spines.
    ax_hi.spines["bottom"].set_visible(False)
    ax_lo.spines["top"].set_visible(False)
    ax_hi.tick_params(bottom=False, labelbottom=False)
    kw = dict(marker=[(-1, -0.4), (1, 0.4)], markersize=7, linestyle="none",
              color="0.35", mec="0.35", mew=1, clip_on=False)
    ax_hi.plot([0, 1], [0, 0], transform=ax_hi.transAxes, **kw)
    ax_lo.plot([0, 1], [1, 1], transform=ax_lo.transAxes, **kw)

    rate = n_mem / n_wt_ok if n_wt_ok else float("nan")
    ax_hi.text(0.03, 0.92, f"memorized {n_mem}/{n_wt_ok} = {rate:.3f}",
               transform=ax_hi.transAxes, va="top", ha="left", fontsize=8,
               bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="0.7", lw=0.6))
    ax_hi.set_title(method, fontsize=11, pad=8)


def build_figure(data, order, pretty, variants: tuple[str, ...], subtitle: str):
    ncol = len(order)
    fig, axes = plt.subplots(2, ncol, figsize=(2.9 * ncol, 4.6), squeeze=False,
                             sharex="col", gridspec_kw=dict(height_ratios=[1, 3.2],
                                                            hspace=0.06))
    for i, m in enumerate(order):
        panel(axes[0][i], axes[1][i], pretty.get(m, m), data[m], variants)
    for i in range(1, ncol):
        axes[0][i].tick_params(labelleft=False)
        axes[1][i].tick_params(labelleft=False)
    fig.supylabel("mutant-pocket ligand RMSD (Å)", fontsize=10, x=0.055)
    for ax in axes[1]:
        ax.set_xlabel("WT ligand RMSD (Å)")

    fig.suptitle(subtitle, fontsize=13)
    fig.tight_layout()
    # A single-variant figure has one series per panel, so it needs no legend --
    # the title names it (see the dataviz skill).
    y_txt = -0.155
    if len(variants) > 1:
        h, l = axes[1][0].get_legend_handles_labels()
        seen, hh, ll = set(), [], []
        for a, b in zip(h, l):
            if b not in seen:
                seen.add(b); hh.append(a); ll.append(b)
        fig.legend(hh, ll, loc="lower center", ncol=3, frameon=False,
                   fontsize=9, bbox_to_anchor=(0.5, -0.07))
    else:
        y_txt = -0.10
    fig.text(0.5, y_txt,
             "Each point is one system.  Under the horizontal line = still native on a "
             "destroyed pocket (MEMORIZED, the bad outcome); above it = ligand moved "
             "(desired).\nGrey band = the method failed on WT, so its mutant cell is "
             "uninformative — that is exactly what WT-conditioning removes.\n"
             "y-axis is broken at 20 Å (linear in both segments) so the 10% of points "
             "out to 92 Å are shown without squashing the 2 Å threshold; nothing in y "
             "is clipped.  x is pinned at 15 Å (3% of WT values, all inside the "
             "discarded grey band).",
             ha="center", va="bottom", fontsize=8.5, color="0.3")
    return fig


def main() -> int:
    data = load_pairs()
    # Docking first, then co-folding — so the family boundary is one vertical cut.
    order = [m for m in ("surfdock", "unidock2", "icm", "gnina", "AF3+MSA", "Boltz2")
             if m in data]
    pretty = {"surfdock": "SurfDock", "gnina": "GNINA", "unidock2": "UniDock2",
              "icm": "ICM", "AF3+MSA": "AF3+MSA", "Boltz2": "Boltz-2"}
    if not order:
        print("no data found", file=sys.stderr)
        return 1

    jobs = [(VARIANTS, "paired_rmsd_wt_vs_mutant.png",
             "Per-system ligand RMSD, wild-type vs mutated pocket  (top-1 pose)")]
    jobs += [((v,), f"paired_rmsd_{v}.png",
              f"Per-system ligand RMSD, wild-type vs {VARIANT_LABEL[v]}  (top-1 pose)")
             for v in VARIANTS]

    for variants, fname, subtitle in jobs:
        fig = build_figure(data, order, pretty, variants, subtitle)
        out = FIG_DIR / fname
        fig.savefig(out, dpi=170, bbox_inches="tight")
        plt.close(fig)
        n = sum(len(data[m].get(v, [])) for m in order for v in variants)
        print(f"Wrote {out}   ({n} points, variants={','.join(variants)})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
