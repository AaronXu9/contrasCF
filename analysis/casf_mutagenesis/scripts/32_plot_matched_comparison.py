"""Like-for-like cross-method comparison on a COMMON system set.

Answers the fairness objection the per-method table cannot: every method in
that table has its own denominator (SurfDock 249, Boltz-2 229, ...), and
co-folding is additionally credited on systems where it never got the fold
right. Docking is *handed* a receptor; co-folding must *predict* one, so a
Ca filter is only definable for co-folding — meaning a Ca-conditioned table
is still not symmetric across arms.

The symmetric construction is a common-system intersection:

  A  systems every one of the 6 methods produced for all 4 variants   n=209
  B  A, plus BOTH co-folding models solved WT (ligand < 2 A AND Ca < 2 A)  n=96

In B every method is scored on the identical systems and the co-folding
models are guaranteed to have solved WT, so neither "different denominators"
nor "credited on systems it cannot solve" applies. Retention is then directly
comparable across all six.

Output: figures/matched_comparison.png
"""
from __future__ import annotations
import csv
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
O = REPO / "analysis/casf_mutagenesis/outputs"
FIG = REPO / "analysis/casf_mutagenesis/figures"
CASF = Path("/home/aoxu/projects/VLS-Benchmark-Dataset/data/pdbbind_cleansplit")
V = ("rem", "pack", "inv")
TEAL, CORAL, INK, MUTED = "#0E9C93", "#E4572E", "#141B24", "#697687"

METH = [("SurfDock", "d", "surfdock"), ("UniDock2", "d", "unidock2"),
        ("ICM", "d", "icm"), ("GNINA", "d", "gnina"),
        ("Boltz-2", "c", "Boltz2"), ("AF3+MSA", "c", "AF3+MSA")]


def load():
    built = {s for s in os.listdir(CASF / "raw")
             if all((O / s / v / "af3.json").exists() for v in ("wt",) + V)}
    cof = {}
    for r in csv.DictReader(open(O / "results_full.csv")):
        if r["pose_idx"] != "0" or r["status"] != "ok" or r["pdbid"] not in built:
            continue
        try:
            lig = float(r["ligand_rmsd_a"])
        except (TypeError, ValueError):
            continue
        try:
            ca = float(r["ca_rmsd_a"])
        except (TypeError, ValueError):
            ca = None
        cof[(r["model"], r["pdbid"], r["variant"])] = (lig, ca)
    dock = {}
    rows = [r for r in csv.DictReader(open(O / "docking_results.csv"))
            if r.get("module") == "casf"]
    rows += list(csv.DictReader(open(O / "icm_results.csv")))
    for r in rows:
        if r["status"] != "ok" or r["system"] not in built:
            continue
        try:
            dock[(r["engine"], r["system"], r["variant"])] = float(r["rmsd_a"])
        except (TypeError, ValueError):
            pass
    return built, cof, dock


def main() -> int:
    built, cof, dock = load()
    cov = {}
    for lbl, k, key in METH:
        d = dock if k == "d" else cof
        cov[lbl] = {s for s in built if all((key, s, v) in d for v in ("wt",) + V)}
    A = set.intersection(*cov.values())

    def cof_ok(m, s):
        r = cof.get((m, s, "wt"))
        return r is not None and r[0] < 2 and r[1] is not None and r[1] < 2
    B = {s for s in A if cof_ok("Boltz2", s) and cof_ok("AF3+MSA", s)}

    def retention(key, kind, S):
        d = dock if kind == "d" else cof
        out = []
        for v in V:
            xs = [(d[(key, s, v)][0] if kind == "c" else d[(key, s, v)])
                  for s in S if (key, s, v) in d]
            out.append((np.array(xs) < 2).mean() if xs else np.nan)
        return float(np.mean(out))

    labels = [m[0] for m in METH]
    rA = [retention(m[2], m[1], A) for m in METH]
    rB = [retention(m[2], m[1], B) for m in METH]
    order = np.argsort(rB)
    labels = [labels[i] for i in order]
    rA = [rA[i] for i in order]
    rB = [rB[i] for i in order]
    kinds = [METH[i][1] for i in order]

    y = np.arange(len(labels))[::-1]
    h = 0.35
    fig, ax = plt.subplots(figsize=(9.2, 4.5))
    fig.subplots_adjust(left=0.155, right=0.97, top=0.755, bottom=0.20)

    ax.barh(y + h / 1.85, rA, height=h, color="#C9D2DA", edgecolor=MUTED,
            linewidth=0.5, label="A · all-methods common set (n=209)", zorder=3)
    cols = [CORAL if k == "c" else TEAL for k in kinds]
    ax.barh(y - h / 1.85, rB, height=h, color=cols, zorder=3,
            label="B · + both co-folding models solved WT (n=96)")
    for yy, a, b in zip(y, rA, rB):
        ax.text(a + .006, yy + h / 1.85, f"{a:.3f}", va="center", fontsize=8.3,
                color=MUTED, family="monospace")
        ax.text(b + .006, yy - h / 1.85, f"{b:.3f}", va="center", fontsize=8.7,
                color=INK, family="monospace", fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9.5)
    ax.set_xlim(0, 0.50)
    ax.set_ylim(-0.6, len(labels) + 0.3)
    ax.set_xlabel("adversarial retention  (fraction still < 2 Å after the pocket is destroyed)",
                  labelpad=7)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.legend(loc="upper right", fontsize=8.4, frameon=False, bbox_to_anchor=(1.0, 1.03))

    fig.text(0.006, 0.965, "Like-for-like comparison on a common system set",
             fontsize=12.5, fontweight="bold", color=INK)
    fig.text(0.006, 0.905,
             "Every method scored on the IDENTICAL systems — removing both the different-denominator "
             "and the credited-on-unsolvable objections.",
             fontsize=8.6, color=MUTED)
    fig.text(0.006, 0.855, "teal = docking · coral = co-folding · lower = less memorisation.",
             fontsize=8.6, color=MUTED)
    fig.text(0.006, 0.018,
             "In set B the co-folding models are WT-correct by construction (ligand < 2 Å and Cα < 2 Å), so their "
             "WT rate is 1.000.\nTightening from A to B RAISES co-folding retention and leaves docking flat — "
             "the separation widens under the stricter, fairer test.",
             fontsize=7.4, color=MUTED, linespacing=1.5)

    out = FIG / "matched_comparison.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")
    print(f"  set A n={len(A)}   set B n={len(B)}")
    for l, a, b in zip(labels, rA, rB):
        print(f"    {l:<9} A {a:.3f}   B {b:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
