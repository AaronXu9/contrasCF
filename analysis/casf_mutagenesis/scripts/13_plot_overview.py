"""Cross-method memorization plots for the CASF-mutagenesis study.

Reads aggregate CSVs from analysis/casf_mutagenesis/outputs/:
  - memorization_full.csv         (co-folding: Boltz-2, AF3, AF3+MSA — pocket mutation)
  - docking_memorization.csv      (GNINA, UniDock2 — pocket mutation AND ligand mutation)
  - paired_affinity_full.csv      (Boltz-2 affinity Δ — pocket mutation)

Produces a 4-panel overview:
  (a) Pocket-mutation: top-1 ligand RMSD < 2 Å memorization rate, by method × variant.
      WT bar = "this method places the native ligand correctly" (success ceiling).
      Adv bars = "this method places the ligand correctly even after the pocket is broken"
                 (lower = better physics; high = memorization).
  (b) Ligand-mutation: same axes, ligand-side variants (halogenation, charge swap, methylation).
      Only docking engines have results today (GNINA, UniDock2); co-folding rows absent.
  (c) Boltz-2 affinity Δaff distribution per pocket-mutation variant (full CASF).
      Recognized perturbations → right-shifted (positive); memorized → centered at 0.
  (d) Boltz-2 P(binder) Δ per pocket-mutation variant.
      Recognized → left-shifted (negative); memorized → centered at 0.

Output: analysis/casf_mutagenesis/figures/overview_full.png
"""
from __future__ import annotations
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
OUT = REPO_ROOT / "analysis/casf_mutagenesis/outputs"
FIG_DIR = REPO_ROOT / "analysis/casf_mutagenesis/figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ---- consistent palette for variants ----
POCKET_VARIANTS = ("wt", "rem", "pack", "inv")
POCKET_COLORS = {
    "wt":   "#888888",
    "rem":  "#4C72B0",
    "pack": "#DD8452",
    "inv":  "#55A868",
}
POCKET_LABELS = {
    "wt":   "wt",
    "rem":  "rem (→G)",
    "pack": "pack (→F)",
    "inv":  "inv (Miyata)",
}

# Ligand-side variant grouping (collapse halo_Br/Cl/F → "halo", chrg_* → "chrg",
# meth_1..5 → "meth"). Keeps the bar plot readable.
LIG_GROUPS = {
    "halo":  ("halo_F_1", "halo_Cl_1", "halo_Br_1"),
    "chrg-": ("chrg_neu_methyl", "chrg_neu_ethyl", "chrg_neu_propyl"),
    "chrg+": ("chrg_pos_1", "chrg_pos_2", "chrg_pos_3"),
    "meth":  ("meth_1", "meth_2", "meth_3", "meth_4", "meth_5"),
}
LIG_GROUP_ORDER = ("wt", "halo", "chrg-", "chrg+", "meth")
LIG_COLORS = {
    "wt":    "#888888",
    "halo":  "#4C72B0",
    "chrg-": "#DD8452",
    "chrg+": "#C44E52",
    "meth":  "#55A868",
}


def load_memorization_full() -> dict[tuple[str, str], dict]:
    """Return {(method, variant): {'n': int, 'rate_2A': float}}.

    Co-folding side: includes AF3 (subset20 only, low n), AF3+MSA, Boltz2.
    Skip "wt" rows because memorization_full only reports adversarial variants.
    We'll derive WT separately from results_full.csv.
    """
    path = OUT / "memorization_full.csv"
    out: dict[tuple[str, str], dict] = {}
    with path.open() as f:
        for r in csv.DictReader(f):
            out[(r["model"], r["variant"])] = {
                "n": int(r["n_total"]),
                "rate_2A": float(r["memorization_rate_2A"]),
            }
    return out


def load_wt_rates_from_results(path_csv: Path) -> dict[str, dict]:
    """Compute per-method WT 'memorization rate' (= success rate on WT) from
    results_full.csv. WT < 2 Å = model placed native ligand correctly."""
    by_model = defaultdict(list)
    with path_csv.open() as f:
        for r in csv.DictReader(f):
            if r["variant"] != "wt" or r["status"] != "ok":
                continue
            if int(r["pose_idx"]) != 0:
                continue   # rank-0 / top-1-by-confidence
            try:
                rmsd = float(r["ligand_rmsd_a"])
            except (ValueError, KeyError):
                continue
            by_model[r["model"]].append(rmsd)
    return {
        m: {"n": len(rmsds), "rate_2A": sum(1 for x in rmsds if x < 2.0) / len(rmsds)}
        for m, rmsds in by_model.items()
    }


def load_docking_memorization() -> dict[tuple[str, str, str], dict]:
    """Return {(module, engine, variant): {'n', 'rate_2A'}}."""
    path = OUT / "docking_memorization.csv"
    out: dict[tuple[str, str, str], dict] = {}
    with path.open() as f:
        for r in csv.DictReader(f):
            out[(r["module"], r["engine"], r["variant"])] = {
                "n": int(r["n"]),
                "rate_2A": float(r["<2_A"]),
            }
    return out


def load_paired_affinity(scope: str = "full") -> dict[str, list[dict]]:
    path = OUT / f"paired_affinity_{scope}.csv"
    by_v: dict[str, list[dict]] = {"rem": [], "pack": [], "inv": []}
    with path.open() as f:
        for r in csv.DictReader(f):
            if not r["delta_affinity"]:
                continue
            by_v[r["variant"]].append({
                "d_aff":  float(r["delta_affinity"]),
                "d_prob": float(r["delta_probability"]),
            })
    return by_v


# ---------------------------------------------------------------------------

def panel_a_pocket(ax) -> None:
    """Memorization rate (<2 Å), pocket mutation, all methods."""
    cofold = load_memorization_full()
    wts = load_wt_rates_from_results(OUT / "results_full.csv")
    dock = load_docking_memorization()

    # Methods in plot order; only include if data exists
    methods: list[tuple[str, str]] = []  # (label, key into rate-getter)
    if ("Boltz2", "rem") in cofold:
        methods.append(("Boltz-2", "cofold:Boltz2"))
    if ("AF3+MSA", "rem") in cofold:
        methods.append(("AF3+MSA", "cofold:AF3+MSA"))
    if ("AF3", "rem") in cofold and cofold[("AF3", "rem")]["n"] >= 50:
        methods.append(("AF3", "cofold:AF3"))
    if ("casf", "gnina", "rem") in dock:
        methods.append(("GNINA", "dock:gnina"))
    if ("casf", "unidock2", "rem") in dock:
        methods.append(("UniDock2", "dock:unidock2"))
    if ("casf", "surfdock", "rem") in dock:
        methods.append(("SurfDock", "dock:surfdock"))

    width = 0.18
    x = np.arange(len(methods))
    for i, variant in enumerate(POCKET_VARIANTS):
        rates = []
        ns = []
        for _, key in methods:
            if key.startswith("cofold:"):
                model = key.split(":", 1)[1]
                if variant == "wt":
                    r = wts.get(model)
                    rates.append(r["rate_2A"] if r else 0.0)
                    ns.append(r["n"] if r else 0)
                else:
                    r = cofold.get((model, variant))
                    rates.append(r["rate_2A"] if r else 0.0)
                    ns.append(r["n"] if r else 0)
            else:
                _, engine = key.split(":", 1)
                r = dock.get(("casf", engine, variant))
                rates.append(r["rate_2A"] if r else 0.0)
                ns.append(r["n"] if r else 0)
        bars = ax.bar(x + (i - 1.5) * width, rates, width,
                      label=POCKET_LABELS[variant], color=POCKET_COLORS[variant],
                      edgecolor="white", linewidth=0.5)
        # annotate WT bars with n
        if variant == "wt":
            for bar, n in zip(bars, ns):
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.015,
                        f"n={n}", ha="center", va="bottom", fontsize=7, color="grey")

    ax.set_xticks(x)
    ax.set_xticklabels([m[0] for m in methods], fontsize=9)
    ax.set_ylabel("Top-1 ligand RMSD < 2 Å rate")
    ax.set_title("(a) Pocket mutation — rate of placing ligand near native\n"
                 "WT bar = success ceiling; adversarial bars: low = recognized, high = memorized",
                 fontsize=10)
    ax.legend(fontsize=8, loc="upper right", ncol=2, framealpha=0.95)
    ax.set_ylim(0, 1)
    ax.grid(axis="y", alpha=0.3)


def panel_b_ligand(ax) -> None:
    """Same as panel a but for ligand-side variants (only docking engines have data)."""
    dock = load_docking_memorization()
    engines = ("gnina", "unidock2", "surfdock")
    LABELS = {"gnina": "GNINA", "unidock2": "UniDock2", "surfdock": "SurfDock"}
    methods = []
    for eng in engines:
        if ("ligand", eng, "wt") in dock:
            methods.append((LABELS[eng], eng))

    # collapse per ligand group: average rate, weighted by n
    def grouped_rate(eng: str, group_key: str) -> tuple[float, int]:
        if group_key == "wt":
            r = dock.get(("ligand", eng, "wt"))
            return (r["rate_2A"], r["n"]) if r else (0.0, 0)
        ns = []
        rates = []
        for v in LIG_GROUPS[group_key]:
            r = dock.get(("ligand", eng, v))
            if not r:
                continue
            ns.append(r["n"])
            rates.append(r["rate_2A"])
        if not ns:
            return (0.0, 0)
        weighted = sum(rates[i] * ns[i] for i in range(len(ns))) / sum(ns)
        return (weighted, sum(ns))

    width = 0.16
    x = np.arange(len(methods))
    for i, gkey in enumerate(LIG_GROUP_ORDER):
        rates = []
        ns = []
        for _, eng in methods:
            r, n = grouped_rate(eng, gkey)
            rates.append(r); ns.append(n)
        bars = ax.bar(x + (i - 2) * width, rates, width,
                      label=gkey, color=LIG_COLORS[gkey],
                      edgecolor="white", linewidth=0.5)
        if gkey == "wt":
            for bar, n in zip(bars, ns):
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.015,
                        f"n={n}", ha="center", va="bottom", fontsize=7, color="grey")

    ax.set_xticks(x)
    ax.set_xticklabels([m[0] for m in methods], fontsize=9)
    ax.set_ylabel("Top-1 ligand RMSD < 2 Å rate")
    ax.set_title("(b) Ligand mutation — rate of placing ligand near native\n"
                 "(co-folding rows absent — no Boltz-2/AF3 runs on ligand_mutagenesis yet)",
                 fontsize=10)
    ax.legend(fontsize=8, loc="upper right", ncol=2, framealpha=0.95)
    ax.set_ylim(0, 1)
    ax.grid(axis="y", alpha=0.3)


def panel_c_affinity_delta(ax) -> None:
    by_v = load_paired_affinity("full")
    bins = np.linspace(-2.5, 5, 30)
    for v in ("rem", "pack", "inv"):
        d = [x["d_aff"] for x in by_v[v]]
        ax.hist(d, bins=bins, alpha=0.55, color=POCKET_COLORS[v],
                label=f"{POCKET_LABELS[v]}  (n={len(d)}, median={np.median(d):+.2f})",
                edgecolor="white", linewidth=0.4)
    ax.axvline(0, color="k", lw=0.9, ls="--", label="memorized (Δ=0)")
    ax.axvline(1, color="r", lw=0.9, ls=":", alpha=0.7, label="10× weaker (Δ=+1)")
    ax.set_xlabel("Δ log[IC50]  (adversarial − WT)")
    ax.set_ylabel("count")
    ax.set_title("(c) Boltz-2 affinity Δ (full CASF, n=229 per variant)\n"
                 "← memorized          recognized →", fontsize=10)
    ax.legend(fontsize=7.5, loc="upper right")
    ax.grid(alpha=0.3)


def panel_d_probability_delta(ax) -> None:
    by_v = load_paired_affinity("full")
    bins = np.linspace(-1, 1, 21)
    for v in ("rem", "pack", "inv"):
        d = [x["d_prob"] for x in by_v[v]]
        ax.hist(d, bins=bins, alpha=0.55, color=POCKET_COLORS[v],
                label=f"{POCKET_LABELS[v]}  (median={np.median(d):+.3f})",
                edgecolor="white", linewidth=0.4)
    ax.axvline(0, color="k", lw=0.9, ls="--", label="memorized (Δ=0)")
    ax.axvline(-0.3, color="r", lw=0.9, ls=":", alpha=0.7, label="−0.3 (clear drop)")
    ax.set_xlabel("Δ P(binder)  (adversarial − WT)")
    ax.set_ylabel("count")
    ax.set_title("(d) Boltz-2 binding-probability Δ (full CASF)\n"
                 "recognized ←          → memorized", fontsize=10)
    ax.legend(fontsize=7.5, loc="upper left")
    ax.grid(alpha=0.3)


def main() -> int:
    fig, axes = plt.subplots(2, 2, figsize=(13, 11))
    panel_a_pocket(axes[0, 0])
    panel_b_ligand(axes[0, 1])
    panel_c_affinity_delta(axes[1, 0])
    panel_d_probability_delta(axes[1, 1])
    fig.suptitle("CASF-mutagenesis: cross-method memorization overview",
                 fontsize=13, y=0.995)
    fig.tight_layout()
    out = FIG_DIR / "overview_full.png"
    fig.savefig(out, dpi=170, bbox_inches="tight")
    print(f"Wrote {out}")
    fig.savefig(FIG_DIR / "overview_full.pdf", bbox_inches="tight")
    return 0


if __name__ == "__main__":
    sys.exit(main())
