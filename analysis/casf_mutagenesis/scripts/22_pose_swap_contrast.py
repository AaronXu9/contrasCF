#!/usr/bin/env python3
"""Steep-vs-flat contrast: Boltz-2's affinity head (flat under ligand ejection)
vs GNINA (genuine-physics scorer whose Vina term collapses). Reads the Boltz
panel (19/20 output) + the GNINA panel (21 output), joins by system, renders the
two-panel figure, and prints the summary that lands the conclusion.

Run (protenix env): .../protenix/bin/python 22_pose_swap_contrast.py \
    --boltz /tmp/poseswap_panel2 --gnina /tmp/gnina_ps --out /tmp/poseswap_panel2
"""
from __future__ import annotations
import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd

CLEAR = [0.0, 5.0, 15.0, 30.0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--boltz", default="/tmp/poseswap_panel2")
    ap.add_argument("--gnina", default="/tmp/gnina_ps")
    ap.add_argument("--out", default="/tmp/poseswap_panel2")
    args = ap.parse_args()

    boltz = {}
    for f in glob.glob(f"{args.boltz}/*/poseswap_*.csv"):
        df = pd.read_csv(f)
        lad = df[df.axis.isin(["none", "radial"])].sort_values("disp_A")
        boltz[df.tag.iloc[0]] = dict(zip(lad.disp_A, lad.boltz_aff))
    gnina = {}
    for f in glob.glob(f"{args.gnina}/*/gnina_*.csv"):
        df = pd.read_csv(f).sort_values("clearance")
        gnina[df.system.iloc[0]] = df.set_index("clearance")

    systems = sorted(set(boltz) & set(gnina))
    rows = []
    for s in systems:
        b, gz = boltz[s], gnina[s]
        if not all(c in b for c in CLEAR) or not all(c in gz.index for c in CLEAR):
            continue
        rows.append(dict(
            system=s,
            boltz_gap=b[30.0] - b[0.0],
            vina_native=gz.loc[0.0, "vina"], vina_eject=gz.loc[30.0, "vina"],
            vina_loss=gz.loc[0.0, "vina"] - gz.loc[30.0, "vina"],   # native is negative; loss = |binding| lost
            cnnaff_drop=gz.loc[0.0, "cnnaffinity"] - gz.loc[30.0, "cnnaffinity"],
            far_mindist=gz.loc[30.0, "min_lig_prot"],
        ))
    C = pd.DataFrame(rows).set_index("system")
    print(f"\n=== POSE-SWAP CONTRAST  (n={len(C)} systems with both Boltz + GNINA) ===")
    print(C.round(3).to_string())
    print("\n--- summary (ligand ejected ~30+ A into solvent) ---")
    print(f"  Boltz-2 affinity head: median Δlog[IC50] = {C.boltz_gap.median():+.3f}  (flat; physics wants +3..+6)")
    print(f"  GNINA Vina (physics):  median native = {C.vina_native.median():.2f} kcal/mol -> median ejected = "
          f"{C.vina_eject.median():.2f}  (collapses; median binding-energy lost = {C.vina_loss.median():.2f} kcal/mol)")
    print(f"  GNINA CNNaffinity:     median drop = {C.cnnaff_drop.median():+.2f} pK  (CNN has its own memorization floor)")
    print(f"  fraction GNINA Vina -> ~0 (|eject|<0.5): {(C.vina_eject.abs() < 0.5).mean():.0%}")

    C.round(4).to_csv(Path(args.out) / "pose_swap_contrast.csv")

    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    grid = np.array(CLEAR)
    fig, ax = plt.subplots(1, 2, figsize=(12, 5), sharex=True)
    # (a) Boltz affinity Δ
    for s in C.index:
        y = [boltz[s][c] - boltz[s][0.0] for c in CLEAR]
        ax[0].plot(grid, y, color="#888", alpha=.5, lw=1)
    medb = [np.median([boltz[s][c] - boltz[s][0.0] for s in C.index]) for c in CLEAR]
    ax[0].plot(grid, medb, color="#36c", lw=3, marker="o", label="median")
    ax[0].axhline(0, color="k", lw=.6); ax[0].axhline(3, color="green", ls=":", lw=1, label="physics floor (+3)")
    ax[0].set_ylim(-0.5, 3.4)
    ax[0].set_xlabel("ligand ejection — clearance beyond protein surface (Å)")
    ax[0].set_ylabel("Δ predicted log[IC50] (+ = weaker)")
    ax[0].set_title(f"Boltz-2 affinity head — FLAT\nmedian gap = {C.boltz_gap.median():+.2f} log units"); ax[0].legend(fontsize=8)
    # (b) GNINA Vina
    for s in C.index:
        y = [gnina[s].loc[c, "vina"] for c in CLEAR]
        ax[1].plot(grid, y, color="#888", alpha=.5, lw=1)
    medg = [np.median([gnina[s].loc[c, "vina"] for s in C.index]) for c in CLEAR]
    ax[1].plot(grid, medg, color="#c00", lw=3, marker="o", label="median")
    ax[1].axhline(0, color="k", lw=.6)
    ax[1].set_xlabel("ligand ejection — clearance beyond protein surface (Å)")
    ax[1].set_ylabel("GNINA Vina binding energy (kcal/mol)")
    ax[1].set_title(f"GNINA physics term — COLLAPSES\nmedian {C.vina_native.median():.1f} → {C.vina_eject.median():.1f} kcal/mol"); ax[1].legend(fontsize=8)
    fig.suptitle("Pose-swap: a real physics scorer (GNINA) loses all binding energy when the ligand is ejected — Boltz-2's affinity head does not notice", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, .95])
    out = Path(args.out) / "pose_swap_contrast.png"
    fig.savefig(out, dpi=130)
    print(f"[wrote] {out}\n[wrote] {Path(args.out)/'pose_swap_contrast.csv'}")


if __name__ == "__main__":
    main()
