#!/usr/bin/env python
"""Iterate 4 models × 16 cases, compute RMSD/confidence/clashes, write results.csv."""
from __future__ import annotations
import os, sys
from pathlib import Path

# Ensure conda env's libstdc++ is loaded (system lib is older than RDKit requires).
_ENV_LIB = "/home/aoxu/miniconda3/envs/rdkit_env/lib"
if _ENV_LIB not in os.environ.get("LD_LIBRARY_PATH", ""):
    os.environ["LD_LIBRARY_PATH"] = _ENV_LIB + os.pathsep + os.environ.get("LD_LIBRARY_PATH", "")

import pandas as pd
from rdkit import RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
RDLogger.DisableLog("rdApp.*")

from config import CASES, MODELS, RESULTS_DIR  # noqa: E402
from pipeline import run_one  # noqa: E402


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    models = list(MODELS.keys())
    cases = list(CASES.keys())
    print(f"Running pipeline on {len(models)} models × {len(cases)} cases = {len(models)*len(cases)} cells")
    for model in models:
        for case in cases:
            row = run_one(model, case)
            rows.append(row)
            tag = "OK" if row.get("error") is None else "ERR"
            print(
                f"  [{tag}] {model:8s}/{case:22s}"
                f"  all={row['ligand_rmsd_all']:.3f}"
                f"  common={row['ligand_rmsd_common']:.3f}"
                f"  bestfit={row['ligand_rmsd_bestfit']:.3f}"
                f"  ca={row['ca_rmsd']:.3f}"
                f"  clashes={row['n_clashes']}"
                + (f"   err={row['error']}" if row.get("error") else "")
            )
    df = pd.DataFrame(rows)
    out = RESULTS_DIR / "results.csv"
    df.to_csv(out, index=False)
    print(f"\nWrote {out} ({len(df)} rows, {len(df.columns)} cols)")

    # --- verification ----------------------------------------------------
    ok_per_model = df[df["error"].isna()].groupby("model").size()
    print("\nCoverage:")
    for m, n in ok_per_model.items():
        print(f"  {m}: {n}/{len(cases)} OK")
    bad = df[df["error"].notna()]
    if len(bad):
        print("\nFailed cells:")
        for _, r in bad.iterrows():
            print(f"  {r['model']}/{r['case']}: {r['error']}")

    # WT sanity check: for bindingsite_wt, ligand_rmsd_all should equal ligand_rmsd_common
    # (full ATP vs full ATP). Also cross-check ligand_rmsd_bestfit against the
    # paper's quoted WT numbers (Masters et al. 2025, p. 2: AF3=0.2, RFAA=2.2,
    # Boltz and Chai intermediate).
    wt = df[df["case"] == "bindingsite_wt"]
    paper_wt = {"AF3": 0.2, "RFAA": 2.2, "Boltz": 1.0, "Chai": 1.8}   # approximate
    print("\nWT sanity (bindingsite_wt):")
    print("  pose RMSD (all vs common) should match within 1e-3 for WT:")
    for _, r in wt.iterrows():
        diff = abs(r["ligand_rmsd_all"] - r["ligand_rmsd_common"])
        tag = "✓" if diff < 1e-3 else "!"
        print(
            f"    {tag} {r['model']:8s}: all={r['ligand_rmsd_all']:.4f}  common={r['ligand_rmsd_common']:.4f}"
        )
    print("\n  bestfit (paper-style, optimal rigid-body fit) vs paper's quoted WT values:")
    for _, r in wt.iterrows():
        m = r["model"]
        quoted = paper_wt.get(m)
        bf = r["ligand_rmsd_bestfit"]
        if quoted is None:
            print(f"    - {m:8s}: bestfit={bf:.3f}  (no paper reference)")
            continue
        tag = "✓" if abs(bf - quoted) < 0.3 else "!"
        print(
            f"    {tag} {m:8s}: bestfit={bf:.3f} Å   paper≈{quoted} Å   (Δ={bf - quoted:+.3f})"
        )


if __name__ == "__main__":
    main()
