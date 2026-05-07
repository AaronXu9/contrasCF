"""Analyze Boltz-2 + AF3 predictions on the CASF subset20.

For every (pdbid, variant, model) cell with a predicted mmCIF, compute:
  - Cα RMSD vs the crystal protein (sanity: should be small for WT, larger
    for adversarial variants where AF3 may misfold the mutated pocket).
  - Ligand heavy-atom RMSD vs the crystal ligand pose, after Cα superposition.

Aggregate memorization rate per (model, variant) — fraction of *adversarial*
cases (rem / pack / inv) with ligand RMSD < 2 Å (paper's threshold).

Run:
    LD_LIBRARY_PATH=/home/aoxu/miniconda3/envs/rdkit_env/lib:$LD_LIBRARY_PATH \\
        /home/aoxu/miniconda3/envs/rdkit_env/bin/python \\
        analysis/casf_mutagenesis/scripts/05_analyze_subset20.py
"""
from __future__ import annotations
import csv
import json
import sys
from dataclasses import asdict
from pathlib import Path

REPO_ROOT = Path("/mnt/katritch_lab2/aoxu/contrasCF")
sys.path.insert(0, str(REPO_ROOT / "analysis"))

from casf_mutagenesis.analysis import (
    MEMORIZATION_THRESHOLDS_A, MODEL_FILES, analyze_prediction,
    memorization_stats,
)
from casf_mutagenesis.config import OUTPUT_ROOT, SUBSET20_JSON, VARIANTS


def main() -> int:
    ids = json.loads(SUBSET20_JSON.read_text())["casf2016"]
    rows = []
    print(f"subset20: {len(ids)} systems × {len(VARIANTS)} variants × "
          f"{len(MODEL_FILES)} models = "
          f"{len(ids) * len(VARIANTS) * len(MODEL_FILES)} cells\n")
    for pdbid in ids:
        for variant in VARIANTS:
            for model in MODEL_FILES:
                rec = analyze_prediction(pdbid, variant, model)
                rows.append(rec)
                if rec.status == "ok":
                    print(
                        f"  {pdbid:6s} {variant:5s} {model:7s} "
                        f"lig_rmsd={rec.ligand_rmsd_a:6.2f} Å "
                        f"ca_rmsd={rec.ca_rmsd_a:5.2f} Å  "
                        f"({rec.n_ca_paired} Cα; "
                        f"{rec.n_heavy_matched}/{rec.n_heavy_native} heavy)"
                    )
                else:
                    print(f"  {pdbid:6s} {variant:5s} {model:7s} {rec.status}"
                          f"{(': ' + rec.error) if rec.error else ''}")

    # CSV
    csv_path = OUTPUT_ROOT / "results_subset20.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(asdict(rows[0]).keys()))
        w.writeheader()
        for r in rows:
            w.writerow(asdict(r))
    print(f"\nPer-cell results: {csv_path}")

    # Memorization aggregates
    stats = memorization_stats(rows)
    summary_path = OUTPUT_ROOT / "memorization_subset20.csv"
    with summary_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "model", "variant", "n_total",
            "memorization_rate_2A", "memorization_rate_4A", "median_rmsd_A",
        ])
        for (model, variant), s in sorted(stats.items()):
            w.writerow([
                s.model, s.variant, s.n_total,
                f"{s.rate(2.0):.3f}", f"{s.rate(4.0):.3f}",
                f"{s.median_rmsd_a:.2f}" if s.median_rmsd_a is not None else "",
            ])
    print(f"Memorization summary: {summary_path}")
    print()
    print(f"  {'model':<7s} {'variant':<6s} {'n':>3s} {'<2 Å':>6s} {'<4 Å':>6s} {'median':>7s}")
    for (model, variant), s in sorted(stats.items()):
        print(f"  {s.model:<7s} {s.variant:<6s} {s.n_total:>3d} "
              f"{s.rate(2.0):>6.2f} {s.rate(4.0):>6.2f} "
              f"{s.median_rmsd_a:>7.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
