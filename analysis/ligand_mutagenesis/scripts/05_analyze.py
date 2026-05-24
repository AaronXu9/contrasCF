"""Analyze Boltz-2 (with affinity) predictions on the ligand_mutagenesis module.

Sibling of `casf_mutagenesis/scripts/05_analyze_subset20.py`. Walks
`manifest_full.json` for each (pdbid, variant) with predicted CIFs, then:
  - Loads model_0..N CIFs (best-of-5 if present), computes ligand RMSD vs the
    crystal ligand via the same substructure-matching machinery used for
    casf_mutagenesis (handles halo/methyl variants where pred has extra atoms
    relative to crystal — RDKit finds the shared scaffold).
  - Reads the affinity sidecar (`affinity_<prefix>.json`) — broadcast to every
    pose record for that cell since affinity is per-system in Boltz-2.
  - Computes paired Δaffinity vs WT per system for halo / charge / methyl
    variants.

The ligand_mutagenesis module only has Boltz-2 results (no AF3 runs), so the
output CSVs are Boltz-2-only.

Outputs (in OUTPUT_ROOT = analysis/ligand_mutagenesis/outputs/):
  - results_ligand.csv          — per-pose data
  - paired_affinity_ligand.csv  — WT↔adversarial Δaff / Δprob per (pdbid, variant)

Run:
    LD_LIBRARY_PATH=/home/aoxu/miniconda3/envs/rdkit_env/lib:$LD_LIBRARY_PATH \\
        /home/aoxu/miniconda3/envs/rdkit_env/bin/python \\
        analysis/ligand_mutagenesis/scripts/05_analyze.py
"""
from __future__ import annotations
import csv
import json
import os
import sys
from dataclasses import asdict
from pathlib import Path

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
sys.path.insert(0, str(REPO_ROOT / "analysis"))

from casf_mutagenesis.analysis import (  # noqa: E402
    MODEL_FILES, PredictionRecord, _analyze_single_pose,
    affinity_paired_stats,
)
from ligand_mutagenesis.config import OUTPUT_ROOT  # noqa: E402


def main() -> int:
    manifest_path = OUTPUT_ROOT / "manifest_full.json"
    if not manifest_path.exists():
        print(f"manifest not found: {manifest_path}", file=sys.stderr)
        return 1
    manifest = json.loads(manifest_path.read_text())
    systems = [s for s in manifest["systems"] if s["status"] == "ok"]
    print(f"manifest: {len(systems)} ok systems")

    spec = MODEL_FILES["Boltz2"]
    rows: list[PredictionRecord] = []
    n_cells = n_missing = n_ok = 0

    for sys_entry in systems:
        pdbid = sys_entry["pdbid"]
        sys_dir_base = OUTPUT_ROOT / pdbid
        for v_name, v_info in sys_entry["variants"].items():
            variant_smiles = v_info.get("smiles")
            if not variant_smiles:
                continue
            v_dir = sys_dir_base / v_name
            prefix = f"{pdbid}_{v_name}"
            cifs = sorted(v_dir.glob(spec["cif_glob"].format(prefix=prefix)))
            n_cells += 1
            if not cifs:
                rows.append(PredictionRecord(
                    pdbid=pdbid, variant=v_name, model="Boltz2",
                    pose_idx=0, status="missing_cif",
                ))
                n_missing += 1
                continue
            for rank, cif in enumerate(cifs):
                rec = _analyze_single_pose(
                    pdbid, v_name, "Boltz2", rank, cif, spec, v_dir,
                    smiles_override=variant_smiles,
                )
                rows.append(rec)
                if rec.status == "ok":
                    n_ok += 1

    print(f"cells: {n_cells} ({n_missing} missing CIFs); "
          f"poses analyzed: {n_ok} ok / {len(rows)} total")

    # Per-pose CSV
    results_path = OUTPUT_ROOT / "results_ligand.csv"
    with results_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(asdict(rows[0]).keys()))
        w.writeheader()
        for r in rows:
            w.writerow(asdict(r))
    print(f"Per-pose results: {results_path}")

    # Paired affinity Δ vs WT — uses rank-0 (top-1-by-confidence) per cell
    aff_paired = affinity_paired_stats(rows)
    if aff_paired:
        paired_path = OUTPUT_ROOT / "paired_affinity_ligand.csv"
        with paired_path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(asdict(aff_paired[0]).keys()))
            w.writeheader()
            for r in aff_paired:
                w.writerow(asdict(r))
        print(f"Paired affinity (Boltz-2, ligand-side): {paired_path}")

        # Stdout summary: median Δaff / Δprob per variant group
        import statistics
        from collections import defaultdict
        by_v: dict[str, list[tuple[float, float]]] = defaultdict(list)
        for r in aff_paired:
            if r.delta_affinity is None or r.delta_probability is None:
                continue
            by_v[r.variant].append((r.delta_affinity, r.delta_probability))
        print()
        print(f"  Affinity Δ (adv − wt; positive ⇒ model recognized perturbation):")
        print(f"  {'variant':<20s} {'n':>3s} {'median Δaff':>13s} {'median Δprob':>14s}")
        for variant in sorted(by_v):
            pairs = by_v[variant]
            d_aff = statistics.median(p[0] for p in pairs)
            d_prob = statistics.median(p[1] for p in pairs)
            print(f"  {variant:<20s} {len(pairs):>3d} {d_aff:>13.3f} {d_prob:>14.3f}")
    else:
        print("No affinity sidecars found.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
