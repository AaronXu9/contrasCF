"""Run GNINA-pose RMSD analysis on both casf_mutagenesis WT cells and
ligand_mutagenesis variant cells.

Walks two output trees:
  - casf_mutagenesis/outputs/<sys>/wt/gnina/poses.sdf   (binding-site WT)
  - ligand_mutagenesis/outputs/<sys>/<variant>/gnina/poses.sdf

For each cell:
  - Top-1 GNINA pose is loaded (first SDF record).
  - Crystal ligand loaded from data/casf2016/crystal_ligands/<sys>_ligand.sdf.
  - Heavy-atom RMSD computed via MCS (handles modified ligands).
  - Pose is already in crystal frame (GNINA docks against crystal protein
    with box on crystal-lig centroid), so no superposition is needed.

Writes:
  outputs/gnina_results.csv  — per-cell rows
  outputs/gnina_memorization.csv — aggregate per (module, variant)

Lower-is-better is the standard pose-prediction direction for the WT
binding-site case. For ligand_mutagenesis variants, the paper's framing
(higher adversarial RMSD = more physics-aware) DOES apply: a model that
keeps the modified ligand at the crystal pose despite chemistry changes
is memorising; one that moves it away is responding to physics.
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

from casf_mutagenesis.gnina_analysis import GninaRecord, analyze_gnina
from casf_mutagenesis.config import CASF_LIGANDS

CASF_OUTPUTS = REPO_ROOT / "analysis" / "casf_mutagenesis" / "outputs"
LIGAND_OUTPUTS = REPO_ROOT / "analysis" / "ligand_mutagenesis" / "outputs"


def discover_cells() -> list[tuple[str, str, str, Path]]:
    """Return [(module, system, variant, gnina_sdf), ...] for every cell."""
    cells: list[tuple[str, str, str, Path]] = []
    # casf binding-site: walk all variants (wt now joined by rem/pack/inv)
    for module, root in (("casf", CASF_OUTPUTS), ("ligand", LIGAND_OUTPUTS)):
        if not root.is_dir():
            continue
        for sys_dir in sorted(root.iterdir()):
            if not sys_dir.is_dir() or len(sys_dir.name) != 4:
                continue
            for var_dir in sorted(sys_dir.iterdir()):
                if not var_dir.is_dir():
                    continue
                sdf = var_dir / "gnina" / "poses.sdf"
                if sdf.exists():
                    cells.append((module, sys_dir.name, var_dir.name, sdf))
    return cells


def main() -> int:
    cells = discover_cells()
    print(f"Discovered {len(cells)} GNINA poses to analyze")
    rows: list[GninaRecord] = []
    for i, (module, system, variant, sdf) in enumerate(cells, 1):
        crystal = CASF_LIGANDS / f"{system}_ligand.sdf"
        # receptor.pdb is alongside the docking inputs — sdf.parent ==
        # <v_dir>/gnina, so <v_dir>/docking/receptor.pdb.
        receptor_pdb = sdf.parent.parent / "docking" / "receptor.pdb"
        rec = analyze_gnina(
            system, variant, module, sdf, crystal,
            receptor_pdb=receptor_pdb if receptor_pdb.exists() else None,
        )
        rows.append(rec)
        if rec.status == "ok":
            print(
                f"  [{i:4d}/{len(cells)}] {module:6s} {system:6s} {variant:<18s} "
                f"rmsd={rec.rmsd_a:6.2f} Å  matched={rec.n_matched_heavy}/"
                f"{rec.n_heavy_crystal}"
            )
        else:
            print(f"  [{i:4d}/{len(cells)}] {module:6s} {system:6s} {variant:<18s} "
                  f"{rec.status}{(': ' + rec.error) if rec.error else ''}")

    out_csv = REPO_ROOT / "analysis" / "casf_mutagenesis" / "outputs" / "gnina_results.csv"
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(asdict(rows[0]).keys()))
        w.writeheader()
        for r in rows: w.writerow(asdict(r))
    print(f"\nPer-cell GNINA results: {out_csv}")

    # Aggregates: <2 Å, <4 Å, median by (module, variant)
    from collections import defaultdict
    import statistics as st
    by = defaultdict(list)
    for r in rows:
        if r.status == "ok" and r.rmsd_a is not None:
            by[(r.module, r.variant)].append(r.rmsd_a)

    agg_csv = REPO_ROOT / "analysis" / "casf_mutagenesis" / "outputs" / "gnina_memorization.csv"
    with agg_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["module", "variant", "n", "<2_A", "<4_A", "median_rmsd_A"])
        print(f"\n  {'module':<7s} {'variant':<18s} {'n':>3s} {'<2Å':>5s} {'<4Å':>5s} {'median':>7s}")
        for (module, variant), vals in sorted(by.items()):
            n = len(vals)
            b2 = sum(1 for x in vals if x < 2)
            b4 = sum(1 for x in vals if x < 4)
            med = st.median(vals)
            w.writerow([module, variant, n, f"{b2/n:.3f}", f"{b4/n:.3f}", f"{med:.2f}"])
            print(f"  {module:<7s} {variant:<18s} {n:>3d} {b2/n:>5.2f} {b4/n:>5.2f} {med:>7.2f}")
    print(f"\nAggregate: {agg_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
