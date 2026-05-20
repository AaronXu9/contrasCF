"""Pose RMSD analysis for both GNINA and UniDock2 across both modules.

Generalised version of 09_analyze_gnina.py. Walks
  <module_outputs>/<system>/<variant>/{gnina,unidock2}/poses.sdf
for module in {casf_mutagenesis, ligand_mutagenesis} and analyses every
cell that has a poses.sdf.

Outputs:
  outputs/docking_results.csv  — per-cell (with `engine` column)
  outputs/docking_memorization.csv — aggregate per (module, engine, variant)
"""
from __future__ import annotations
import csv
import os
import sys
from collections import defaultdict
from dataclasses import asdict
from pathlib import Path

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
sys.path.insert(0, str(REPO_ROOT / "analysis"))

from casf_mutagenesis.gnina_analysis import GninaRecord, analyze_gnina
from casf_mutagenesis.config import CASF_LIGANDS

CASF_OUTPUTS = REPO_ROOT / "analysis" / "casf_mutagenesis" / "outputs"
LIGAND_OUTPUTS = REPO_ROOT / "analysis" / "ligand_mutagenesis" / "outputs"
ENGINES = ("gnina", "unidock2")


def discover_cells() -> list[tuple[str, str, str, str, Path]]:
    """Return [(module, system, variant, engine, poses_sdf), ...]"""
    cells: list[tuple[str, str, str, str, Path]] = []
    for module, root in (("casf", CASF_OUTPUTS), ("ligand", LIGAND_OUTPUTS)):
        if not root.is_dir():
            continue
        for sys_dir in sorted(root.iterdir()):
            if not sys_dir.is_dir() or len(sys_dir.name) != 4:
                continue
            for var_dir in sorted(sys_dir.iterdir()):
                if not var_dir.is_dir():
                    continue
                for engine in ENGINES:
                    sdf = var_dir / engine / "poses.sdf"
                    if sdf.exists():
                        cells.append((module, sys_dir.name, var_dir.name, engine, sdf))
    return cells


def main() -> int:
    cells = discover_cells()
    by_engine = defaultdict(int)
    for c in cells:
        by_engine[c[3]] += 1
    print(f"Discovered {len(cells)} pose cells: {dict(by_engine)}")
    rows: list[GninaRecord] = []
    for i, (module, system, variant, engine, sdf) in enumerate(cells, 1):
        crystal = CASF_LIGANDS / f"{system}_ligand.sdf"
        receptor_pdb = sdf.parent.parent / "docking" / "receptor.pdb"
        rec = analyze_gnina(
            system, variant, module, sdf, crystal,
            receptor_pdb=receptor_pdb if receptor_pdb.exists() else None,
            engine=engine,
        )
        rows.append(rec)
        if rec.status == "ok":
            print(
                f"  [{i:4d}/{len(cells)}] {engine:<8s} {module:6s} {system:6s} "
                f"{variant:<18s} rmsd={rec.rmsd_a:6.2f} Å  matched={rec.n_matched_heavy}"
            )
        else:
            print(f"  [{i:4d}/{len(cells)}] {engine:<8s} {module:6s} {system:6s} "
                  f"{variant:<18s} {rec.status}"
                  f"{(': ' + rec.error) if rec.error else ''}")

    out_csv = REPO_ROOT / "analysis" / "casf_mutagenesis" / "outputs" / "docking_results.csv"
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(asdict(rows[0]).keys()))
        w.writeheader()
        for r in rows: w.writerow(asdict(r))
    print(f"\nPer-cell results: {out_csv}")

    # Aggregates per (module, engine, variant)
    import statistics as st
    by = defaultdict(list)
    for r in rows:
        if r.status == "ok" and r.rmsd_a is not None:
            by[(r.module, r.engine, r.variant)].append(r.rmsd_a)

    agg_csv = REPO_ROOT / "analysis" / "casf_mutagenesis" / "outputs" / "docking_memorization.csv"
    with agg_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["module", "engine", "variant", "n", "<2_A", "<4_A", "median_rmsd_A"])
        print(f"\n  {'module':<7s} {'engine':<9s} {'variant':<18s} {'n':>3s} "
              f"{'<2Å':>5s} {'<4Å':>5s} {'median':>7s}")
        for (module, engine, variant), vals in sorted(by.items()):
            n = len(vals)
            b2 = sum(1 for x in vals if x < 2)
            b4 = sum(1 for x in vals if x < 4)
            med = st.median(vals)
            w.writerow([module, engine, variant, n,
                        f"{b2/n:.3f}", f"{b4/n:.3f}", f"{med:.2f}"])
            print(f"  {module:<7s} {engine:<9s} {variant:<18s} {n:>3d} "
                  f"{b2/n:>5.2f} {b4/n:>5.2f} {med:>7.2f}")
    print(f"\nAggregate: {agg_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
