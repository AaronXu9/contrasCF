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
ENGINES = ("gnina", "unidock2", "surfdock")


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

    # Aggregates per (module, engine, variant).
    #
    # Two rates are emitted per row:
    #   * `*_wtok`  -- WT-CONDITIONED (primary). Restricted to systems where THIS
    #     engine placed the wild-type ligand correctly (< 2 A). Without this,
    #     an engine is credited for "responding" on systems it simply cannot
    #     solve, which flatters low-WT-accuracy engines. E.g. UniDock2 solves
    #     only 57.8% of WT, and conditioning moves its adversarial rate
    #     0.071-0.092 -> 0.104-0.141.
    #   * the bare columns -- UNCONDITIONED companion, kept so existing
    #     consumers and previously published numbers stay reproducible.
    import statistics as st
    by = defaultdict(list)
    per_system = defaultdict(dict)   # (module, engine) -> {system: {variant: rmsd}}
    for r in rows:
        if r.status == "ok" and r.rmsd_a is not None:
            by[(r.module, r.engine, r.variant)].append(r.rmsd_a)
            per_system[(r.module, r.engine)][r.system] = \
                {**per_system[(r.module, r.engine)].get(r.system, {}), r.variant: r.rmsd_a}

    # Systems whose WT cell this engine got right.
    wt_ok = {}
    for key, sysd in per_system.items():
        wt_ok[key] = ({s for s, d in sysd.items() if d.get("wt") is not None and d["wt"] < 2.0},
                      {s for s, d in sysd.items() if d.get("wt") is not None})

    agg_csv = REPO_ROOT / "analysis" / "casf_mutagenesis" / "outputs" / "docking_memorization.csv"
    with agg_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["module", "engine", "variant",
                    "n_wtok", "<2_A_wtok", "<4_A_wtok", "median_rmsd_A_wtok",
                    "n", "<2_A", "<4_A", "median_rmsd_A",
                    "wt_correct_systems", "wt_total_systems"])
        print(f"\n  {'module':<7s} {'engine':<9s} {'variant':<18s} "
              f"{'n_wtok':>6s} {'<2Å*':>5s} | {'n':>4s} {'<2Å':>5s} {'median':>7s}")
        for (module, engine, variant), vals in sorted(by.items()):
            n = len(vals)
            b2 = sum(1 for x in vals if x < 2)
            b4 = sum(1 for x in vals if x < 4)
            med = st.median(vals)

            ok_set, all_set = wt_ok.get((module, engine), (set(), set()))
            cvals = [d[variant] for s, d in per_system[(module, engine)].items()
                     if s in ok_set and variant in d]
            cn = len(cvals)
            if cn:
                c2, c4 = sum(1 for x in cvals if x < 2), sum(1 for x in cvals if x < 4)
                cw = [f"{cn}", f"{c2/cn:.3f}", f"{c4/cn:.3f}", f"{st.median(cvals):.2f}"]
                cshow = f"{cn:>6d} {c2/cn:>5.2f}"
            else:
                # No WT-correct systems -> the conditional is undefined, not zero.
                cw = ["0", "", "", ""]
                cshow = f"{0:>6d} {'n/a':>5s}"
            w.writerow([module, engine, variant] + cw +
                       [n, f"{b2/n:.3f}", f"{b4/n:.3f}", f"{med:.2f}",
                        len(ok_set), len(all_set)])
            print(f"  {module:<7s} {engine:<9s} {variant:<18s} {cshow} | "
                  f"{n:>4d} {b2/n:>5.2f} {med:>7.2f}")
    print(f"\nAggregate: {agg_csv}   (*_wtok = WT-conditioned, primary)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
