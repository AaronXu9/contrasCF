"""Analyze Boltz-2 + AF3 predictions on a CASF scope.

For every (pdbid, variant, model) cell with a predicted mmCIF, compute:
  - Cα RMSD vs the crystal protein (sanity: should be small for WT, larger
    for adversarial variants where AF3 may misfold the mutated pocket).
  - Ligand heavy-atom RMSD vs the crystal ligand pose, after Cα superposition.

Aggregate memorization rate per (model, variant) — fraction of *adversarial*
cases (rem / pack / inv) with ligand RMSD < 2 Å (paper's threshold).

Scope (env-var controlled, defaults to subset20):
  CONTRASCF_SCOPE=subset20  → 20 PDBs from PDBbind_casf2016_subset20.json
  CONTRASCF_SCOPE=full      → 285 PDBs from PDBbind_data_split_cleansplit.json["casf2016"]

Output filenames carry the scope suffix:
  results_<scope>.csv, memorization_<scope>.csv

Run:
    source env/lab.sh         # or env/carc.sh
    CONTRASCF_SCOPE=full $CONTRASCF_PY analysis/casf_mutagenesis/scripts/05_analyze_subset20.py
"""
from __future__ import annotations
import csv
import json
import sys
from dataclasses import asdict
from pathlib import Path

import os
REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
sys.path.insert(0, str(REPO_ROOT / "analysis"))

from casf_mutagenesis.analysis import (
    MEMORIZATION_THRESHOLDS_A, MODEL_FILES, analyze_predictions,
    memorization_stats,
)
from casf_mutagenesis.config import OUTPUT_ROOT, SPLIT_JSON, SUBSET20_JSON, VARIANTS


def _resolve_ids() -> tuple[list[str], str]:
    scope = os.environ.get("CONTRASCF_SCOPE", "subset20")
    if scope == "full":
        ids = json.loads(SPLIT_JSON.read_text())["casf2016"]
    else:
        ids = json.loads(SUBSET20_JSON.read_text())["casf2016"]
    return ids, scope


def main() -> int:
    ids, scope = _resolve_ids()
    rows = []
    print(f"{scope}: {len(ids)} systems × {len(VARIANTS)} variants × "
          f"{len(MODEL_FILES)} models = "
          f"{len(ids) * len(VARIANTS) * len(MODEL_FILES)} cells\n")
    for pdbid in ids:
        for variant in VARIANTS:
            for model in MODEL_FILES:
                recs = analyze_predictions(pdbid, variant, model)
                rows.extend(recs)
                for rec in recs:
                    if rec.status == "ok":
                        conf = (rec.confidence_score
                                if rec.confidence_score is not None
                                else (rec.ranking_score or 0.0))
                        print(
                            f"  {pdbid:6s} {variant:5s} {model:7s} "
                            f"pose={rec.pose_idx} "
                            f"lig_rmsd={rec.ligand_rmsd_a:6.2f} Å "
                            f"ca_rmsd={rec.ca_rmsd_a:5.2f} Å  "
                            f"conf={conf:.3f}"
                        )
                    else:
                        print(f"  {pdbid:6s} {variant:5s} {model:7s} "
                              f"pose={rec.pose_idx} {rec.status}"
                              f"{(': ' + rec.error) if rec.error else ''}")

    # CSV
    csv_path = OUTPUT_ROOT / f"results_{scope}.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(asdict(rows[0]).keys()))
        w.writeheader()
        for r in rows:
            w.writerow(asdict(r))
    print(f"\nPer-cell results: {csv_path}")

    # Memorization aggregates (top-1-by-confidence semantics: rank-0 only).
    from casf_mutagenesis.analysis import bootstrap_memorization_ci
    rank0 = [r for r in rows if r.pose_idx == 0]
    stats = memorization_stats(rank0)
    ci_2a = bootstrap_memorization_ci(rows, threshold_a=2.0)
    ci_4a = bootstrap_memorization_ci(rows, threshold_a=4.0)
    summary_path = OUTPUT_ROOT / f"memorization_{scope}.csv"
    with summary_path.open("w", newline="") as f:
        w = csv.writer(f)
        # WT-CONDITIONED rate (primary) alongside the unconditioned one.
        # Restricting to systems where THIS model placed the WT ligand correctly
        # stops a model being credited for "responding" on systems it cannot
        # solve at all. The effect is largest where WT accuracy is lowest --
        # Boltz-2 solves 59% of WT, and conditioning moves its memorization
        # rate 0.166-0.240 -> 0.257-0.368. A model with zero WT-correct systems
        # (e.g. AF3 without MSA) has an UNDEFINED conditional, written as blank,
        # never as 0.
        wt_ok = {}
        for r in rank0:
            if r.variant == "wt" and r.status == "ok" and r.ligand_rmsd_a is not None:
                wt_ok.setdefault(r.model, set())
                if r.ligand_rmsd_a < 2.0:
                    wt_ok[r.model].add(r.pdbid)
        cond = {}
        for r in rank0:
            if r.variant == "wt" or r.status != "ok" or r.ligand_rmsd_a is None:
                continue
            if r.pdbid in wt_ok.get(r.model, set()):
                cond.setdefault((r.model, r.variant), []).append(r.ligand_rmsd_a)

        w.writerow([
            "model", "variant",
            "n_wtok", "memorization_rate_2A_wtok", "memorization_rate_4A_wtok",
            "median_rmsd_A_wtok", "wt_correct_systems",
            "n_total",
            "memorization_rate_2A", "ci_lo_2A", "ci_hi_2A",
            "memorization_rate_4A", "ci_lo_4A", "ci_hi_4A",
            "median_rmsd_A",
        ])
        import statistics as _st
        for (model, variant), s in sorted(stats.items()):
            r2, lo2, hi2 = ci_2a.get((model, variant), (s.rate(2.0), 0.0, 0.0))
            r4, lo4, hi4 = ci_4a.get((model, variant), (s.rate(4.0), 0.0, 0.0))
            cv = cond.get((model, variant), [])
            if cv:
                cw = [len(cv),
                      f"{sum(1 for x in cv if x < 2.0) / len(cv):.3f}",
                      f"{sum(1 for x in cv if x < 4.0) / len(cv):.3f}",
                      f"{_st.median(cv):.2f}"]
            else:
                cw = [0, "", "", ""]
            w.writerow([
                s.model, s.variant,
                *cw, len(wt_ok.get(model, set())),
                s.n_total,
                f"{r2:.3f}", f"{lo2:.3f}", f"{hi2:.3f}",
                f"{r4:.3f}", f"{lo4:.3f}", f"{hi4:.3f}",
                f"{s.median_rmsd_a:.2f}" if s.median_rmsd_a is not None else "",
            ])
    print(f"Memorization summary: {summary_path}")
    print()
    print(f"  {'model':<7s} {'variant':<6s} {'n':>3s} {'<2 Å':>6s} {'<4 Å':>6s} {'median':>7s}")
    for (model, variant), s in sorted(stats.items()):
        print(f"  {s.model:<7s} {s.variant:<6s} {s.n_total:>3d} "
              f"{s.rate(2.0):>6.2f} {s.rate(4.0):>6.2f} "
              f"{s.median_rmsd_a:>7.2f}")

    # Paired-WT framing: per (pdbid, model), Δ RMSD adversarial − WT
    from casf_mutagenesis.analysis import affinity_paired_stats, paired_stats
    for selector, suffix in [("top1", ""), ("oracle", "_oracle")]:
        paired = paired_stats(rows, pose_selector=selector)
        if not paired:
            continue
        path = OUTPUT_ROOT / f"paired{suffix}_{scope}.csv"
        with path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(asdict(paired[0]).keys()))
            w.writeheader()
            for r in paired:
                w.writerow(asdict(r))
        print(f"Paired ({selector}): {path}")

    # Paired affinity (Boltz-2 only): WT vs adversarial Δ in predicted log
    # affinity and binding probability. Empty list if no affinity sidecars
    # were copied (Boltz YAML didn't request affinity before this commit).
    aff_paired = affinity_paired_stats(rows)
    if aff_paired:
        path = OUTPUT_ROOT / f"paired_affinity_{scope}.csv"
        with path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(asdict(aff_paired[0]).keys()))
            w.writeheader()
            for r in aff_paired:
                w.writerow(asdict(r))
        print(f"Paired affinity (Boltz-2): {path}")
        # Brief on-screen summary: per-model per-variant median deltas
        print()
        print(f"  Affinity Δ (adv − wt; positive ⇒ model recognized perturbation):")
        print(f"  {'model':<7s} {'variant':<6s} {'n':>3s} {'median Δaff':>12s} {'median Δprob':>13s}")
        by_key: dict[tuple[str, str], list[tuple[float, float]]] = {}
        for r in aff_paired:
            if r.delta_affinity is None or r.delta_probability is None:
                continue
            by_key.setdefault((r.model, r.variant), []).append(
                (r.delta_affinity, r.delta_probability)
            )
        import statistics
        for (model, variant), pairs in sorted(by_key.items()):
            d_aff = statistics.median(p[0] for p in pairs)
            d_prob = statistics.median(p[1] for p in pairs)
            print(f"  {model:<7s} {variant:<6s} {len(pairs):>3d} "
                  f"{d_aff:>12.3f} {d_prob:>13.3f}")
    else:
        print("No Boltz-2 affinity records found — re-run Boltz-2 with "
              "affinity-enabled YAMLs to populate paired_affinity_<scope>.csv.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
