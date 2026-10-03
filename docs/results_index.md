# Results index

Every headline number in this project, mapped to the file that holds it and the
script that writes it. Values were read from the files on **2026-10-03**, after the
co-folding RMSD fix (TODO item 16) and the S3 restructure. If a number here
disagrees with a doc, trust the file and re-check the doc.

Paths are relative to `analysis/`. Run everything after `source env/lab.sh`;
`$PY` below means `$CONTRASCF_PY`. Old script names resolve via each arm's
`scripts/RENAMES.md`.

## Regenerate all headline tables

```bash
source env/lab.sh; PY=$CONTRASCF_PY; C=analysis/casf_mutagenesis/scripts; G=analysis/ligand_mutagenesis/scripts
$PY $C/analyze/01_analyze_cofold.py           # co-folding, full CASF (default scope)
$PY $C/analyze/02_analyze_gnina.py            # GNINA, both arms
$PY $C/analyze/03_analyze_docking_engines.py  # GNINA + UniDock2 + SurfDock, both arms
$PY $C/analyze/04_analyze_icm.py              # ICM
$PY $C/analyze/06_confidence_response.py      # pocket-axis confidence
$PY $C/analyze/07_conf_aff_rmsd.py            # pocket-axis three heads
$PY $G/analyze/01_analyze_cofold.py           # ligand arm co-folding (Boltz-2, AF3+MSA)
$PY $G/analyze/02_confidence_response.py      # ligand-axis confidence (Boltz-2)
$PY $G/analyze/03_retention_table.py          # ligand arm, five methods
for s in all cdk2 gdh mek1; do $PY analysis/paper_repro/scripts/02_run_analysis.py --scope $s; done
$PY -m analysis.ramd_pilot.tests.test_smoke   # tau-RAMD analysis stack self-test
```

These are CPU-only readers of predictions already on disk; the whole set takes
about 15 minutes. All are deterministic — the S3 gate reproduced every file
byte-for-byte.

## 1. Pocket mutagenesis (CASF-2016, wt / rem / pack / inv)

Primary metric: top-1 ligand RMSD < 2 Å, **conditioned on the method solving that
system's WT** (`*_wtok`). Low adversarial = good. WT shown unconditioned.

| claim | value | file (`casf_mutagenesis/outputs/`) | written by (`casf_mutagenesis/scripts/`) |
|---|---|---|---|
| AF3+MSA WT accuracy | 0.824 (n=238) | `results_full.csv` | `analyze/01_analyze_cofold.py` |
| Boltz-2 WT accuracy | 0.594 (n=229) | `results_full.csv` | `analyze/01_analyze_cofold.py` |
| AF3 without MSA | 0/19 WT — conditional undefined | `memorization_full.csv` | `analyze/01_analyze_cofold.py` |
| AF3+MSA memorization rem / pack / inv | 0.439 / 0.367 / 0.301 | `memorization_full.csv` | `analyze/01_analyze_cofold.py` |
| Boltz-2 memorization rem / pack / inv | 0.368 / 0.360 / 0.257 | `memorization_full.csv` | `analyze/01_analyze_cofold.py` |
| GNINA WT; rem / pack / inv | 0.729; 0.180 / 0.157 / 0.151 | `docking_memorization.csv` | `analyze/03_analyze_docking_engines.py` |
| UniDock2 WT; rem / pack / inv | 0.578; 0.141 / 0.141 / 0.104 | `docking_memorization.csv` | `analyze/03_analyze_docking_engines.py` |
| SurfDock WT; rem / pack / inv | 0.876; 0.024 / 0.043 / 0.030 | `docking_memorization.csv` | `analyze/03_analyze_docking_engines.py` |
| ICM WT; rem / pack / inv | 0.685; 0.133 / 0.170 / 0.114 | `icm_results.csv` | `analyze/04_analyze_icm.py` |
| best-of-5 (oracle) companion | per cell | `paired_oracle_full.csv` | `analyze/01_analyze_cofold.py` |
| per-system WT vs mutant pairs | per cell | `paired_full.csv` | `analyze/01_analyze_cofold.py` |
| headline figure | `figures/overview_full.png` | — | `plot/02_plot_overview.py` (**panel b stale, TODO 19**) |
| like-for-like common system set | `figures/matched_comparison.png` | `docking_results.csv`, `results_full.csv` | `plot/05_plot_matched_comparison.py` |
| per-system paired RMSD scatters | `figures/paired_rmsd_{wt_vs_mutant,rem,pack,inv}.png` | `docking_results.csv`, `icm_results.csv` | `plot/03_plot_paired_rmsd.py` |
| dose dependence on mutation count | journal `2026-09-04-memorization-is-dose-dependent` | `results_full.csv` | see journal entry |

**Three heads dissociate (Lark §3.5):**

| claim | file | written by |
|---|---|---|
| confidence registers pocket damage | `confidence_stats_full.csv`, `confidence_conditional_full.csv`, `paired_confidence_full.csv` | `analyze/06_confidence_response.py` |
| affinity head vs pose (Q1–Q3) | `q1_affinity_strata.csv`, `q2_within_case_rho.csv`, `q3_corr_matrix.csv`, `paired_conf_aff_rmsd.csv` | `analyze/07_conf_aff_rmsd.py` |
| claim figure | `figures/three_heads_dissociate.png` | `plot/06_plot_three_heads.py` |
| pose-swap capstone | `figures/pose_swap_{affinity,contrast}.png` | `analyze/09`–`12_pose_swap_*.py` (**tables not persisted, see Gaps**) |

**Physical validity of retained poses:**

| claim | file | written by |
|---|---|---|
| retained poses are not strained (PoseBusters, rotamers) | `cofold_validity.csv`, `cofold_rotamers.csv`, `crystal_rotamer_baseline.csv` | `analyze/13_validate_cofold_pockets.py` |
| crystal pose rescored in the mutated pocket | `crystal_pose_rescore.csv` | `analyze/14_rescore_crystal_pose.py` |

## 2. Ligand mutagenesis (CASF-2016, halo / meth / charge)

Metric: **retention** — top-1 RMSD < 2 Å given the method solved WT. Not
memorization: keeping the pose after a one-atom edit is plausibly correct. Always
quote with Δheavy (heavy atoms added).

| claim | value (Δheavy 1 / 2 / ≥3) | file (`ligand_mutagenesis/outputs/`) | written by (`ligand_mutagenesis/scripts/`) |
|---|---|---|---|
| AF3+MSA retention | 0.899 / 0.400 / 0.318 | `retention_by_dheavy.csv` | `analyze/03_retention_table.py` |
| Boltz-2 retention | 0.820 / 0.424 / 0.333 | `retention_by_dheavy.csv` | `analyze/03_retention_table.py` |
| GNINA retention | 0.709 / 0.278 / 0.223 | `retention_by_dheavy.csv` | `analyze/03_retention_table.py` |
| UniDock2 retention | 0.653 / 0.219 / 0.104 | `retention_by_dheavy.csv` | `analyze/03_retention_table.py` |
| SurfDock | **invalid — do not quote** (TODO 15) | `retention_by_dheavy.csv` | — |
| per method × variant, Wilson CIs | 75 rows | `memorization_ligand.csv` | `analyze/03_retention_table.py` |
| co-folding per pose | 12324 rows | `results_ligand.csv` | `analyze/01_analyze_cofold.py` |
| docking per cell / aggregate | — | `docking_results_ligand.csv`, `docking_memorization_ligand.csv`, `gnina_*_ligand.csv` | `casf_mutagenesis/scripts/analyze/02`, `03_*` |
| Boltz-2 retained / moved under ligand edits | 504 / 224 adversarial cells | `paired_conf_aff_rmsd_ligand.csv` | `analyze/02_confidence_response.py` |
| Boltz-2 confidence ranks poses on variants | AUROC 0.87 | `s2c_pooled_pose_auroc_ligand.csv` | `analyze/02_confidence_response.py` |
| Boltz-2 affinity Δ vs WT | per variant | `paired_affinity_ligand.csv` | `analyze/01_analyze_cofold.py` |

Gated by `experiments/2026-10-01-ligand-arm-five-method-memorization/`. Any
ligand-arm co-folding number produced before 2026-10-01 is wrong (TODO 16).

## 3. Paper reproduction (20 hand-built cases: CDK2, GDH, MEK1)

| claim | file (`results/`) | written by (`paper_repro/scripts/`) |
|---|---|---|
| per model × case RMSD, confidence, clashes | `all/results.csv` (+ `cdk2/`, `gdh/`, `mek1/`) | `02_run_analysis.py --scope <s>` |
| figures | `<scope>/figures/` | `03_make_plots.py`, `04_render_figures.py` (PyMOL env) |
| physics analysis | `<scope>/` | `05_physics_analysis.py` |

`results/gdh/results.csv` was stale until 2026-10-03 (TODO 17); figures under
`results/gdh/figures/` have not been regenerated.

## 4. τ-RAMD physics validation (CDK2 WT vs all-Phe pocket)

| claim | value | file (`ramd_pilot/outputs/production_force6/`) | written by (`ramd_pilot/scripts/`) |
|---|---|---|---|
| verdict at force 6 | **MARGINAL** (gate R ≥ 5) | `pilot_decision.json` | `04_analyze.py` |
| residence-time ratio R | 4.25 | `pilot_decision.json` | `04_analyze.py` |
| τ_KM WT / pack | 7.04 ns / 1.66 ns | `pilot_decision.json`, `system_stats.json` | `04_analyze.py` |
| per-replica exits | n = 15 + 15 | `master_table.csv` | `03_collect_results.py` |

Details: `docs/ramd_pilot_results.md`. MD runs on CARC only.

## 5. CounterFold (Lark §4)

Lives on branch `worktree-counterfold` (`counterfold/`), not on `main`. Not indexed here.

## Gaps — headline material without reproducible provenance

| item | problem | TODO |
|---|---|---|
| `casf_mutagenesis/outputs/mutant_receptor_alignment.csv` | no script in the repo writes it (dated 2026-09-03); cited as the frame-error gate in the contrascf-casf skill | 20 |
| pose-swap tables `pose_swap_summary.csv`, `pose_swap_contrast.csv` | `analyze/10`, `12` default `--out` to `/tmp/poseswap_panel2`; the tables are gone, only the figures survive | 20 |
| ligand-arm SurfDock | invalid input footprint | 15 |
| overview panel (b) | drawn from pre-fix ligand data | 19 |
