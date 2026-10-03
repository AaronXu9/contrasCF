# Organize Codes Results

**Type:** plan
**Date:** 2026-10-01
**Status:** in_progress
**Project:** contrasCF

## Goal
Turn the repo from an accreted history into three clearly separated arms over one
shared core — **paper reproduction** (20 hand-built Masters et al. cases),
**CASF pocket mutagenesis**, **CASF ligand mutagenesis** — with every headline
number traceable to the file and script that produce it. Success criterion: a
fresh reader can regenerate any arm's headline table from `docs/results_index.md`
alone, on lab or CARC, and no module reaches into another arm's `scripts/` by path.

## Background / context
Two ligand-arm runs finished on 2026-09-04 but were never ingested: SurfDock on
lab (1265/1300 ok, 31 failed across 8 systems) and AF3+MSA on CARC (26/26 array
tasks COMPLETED, 1218/1300 structures, results still only on CARC). The repo has
43 unmerged commits on `fix/af3-multichain-and-data-audit`, 49 dirty files,
colliding script numbers (10/11/12/34 each used twice in casf `scripts/`), a
full-set analyzer named `05_analyze_subset20.py` that defaults to 20 systems,
ligand docking results written into the casf module's tables, and two notebook
locations. User decision 2026-10-01: the paper-reproduction arm must stay fully
regenerable (memory: paper-reproduction-arm). The CASF modules depend on it in
three places: `analysis/src/loaders.py`, `analysis/scripts/13_run_surfdock.py`
(loaded by path as the SurfDock engine), and the paper reference SMILES in
`analysis/src/config.py`.

## Steps

### S1 — Freeze the finished ligand-arm results
- **Pre-condition:** both runs terminal (verified: SurfDock log `finished … exit=0`; sacct 26 COMPLETED).
- **Action:** rsync AF3+MSA ligand outputs CARC→lab; extend `ligand_mutagenesis/scripts/05_analyze.py` to score AF3+MSA alongside Boltz-2; rerun `12_analyze_docking_engines.py` so ligand SurfDock enters `docking_results.csv`; write a WT-conditioned, Δheavy-stratified `memorization_ligand.csv`; record the failure inventory.
- **Post-condition:** ligand arm has scored rows for Boltz-2, AF3+MSA, GNINA, UniDock2, SurfDock; `memorization_ligand.csv` exists; every short `n` explained.
- **Validation entries:**
  - [x] ligand-af3msa-scoring-gate (pass; blind spot: WT-only)
  - [x] cofold-variant-rmsd-fix-gate (pass)
- **Implementation entries:**
  - [x] postmortem cofold-ligand-variant-rmsd-mapping (fixed)
  - [x] postmortem surfdock-ligand-arm-input-footprint (open → TODO 15)
- **Experiment entries:**
  - [x] ligand-arm-five-method-memorization (done; SurfDock excluded as invalid)
- **Status:** done (2026-10-01). Post-condition met for 4 of 5 methods; SurfDock rows ingested but invalid.

### S2 — Commit loose work and merge to main
- **Pre-condition:** S1 committed.
- **Action:** commit the 4 untracked scripts, 12 journal entries, `experiments/`; open a PR `fix/af3-multichain-and-data-audit` → main. `worktree-counterfold` stays separate.
- **Post-condition:** working tree clean; PR open with a reviewable summary.
- **Status:** done (2026-10-01) — tree clean; PR https://github.com/AaronXu9/contrasCF/pull/1 open, **awaiting user review/merge**.

### S3 — Restructure into three arms + shared core
- **Pre-condition:** S2 merged (or restructure branch cut from the S2 head).
- **Action:** extract the 3 cross-arm dependencies into a shared core; renumber scripts into stages; split per-arm result tables; one notebook location; retire `SURFDOCK_FIX.md` into docs. Update CARC sbatch scripts and the `dockstrat`/`contrascf-casf` skills in the same commit as each rename.
- **Post-condition:** every arm's headline regenerates from its new path and matches the pre-restructure CSVs byte-for-byte or numerically; one CARC dry run passes.
- **Status:** done on branch `restructure/three-arms` (2026-10-03). Gate: baseline regenerated from the OLD layout first (so staleness could not masquerade as a change), then every entry point re-run from the NEW layout — 6 pocket-arm analyses, 3 ligand-arm analyses, paper arm × 4 scopes, τ-RAMD smoke test: all exit 0; **42 CSVs byte-identical; the 4 split docking tables row-identical** (old combined = new pocket + new ligand). Import smoke: 76/94 load standalone, failure set identical to the old layout (package-only library modules + PyMOL). Found en route: stale `results/gdh/results.csv` (TODO 17), SurfDock box symlink bridge (TODO 18), overview panel (b) needs regen (TODO 19), τ-RAMD `test_smoke.py` is a script not a pytest file.

### S3 design (written 2026-10-01, awaiting go-ahead)

**Keep** the package names `casf_mutagenesis` and `ligand_mutagenesis`: renaming them breaks every import, doc and skill for little gain. Restructure inside them.

**1. Shared core — `analysis/core/`** (removes every cross-arm reach-in):

| new module | from | why |
|---|---|---|
| `core/ligand_rmsd.py` | `casf_mutagenesis/analysis.py` `_atom_correspondences`, `_matched_rmsd`, `_bestfit_rmsd` + `gnina_analysis.py` `_mcs_match_indices` | both atom matchers in ONE file, moved verbatim. **Not unified in S3**: the docking matcher takes a single match without symmetry enumeration, so unifying it changes docking RMSDs (that is TODO item 10, rate understated ~2.4 pts) and would break S3's numbers-identical gate. Unify later as its own gated change |
| `core/loaders.py` | `analysis/src/loaders.py` | imported by the pocket arm |
| `core/surfdock_engine.py` | helpers in `analysis/scripts/13_run_surfdock.py` | loaded BY FILE PATH by `14_run_surfdock_variants.py` |
| `core/reference_smiles.py` | paper SMILES in `analysis/src/config.py` | used by the ligand arm's verify gate |

Old locations re-export from `core/` so the paper-reproduction arm keeps running unchanged.

**2. Paper-reproduction arm** — move `analysis/src/` + `analysis/scripts/` → `analysis/paper_repro/{lib,scripts}/`. Stays fully regenerable (memory: paper-reproduction-arm). `analysis/native/` and `analysis/results/` stay in place (shared data / existing result paths).

**3. Scripts grouped by stage, renumbered inside each stage** (pocket arm, 43 scripts; collisions at 10/11/12/34 disappear):

| stage | pocket-arm scripts (current numbers) |
|---|---|
| `build/` | 00, 01, 02, 10_build_mutant_docking, 11_build_mek1, 25, 26, 37, audit_mutation_presence |
| `run/` | 03, 04, 06, 07, 08, 11_run_unidock2, 12_run_mek1, 14, 27, 28, 36 |
| `analyze/` | 05 → `analyze_cofold.py` (**default scope becomes `full`**; sbatch set `subset20` explicitly), 09, 12_analyze_docking, 15–22, 30, 34_validate, 35 |
| `figures/` | 10_plot_affinity, 13, 29, 31, 32, 33, 34_plot_three_heads |
| `export/` | 23, 24 |

Ligand arm (8 scripts) gets the same four stages. Each `scripts/` gets `RENAMES.md` (old → new) so historical journal references still resolve; journal entries themselves stay immutable.

**4. Per-arm result tables.** `12_analyze_docking_engines.py` writes ligand rows to `ligand_mutagenesis/outputs/docking_results_ligand.csv` instead of the pocket arm's table with a `module` column. Consumers updated in the same commit.

**5. Physics-validation track — τ-RAMD (added 2026-10-03 after user review).** `analysis/ramd_pilot/` stays where it is: it is already a clean package (`prep/`, `ramd/`, `analysis/`, `scripts/00-04`, `tests/`, `carc_setup/`), imports nothing from the three arms, and runs on CARC through `CONTRASCF_RAMD_OUT` / `CONTRASCF_GMX_RAMD` in `env/carc.sh`. It is a fourth component beside the three arms — a physics oracle (exit-time ranking, Kokh & Wade 2018) validating paper cases, verdict MARGINAL at force 6 (R = 4.25 vs gate 5, p < 0.001). Its only input dependency is `analysis/native/` (paper native CIFs), which therefore **stays at `analysis/native/` as shared data** rather than moving under `paper_repro/`. S3 gate adds: `pytest analysis/ramd_pilot/tests` passes.

**6. Housekeeping.** `docs/lab_notebook/` (1 entry) → `journal/`; top-level `SURFDOCK_FIX.md` (contains retracted advice) → `docs/archive/` behind a banner.

**Same-commit updates for every rename:** 8 `slurm/*.sh`, `env/*.sh`, the `dockstrat` and `contrascf-casf` skills, and current docs (`casf_overview`, `casf_mutagenesis`, `ligand_mutagenesis`, `data_store_map`, `data_prep_todo`, `project_notes`). Historical specs/plans under `docs/superpowers/` are left as written.

**Verification gate (S3 post-condition):** from the new paths, regenerate every arm's analysis and diff against the S1 frozen CSVs (`results_full`, `memorization_full`, `paired_*_full`, `docking_results`, `results_ligand`, `memorization_ligand`, paper arm `analysis/results/all/results.csv`): numerically identical. Then on CARC: `git pull` in the worktree, import check, and one `--array=0` dry run of each touched sbatch.

**Branching:** cut `restructure/three-arms` from the S2 head. If PR #1 changes in review, rebase before merging.

### S4 — Results index
- **Pre-condition:** S3 done.
- **Action:** write `docs/results_index.md` mapping each headline number → CSV → producing script → command.
- **Post-condition:** every headline in `casf_overview.md` and the Lark draft has a row.
- **Status:** pending

## Risk register
| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Rename breaks CARC jobs or skills that call scripts by path | high | med | update sbatch + skills in the same commit; CARC dry run gate in S3 |
| Restructure silently changes a number | med | high | S3 post-condition: regenerated CSVs diffed against frozen S1 copies |
| AF3 ligand rsync pulls GBs of scratch | med | low | rsync include-filter on `af3msa_*` + `af3_msa.json` only; size check first |
| Merging 43 commits hides a regression | low | med | PR, not direct merge; user reviews |

## Open questions
- [ ] Do the 8 SurfDock failure systems overlap the protein arm's 13c exclusions? Moot until the ligand-arm SurfDock re-run (TODO 15).
- [ ] Re-read the June ligand-axis "severity-gated" claim against the regenerated confidence tables.

## Decisions log
- 2026-10-01 — Freeze results before restructuring — reorganizing first means re-pointing analysis at moving paths and loses the byte-diff baseline.
- 2026-10-01 — Keep paper-reproduction arm regenerable — user instruction; it is the direct comparison to Masters et al.
- 2026-10-01 — Merge via PR, not direct merge — 43 commits deserve review.
- 2026-10-01 — Fix the shared co-folding scorer inside S1 rather than defer — the five-method table would otherwise freeze garbage; verified the protein arm has zero rate/flag changes across all 7896 rows.
- 2026-10-03 — S3 implemented and gated; CARC dry run is the last check before the PR.
- 2026-10-03 — User approved PR #1 and the S3 layout; τ-RAMD track added to the layout (was omitted); matcher unification split out of S3 because it changes docking numbers.
- 2026-10-01 — Do NOT re-run ligand-arm SurfDock inside this plan — ~13 GPU-h and a methodological choice (what defines the pocket); tracked as TODO 15 for the user.
