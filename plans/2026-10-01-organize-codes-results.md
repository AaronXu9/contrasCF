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
- **Status:** pending

### S3 — Restructure into three arms + shared core
- **Pre-condition:** S2 merged (or restructure branch cut from the S2 head).
- **Action:** extract the 3 cross-arm dependencies into a shared core; renumber scripts into stages; split per-arm result tables; one notebook location; retire `SURFDOCK_FIX.md` into docs. Update CARC sbatch scripts and the `dockstrat`/`contrascf-casf` skills in the same commit as each rename.
- **Post-condition:** every arm's headline regenerates from its new path and matches the pre-restructure CSVs byte-for-byte or numerically; one CARC dry run passes.
- **Status:** pending

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
- 2026-10-01 — Do NOT re-run ligand-arm SurfDock inside this plan — ~13 GPU-h and a methodological choice (what defines the pocket); tracked as TODO 15 for the user.
