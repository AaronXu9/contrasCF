# Audit Followups

**Type:** plan
**Date:** 2026-11-03
**Status:** done (S1 here; S2-S5 owned by another session)
**Project:** contrasCF

## Goal
Execute the follow-up sequence of the 2026-10-21 pipeline audit
(`experiments/2026-10-21-pipeline-audit/AUDIT.md`, §3) so the pocket-arm method
table can be reissued. Success: a corrected six-method pocket table on the full
cohort, scored with mutation-preserving correspondence and no-op exclusion, that
reproduces the audit's Set B numbers where they overlap; plus the receptor-preparation
and Boltz-2 instrumentation evidence the audit asks for before interpreting it.

## Background / context
The audit (another session, 2026-10-21) found: legacy chain selection can align a
multimer onto the wrong, differently mutated site (6 retained rows were alternate
sites); GNINA mutant inputs had uppercase CL/BR → dummy atoms (129 cells rerun);
4EO8 is a no-op; UniDock2 actually used a 30 Å box (`size` vs `box_size`); Boltz-2's
saved structure is not the affinity-scored sample. Its scorer
(`pocket/audit_pockets.py`) covers only Set B (96 systems), pack/inv, GNINA/UniDock2/ICM
(+SurfDock supplement). Boltz-2 and AF3+MSA were not re-audited, so the six-method
ranking is not revalidated. User instruction 2026-11-03: "merge PR#2, and audit
follow-ups, and update the table."

## Steps

### S1 — Corrected six-method pocket table, full cohort (audit follow-up 1)
- **Pre-condition:** corrected GNINA replays available (`gnina_replay/replacement_manifest.csv`, 129 cells); ICM aligned receptors for all 717 mutant cells (copied from CARC to `experiments/2026-11-03-icm-aligned-receptors/`).
- **Action:** one scorer for GNINA (corrected), UniDock2, ICM, SurfDock, Boltz-2, AF3+MSA over wt/rem/pack/inv: exact input-sequence chain mapping, native-chain 5 Å shell Cα fit, symmetric ligand RMSD, no-op flag, changed-residue contact and native-in-box diagnostics.
- **Post-condition:** gate reproduces the audit's Set B pack/inv counts; table reissued in `docs/results_index.md` with own-denominator and matched-set views.
- **Status:** done (2026-11-03). Gate PASS (99-100% per-row agreement). Corrected table in `docs/results_index.md` §1a. Finding: 106 of 117 co-folding retentions lost under the corrected protocol are relocations to a homologous intact site in multimers (AF3+MSA rem 0.439→0.326, Boltz-2 0.368→0.286); docking nearly unchanged. Two defects fixed en route: core RMSD coordinate indexing with retained explicit H; ICM bond orders.

### S2 — Separate pocket damage from receptor preparation (audit follow-up 2)
- **Pre-condition:** S1 done.
- **Action:** pilot on a predeclared stratified subset (no-op / low / high dose, retained and lost): build in-place side-chain-substituted receptors on the crystal with restrained relaxation; dock GNINA/UniDock2 with a WT-defined native-site box; compare with the AF3-refolded receptor.
- **Post-condition:** per-case contact change, backbone movement, clash and native-pose accommodation, plus retention under both receptor sources.
- **Status:** ALREADY DONE by another session (`experiments/2026-10-21-matched-receptor-*`, `2026-10-21-matched-receptors/RESULTS.md`): with in-place prepared receptors, 2YKI retains 0/9 mutant vs 9/9 WT (AF3 receptor retained 3/3); retention is strongly receptor-preparation-dependent. Not re-run here.

### S3 — Instrument Boltz-2 on a predeclared panel (audit follow-up 3)
- **Pre-condition:** the audit's single-system direct-channel gate (`boltz/direct_pose_channel_3u5j_v2`) passes.
- **Action:** extend it to a small predeclared panel × 3 seeds, saving both inference passes and the affinity-scored sample, with native / ejection / translation controls.
- **Post-condition:** per-system coordinate-sensitivity estimates with seed spread.
- **Status:** ALREADY DONE / extended by another session (`2026-10-27-boltz-*`, `2026-10-30-boltz-head-replay-*`, `2026-11-02-boltz-full-cohort-*`): full cohort, 687 comparisons — pose motion tracks ligand-ipTM change (rho -0.61) far more than affinity change (rho -0.20). Not re-run here.

### S4 — Experimental protein variants; S5 — interaction-aware model
- IN PROGRESS in another session as of 2026-11-03 (`2026-11-02-experimental-pocket-*`, `2026-11-02/03-pocket-ranker-*`, entries marked running). Not touched here to avoid colliding with active work.

## Risk register
| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| New scorer disagrees with the audit's Set B numbers | med | high | S1 gate; explain every disagreeing row before reissuing |
| Co-folded chains do not map by exact sequence (stripped residues) | med | med | report as `receptor_error`, never fall back to nearest chain |
| S2 relaxation protocol choices bias retention | med | med | predeclare subset and protocol before looking at outcomes |

## Decisions log
- 2026-11-03 — Interpret "update the table" as reissuing the pocket-arm six-method table once S1's gate passes; the audit forbids reissuing it before the correspondence checks.
