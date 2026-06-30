# Three Head Dissociation

**Type:** milestone
**Date:** 2026-06-30
**Status:** done
**Project:** contrasCF

## Headline
On a broken binding pocket, Boltz-2's three output heads **dissociate**: the **structure** head half-memorizes the native pose (~25–30% of cases, RMSD < 2 Å despite a destroyed pocket), the **confidence** head *registers* the break (interface ipTM / P(binder) drop significantly — "the model knows it memorized"), and the **affinity** head is **intrinsically pose-blind** — confirmed by a direct pose-swap intervention to be flat even when the ligand is ejected 35 Å into solvent.

## Why it matters
This is the project's central diagnosis and it is now closed on three independent lines of evidence (correlational confidence analysis on full CASF, an AF3+MSA cross-model replication, and a *construction-based* pose-swap that removes the trunk confound). It reframes the Masters et al. 2025 "affinity is unchanged after mutation" observation: the affinity invariance is **not** evidence that the structure was preserved — structure, confidence, and affinity are decoupled, and the affinity head's learned function is *flat in pose* by construction, not because the mutations were too gentle. A genuine-physics reference (GNINA Vina) collapses 100% on the same decoys, placing Boltz-2's affinity head at the extreme end of a memorization hierarchy. This is the load-bearing result the write-up (`docs/casf_confidence.md`) and any paper claim rest on.

## The three heads (summary of evidence)

| Head | Behavior on a broken pocket | Evidence | Metric |
|---|---|---|---|
| **Structure** | half-memorizes the native pose | ~25–30% of mutated systems still RMSD < 2 Å | per-pose ligand RMSD, best-of-5 |
| **Confidence** | *registers* the break — "knows it memorized" | interface ipTM / P(binder) drop significantly even when structure memorizes | Wilcoxon signed-rank, AUROC (Mann–Whitney), bootstrap CI, Spearman vs RMSD |
| **Affinity** | intrinsically blind to ligand geometry | flat under pocket mutation AND under direct 35 Å ligand ejection (trunk held byte-identical) | pose-swap Δlog[IC50] = **−0.004**; paired Wilcoxon p = 0.27 (n.s.) |

## The arc (what was done, in order)

1. **Confidence-drop test (Q: does confidence drop significantly after pocket mutation?).** Built `16_confidence_response.py` (pure-stdlib: hand-rolled Wilcoxon signed-rank, AUROC via Mann–Whitney, bootstrap CI, Spearman; validated == scipy/sklearn) on `results_full.csv` + `paired_full.csv`. Finding: **yes** — confidence drops are statistically significant and not noise; the model's confidence head registers the destroyed pocket even on cases where the *structure* head memorized. → memory `confidence-registers-broken-pocket`.

2. **Joint confidence × affinity × RMSD analysis (Q1–Q3).** `17_conf_aff_rmsd.py` (protenix env; best-of-5 pose unit; partial Spearman for regression-to-the-mean control; figure). Established the **three-head dissociation** as the headline: confidence tracks RMSD, affinity does not. The faint +0.12 ΔRMSD residual on affinity looked like weak pose-reading but was later shown (pose-swap) to be noise, not signal. → memory `affinity-head-pose-insensitive`.

3. **AF3 + MSA cross-model replication.** `18_ingest_af3msa_confidence.py` ingests AF3 confidence JSONs from CARC; submitted/collected CARC jobs to get full-CASF RMSD (n = 239, up from a subset-20). Confirms the confidence-registers-the-break pattern is not Boltz-specific.

4. **Boltz-2 affinity architecture read.** Traced `AffinityModule.forward(s_inputs, z, x_pred, feats)`: pose enters **only** via a distogram of `x_pred` over protein–ligand cross-pairs (cdist of representative atoms, 64 bins, max 22 Å); the trunk (`s_inputs`, `z`) is pose-independent; `coords_affinity = sample_atom_coords[argsort(iptm)[0]]`; ensemble of 2 modules + MW correction (`model_coef=1.03525938, mw_coef=-0.59992683, bias=2.83288489`). This *architectural* read is what made the pose-swap confound removable. → memory `boltz2-affinity-architecture`.

5. **Pose-swap intervention (the capstone).** `19_pose_swap_affinity.py` monkeypatches `AffinityModule.forward` so the production call *also* scores a ladder of decoy `x_pred` (ligand ejected to clearance 5/15/30 Å beyond the protein bounding sphere) with the trunk held byte-identical; native call **is** production (identity automatic); whole-complex translation control → Δ = 0. Aggregator `20_pose_swap_aggregate.py`. **Result (n=29):** median native→eject30 gap **−0.004** log units, **0/29** weaken by ≥+1, per-system slope ~100× below physics, paired Wilcoxon p = **0.27**. → memory `pose-swap-result`. Detailed ELN: `docs/lab_notebook/2026-06-06_pose-swap-test.md`.

6. **GNINA physics reference.** `21_pose_swap_gnina.py` (`gnina --score_only`) + `22_pose_swap_contrast.py` (3-panel contrast). On the identical decoy ladder: **Vina collapses −8.7 → 0.0 kcal/mol (100% of systems)**, CNNaffinity drops +2.5 pK, CNNscore 0.95 → 0.54 — while Boltz-2's affinity head stays flat. **Hierarchy of memorization:** physics (Vina) clean → gnina's learned CNN partial → Boltz-2 affinity head most extreme (flat).

7. **Pose export for visual inspection.** `23_export_poses.py` (crystal-frame: target.pdb + multi-state ligand SDF + .pml) and `24_export_predicted_poses.py` (Boltz-predicted-frame, reads `*_model_0.cif` via gemmi) — to eyeball that the ejected ligand is genuinely in solvent.

8. **Ligand-axis replication.** `analysis/ligand_mutagenesis/scripts/07_confidence_response.py` (subagent): the confidence-registers pattern holds on the ligand axis too, but **severity-gated** (registers only for sufficiently disruptive ligand changes).

## Critical confound found & fixed (recorded for methods)
The first pose-swap decoy scheme translated the ligand along a *local* pocket-exit vector. For **buried** pockets that vector is near-arbitrary, so a 20 Å translation just **plowed the ligand through the protein** — min ligand–protein distance stayed ~3 Å (not ejected). Caught by an added min-distance diagnostic (prompted by a "which direction?" review question). **Fix:** eject to `centroid + outward·(R_protein + clearance)` and **verify** via recorded min-distance (native ~3 Å → eject30 30–45 Å). *Lesson: verify the achieved separation, not the nominal translation.* A second fix: production `affinity_pred_value` is the **raw ensemble** (pre-MW); MW is a constant offset that cancels in a fixed-ligand Δ ladder — switched to raw to match production exactly.

## Where everything lives
- **Write-up (topic doc):** `docs/casf_confidence.md` — full Q1–Q3 analysis + the pose-swap capstone section. Headline conclusion folded in (commit `6923ebd`).
- **Pose-swap ELN:** `docs/lab_notebook/2026-06-06_pose-swap-test.md` (+ GNINA addendum).
- **Spec:** `docs/superpowers/specs/2026-06-06-pose-swap-test-design.md`.
- **Scripts:** `analysis/casf_mutagenesis/scripts/16`–`24`; `analysis/ligand_mutagenesis/scripts/07`.
- **Figures:** `analysis/casf_mutagenesis/figures/pose_swap_affinity.png`, `pose_swap_contrast.png`.
- **Memories:** `confidence-registers-broken-pocket`, `affinity-head-pose-insensitive`, `boltz2-affinity-architecture`, `pose-swap-result`.
- **Git:** `pose-swap` branch merged to `main`; all commits through `6923ebd` on `origin/main`.

## Provenance (models / data / env)
| Item | Exact id / value | Source |
|---|---|---|
| Boltz-2 | v2.2.1; checkpoints `boltz2_conf.ckpt` + `boltz2_aff.ckpt` + `ccd.pkl`/mols | rsync'd CARC `~/.boltz` → local `/home/aoxu/.boltz`; env `boltzina_env` (torch 2.7.0+cu126) |
| GNINA | v1.3 | `/mnt/katritch_lab2/aoxu/envs/gnina/bin/gnina` |
| Dataset | CASF mutagenesis panel — 252 affinity-enabled systems; pose-swap on n=29 diverse subset; AF3 confidence n=239 | `analysis/casf_mutagenesis/` (secondary analysis of Masters et al. 2025) |
| Analysis env | `protenix` (pandas/scipy/matplotlib/sklearn); `rdkit_env` (gemmi/rdkit/Bio) | `/mnt/katritch_lab2/aoxu/envs/` |
| Compute | lab RTX 4090 (24 GB), ~1 min/system; CARC `katritch_223` for AF3 RMSD | — |

## References
- Masters et al. 2025 — co-folding-models physics-vs-memorization study (the primary work this project does secondary analysis of). [cite DOI/arXiv in final write-up]
- Boltz-2 — Wohlwend et al. (affinity-head architecture; `boltz.model.modules.affinity.AffinityModule`).
- GNINA v1.3 — McNutt et al. (Vina + CNNaffinity/CNNscore scorer used as the physics reference).

## Code ↔ theory alignment
- Pose enters affinity only via `x_pred` distogram → `boltz/model/modules/affinity.py` (cdist of representative atoms, 64 bins, max 22 Å); trunk `s_inputs`/`z` pose-independent. Pose-swap hook in `19_pose_swap_affinity.py` (monkeypatch `AffinityModule.forward`).
- "Best-ipTM sample scored" → `coords_affinity = sample_atom_coords[argsort(iptm)[0]]` in the affinity path; pose-swap holds this fixed and varies only the decoy `x_pred`.
- Confidence metric = interface ipTM / P(binder) → `16_confidence_response.py` reads `results_full.csv`/`paired_full.csv`; stats hand-rolled and validated against scipy/sklearn.

## Open threads
- [suggested] Run a *confidence* pose-swap (matched intervention on the confidence head) to show it reacts where the affinity head does not — symmetric to the affinity pose-swap.
- [suggested] Scale the affinity pose-swap from n=29 to all 252 systems overnight — verdict won't change, only tightens CIs.
- [open] Delete the merged `pose-swap` branch (local + remote) once the merge is confirmed stable on `main`.
- [tentative] The ligand-axis severity-gating (`07_confidence_response.py`) deserves the same RTM-hardened statistical treatment as the pocket axis before it goes in the write-up.
