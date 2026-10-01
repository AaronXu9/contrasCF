# cofold-variant-rmsd-fix-gate

**Type:** validation
**Kind:** validation
**Date:** 2026-10-01
**Status:** done

## What this validates
**Plan step:** S1
**Validates iteration:** ligand-arm-five-method-memorization
**Validates upcoming experiment:** ligand-arm-five-method-memorization

## Validation kind
reference-reproduction (synthetic known answer) + equivalence (full protein-arm regression)

## Assertion
The new `_atom_correspondences` in `analysis/casf_mutagenesis/analysis.py` (F1) scores a ligand against itself-plus-one-atom, itself-minus-one-atom and itself-with-an-element-swap at 0 Å in-place and best-fit, through the superset, subset and MCS paths respectively, while the OLD code is shown to fail the same case once atoms are reordered; (F2) leaves the protein arm's rates and flags unchanged; (F3) brings a fluorine variant's pocket-blind best-fit RMSD back to its wild type's level.

## Setup
- **Code:** `8c756b400ee432a50dd11d28b8f1e70e5fbe6a32` on `fix/af3-multichain-and-data-audit` (dirty)
- **Env:** `unidock2` (Python 3.10.12)
- **Host:** katlab

## Model & data provenance
| artifact | exact id / name | size | source (repo@sha · HF id · DOI · URL) | version / sha256 | role |
|---|---|---|---|---|---|
| crystal ligand 1bcu | `crystal_ligands/1bcu_ligand.sdf` | 16 heavy atoms | PDBbind clean-split, `/home/aoxu/projects/VLS-Benchmark-Dataset/data/pdbbind_cleansplit` | CASF-2016 core | F1 synthetic base |
| protein-arm results | `analysis/casf_mutagenesis/outputs/results_full.csv` | 7896 rows | contrasCF, frozen copy pre-fix | 2026-09-04 | F2 reference |
| ligand-arm predictions | Boltz-2 2.2.1, AF3 v3.0.1 + MSA | 12324 pose rows | `analysis/ligand_mutagenesis/outputs/` | 2026-05-23 / 2026-09-04 | F3 |

## Validation code
Snapshot: `validation-snapshots/fix_gate.py` (F1, F2-sample); the pre-fix module is preserved as `validation-snapshots/analysis_pre_fix.py`. F2-full and F3 are the diff blocks recorded in Actual outcome, run against `frozen_pre_s1/` copies.

## Commands & run log
Captured by `lab-notebook run <slug> -- <command>` → see `commands.log` in this entry's folder.

## Expected outcome
F1: 0 Å (< 1e-6) on all three synthetic cases, OLD > 0.5 Å on the reordered case. F2: no rate or boolean flag changes in `memorization_full.csv` or any `paired_*_full.csv`. F3: halo_F_1 median bestfit within 0.3 Å of WT for both models.

## Actual outcome
```
[PASS] F1 new +F: in-place 0.00e+00  bestfit 1.09e-15  matched 16  mode pred_superset
[PASS] F1 new -1 atom: in-place 0.00e+00  bestfit 7.96e-16  matched 15  mode crystal_superset
[PASS] F1 new element swap: in-place 0.00e+00  bestfit 7.91e-16  matched 15  mode mcs
[PASS] F1 test bites: OLD code wrong on reordered +F: old 3.47 A vs new 0.00e+00 A
[PASS] F2 protein arm unchanged (sample): 276 cells re-scored, max |Δ| 0.0000 A
[PASS] F3 Boltz2: halo_F_1 bestfit 0.77 vs wt 0.76
[PASS] F3 AF3+MSA: halo_F_1 bestfit 0.46 vs wt 0.47
```
F2 full re-score (`CONTRASCF_SCOPE=full 05_analyze_subset20.py`), all 7896 rows: 5 `ligand_rmsd_a` values changed (1eby pack Boltz-2 pose 4, 1y6r wt Boltz-2 pose 0, 2w66 inv Boltz-2 pose 4, 4tmn rem AF3+MSA pose 0, 5c28 pack Boltz-2 pose 4) — molecules failing RDKit sanitization, which the old code silently paired by file order and the new code maps by MCS. **`memorization_full.csv`: no cell changed. `paired_full` / `paired_oracle_full`: only RMSD values moved, no `wt_correct_2A` / `memorized_given_wt` flag changed. `paired_confidence_full`, `paired_affinity_full`: byte-identical.** A topology-keyed correspondence cache was added after the first F1/F2 pass and both were re-run green.

## Pass / fail
pass

## Code ↔ theory alignment
| claim / property checked | code (file:line / function) | justified by (ref ↓) | match? |
|---|---|---|---|
| variant scored on shared scaffold via superset / subset / MCS correspondence | `analysis/casf_mutagenesis/analysis.py` `_atom_correspondences_uncached` | RDKit FMCS (Dalke & Hastings 2013) | exact |
| in-place RMSD minimised over symmetry-equivalent mappings | `analysis.py` `_matched_rmsd` | Masters et al. 2025 RMSD protocol | exact |
| pocket-blind Kabsch RMSD over the same mappings | `analysis.py` `_bestfit_rmsd`, `_kabsch_rmsd` | Kabsch 1976 | exact |
| cache is behaviour-preserving (index-exact topology key) | `analysis.py` `_topology_key` | F1/F2 re-run after adding it | exact |

## References
- Masters, M. R. et al. 2025. *Investigating whether deep learning models for co-folding learn the physics of protein–ligand interactions.* Nature Communications. (PDF in `docs/`.)
- Dalke, A. & Hastings, J. 2013. *FMCS: a novel algorithm for the multiple MCS problem.* J. Cheminform. 5(Suppl 1), O6. doi:10.1186/1758-2946-5-S1-O6
- Kabsch, W. 1976. *A solution for the best rotation to relate two sets of vectors.* Acta Cryst. A32, 922–923. doi:10.1107/S0567739476001873
- RDKit: Open-source cheminformatics. https://www.rdkit.org

## Optional: If failed — what to fix
n/a — passed.
