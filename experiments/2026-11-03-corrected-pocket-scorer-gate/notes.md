# corrected-pocket-scorer-gate

**Type:** validation
**Kind:** validation
**Date:** 2026-11-03
**Status:** done

## What this validates
**Plan step:** S1
**Validates iteration:** corrected-six-method-pocket-table
**Validates upcoming experiment:** corrected-six-method-pocket-table

## Validation kind
reference-reproduction (against the 2026-10-21 audit) + regression

## Assertion
The new full-cohort scorer `15_corrected_pocket_scores.py` reproduces the audit's per-row retained/not classification on Set B (96 systems, pack/inv) for GNINA (corrected replays), UniDock2, ICM and SurfDock at ≥95% agreement, with every disagreement explained by the symmetry-aware matcher lowering an RMSD below 2 Å. Regression: the `core.ligand_rmsd` coordinate-indexing fix and the gapped-exact mapping mode change no previously scored row.

## Setup
- **Code:** `27034932cd150944860948f1e033ad4d8d834144` on `audit-followups` (dirty)
- **Env:** `unidock2` (Python 3.10.12)
- **Host:** katlab

## Model & data provenance
| artifact | exact id / name | size | source (repo@sha · HF id · DOI · URL) | version / sha256 | role |
|---|---|---|---|---|---|
| audit reference | `pose_diagnostics.csv`, `surfdock_pose_diagnostics.csv` | 576 + 192 rows | `experiments/2026-10-21-pipeline-audit/pocket/` | 2026-10-21 | UniDock2/ICM/SurfDock reference |
| GNINA corrected | `merged_gnina192.csv`, `replacement_manifest.csv` | 192 / 129 rows | `experiments/2026-10-21-pipeline-audit/` | 2026-10-21 | GNINA reference + corrected poses |
| ICM aligned receptors | `receptor_aligned.pdb` | 717 files | CARC main clone outputs | copied 2026-11-03 | ICM frame |

## Validation code
`validation-snapshots/gate.py` in this folder; scorer at `analysis/casf_mutagenesis/scripts/analyze/15_corrected_pocket_scores.py`.

## Commands & run log
Captured by `lab-notebook run <slug> -- <command>` → see `commands.log` in this entry's folder.

## Expected outcome
Each method passes (agreement ≥95%, no unexplained disagreement); `GATE PASS`. Previously scored rows unchanged after each scorer/core change.

## Actual outcome
Final run: GNINA 192/192 scored, agreement 0.995, retained pack/inv audit 15/15 vs ours 16/15; UniDock2 1.000, 11/7 vs 11/7; ICM 0.990, 15/11 vs 16/12; SurfDock 0.995, 4/3 vs 4/4; no unexplained disagreement. `GATE PASS`.
Regressions: core coordinate-indexing fix changed 0 of 7928 protein-arm and 0 of 12356 ligand-arm co-folding rows (summaries byte-identical); it recovered docking rows previously lost to IndexError. Gapped-exact mapping left all 5392 previously scored rows unchanged and added coverage (errors ~300 → 43 rows).
Defects found on the way: (1) `core.ligand_rmsd` took heavy-atom coordinates from a filtered array while matches indexed a molecule that can retain explicit H — fixed; (2) ICM SDF bond orders fail RDKit sanitization — bond orders taken from the crystal ligand.

## Pass / fail
pass

## Code ↔ theory alignment
| claim / property checked | code (file:line / function) | justified by (ref ↓) | match? |
|---|---|---|---|
| mutation-preserving exact chain correspondence | `15_corrected_pocket_scores.py` `correspondence` | 2026-10-21 audit §1, `audit_pockets.py:one` | exact (+ labeled gapped mode) |
| native-chain 5 Å shell Cα frame, ligand never fitted | `15_corrected_pocket_scores.py` `frame` | audit `audit_pockets.py`; Kabsch 1976 | exact |
| symmetry-aware ligand RMSD | `core/ligand_rmsd.py` `_matched_rmsd` | Masters et al. 2025 protocol | exact |

## References
- 2026-10-21 pipeline audit: `experiments/2026-10-21-pipeline-audit/AUDIT.md` (this repo).
- Masters, M. R. et al. 2025. Nature Communications. doi:10.1038/s41467-025-63947-5
- Kabsch, W. 1976. Acta Cryst. A32, 922–923. doi:10.1107/S0567739476001873

## Optional: If failed — what to fix
n/a
