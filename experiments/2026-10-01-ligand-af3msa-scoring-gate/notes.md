# ligand-af3msa-scoring-gate

**Type:** validation
**Kind:** validation
**Date:** 2026-10-01
**Status:** done

## What this validates
**Plan step:** S1
**Validates iteration:** ligand-arm-five-method-memorization
**Validates upcoming experiment:** ligand-arm-five-method-memorization

## Validation kind
reference-reproduction + equivalence (regression)

## Assertion
Extending `ligand_mutagenesis/scripts/05_analyze.py` to score AF3+MSA, and re-running `12_analyze_docking_engines.py` to ingest ligand-arm SurfDock, (G1) leaves every Boltz-2 row unchanged, (G2) reproduces the protein arm's AF3+MSA wild-type top-1 result — rate within 0.08 and per-system pass/fail concordance >= 0.85 — because the ligand arm's `wt` cells use the same protein chains, the same native SMILES (`ligand_io.resolve_wt_smiles` wraps the casf resolver) and the same MSA cache, and (G3) leaves every pre-existing docking row unchanged while adding the ligand SurfDock rows.

## Setup
- **Code:** `8c756b400ee432a50dd11d28b8f1e70e5fbe6a32` on `fix/af3-multichain-and-data-audit` (dirty)
- **Env:** `unidock2` (Python 3.10.12)
- **Host:** katlab

## Model & data provenance
| artifact | exact id / name | size | source (repo@sha · HF id · DOI · URL) | version / sha256 | role |
|---|---|---|---|---|---|
| AlphaFold3 | AF3 v3.0.1, `af3.bin` weights | ~1.1 GB weights | github.com/google-deepmind/alphafold3; CARC `/project2/katritch_223/aoxu/dockStrat/forks/alphafold3/models/af3.bin` | v3.0.1 | AF3+MSA ligand-arm predictions (CARC array 11750008/11750009) |
| Boltz-2 | boltz 2.2.1 | n/a | github.com/jwohlwend/boltz | 2.2.1 (`boltzina_env`) | Boltz-2 rows (regression baseline) |
| SurfDock | SurfDock, `model_weights/{docking,posepredict}` | 140 MB weights | github.com/CAODH/SurfDock; lab `/home/aoxu/projects/SurfDock` | dockStrat CogLigandBench `2107de18e9` crop fix | ligand-arm SurfDock rows |
| CASF-2016 core, ligand arm | `ligand_mutagenesis/outputs/manifest_full.json` | 251 systems / 1300 cells | contrasCF@`8c756b4` | manifest 2026-05-07 | evaluated cells |
| CASF-2016 core, protein arm | `casf_mutagenesis/outputs/results_full.csv` | AF3+MSA wt n=238 | contrasCF@`8c756b4` | frozen 2026-09-04 | G2 reference |

## Validation code
Snapshot: `validation-snapshots/gate.py` in this folder (argument: scratchpad dir holding the pre-change `results_ligand.pre_af3.csv` and `frozen_pre_s1/` docking tables).

## Commands & run log
Captured by `lab-notebook run <slug> -- <command>` → see `commands.log` in this entry's folder.

## Expected outcome
All five checks print PASS; `GATE PASS`.

## Actual outcome
```
[PASS] G1 Boltz-2 rows unchanged: new (6152, 24) vs old (6152, 24)
[PASS] G2a WT rate reproduces: overlap n=236  protein-arm 0.822  ligand-arm 0.835  (full ligand-arm n=237 rate 0.835)
[PASS] G2b per-system concordance: 0.987  median |dRMSD|=0.00 A
[PASS] G3a pre-existing docking rows unchanged: (5473, 10) vs (5473, 10)
[PASS] G3b ligand SurfDock ingested: 1269 rows, 249 systems, ok=1268
GATE PASS
```
G2's median |ΔRMSD| of 0.00 Å means identical inputs produce effectively identical AF3 predictions across the two arms; the 3 of 236 discordant systems sit near the 2 Å threshold.

## Pass / fail
pass

## Code ↔ theory alignment
| claim / property checked | code (file:line / function) | justified by (ref ↓) | match? |
|---|---|---|---|
| ligand RMSD is in-place (Cα-superpose protein, apply R/t to ligand), symmetry-corrected, MCS atom map | `analysis/casf_mutagenesis/analysis.py` `_analyze_single_pose` | Masters et al. 2025 RMSD protocol; contrascf-casf skill | exact |
| top-1 = `pose_idx == 0`, ranked by model confidence | `analysis/ligand_mutagenesis/scripts/05_analyze.py` (rank enumerate over sorted CIFs) | Masters et al. 2025 | exact |
| AF3+MSA file layout `af3msa_{prefix}_model_*.cif` | `analysis/casf_mutagenesis/analysis.py:90` `MODEL_FILES` | runner `06_run_af3_msa_subset20.py` `_copy_all_samples` | exact |

## References
- Masters, M. R. et al. 2025. *Investigating whether deep learning models for co-folding learn the physics of protein–ligand interactions.* Nature Communications. (PDF shipped in `docs/`.)
- Abramson, J. et al. 2024. *Accurate structure prediction of biomolecular interactions with AlphaFold 3.* Nature 630, 493–500. doi:10.1038/s41586-024-07487-w
- Passaro, S. et al. 2025. *Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction.* bioRxiv. doi:10.1101/2025.06.14.659707
- Cao, D. et al. 2025. *SurfDock is a surface-informed diffusion generative model for reliable and accurate protein–ligand complex prediction.* Nature Methods 22, 310–322. doi:10.1038/s41592-024-02516-y
- Code: AaronXu9/contrasCF @ fix/af3-multichain-and-data-audit

## Optional: If failed — what to fix
[only fill if fail]
