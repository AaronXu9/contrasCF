# surfdock-ligand-arm-input-footprint

**Type:** postmortem
**Date:** 2026-10-01
**Severity:** [data-loss | experiment-invalidated | time-lost-hours | cosmetic]

## What broke
The ligand-arm SurfDock sweep (lab, 2026-09-04, 1265/1300 cells ok) is not a valid measurement of SurfDock. Every cell's surface patch was defined by a wrongly oriented input ligand.

## What I saw
WT-correct rate (top-1 < 2 Å) on cells that are physically identical across arms — same crystal receptor, same native ligand:

| engine | protein arm WT | ligand arm WT |
|---|---|---|
| SurfDock | 218/249 = 0.876 | **43/246 = 0.175** |
| GNINA | 183/251 = 0.729 | 157/251 = 0.625 |
| UniDock2 | 145/251 = 0.578 | 111/249 = 0.446 |

The ligand arm's `docking/ligand.sdf` is an RDKit embedding whose centroid sits 22–173 Å from the crystal ligand (12/12 spot-checked); the protein arm's is the crystal pose itself (centroid offset 0.00 Å). The 2026-09-04 smoke test compared only `1bcu`, a system SurfDock fails in both arms (4.06 vs 5.65 Å), so it could not detect this.

## What I expected
Ligand-arm SurfDock WT ≈ protein-arm WT (~0.88), since the WT cells are the same complex.

## Hypotheses for why
1. SurfDock's pocket/interface definition depends on input ligand coordinates — **confirmed**: `_translate_ligand_to_pocket` (`analysis/scripts/13_run_surfdock.py:149`) moves the centroid to the box center but keeps the embedding's arbitrary orientation, and the interface crop keeps mesh faces within 3 Å of those atoms.
2. Interface-crop regression (the 2026-08-25 bug) — not the cause: the crop code is present and identical to the one producing the healthy protein-arm run.

## Root cause
`ligand_mutagenesis/build.py` writes a freshly RDKit-embedded `ligand.sdf` per variant (correct for GNINA/UniDock2, which take only `box.json` as the site and use the ligand as a starting conformer). SurfDock additionally uses the ligand's atom positions to choose which surface to show the model. With a randomly oriented conformer, the 3 Å interface shell covers a different, partly wrong patch than the crystal footprint the model was evaluated with — a train/serve mismatch of the same kind as the August crop bug, introduced by input provenance rather than code.

## Fix
Not yet applied — needs a decision (TODO item 15). Proposed: give SurfDock the **crystal ligand pose** as its pocket-defining reference for every variant (the variant shares the scaffold, so the pocket is the same), and the variant molecule as the thing to dock. Requires separating "pocket reference ligand" from "ligand to dock" in `14_run_surfdock_variants.py` and re-running ~1300 cells (~13 GPU-h on the lab 4090; CARC cannot do surface prep, TODO item 14).

## Prevention
- [confirmed] A cross-arm smoke test must use systems the engine SOLVES in the reference arm; a pair of failures cannot show a difference.
- [confirmed] Before any new sweep, compare the WT rate against the other arm's WT on identical complexes — a > 0.1 gap is a stop signal, not a finding.
- [tentative] Record, per engine, which input fields it actually consumes (box only vs ligand coordinates); the arms differed in an input only one engine reads.

## Blast radius
- [confirmed] All 1265 ligand-arm SurfDock cells are INVALID as a SurfDock measurement. Exclude SurfDock from every ligand-arm comparison until re-run; `memorization_ligand.csv` keeps the rows but they must not be quoted.
- [confirmed] GNINA and UniDock2 ligand-arm results remain valid: both take the site from `box.json`, which is identical across arms (verified on 1e66).
- [open] Secondary protocol difference, not a bug: protein-arm WT docking starts from the crystal conformer (re-docking), the ligand arm from an RDKit conformer. That plausibly explains GNINA/UniDock2's lower ligand-arm WT rates (0.625 vs 0.729; 0.446 vs 0.578). Cross-arm docking comparisons must state it; within-arm comparisons are unaffected.
