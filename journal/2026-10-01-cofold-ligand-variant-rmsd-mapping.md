# cofold-ligand-variant-rmsd-mapping

**Type:** postmortem
**Date:** 2026-10-01
**Severity:** [data-loss | experiment-invalidated | time-lost-hours | cosmetic]

## What broke
Every co-folding (Boltz-2, AF3+MSA) ligand RMSD on a ligand-mutagenesis **variant** cell was computed with a wrong atom correspondence, in the shared scorer `analysis/casf_mutagenesis/analysis.py` (`_matched_rmsd`, `_bestfit_rmsd`). Since 2026-05-23 for Boltz-2; caught 2026-10-01 while building the five-method retention table.

## What I saw
Five-method retention at Δheavy = 1 (one added atom): GNINA 0.71, UniDock2 0.65, SurfDock 0.47, but **Boltz-2 0.037 and AF3+MSA 0.034**, while both co-folding models solve ~83% / ~59% of wild types. Smoking gun: the *self-superposed* (pocket-blind) RMSD jumped from 0.76 Å on wt to **3.73 Å on halo_F_1** for Boltz-2 (AF3+MSA 0.47 → 3.71). A single fluorine cannot change scaffold geometry by Å; a wrong atom map can. `n_heavy_matched / n_heavy_native` read 1.00, so the rows looked fully matched.

## What I expected
Co-folding to retain the pose under a one-atom edit at least as often as docking — Masters et al.'s central finding is that co-folding models largely ignore ligand perturbations.

## Hypotheses for why
1. Atom correspondence fails when crystal and variant differ in heavy-atom count (consistent: wt fine, every variant broken). **Confirmed.**
2. Frame error from the Cα superposition — rejected: `ca_rmsd_a` medians are unchanged (0.67 vs 0.75 Å), and the bestfit metric, which ignores frames, is also broken.
3. Real model behaviour — rejected by the bestfit argument above.

## Root cause
`_matched_rmsd` correctly enumerated `pred.GetSubstructMatches(crystal)` (crystal is a substructure of a halogenated/methylated variant), then a guard

```python
if pred_pts.shape[0] != n_match_atoms or crystal_pts.shape[0] != n_match_atoms:
    d = pred_xyz - crystal_xyz          # after truncating both to min length
```

threw every valid match away whenever the counts differed and fell back to pairing atoms by **file order**, reporting `n_matched = n_crystal`. `_bestfit_rmsd` had the same guard. The reverse-direction branch (`crystal ⊃ pred`) was separately wrong: it replaced matches with the identity mapping. Charge-swap variants, where neither ligand contains the other, had no correct path at all. The protein arm never hit any of this because its ligand is identical to the crystal, and the 2026-10-01 validation gate (`ligand-af3msa-scoring-gate`) only checked wt cells, where counts are equal — the gate had a blind spot exactly where the bug lived.

## Fix
New `_atom_correspondences(crystal_heavy, pred_heavy)` returns all symmetry-equivalent (crystal_idx, pred_idx) pairs in one of three modes — `pred_superset` (variant adds atoms), `crystal_superset` (removes atoms), `mcs` (charge swaps; elements compared, any bond) — and both scorers minimise over them on the matched atoms only. No correspondence now **raises** (cell → `status=error`) instead of returning a number. Gated by `cofold-variant-rmsd-fix-gate`: synthetic known-answer cases (+F, −1 atom, element swap) give 0.00 Å in all three modes; the OLD code gives 3.47 Å on the +F case once atoms are reordered as an independent SDF writer would; 276 protein-arm cells re-scored match `results_full.csv` to 0.0000 Å.

## Prevention
- [confirmed] Any validation gate for a scorer must include at least one cell from every input shape it handles — here wt (equal size), additive, subtractive and MCS variants — not only the stratum with a reference value.
- [confirmed] Report the atom-correspondence mode per cell and never let a fallback report a full match count; a silent fallback that looks like success is the failure pattern here and in the SurfDock crop bug.
- [tentative] Add a standing sanity check to the ligand analyzer: flag any variant whose pocket-blind bestfit RMSD exceeds its wt's by > 1.5 Å, since a small ligand edit cannot do that.

## Blast radius
- [confirmed] `results_ligand.csv` variant rows (Boltz-2 since 2026-05-23; AF3+MSA added today) — wrong `ligand_rmsd_a`, `ligand_rmsd_fullca_a`, `bestfit_rmsd_a`. Regenerated with the fix.
- [confirmed] Every ligand-arm product of `07_confidence_response.py` (`s2a_*`, `s2b_*`, `s2c_*`, `s3_*`, `s4_*`, `paired_conf_aff_rmsd_ligand.csv`, `confidence_*_ligand.csv`, `figures/conf_aff_rmsd_ligand.*`) stratified memorized/responded by that RMSD. Must be regenerated.
- [confirmed] The claim in `journal/2026-06-30-three-head-dissociation.md` item 8 — "the confidence-registers pattern holds on the ligand axis too, but severity-gated" — rests on those strata and is **unverified** until re-run.
- [confirmed] NOT affected: the protein arm (identical ligands, verified 276/276 unchanged), all docking numbers (separate MCS matcher in `gnina_analysis.py`), the Lark draft and `casf_confidence.md` (both scoped to the pocket axis, ligand axis listed as TODO).
- [confirmed] NOT affected: the paper-reproduction arm. Its in-place RMSD (`analysis/src/pipeline.py:74,82`) uses explicit MCS pair lists from `ligand_match.match_all` / `match_common`, never a file-order fallback. Its separate `ligand_rmsd_bestfit` failure on modified ligands (GetBestRMS) is already TODO item 11.
