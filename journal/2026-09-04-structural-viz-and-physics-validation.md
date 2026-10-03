# structural-viz-and-physics-validation

**Type:** iteration
**Date:** 2026-09-04
**Code:** `9160bbc87e6c8f83a0af76049959d80745b53aed` on `fix/af3-multichain-and-data-audit` (dirty)

## Plan reference
**Plan:** ad-hoc
**Step:** n/a

## Change

Three new scripts. `33_render_mutation_views.py` (667 lines) renders the pocket contrast as a **4×4 PyMOL grid** — ICM / SurfDock / Boltz-2 / AF3+MSA × wt/rem/pack/inv — with every object placed in the crystal frame by the *same* transforms the RMSD analysis uses. `34_validate_cofold_pockets.py` (439) runs PoseBusters plus a side-chain rotamer/orientation analysis over the co-folded complexes. `35_rescore_crystal_pose.py` (295) rescores fixed poses with GNINA `--score_only` across four arms of decreasing model contamination.

## Why

The module had aggregate plots but no structural viz at all: the only pose exporters (`23_`, `24_`) took one `--system` at a time, wrote to `/tmp`, and covered the pose-swap ejection geometry rather than the mutation contrast. Nothing let anyone *look* at what "memorized in a destroyed pocket" means. Once the pictures existed they raised a physics question the RMSD cannot answer — whether the retained pose is physically legal — which is what `34_` and `35_` were built to settle. Findings in [[2026-09-04-memorization-is-dose-dependent]].

## Why this design over alternatives

Exporting the *analysis-time* transform to disk rather than re-deriving it in PyMOL: the stored mutant SDFs live in the AF3-predicted receptor frame, and the Cα alignment is applied inside `analyze_gnina` and never persisted, so loading them directly puts the ligand 80–90 Å from the pocket — I made exactly that error and had to retract a table before catching it. Rejected a "load everything and align in PyMOL" design for that reason: the alignment must be the pipeline's own, or the figure and the number disagree. For `35_` I rejected an all-in-one arm and built four, because the obvious arm (crystal pose into the AF3 mutant receptor) turned out to be dominated by Cα-fit error — a 1.3 Å fit produced a spurious +111 kcal/mol on `3zso` — so the alignment-free arms (`*_self`, `crystal_truncated_in_place`) are the only quotable ones.

## Files touched

- `analysis/casf_mutagenesis/scripts/33_render_mutation_views.py` — new. `docking_receptor()` writes one shared aligned receptor per cell (both engines were handed the same one); `export_icm()` deliberately does **not** transform (ICM docked into a pre-aligned receptor, so its poses are already in the crystal frame); `export_cofold()` handles both co-folding models; `_write_pdb()`/`_write_sdf()` strip hydrogens so the protonated-crystal WT receptor does not render ~2× the sticks of its own mutant panels; `pick_systems()` gained `--min-mutations` and `low-mutation`/`high-mutation`/`extremes` modes.
- `analysis/casf_mutagenesis/scripts/34_validate_cofold_pockets.py` — new. `CORE_CHECKS` restricts PoseBusters to the 11 non-vacuous checks; `residue_chis()`/`chi1_offset()` compute rotamer strain; `crystal_baseline()` supplies the experimental control.
- `analysis/casf_mutagenesis/scripts/35_rescore_crystal_pose.py` — new. `crystal_rem_receptor()` builds the model-free `rem` control by truncating crystal residues to backbone; `score_only(..., cleanup=)` deletes per-cell temp files (an earlier run leaked ~6 000 of them and died at 240/251 with ENOSPC); results stream to CSV per row so a multi-hour run stays readable.
- `.gitignore` +8 — keeps `*_contact.png` and `manifest.json` under `figures/structures/`, drops the regenerable PDB/SDF/per-scene renders (88 MB on disk, 19 MB tracked).

## Diff snapshot

The CLI's auto-captured diff is against `HEAD~1` and shows unrelated AF3/CARC work, not this change. The relevant delta is three new untracked scripts plus the `.gitignore` rule:

```diff
+ analysis/casf_mutagenesis/scripts/33_render_mutation_views.py   | 667 ++++++++++
+ analysis/casf_mutagenesis/scripts/34_validate_cofold_pockets.py | 439 +++++++
+ analysis/casf_mutagenesis/scripts/35_rescore_crystal_pose.py    | 295 ++++++
  .gitignore                                                      |   8 ++++
--- a/.gitignore
+++ b/.gitignore
+analysis/casf_mutagenesis/figures/structures/**
+!analysis/casf_mutagenesis/figures/structures/**/
+!analysis/casf_mutagenesis/figures/structures/**/*_contact.png
+!analysis/casf_mutagenesis/figures/structures/**/manifest.json
```

## Validation evidence

Inline gate on `1e66`, re-reading the *exported* files with the pipeline's own RMSD code — this is the check that the persisted transform is correct, not just the in-memory one:

```
             reported   recomputed from exported files
surfdock wt   0.33        0.32     (0.01 = symmetry-equivalent MCS mapping on re-read)
surfdock rem  3.49        3.49
surfdock pack 2.39        2.39
surfdock inv 12.42       12.42
boltz2   wt   0.62        0.62
boltz2   rem  1.47        1.47
boltz2   pack 1.14        1.14
boltz2   inv  3.08        3.08
```

Second gate, `34_`: the crystal rotamer baseline returns a 1.1 % χ1 outlier rate over 6 726 residues, matching the ~1 % expected from the Lovell library — so the dihedral code and the 40° tolerance are calibrated before being applied to predictions.

## Code ↔ theory alignment

| changed behaviour | code (file:line) | justified by (ref ↓) | match / divergence |
|---|---|---|---|
| mutant docking poses Cα-transformed to crystal frame; ICM poses are not | `33_render_mutation_views.py:docking_receptor,export_icm` | mirrors `gnina_analysis.py:_ca_transform`; ICM used a pre-aligned receptor per `30_analyze_icm.py` docstring | exact |
| co-folding aligned by ligand-near-chain Cα, not largest chain | `33_render_mutation_views.py:export_cofold` | mirrors `analysis.py:_analyze_single_pose`; needed for homo-multimers (4w9l) | exact |
| PoseBusters run without cofactor/water checks | `34_validate_cofold_pockets.py:CORE_CHECKS` | Buttenschoen 2024 `dock` config | **divergence**: co-folding emits no waters/cofactors, so those checks fail 100 % for every variant incl. wt; including them would report 0 % validity everywhere |
| χ1 outlier = >40° from −60/60/180 | `34_validate_cofold_pockets.py:CANONICAL_CHI1,CHI1_TOL` | Lovell et al. 2000 penultimate rotamer library | approximate — a tolerance band, not a full library lookup; calibrated against the crystal baseline |
| model-free `rem` = delete pocket side chains from crystal | `35_rescore_crystal_pose.py:crystal_rem_receptor` | `rem` is →Gly, i.e. pure deletion (Masters et al.) | exact — no packing needed, so no model and no superposition enters |
| hydrogens stripped from all exported structures | `33_render_mutation_views.py:_write_pdb` | cosmetic only; RMSD is heavy-atom throughout | exact |

## References

- Masters, M. R. et al. (2025). Pocket-mutagenesis counterfactual protocol (`rem`/`pack`/`inv`). Repo PDF; Zenodo 14749304.
- Buttenschoen, M., Morris, G. M., Deane, C. M. (2024). PoseBusters. *Chem. Sci.* 15, 3130. DOI 10.1039/D3SC04185A.
- McNutt, A. T. et al. (2021). GNINA 1.0. *J. Cheminform.* 13, 43. DOI 10.1186/s13321-021-00522-2.
- Lovell, S. C. et al. (2000). The penultimate rotamer library. *Proteins* 40, 389.
- Schrödinger, LLC. PyMOL 3.0.0 — `/home/aoxu/miniconda3/envs/PyMOL-PoseBench/bin/pymol`, headless `-cq`.

## Outcome

Works — 17 systems rendered as 4×4 grids, 916 PoseBusters cells, 3 838 rotamers, ~3 000 rescoring cells.

- [confirmed] `--min-mutations` defaults to 5; unfiltered selection surfaces the weakest evidence, see [[2026-09-04-perturbation-strength-blind-spot]].
- [confirmed] Only the `*_self` and `crystal_truncated_in_place` rescoring arms are alignment-free; `crystal_in_docking_receptor` must not be quoted.
- [open] `35_` has not finished the last systems of the 251-system sweep; the aggregates quoted so far are from 35 systems and must be restated on the full set.
- [open] Scripts are untracked; nothing in this iteration is committed.
- [suggested] Fold `33_`/`34_` into the `contrascf-casf` skill's standard workflow so structural viz runs alongside the aggregate plots rather than on request.
