# Surfdock Restored Physics Signature

**Type:** milestone
**Date:** 2026-08-28
**Status:** logged
**Project:** contrasCF

## Headline
After fixing the missing interface crop, SurfDock goes from a retracted total failure (0/244 WT under 2 Å) to **87.6 % WT under 2 Å at 1.06 Å median** — and now shows the **steepest WT→adversarial collapse of any method in the study** (+0.842), the cleanest physics signature in the cross-method matrix.

## Why it matters

The study's cross-method claim changes shape. SurfDock was previously written up as "not a reliable substitute for Vina-family docking" and excluded from interpretation; it is in fact the **best-separating** engine we have, holding both the highest WT ceiling of any docking method *and* the lowest adversarial rates:

| engine | WT <2 Å | rem | pack | inv | WT − mean(adv) |
|---|---|---|---|---|---|
| **SurfDock** | **0.876** | 0.029 | 0.046 | 0.026 | **+0.842** |
| GNINA | 0.729 | 0.138 | 0.130 | 0.117 | +0.601 |
| UniDock2 | 0.578 | 0.092 | 0.092 | 0.071 | +0.493 |

n = 249 / 238 / 238 / 231 (wt / rem / pack / inv) of 251 / 239 / 239 / 239 attempted; 12 cells excluded (1.24 %). Source: `analysis/casf_mutagenesis/outputs/docking_memorization.csv`, regenerated 2026-08-26.

Against the co-folding arms the contrast sharpens the memorization argument rather than blunting it: SurfDock's WT ceiling (0.876) essentially matches AF3+MSA (0.89), but its adversarial rate (0.026–0.046) is **~10× lower** than AF3+MSA's (0.16–0.31) and ~7× lower than Boltz-2's (0.17–0.24). A method can reach co-folding-grade accuracy on native structures while still collapsing on broken pockets — so the co-folding models' adversarial retention is not the price of accuracy.

This also closes a three-month-old wrong conclusion: the failure was never SurfDock or CASF, it was our own preprocessing handing the model ~10× the mesh it was trained on. Full diagnosis in postmortem `2026-08-28-surfdock-interface-crop`.

## Open threads
- [confirmed] Rewrite the retracted `docs/casf_overview.md` SurfDock section with these numbers, the per-variant n, and the 12 exclusions; the figure is already regenerated.
- [open] 8 of the 12 exclusions fail with "0 graphs" and the cause is unknown. Mesh size does **not** predict it (failed 42–158 vertices vs succeeded 7–166, `analysis/casf_mutagenesis/mesh_census.json`); 7 of 8 are `inv`, so it is variant-linked and nothing further is established.
- [open] Single seed (seed 42, 40 diffusion samples per complex). The WT/adversarial gap is far too large to be seed noise, but the per-variant adversarial rates (0.026 vs 0.046) should not be compared to each other without ≥3 seeds.
- [tentative] SurfDock being the sharpest discriminator makes it the natural physics reference for the pose-swap contrast, where GNINA's Vina term currently plays that role.
- [suggested] The ligand-mutation axis has no SurfDock arm; the runner would need to walk `analysis/ligand_mutagenesis/outputs/` the way `14_` walks the pocket module.
