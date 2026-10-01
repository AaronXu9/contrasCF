# surfdock-receptor-provenance-artifact

**Type:** postmortem
**Date:** 2026-09-04
**Severity:** experiment-invalidated

## What broke

The two **zero-mutation** systems (`3mss`, `4eo8` — the §4d generation no-ops) are unintentional but perfect controls for the receptor-provenance confound: their `rem`/`pack`/`inv` sequences are *identical to WT*, so any wt→variant change must be pure artifact. SurfDock moves **5.2–5.5 Å** on both. That casts doubt on how much of SurfDock's headline "+0.842, steepest WT→adversarial collapse of any method" ([[2026-08-28-surfdock-restored-physics-signature]]) is mutation response versus a crystal→predicted receptor swap.

## What I saw

Top-1 ligand RMSD across the four variants of the two zero-mutation systems:

```
system  method      wt    rem   pack    inv   max drift from wt
3mss    boltz2    2.15   2.15   2.15   2.15   0.00  (bit-identical)
3mss    af3msa    2.33   2.33   2.33   2.34   0.01
3mss    icm      10.54  10.54  10.53  10.54   0.01
3mss    gnina    10.57   9.82  10.09  10.07   0.75
3mss    unidock2  2.52   2.56   3.47   2.62   0.95
3mss    surfdock  1.61   7.04   7.02   7.09   5.48   <--
4eo8    boltz2     0.67   0.67   0.67   0.67   0.00  (bit-identical)
4eo8    af3msa     0.93   0.93   0.93   0.93   0.00
4eo8    icm        1.10   0.93   0.92   0.95   0.17
4eo8    gnina      1.03   1.14   1.20   1.14   0.17
4eo8    unidock2   1.13   2.52   1.11   2.21   1.39
4eo8    surfdock   0.99   6.22   6.22   4.53   5.23   <--
```

The co-folding rows returning bit-identical values is the positive control: same sequence in, same structure out, so the pipeline is deterministic and the comparison machinery is sound. The SurfDock rows are the anomaly.

Direct provenance check on the receptor files themselves — the WT docking receptor carries `1.00 0.00` B-factors and explicit hydrogens (a protonated crystal structure), while the mutant receptor carries pLDDT-like values (`78.16`, `87.58`) and no hydrogens:

```
outputs/1e66/wt/docking/receptor.pdb   ATOM 2  H  SER A 4  ... 1.00  0.00
outputs/1e66/rem/docking/receptor.pdb  ATOM 2  CA SER A 1  ... 1.00 78.16
```

## What I expected

That with zero mutations every engine would return its WT number for all four variants, as the co-folding arms and ICM do.

## Hypotheses for why

- **Stochastic seed noise** — rejected. SurfDock is run at fixed seed 42, and the effect is 5 Å on both systems in the same direction, not scatter.
- **Bad Cα alignment inflating the RMSD** — rejected for these cells: `mutant_receptor_alignment.csv` gives fits well under the 5 Å concern threshold for both systems.
- **SurfDock is unusually sensitive to receptor provenance** — **confirmed, mechanistically plausible.** SurfDock consumes an MSMS molecular-*surface* mesh, and a surface is a second-order function of atomic coordinates: sub-Å backbone deviations between a crystal and an AF3 prediction move the surface far more than they move Cα positions. The Vina-family engines (GNINA, UniDock2) read atoms directly and drift ≤ 1.4 Å; ICM does not drift at all.

## Root cause

`scripts/10_build_mutant_docking.py:135` builds every mutant docking receptor from `af3msa_<prefix>_model_0.cif` — an **AF3+MSA co-folded holo prediction** with the ligand stripped afterwards — while the WT cell docks into the **protonated crystal** receptor. So the WT→mutant contrast confounds two changes at once: the mutation, and crystal→predicted receptor provenance. This is TODO §4a, previously logged as "mild"; the zero-mutation controls show it is not mild for SurfDock.

## Fix

None yet — this is a measurement, not a code defect. The controls were surfaced by reading `n_mutations` out of the manifest while building the representative-system selector ([[2026-09-04-structural-viz-and-physics-validation]]).

## Prevention

- [confirmed] `3mss` and `4eo8` are now documented as standing zero-mutation controls: any engine that moves on them is showing provenance sensitivity, not mutation response.
- [open] Run the decisive experiment: dock WT into the **AF3-predicted WT** receptor for all 251 systems, per engine. That isolates provenance from mutation with n = 251 instead of n = 2. Inputs already exist (`outputs/<sys>/wt/af3msa_<sys>_wt_model_0.cif`); it needs a variant of `10_build_mutant_docking.py` pointed at the wt cell plus one sweep per engine.
- [open] Until that runs, SurfDock's WT→adversarial gap must be quoted with the caveat that an unknown fraction is receptor provenance.
- [suggested] Report a surface-mesh similarity (vertex count / Hausdorff distance, crystal vs predicted) alongside Cα RMSD; Cα agreement demonstrably does not imply surface agreement, which is the quantity SurfDock actually consumes.

## Blast radius

- [tentative] `2026-08-28-surfdock-restored-physics-signature` claims SurfDock shows "the steepest WT→adversarial collapse of any method in the study (+0.842)" and "the cleanest physics signature in the cross-method matrix". The direction of that claim likely survives — SurfDock's adversarial rates are 0.024–0.046 and the effect is far larger than 5 Å on most systems — but the *magnitude* is inflated by an unquantified provenance component and should not go into a paper as-is.
- [confirmed] GNINA (≤ 0.75 Å) and UniDock2 (≤ 1.39 Å) drift is small; their numbers are not materially affected.
- [confirmed] Co-folding arms are entirely unaffected — they never consume a docking receptor.
- [confirmed] No CSV is invalid; the cells are correct measurements of a confounded comparison.
- [open] Whether the SurfDock effect generalises beyond n = 2 is unknown until the WT-into-predicted-WT sweep runs.
