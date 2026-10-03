# perturbation-strength-blind-spot

**Type:** postmortem
**Date:** 2026-09-04
**Severity:** experiment-invalidated

## What broke

The first batch of structural "memorization exhibits" was selected by outcome alone (WT solved, adversarial still < 2 Å) with **no filter on how many residues were actually mutated**. Because retention is strongly dose-dependent, that criterion preferentially surfaces the *weakest possible evidence*. Two of the twelve rendered systems, including the one I presented as the headline case, are not evidence of anything:

| system | mutations applied | what I claimed |
|---|---|---|
| `3u5j` | **1** (`N140F`) | "reproduces the crystal pose to 0.25/0.63/0.77 Å in a pocket mutated to all-Gly, all-Phe and Miyata-inverted" |
| `4eo8` | **0** (generation no-op, TODO §4d) | included silently in the "memorized" set |

A separate, second defect surfaced in the same pass: the receptor→crystal Cα fit used to place mutant poses is never quality-checked, and a poor fit inflates the reported ligand RMSD.

## What I saw

`outputs/manifest_full_casf.json` records `variants.<v>.n_mutations` per system. Reading it across the 251-system set:

```
rem  n_mutations: median 5   mean 5.8   min 0   max 13
     <=2 mutations: 11.6 %   <=5: 50.6 %
     two systems (3mss, 4eo8) have ZERO
```

So "every pocket residue is substituted" — the phrasing carried in the module docstrings and in my own summaries — describes a median of **five** residues, and `3u5j`'s "destroyed pocket" is a single point mutation.

Alignment audit over all 717 mutant docking cells (`outputs/mutant_receptor_alignment.csv`):

```
median Ca-fit 0.95 A
Ca-fit >  2 A : 150  (20.9 %)
Ca-fit >  5 A :  47  ( 6.6 %)
Ca-fit > 10 A :  30  ( 4.2 %)
worst: 5c2h/inv 18.63, 2wn9/rem 18.39, 4bkt/inv 15.01, 4w9h/rem 14.52
```

`4w9h` was in the rendered set; its "56 Å ejection" is frame error, not displacement.

## What I expected

That an outcome-based filter (`memorized`) would surface the systems that best carry the argument. It does the opposite: conditioning on retention while retention is dose-dependent is conditioning on *low perturbation*.

## Hypotheses for why

- **The selector was buggy** — rejected. `pick_systems` did exactly what it was written to do; the criterion itself was wrong.
- **Pocket size is roughly constant across CASF so the filter would not matter** — rejected by the manifest: range 0–13, and retention runs 78 % → 15 % across that range ([[2026-09-04-memorization-is-dose-dependent]]).
- **The docstring wording was merely loose** — partly true and worse than loose: "every pocket residue is substituted" is accurate about the *rule* (every residue within 3.5 Å) but reads as a claim about *magnitude*, and I propagated it as magnitude.

## Root cause

`POCKET_CUTOFF_A = 3.5` (side-chain heavy atom → ligand heavy atom) selects a small contact shell, so the perturbation is much weaker than the variant names suggest. **Perturbation strength was never carried alongside the RMSD** — not in `docking_results.csv`, not in `results_full.csv`, not on any figure — so nothing in the pipeline could flag that a "memorized" cell had one mutated residue.

## Fix

- `scripts/33_render_mutation_views.py`: `--min-mutations` (default 5, the dataset median) filters the pick; candidates are ranked **most-mutated first**; the mutation count is printed in every figure title; new `--pick low-mutation | high-mutation | extremes` modes select representative perturbation extremes with a complete 4×4 grid, Cα fit < 3 Å and all four methods solving WT.
- `outputs/mutant_receptor_alignment.csv` (new): Cα-fit quality for all 717 mutant cells, so a bad frame is detectable.
- Code carried by [[2026-09-04-structural-viz-and-physics-validation]].

## Prevention

- [confirmed] Mutation count is printed on every structural figure and is a required companion to any quoted RMSD.
- [confirmed] `--min-mutations` defaults to 5 rather than 0, so the failure mode is now opt-in rather than default.
- [confirmed] `manifest.json` in each figure directory records the per-cell Cα fit; > 5 Å means the RMSD is frame error.
- [confirmed] `3mss` and `4eo8` are documented as zero-mutation controls, not as data — see [[2026-09-04-surfdock-receptor-provenance-artifact]] for the use they *are* good for.
- [open] `docking_results.csv` and `results_full.csv` still do not carry `n_mutations` or `ca_fit_rmsd_a`; joining them by hand is the remaining foot-gun. Add both columns in `12_analyze_docking_engines.py` and `05_analyze_subset20.py`.
- [suggested] Add a generator guard that refuses to emit a variant whose mutation spec equals WT (the §4d no-op), rather than measuring the damage afterwards.

## Blast radius

- [confirmed] No published CSV or aggregate rate is wrong — the defect is in *selection for display* and in my prose, not in the computation. `memorization_full.csv` and `docking_memorization.csv` stand.
- [confirmed] The first 12 rendered systems were re-selected and all 17 re-rendered; `3u5j` and `4eo8` are retained but labelled with their true counts (1 and 0).
- [confirmed] Aggregate impact of the bad Cα fits is negligible: filtering to fit < 2 Å moves SurfDock 2.4 → 2.4 %, GNINA 18.0 → 19.4 %, UniDock2 14.1 → 14.2 %. Individual cells (`4w9h`, `5c2h`, `2wn9`, `4bkt`) are untrustworthy.
- [tentative] Any earlier verbal or written description of the perturbation as "every pocket residue" / "the pocket that held it destroyed" overstates magnitude and should be restated as "every residue within 3.5 Å — a median of 5".
- [open] `docs/casf_overview.md`, `docs/casf_confidence.md` and the two Feishu docs still carry pooled rates with no dose stratification.
