# ligand-arm-five-method-memorization

**Type:** experiment
**Kind:** experiment
**Date:** 2026-10-01
**Status:** done

## Plan reference
**Plan:** plans/2026-10-01-organize-codes-results.md
**Step:** S1

## Goal
How does top-1 pose retention under a ligand edit vary with edit size (Δ heavy atoms) and chemistry type across all five methods on the CASF ligand-mutagenesis arm — Boltz-2, AF3+MSA, GNINA, UniDock2, SurfDock?

## Hypothesis

### Quantitative prediction
Metric: WT-conditioned top-1 retention, P(ligand RMSD < 2 Å | method solved this system's WT < 2 Å). **This is retention, not memorization**: keeping the pose after a one-atom halogen swap is plausibly the physically correct answer, so the protein arm's "low is good" reading does not transfer.

Disclosure: GNINA and UniDock2 per-variant ligand rates were already inspected on 2026-09-04 (GNINA halo 0.63-0.78, meth ladder 0.69→0.00, chrg_pos_1 0.17). AF3+MSA and SurfDock ligand rates, and Boltz-2's WT-conditioned ligand table, have not been looked at.

- **H1 (dose-response, all five methods):** retention pooled over Δheavy = 1 (halo_* + meth_1) exceeds retention pooled over Δheavy ≥ 3 (meth_3..5, chrg_neu_*, chrg_pos_2/3) by ≥ 0.15 for every method.
- **H2 (co-folding retains more at small edits):** at Δheavy = 1, both co-folding models retain ≥ 0.80, and each exceeds every docking engine.
- **H3 (charge costs more per atom, tentative):** `chrg_pos_1` (Δheavy 2) retains less than `meth_2` (Δheavy 2) for each docking engine.

### Falsification criterion
H1 is refuted for a method if its Δheavy≥3 retention is ≥ its Δheavy=1 retention minus 0.15. H2 is refuted if either co-folding model retains < 0.80 at Δheavy=1 or falls below any docking engine there. H3 is refuted for an engine if chrg_pos_1 ≥ meth_2.

### Pre-registered analysis
Top-1 everywhere (co-folding `pose_idx==0`; docking rank-1 SDF). Conditioning set = systems where that method's WT top-1 < 2 Å. Rates per (method, variant) with Wilson 95% CI; pooled Δheavy strata per method. Δheavy per cell = heavy atoms(variant SMILES) − heavy atoms(WT SMILES) from `manifest_full.json`. Comparisons are descriptive with CIs; cells with n_wtok < 10 are reported but not used to decide H1–H3. Failure cells (missing_cif / docking error) are counted, never silently dropped.

### Decision criteria
| Outcome | Action |
|---|---|
| H1 holds for all 5 | [confirmed] state "every method reads ligand chemistry dose-dependently" in the draft with the Δheavy table |
| H2 holds | [tentative] state co-folding retains more under small edits, explicitly framed as retention not memorization |
| H2 refuted | drop any co-folding-vs-docking claim for the ligand arm |
| H3 holds for ≥2 of 3 engines | report charge-specific cost as a secondary finding |

## Setup
- **Code:** `8c756b400ee432a50dd11d28b8f1e70e5fbe6a32` on branch `fix/af3-multichain-and-data-audit` (dirty)
- **Env:** `unidock2` (Python 3.10.12)
- **Host:** katlab
- **GPU:** GPU 0: NVIDIA GeForce RTX 4090 (UUID: GPU-1c301dd0-32b4-0ba7-bea0-e4f51afea09d)
- **pip freeze hash:** `n/a`
- **Data:** `analysis/ligand_mutagenesis/outputs/manifest_full.json` (251 systems, 1300 cells); co-folding `results_ligand.csv`; docking `analysis/casf_mutagenesis/outputs/docking_results.csv` rows with `module == ligand`
- **Config:** `analysis/ligand_mutagenesis/scripts/08_retention_table.py` (THRESH = 2.0 Å, top-1, WT-conditioned)
- **Determinism:** analysis is deterministic; predictions are single-seed (see Seeds)

## Model & data provenance
| artifact | exact id / name | size | source (repo@sha · HF id · DOI · URL) | version / sha256 | role |
|---|---|---|---|---|---|
| AlphaFold3 + MSA | AF3 v3.0.1, `af3.bin` | ~1.1 GB weights | github.com/google-deepmind/alphafold3; CARC dockStrat fork | v3.0.1 | co-folding, ligand arm (CARC 11750008/9) |
| Boltz-2 | boltz 2.2.1 | n/a | github.com/jwohlwend/boltz | 2.2.1 | co-folding |
| GNINA | gnina binary | n/a | github.com/gnina/gnina; `/home/aoxu/projects/PoseBench/forks/GNINA/gnina` | seed 42, CNN rescore | docking |
| UniDock2 | unidock2 | n/a | github.com/dptech-corp/Uni-Dock; conda env `unidock2` | seed 42 | docking |
| SurfDock | SurfDock weights `docking`/`posepredict` | 140 MB | github.com/CAODH/SurfDock | CogLigandBench `2107de18e9` | docking — **INVALID on this arm**, see postmortem |
| CASF-2016 ligand arm | `manifest_full.json` | 251 systems / 1300 cells | contrasCF `analysis/ligand_mutagenesis` | built 2026-05-07 | evaluation set |

## Validation gates
- `ligand-af3msa-scoring-gate` — AF3+MSA scoring + SurfDock ingestion leave prior rows unchanged; WT reproduces protein arm (0.835 vs 0.822, concordance 0.987) — **pass**. Blind spot: only exercised WT cells.
- `cofold-variant-rmsd-fix-gate` — atom-correspondence fix: synthetic known-answer cases 0.00 Å in all three modes, OLD code shown wrong (3.47 Å), protein-arm sample unchanged, halo_F bestfit returns to WT level (0.77 vs 0.76; 0.46 vs 0.47) — **pass**.

## Commands & run log
Captured mechanically by `lab-notebook run <slug> -- <command>` → see `commands.log`
in this entry's folder; the latest run's timing / exit code / peak memory is in `meta.yaml`.

## Seeds
Single prediction per cell: AF3+MSA fixed model seed; docking seed 42; Boltz-2 5 samples, top-1 by confidence. Variance is across 251 systems (Wilson CIs), not seeds. Seed sensitivity for AF3+MSA was measured on the protein arm (`28_af3_multiseed_check.py`: inv 0.529 ± 0.020 over 3 seeds); treat ligand-arm AF3+MSA rates as single-seed, preliminary at that ±0.02 scale.

## Results
Corrected data (after the atom-mapping fix). WT-conditioned top-1 retention < 2 Å, Wilson 95% CI, pooled by Δheavy stratum (`retention_by_dheavy.csv`):

| method | Δheavy = 1 | Δheavy = 2 | Δheavy ≥ 3 |
|---|---|---|---|
| AF3+MSA | **0.899** [0.872, 0.922] n=567 | 0.400 [0.290, 0.521] n=65 | 0.318 [0.259, 0.383] n=214 |
| Boltz-2 | **0.820** [0.780, 0.855] n=401 | 0.424 [0.272, 0.592] n=33 | 0.333 [0.254, 0.423] n=117 |
| GNINA | 0.709 [0.666, 0.748] n=470 | 0.278 [0.176, 0.409] n=54 | 0.223 [0.166, 0.292] n=166 |
| UniDock2 | 0.653 [0.599, 0.703] n=317 | 0.219 [0.110, 0.388] n=32 | 0.104 [0.059, 0.176] n=106 |
| SurfDock | *invalid* (0.468) | *invalid* | *invalid* |

Charge vs neutral at equal size (Δheavy = 2): chrg_pos_1 vs meth_2 — AF3+MSA 0.342 vs 0.500, Boltz-2 0.364 vs 0.545, GNINA 0.172 vs 0.370, UniDock2 0.167 vs 0.286 (all CIs overlap).

Methylation ladder, AF3+MSA: meth_1 0.849 → meth_2 0.500 → meth_3 0.333; Boltz-2 0.828 → 0.545 → (n=3). Before the fix the co-folding Δheavy=1 cell read 0.034 / 0.037.

### Per-seed
Single prediction per cell (see Seeds); no per-seed breakdown. Uncertainty is across systems.

### Aggregate
**single-seed, preliminary** with respect to model stochasticity (AF3+MSA seed spread measured at ±0.02 on the protein arm). Across-system uncertainty: Wilson 95% CIs in the table above; the Δheavy=1 vs Δheavy≥3 gaps (0.49–0.58) exceed that seed spread by > 20×, so H1 and H2 are robust to it; H3's gaps (0.12–0.20) are within overlapping CIs and stay tentative.

## Confound checks
- **SurfDock input footprint (invalidating):** ligand-arm `ligand.sdf` is an RDKit embedding 22–173 Å from the site; SurfDock's interface crop follows those atoms after only a centroid shift. WT 0.175 vs 0.876 on the protein arm. Excluded from all verdicts. Postmortem `surfdock-ligand-arm-input-footprint`.
- **Atom mapping (fixed):** postmortem `cofold-ligand-variant-rmsd-mapping`; all co-folding numbers here are post-fix.
- **Conditioning sets differ per method** (WT-correct: AF3+MSA 198, GNINA 157, Boltz-2 136, UniDock2 111), so each rate is over a different system subset; co-folding vs docking compares retention given each method's own success.
- **Docking start conformer:** RDKit embedding here vs crystal conformer on the protein arm; affects cross-arm docking comparisons only.
- **Retention ≠ memorization:** the true pose of a variant is unknown. High retention at Δheavy = 1 is plausibly correct physics; high retention at Δheavy ≥ 3 is ambiguous between memorization and genuine robustness.

## Interpretation
- **H1 holds for all four valid methods:** Δheavy=1 minus Δheavy≥3 = 0.58 (AF3+MSA), 0.49 (Boltz-2), 0.49 (GNINA), 0.55 (UniDock2), all ≥ 0.15 with non-overlapping CIs. Every method reads ligand chemistry dose-dependently, co-folding included.
- **H2 holds, Boltz-2 marginally:** at Δheavy=1 AF3+MSA 0.899 and Boltz-2 0.820 both exceed GNINA 0.709 and UniDock2 0.653; Boltz-2's CI lower bound (0.780) sits below the 0.80 threshold.
- **H3 holds directionally for both valid docking engines and both co-folding models**, but every pair's CIs overlap at n = 11–38. Report as tentative.
- Co-folding retains the pose more than docking at every Δheavy, including ≥ 3 (0.32–0.33 vs 0.10–0.22). Whether that excess at large edits is memorization or better modelling cannot be decided from RMSD-to-crystal alone.

## Compute cost
Analysis only: 05_analyze.py 212 s (post-cache), 08_retention_table.py < 1 s, 12_analyze_docking_engines.py 307 s, on katlab CPU. Predictions: AF3+MSA ~26 a40 task-hours on CARC; SurfDock ~12 h on the lab 4090 (invalid).

## Code ↔ theory alignment
| claim / property | code (file:line / function) | justified by (ref ↓) | match? |
|---|---|---|---|
| in-place, symmetry-corrected ligand RMSD on the shared scaffold | `analysis/casf_mutagenesis/analysis.py` `_atom_correspondences`, `_matched_rmsd` | Masters et al. 2025 RMSD protocol | exact (post-fix) |
| WT-conditioned retention, top-1 | `analysis/ligand_mutagenesis/scripts/08_retention_table.py` `main` | contrascf-casf skill; counterfactual-benchmark-validation §1 | exact |
| perturbation strength = heavy-atom delta from manifest SMILES | `08_retention_table.py` `delta_heavy` | counterfactual-benchmark-validation §1 (stratify by strength) | exact |
| Wilson 95% CI | `08_retention_table.py` `wilson` | Wilson 1927 | exact |

## References
- Masters, M. R. et al. 2025. *Investigating whether deep learning models for co-folding learn the physics of protein–ligand interactions.* Nature Communications. (PDF in `docs/`.)
- Abramson, J. et al. 2024. *Accurate structure prediction of biomolecular interactions with AlphaFold 3.* Nature 630, 493–500. doi:10.1038/s41586-024-07487-w
- Passaro, S. et al. 2025. *Boltz-2.* bioRxiv doi:10.1101/2025.06.14.659707
- McNutt, A. T. et al. 2021. *GNINA 1.0: molecular docking with deep learning.* J. Cheminform. 13, 43. doi:10.1186/s13321-021-00522-2
- Cao, D. et al. 2025. *SurfDock.* Nature Methods 22, 310–322. doi:10.1038/s41592-024-02516-y
- Wilson, E. B. 1927. *Probable inference, the law of succession, and statistical inference.* JASA 22, 209–212. doi:10.1080/01621459.1927.10502953

## Next steps
- [confirmed] Do not quote any ligand-arm SurfDock number until it is re-run with the crystal ligand as the pocket reference (TODO item 15).
- [confirmed] Do not quote any ligand-arm co-folding RMSD, retention or confidence result produced before 2026-10-01 (atom-mapping bug).
- [tentative] The Δheavy ≥ 3 co-folding excess is the interesting open question; separating memorization from robustness needs a reference for the variant's true pose (e.g. docking consensus or physics-based refinement), not RMSD-to-crystal.
- [open] Re-read the June ligand-axis "severity-gated" claim against the corrected `s2a/s2c/confidence_conditional_ligand.csv`.
