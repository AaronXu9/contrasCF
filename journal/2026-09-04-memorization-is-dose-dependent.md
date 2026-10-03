# Memorization Is Dose Dependent

**Type:** milestone
**Date:** 2026-09-04
**Status:** done
**Project:** contrasCF

## Headline
Co-folding memorization is not a fixed property of the model — it is a **monotonic function of how many pocket residues were mutated**, and the CASF pocket is far smaller than the write-ups imply (**median 5 residues**, range 0–13). AF3+MSA `rem` retention runs **78.3 % at 1–2 mutations down to 14.7 % at 9+**; at 12 mutations the co-folding models fail as hard as docking.

## Why it matters

Every headline rate in `docs/casf_overview.md` ("Boltz-2 memorizes 36.8 % on `rem`") is a weighted average over a perturbation-strength distribution that happens to be centred at 5 mutations. It is a property of **model × perturbation strength**, not of the model alone. Quoted bare, it invites the reading "this model ignores the pocket", which the dose curve refutes: the models *do* respond, monotonically.

WT-conditioned retention (< 2 Å) vs number of mutated pocket residues:

| n_mut | systems (Boltz-2) | AF3+MSA rem | Boltz-2 rem | SurfDock rem |
|---|---|---|---|---|
| 1–2 | 18 | **78.3 %** | 50.0 % | 4.0 % |
| 3–4 | 40 | 59.2 % | 50.0 % | 5.6 % |
| 5–6 | 31 | 40.9 % | 29.0 % | 0.0 % |
| 7–8 | 26 | 33.3 % | 30.8 % | 0.0 % |
| 9–20 | 20 | **14.7 %** | 15.0 % | 3.0 % |
| point-biserial r | | **−0.412** | −0.278 | −0.102 |

Docking is flat because it is already at the floor. The gap between families is real at *matched* perturbation strength (14.7 % vs 3.0 % at 9+), but it is roughly 5×, not the ~15× the pooled numbers suggest.

Three independent physics checks then reframe *what* the retained poses are (all on the 251-system set):

- **PoseBusters** (`outputs/cofold_validity.csv`, 916 cells): mutant complexes pass at 62–69 % vs 80.8 % for WT, the drop driven by `minimum_distance_to_protein`. But split by outcome, the **retained** poses are the clean ones — `rem` retained 88.0 % pass / 12.0 % min-dist failure, statistically identical to WT's 87.5 % / 12.5 %, while *moved* poses sit at 66.3 % / 32.6 %.
- **Rotamers** (`outputs/cofold_rotamers.csv`, 3 838 introduced side chains): χ1 outlier rate **1.0 %** (pack) and 0.9 % (inv) vs **1.1 %** in the crystal baseline (6 726 residues). The models never strain a side chain to make room.
- **GNINA CNN rescoring** (`outputs/crystal_pose_rescore.csv`): the model-free control — crystal receptor with pocket side chains deleted — costs ΔVina **+2.24**, ΔCNNscore **−0.389**. The models' own `rem` complexes carry ΔVina +1.60 (Boltz-2) / +1.96 (AF3+MSA), i.e. **60–88 % of the true penalty**, not zero.

So the mechanism is **constrained retrieval, not fabrication**: the model recalls the complex, the mutation acts as a feasibility constraint, and when the remembered pose stays accommodatable it is emitted inside a structure that is physically clean by every independent measure. The structure module registers that the site degraded; the ligand-placement decision does not act on that information. Same shape as [[2026-06-30-three-head-dissociation]].

## Model & data provenance

| id | size / N | source | version |
|---|---|---|---|
| Boltz-2 | 229 systems × 4 variants, rank-0 | `boltzina_env` boltz 2.2.1, `--model boltz2` | `outputs/results_full.csv` |
| AF3+MSA | 239 systems × 4 variants, rank-0 | AF3 v3.0.1, ColabFold MSA via `msa_via_boltz` | `outputs/results_full.csv` |
| SurfDock | 249 wt / 238 mutant cells | `/home/aoxu/projects/SurfDock`, post interface-crop fix | `outputs/docking_results.csv` |
| ICM | 251 wt / 232–233 mutant cells | `forks/ICM/`, poses pre-aligned to crystal frame | `outputs/icm_results.csv` |
| CASF-2016 core | 285 ids; 251 with a manifest entry | `pdbbind_cleansplit`, HiQBind + 30 RCSB-recovered | `docs/data_store_map.md` |
| GNINA | v1.3 master:97fa6bc, built 2024-10-03 | `/home/aoxu/projects/PoseBench/forks/GNINA/gnina` | `--score_only --cnn_scoring rescore --seed 42` |
| PoseBusters | 916 cells, `dock` config | `PoseBench` env, posebusters 0.2.13 | `outputs/cofold_validity.csv` |

## References

- Masters, M. R. et al. (2025). Do deep-learning co-folding methods generalise? — the pocket-mutagenesis protocol (`rem`/`pack`/`inv`, Miyata-maximal inversion table) this module reimplements. Repo PDF `contrasCF/`; Zenodo 14749304.
- Buttenschoen, M., Morris, G. M., Deane, C. M. (2024). PoseBusters: AI-based docking methods fail to generate physically valid poses. *Chem. Sci.* 15, 3130. DOI 10.1039/D3SC04185A.
- McNutt, A. T. et al. (2021). GNINA 1.0: molecular docking with deep learning. *J. Cheminform.* 13, 43. DOI 10.1186/s13321-021-00522-2.
- Lovell, S. C. et al. (2000). The penultimate rotamer library. *Proteins* 40, 389. — source of the canonical χ1 wells (−60/60/180) and the ~1 % outlier expectation.
- Miyata, T., Miyazawa, S., Yasunaga, T. (1979). Two types of amino acid substitutions in protein evolution. *J. Mol. Evol.* 12, 219. — the distance used by the `inv` table.

## Code ↔ theory alignment

| claim | code | justification |
|---|---|---|
| pocket = side-chain heavy atom within 3.5 Å of ligand; this sets n_mut | `analysis/casf_mutagenesis/config.py:POCKET_CUTOFF_A` | Masters et al. protocol; the 3.5 Å cutoff is why the median pocket is only 5 residues |
| dose curve binned on the manifest's own `n_mutations` | `outputs/manifest_full_casf.json` `variants.<v>.n_mutations`, written by `scripts/02_build_full_casf.py` | counts what was actually applied, not what was requested — catches the `3mss`/`4eo8` no-ops |
| in-place symmetry-corrected heavy-atom RMSD, never ligand superposition | `analysis/casf_mutagenesis/gnina_analysis.py:_aligned_rmsd` | superposing would measure conformer similarity, not placement |
| PoseBusters run with cofactor/water checks excluded | `scripts/34_validate_cofold_pockets.py:CORE_CHECKS` | **divergence from stock PoseBusters**: co-folding emits no waters/cofactors, so those checks fail 100 % for every variant incl. wt and would swamp the comparison |
| χ1 outlier = >40° from nearest canonical well | `scripts/34_validate_cofold_pockets.py:chi1_offset,CHI1_TOL` | Lovell et al. rotamer wells; validated against a 6 726-residue crystal baseline that reproduces the expected ~1 % |
| model-free `rem` receptor by backbone truncation | `scripts/35_rescore_crystal_pose.py:crystal_rem_receptor` | `rem` is pure deletion, so truncating crystal residues reproduces it exactly with no model and no superposition anywhere in the path |

## Open threads
- [confirmed] Every rate quoted in `docs/casf_overview.md`, `docs/casf_confidence.md` and the Feishu docs must carry the dose curve or a mutation-count stratification; the bare pooled number is misleading on its own.
- [confirmed] Figure selection must filter on mutation count — see [[2026-09-04-perturbation-strength-blind-spot]] for how the unfiltered version went wrong.
- [open] The `crystal_in_docking_receptor` rescoring arm is contaminated by Cα-fit error (a 1.3 Å fit produced a spurious +111 kcal/mol on `3zso`) and must not be quoted; only the `*_self` and `crystal_truncated_in_place` arms are alignment-free.
- [open] `pack`/`inv` have no model-free control, because reproducing them on the crystal needs side-chain packing. Without one, their penalties rest on model-built receptors.
- [tentative] Retention may partly be defensible physics rather than memory in the low-n_mut regime: at 1–2 mutations the site genuinely survives, and a model that kept the pose would be right. The argument should lean on the ≥7-mutation stratum.
- [suggested] Re-run the dose curve per seed; all co-folding numbers here are single-seed (seed 42).
