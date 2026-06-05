# MEK1 binding-site mutagenesis experiment — design spec

**Status:** approved 2026-05-20 (user: Aaron Xu).
**Scope:** reproduce Masters et al. 2025 Fig. 2 (MEK1 / 7XLP binding-site mutagenesis) under the existing 16-case analysis pipeline so the resulting predictions are picked up by `rmsd_heatmap.png`, `grid_*.png`, `per_case/`, and the 6-model physics-vs-cofolding axis.

## Goal

The paper has three protein targets: CDK2 / 1B38, MEK1 / 7XLP, GDH / 2VWH. Our shipped dataset under `contrasCF/data/{model}/<case>/` covers only CDK2 (`bindingsite_*` + `atp_charge_*`) and GDH (`glucose_*`). MEK1 — paper Fig. 2 — was flagged as "not shipped" in [docs/project_notes.md](../../project_notes.md#L51-L52) and is a real gap. This spec fills it for the 6-model set used in the heatmap: AF3, Boltz-1, Boltz-2, UniDock2, GNINA, SurfDock. Chai-1 and RFAA are out of scope for now (the heatmap doesn't use them).

## What ships

Four new cases under the existing layout:

- `mek1_wt` — MEK1 (7XLP chain A) with FZC inhibitor, no mutations.
- `mek1_rem` — same protein, all 7 paper-listed pocket residues → Gly.
- `mek1_pack` — same protein, all 7 → Phe.
- `mek1_inv` — same protein, 7 residues swapped per Miyata table (`I105D, E108W, M110D, S158W, F173D, A40W, A59W` in paper numbering).

The ligand (FZC, 33 heavy atoms, an allosteric MEK1 inhibitor) is constant across all 4 variants.

## Pocket residues

Use the paper's **explicit 7-residue list**: A40, A59, I105, E108, M110, S158, F173 (paper numbering). In 7XLP these are at +36 offset → 76, 95, 141, 144, 146, 194, 209 (auth_seq), confirmed AA-identity-match.

**We bypass the auto-3.5 Å side-chain detector** for MEK1: the paper's list is documented (in `analysis/casf_mutagenesis/scripts/00_verify_reference_systems.py`) to be inconsistent with the strict 3.5 Å rule on this system (3 residues sit at 3.5–3.9 Å sc, only within 3.5 Å when backbone-inclusive). The CDK2 strict-rule gate still passes 11/11, so the rule itself is fine — we just defer to the paper's published list for MEK1 to reproduce Fig. 2 exactly.

## Ligand

- Resname: `FZC` (33 heavy atoms, contains aromatic heterocycles + Cl + chiral cyclohexylamine).
- SMILES source: CCD canonical = `Cc1ccnc(Oc2ccc(c(Cl)c2)c3cc4[nH]nc(C)c4c(O[C@H]5CCC[C@@H](N)C5)c3)n1`. As elsewhere in the project we round-trip via RDKit on the crystal SDF to dodge any kekulization edge cases.
- `COMMON_SUBSETS["FZC_FULL"] = <FZC SMILES>` — variants share an identical ligand, so the full SMILES is the right "common subset" for cross-case RMSD pairing.

## Layout and registry changes

Per-model predictions go in the existing layout:

```
contrasCF/data/AF3/mek1_{wt,rem,pack,inv}/        # *_model_0.cif + summary_confidences
contrasCF/data/Boltz/mek1_*/                       # _model_0.cif (Boltz-1)
contrasCF/data/Boltz2/mek1_*/                      # _model_0.cif (Boltz-2)
contrasCF/data/{UniDock2,GNINA,SurfDock}/mek1_*/   # combined_top1.pdb + poses.sdf
docking/inputs/mek1_*/                             # receptor.pdb + ligand.sdf + box.json
```

Two edits in `analysis/src/config.py`:

1. Add to `NATIVES`:
   ```python
   "mek1_fzc": NativeRef(pdb_id="7XLP", ligand_resname="FZC"),
   ```
2. Add 4 entries to `CASES` with `family="bindingsite_mek1"`, `native_key="mek1_fzc"`, `target_smiles=<FZC>`, `common_smarts_key="FZC_FULL"`.
3. Add `"FZC_FULL": <FZC SMILES>` to `COMMON_SUBSETS`.

Separate `bindingsite_mek1` family (not merged into CDK2's `bindingsite`) keeps the existing `grid_bindingsite.png` unchanged and produces a new `grid_bindingsite_mek1.png`.

## Execution pipeline (in dependency order)

1. **Config edit** — `NATIVES`, `CASES`, `COMMON_SUBSETS`.
2. **Cofolding inputs.** A new script `analysis/casf_mutagenesis/scripts/11_build_mek1.py` calls `build_system` with explicit `pocket_override=[...]`. Renders AF3 JSON (local "alphafold3" dialect, version 2) + Boltz YAML (`msa: empty`). Drops a copy of each input under `contrasCF/data/{AF3,Boltz,Boltz2}/mek1_*/` to slot into the existing analysis pipeline.
3. **MSA fetch.** Reuse `casf_mutagenesis.msa_via_boltz.fetch_msa_via_boltz` — one fetch for the WT MEK1 sequence, then `rewrite_a3m_query` to inject the mutated query line for `rem/pack/inv`. Inject into AF3 inputs via the `chain_msas` arg already supported by `inputs_af3.render_af3`.
4. **AF3+MSA run.** Reuse `casf_mutagenesis/scripts/06_run_af3_msa_subset20.py`-pattern runner, restricted to MEK1's 4 cells. Outputs land in `contrasCF/data/AF3/mek1_*/`.
5. **Boltz-1 + Boltz-2 runs.** Two passes (`--model boltz1` and `--model boltz2`) using `boltzina_env` (v2.2.1 — handles both versions via the `--model` flag).
6. **Docking input prep.** Mirror `analysis/scripts/10_prep_docking_inputs.py` for MEK1: strip AF3 prediction → `receptor.pdb`; copy crystal FZC SDF → `ligand.sdf`; map crystal-ligand centroid into AF3 receptor frame via Cα superpose → `box.json` (25³ Å).
7. **Docking runs.** UniDock2 + GNINA via `11_run_docking.py`; SurfDock via `13_run_surfdock.py` with the `--ligand_to_pocket_center` flag (proven necessary for SurfDock on adversarial ligands per the existing notes).
8. **Analysis.** `02_run_analysis.py` auto-picks up the new cases (since `CASES` was updated) → 24 new rows (4 cases × 6 models) appended to `results.csv`.
9. **Figures.** `03_make_plots.py` + `04_render_figures.py` + `05_physics_analysis.py` regenerate. Heatmap gains 4 columns; new `grid_bindingsite_mek1.png` rendered with the 6 model rows × 4 case columns.
10. **Docs.** Fix the "two proteins" framing in `docs/project_notes.md` and in memory. New short `docs/mek1_experiment.md` with the per-system numbers.

## Risks / non-obvious failure modes (and mitigations)

- **AF3 single-sequence Cα RMSD ~17 Å on novel targets** (proven on CDK2). Mitigation: always run AF3 with MSA via Boltz piggyback. Without this the docking receptor is garbage and the entire physics axis is meaningless.
- **FZC kekulization edge case.** Mitigation: derive SMILES from the crystal SDF via `Chem.MolToSmiles(mol_from_sdf)`, same pattern that fixed 4ih5/4de1/4ivc in the CASF runs.
- **SurfDock diffusion blow-up on novel ligands.** Mitigation: pass `--ligand_to_pocket_center` + pre-shift ligand centroid to pocket centre — the documented workaround.
- **Cα superpose on novel resnum range.** 7XLP starts at resnum 37 (not 1). The original 16-case analysis pipeline scans ±50 offset; this is within range so should work. If it silently fails (zero paired Cα), fall back to `casf_mutagenesis.analysis.superpose_by_index`.
- **`bindingsite_pack` data quirk** (no Mg²⁺ in input) doesn't apply here — MEK1's FZC has no Mg cofactor.

## What's deferred / explicitly out of scope

- Chai-1 and RFAA inputs.
- Funnel-metadynamics simulations (paper Table 1; out of project scope generally).
- Per-residue / per-atom analysis of the FZC pose mapping (would need a custom SMARTS subset analogous to ATP_FULL — using FZC_FULL is fine since all 4 variants share the same ligand graph).

## Acceptance criteria

The work is done when:

1. `results.csv` has 24 new rows for `mek1_{wt,rem,pack,inv}` × 6 models with non-NaN `lig_rmsd_core` values.
2. `rmsd_heatmap.png` shows the 4 new columns alongside the existing 16.
3. `grid_bindingsite_mek1.png` exists with the same layout as `grid_bindingsite.png`.
4. WT median ligand RMSD ≤ 3 Å for AF3+MSA and Boltz-1/2 (sanity gate — the unperturbed system should reproduce reasonably).
5. `docs/project_notes.md` no longer claims "two proteins"; `docs/mek1_experiment.md` summarises the new run.
