# contrasCF — project notes

## What this project is

Secondary analysis of the adversarial co-folding dataset from:

> Masters M.R., Mahmoud A.H., Lill M.A. "Investigating whether deep learning models
> for co-folding learn the physics of protein-ligand interactions."
> *Nat. Commun.* 16:8854 (2025). [PDF](Masters%20et%20al.%20-%202025%20-%20Investigating%20whether%20deep%20learning%20models%20for%20co-folding%20learn%20the%20physics%20of%20protein-ligand%20intera.pdf).

The paper tests four co-folding models — **AlphaFold3**, **RoseTTAFold All-Atom**,
**Chai-1**, **Boltz-1** — on adversarial challenges designed from first principles
to disrupt protein-ligand binding. Finding: all four models largely keep the
ligand near the native pocket despite binding-destroying perturbations, i.e.
they memorize global sequence/structure patterns more than they learn physics.

Data and tool-input templates are shipped under [contrasCF/](../contrasCF/):
- [contrasCF/data/{AF3,Boltz,Chai,RFAA}/](../contrasCF/data/) — predicted structures,
  16 test cases × 4 models = 64 predictions.
- [contrasCF/Cofolding-Tools-main/](../contrasCF/Cofolding-Tools-main/) — per-tool input
  templates (AF3 JSON, Boltz YAML, Chai FASTA, RFAA YAML+SDF) and versions used.
- Source data on Zenodo: [14749304](https://zenodo.org/records/14749304); code on
  [14768184](https://zenodo.org/records/14768184).

## The 16 test cases

Two proteins, three challenge families.

### CDK2 kinase (PDB 1B38) — ATP binding pocket

**Binding-site mutations** — protein changed, ligand = native CCD_ATP.
11 pocket residues targeted: I10, T14, V18, A31, K33, D86, K129, Q131, N132, L134, D145.

| case | mutation | Mg²⁺ in input? |
|---|---|---|
| `bindingsite_wt` | none | yes |
| `bindingsite_rem` | all 11 → Gly | yes |
| `bindingsite_pack` | all 11 → Phe | **no** (input quirk) |
| `bindingsite_inv` | 10 positions → Miyata-dissimilar (I10D, T14W, …) | yes |

**ATP charge modifications** — CDK2 unchanged, ligand replaced via raw SMILES.
Triphosphate substituted by:
- `atp_charge_methyl/ethyl/propyl`: neutral alkyl chains, formal charge 0
- `atp_charge_1/2/3`: quaternary amines, formal charge +1/+2/+3

### Glucose dehydrogenase (PDB 2VWH) — glucose binding site

**Ligand methylation** — GDH unchanged; glucose alcohols methylated 0→5 times.
Input also carries NADP cofactor and Zn²⁺.

### Not shipped
- MEK1 (PDB 7XLP) — described in paper Fig. 2 but no predicted structures in the data.

## Per-model output schema (empirically verified)

| model | canonical file | confidence | seeds | ligand chain |
|---|---|---|---|---|
| AF3 | `*_model_0.cif` (or `*_model.cif` for atp_charge/glucose; `*_model_0_N.cif` for some bindingsite_*) | `*_summary_confidences*.json` + `*_confidences.json` (per-atom pLDDT) | 5 in `seed-1_sample-*/` | B (CDK2), D (GDH glucose — NADP on B) |
| Boltz | `*_model_0.cif` | `confidence_*_model_0.json` + `plddt_*.npz` | none | B (CDK2), C (GDH) |
| Chai | `pred.rank_0.cif` | `scores.rank_0.json` + `pae.rank_*.npy` | 5 ranks | B (CDK2), D (GDH) |
| RFAA | `<case>.pdb` (PDB, not CIF) | B-factor pLDDT only (`.pt` checkpoint otherwise) | none | B (CDK2), C (GDH) |

**Key gotcha**: the target-ligand chain ID is not "B" universally — for glucose cases,
AF3 places NADP on B and glucose on D; Boltz/RFAA on C. The analysis pipeline
(see below) identifies the target ligand by SMILES/heavy-atom match, not chain.

## Analysis pipeline

Everything lives under [../analysis/](../analysis/). Two conda envs used:

- `/home/aoxu/miniconda3/envs/rdkit_env` — analysis (rdkit, gemmi, biopython, pandas, matplotlib, seaborn)
- `/home/aoxu/miniconda3/envs/PyMOL-PoseBench` — headless PyMOL for figure rendering

Both require prepending env's `lib/` to `LD_LIBRARY_PATH` to avoid system libstdc++.

### Structure
- [`analysis/src/config.py`](../analysis/src/config.py) — CASES, MODELS, COMMON_SUBSETS (single source of truth; all 16 case SMILES hard-coded from the AF3 `*_data.json` inputs).
- [`analysis/src/loaders.py`](../analysis/src/loaders.py) — gemmi CIF/PDB loading; `select_target_ligand()` picks the target by MCS-vs-SMILES score, ignoring metals, NADP, waters, and ACE caps.
- [`analysis/src/native.py`](../analysis/src/native.py) — caches fetched 1B38/2VWH references (native ATP resname `ATP`; native glucose resname **`BGC`** = β-D-glucose, not GLC).
- [`analysis/src/align.py`](../analysis/src/align.py) — Bio.SVDSuperimposer Cα superposition over residue-number intersection.
- [`analysis/src/ligand_match.py`](../analysis/src/ligand_match.py) — RDKit MCS pairing with **symmetry-minimizing selection** (enumerates `GetSubstructMatches` and picks the pairing with minimum pre-fitted RMSD — matters for phosphate-O permutations).
- [`analysis/src/confidence.py`](../analysis/src/confidence.py) — per-model pTM/ipTM/ligand-iPTM/pLDDT extractors; RFAA falls back to ligand B-factor mean.
- [`analysis/src/clashes.py`](../analysis/src/clashes.py) — heavy-atom protein-ligand pairs at <2.0 Å.
- [`analysis/src/pipeline.py`](../analysis/src/pipeline.py) — `run_one(model, case) -> row dict`.

### Scripts
- [`analysis/scripts/01_fetch_native.py`](../analysis/scripts/01_fetch_native.py) — downloads 1B38.cif, 2VWH.cif from RCSB into `analysis/native/`.
- [`analysis/scripts/02_run_analysis.py`](../analysis/scripts/02_run_analysis.py) — iterates 4 × 16 = 64 cells, writes `analysis/results/results.csv` (24 columns).
- [`analysis/scripts/03_make_plots.py`](../analysis/scripts/03_make_plots.py) — aggregate heatmap, Fig.3-style bar chart, confidence vs RMSD scatter.
- [`analysis/scripts/04_render_figures.py`](../analysis/scripts/04_render_figures.py) — headless PyMOL renders per-case binding-site PNGs and a grid per family (matches paper Fig. 1/4/5 visual style).

### Run commands
```bash
# analysis (rdkit_env)
export LD_LIBRARY_PATH="/home/aoxu/miniconda3/envs/rdkit_env/lib:$LD_LIBRARY_PATH"
/home/aoxu/miniconda3/envs/rdkit_env/bin/python analysis/scripts/01_fetch_native.py
/home/aoxu/miniconda3/envs/rdkit_env/bin/python analysis/scripts/02_run_analysis.py
/home/aoxu/miniconda3/envs/rdkit_env/bin/python analysis/scripts/03_make_plots.py

# rendering (PyMOL-PoseBench)
export LD_LIBRARY_PATH="/home/aoxu/miniconda3/envs/PyMOL-PoseBench/lib:$LD_LIBRARY_PATH"
/home/aoxu/miniconda3/envs/PyMOL-PoseBench/bin/python analysis/scripts/04_render_figures.py
```

## Results summary

All 64 cells successfully processed (no failures).

Median ligand core-RMSD to native (Å):

| family | AF3 | RFAA | Boltz | Chai |
|---|---|---|---|---|
| atp_charge | 0.89 | 1.00 | 0.75 | 0.98 |
| bindingsite | 4.21 | 3.84 | 6.02 | 2.76 |
| glucose | **13.80** | 2.54 | 2.08 | **24.26** |

- `atp_charge`: all four models stay ~1 Å from native despite drastic charge/chemistry changes → no electrostatics learned (paper's main Fig. 5 finding).
- `bindingsite`: models resist adapting to mutated pockets (Fig. 1).
- `glucose`: AF3 and Chai eject heavily methylated glucose from the pocket (>20 Å); Boltz and RFAA keep it nearby (Fig. 4).

### Artifacts
- [`analysis/results/results.csv`](../analysis/results/results.csv) — 64 rows, 24 columns (long-form).
- [`analysis/results/figures/rmsd_heatmap.png`](../analysis/results/figures/rmsd_heatmap.png)
- [`analysis/results/figures/fig3_bars.png`](../analysis/results/figures/fig3_bars.png)
- [`analysis/results/figures/conf_vs_rmsd.png`](../analysis/results/figures/conf_vs_rmsd.png)
- [`analysis/results/figures/grid_{bindingsite,atp_charge,glucose}.png`](../analysis/results/figures/) — paper-style 4-model × N-case grids.
- [`analysis/results/figures/per_case/{family}/{model}_{case}.png`](../analysis/results/figures/per_case/) — 64 individual PyMOL renders.

## Extensions

### CASF-2016 binding-site mutagenesis sweep

Beyond the 16 hand-built cases, the paper also runs a **systematic CASF-2016
sweep** (Fig. 3, n=285) — automatically detecting pockets by 3.5 Å distance
and generating `wt`/`rem`/`pack`/`inv` variants for every PDB. The original
repo never reproduced this. The
[`analysis/casf_mutagenesis/`](../analysis/casf_mutagenesis/) module
implements it; full details, gotchas, and verification results in
[`docs/casf_mutagenesis.md`](casf_mutagenesis.md). CDK2/1B38 reference
verification matches the paper 11/11 residues.

## Known caveats
- WT ligand RMSD comes out ~0.7 Å higher than the paper's quoted values
  (AF3 WT 0.94 Å vs paper 0.2 Å). Likely because we Cα-superpose the full 290-residue
  intersection; the paper almost certainly uses pocket-only Cα superposition (the standard docking metric).
  Relative ordering across cases is preserved, so trends and rankings match the paper.
- MEK1 (7XLP) cases from paper Fig. 2 are not in the shipped data.
- Funnel-metadynamics (paper Table 1 bound/unbound probabilities) requires OpenMM simulations — out of scope here.
- `bindingsite_pack` uniquely lacks Mg²⁺ in the AF3 input (flagged via `config.CASES["bindingsite_pack"].has_mg_input = False`); other bindingsite cases have it.
- Native glucose in 2VWH is CCD **`BGC`** (β-D-glucose), not GLC.
