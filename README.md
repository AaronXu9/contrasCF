# contrasCF

Secondary analysis of the adversarial co-folding dataset from
**Masters M.R., Mahmoud A.H., Lill M.A.** _Investigating whether deep learning
models for co-folding learn the physics of protein-ligand interactions._
**Nat. Commun.** 16:8854 (2025).

The paper tests four co-folding models (AlphaFold3, RoseTTAFold All-Atom,
Chai-1, Boltz-1) on adversarial protein-ligand challenges designed to break
binding from first principles. Finding: all four models largely keep the
ligand near the native pocket despite binding-destroying perturbations —
i.e., they memorise global sequence/structure patterns more than they learn
physics.

This repo contains:

1. **The 20 hand-built adversarial cases** (CDK2, GDH and MEK1 systems) and the
   analysis pipeline that scores them across four cofolding models plus
   physics-based docking baselines (UniDock2, GNINA, SurfDock, Boltz-2).
   Source under [`analysis/`](analysis/), entry scripts under
   [`analysis/paper_repro/scripts/`](analysis/paper_repro/scripts/).
2. **A new automated CASF-2016 mutagenesis pipeline** that scales the
   paper's binding-site mutagenesis (Fig. 3, n=285) by auto-detecting
   pockets and generating `wt`/`rem`/`pack`/`inv` variants for any
   PDBbind-format complex. Source under
   [`analysis/casf_mutagenesis/`](analysis/casf_mutagenesis/);
   end-to-end docs at [`docs/casf_mutagenesis.md`](docs/casf_mutagenesis.md).
3. **The ligand-side counterpart** — halogenation, methylation and charge-swap
   variants of every CASF ligand, scored across five methods. Source under
   [`analysis/ligand_mutagenesis/`](analysis/ligand_mutagenesis/); docs at
   [`docs/ligand_mutagenesis.md`](docs/ligand_mutagenesis.md).
4. **A τ-RAMD physics-validation pilot** that tests whether biased-MD exit
   times separate a paper binder from a non-binder. Source under
   [`analysis/ramd_pilot/`](analysis/ramd_pilot/); results in
   [`docs/ramd_pilot_results.md`](docs/ramd_pilot_results.md).

## How to read RMSD numbers from this project

The adversarial variants are **designed to break binding**. A physics-aware
model should place the ligand correctly on **WT** (low RMSD) but somewhere
different on `rem`/`pack`/`inv` (high RMSD), since those pockets no longer
accommodate the original ligand. **High ligand RMSD on adversarial cases is
the desired outcome**, the inverse of typical pose-prediction benchmarks.
"Memorisation rate" = fraction of adversarial cases with RMSD < threshold;
**lower is better**. Always report alongside the WT baseline.

## Project layout

Three benchmark **arms** over one shared **core**, plus a physics-validation track.
Each arm's scripts are grouped by stage — `build/` → `run/` → `analyze/` → `plot/`
(→ `export/`) — and numbered within a stage. Old script names resolve through each
arm's `scripts/RENAMES.md` (restructure of 2026-10, `plans/2026-10-01-organize-codes-results.md`).

```
analysis/
├── core/                      # shared by every arm — no arm imports another arm
│   ├── ligand_rmsd.py         # atom correspondence + in-place / best-fit ligand RMSD
│   ├── loaders.py             # structure + ligand loading (gemmi, RDKit)
│   ├── surfdock_engine.py     # SurfDock surface→CSV→ESM→diffusion helpers
│   └── reference_smiles.py    # Masters et al. 2025 reference SMILES
├── paper_repro/               # ARM 1 — the 20 hand-built paper cases (CDK2, GDH, MEK1)
│   ├── lib/                   #   config, align, ligand_match, confidence, clashes, pipeline
│   └── scripts/               #   01_fetch_native → 13_run_surfdock
├── casf_mutagenesis/          # ARM 2 — CASF-2016 pocket mutagenesis (wt/rem/pack/inv)
│   ├── *.py                   #   config, pocket, mutate, inputs_*, analysis, gnina_analysis
│   ├── scripts/{build,run,analyze,plot,export}/
│   ├── outputs/               #   per-cell predictions + result tables (gitignored)
│   └── figures/
├── ligand_mutagenesis/        # ARM 3 — CASF-2016 ligand mutagenesis (halo/meth/charge)
│   ├── rules/                 #   methylation, charge_swap, halogenation
│   ├── scripts/{build,run,analyze,plot}/
│   └── outputs/               #   incl. *_ligand.csv docking tables (gitignored)
├── ramd_pilot/                # PHYSICS VALIDATION — τ-RAMD exit-time oracle (CARC, GROMACS-RAMD)
│   └── {prep,ramd,analysis,scripts,tests,carc_setup}/
├── native/                    # crystal references (1B38, 2VWH, 7XLP) — shared input data
└── results/                   # paper-reproduction result tables + figures (gitignored)
contrasCF/                     # paper-shipped predictions (Zenodo 14749304) + input templates
docking/                       # docking working dirs (SurfDock reads inputs/<case>/box.json)
env/                           # lab.sh / carc.sh — source one before any script
slurm/                         # CARC job scripts
docs/                          # living docs; docs/archive/ holds retracted guidance
journal/  experiments/  plans/ # lab-notebook: postmortems, gated experiments, plans (immutable history)
```

## Setup

```bash
# Conda envs needed:
#   rdkit_env       — analysis (rdkit, gemmi, biopython, pandas, matplotlib)
#   PyMOL-PoseBench — headless PyMOL for figure rendering (optional)
#   unidock2        — UniDock2 + UniDock legacy (16-case docking)
#   boltzina_env    — Boltz-2 v2.2.1 binary
#   alphafold3 (dockstrat-managed) — AF3 v3.0.1 inference

# Source the host-specific env file before any script.
source env/lab.sh    # lab workstation
source env/carc.sh   # USC CARC (discovery.usc.edu)
```

All host-specific paths (project root, Boltz/AF3 binaries, dataset root,
GPU index) flow through `CONTRASCF_*` env vars resolved in
[`analysis/casf_mutagenesis/config.py`](analysis/casf_mutagenesis/config.py).
See the "Running on a different host" section of
[`docs/casf_mutagenesis.md`](docs/casf_mutagenesis.md) for the full
override table.

For the **16-case pipeline**, you also need the paper's predictions data
(~400 MB, gitignored) — download from
https://zenodo.org/records/14749304 and unpack into `contrasCF/data/`.

For the **CASF-2016 sweep**, you need PDBbind-cleansplit data — symlink the
local copy to `data/casf2016`:

```bash
ln -s /path/to/pdbbind_cleansplit data/casf2016
```

## Quick start (CASF mutagenesis sweep on subset20)

```bash
cd /path/to/contrasCF
source env/lab.sh    # or env/carc.sh — sets CONTRASCF_* and $CONTRASCF_PY

# 1. gating test (CDK2 must match paper 11/11)
$CONTRASCF_PY analysis/casf_mutagenesis/scripts/build/00_verify_reference_systems.py

# 2. build inputs for the 20-PDB subset
$CONTRASCF_PY analysis/casf_mutagenesis/scripts/build/01_build_subset20.py

# 3. run Boltz-2 on subset20  (~37 s/job × 76 jobs ≈ 47 min)
$CONTRASCF_PY analysis/casf_mutagenesis/scripts/run/01_run_boltz2.py

# 4. run AF3 with MSA on subset20 (~85 s/job × 76 jobs ≈ 1.8 h)
$CONTRASCF_PY analysis/casf_mutagenesis/scripts/run/04_run_af3_msa.py

# 5. analysis: ligand RMSD + memorisation rates (default scope is the FULL set)
CONTRASCF_SCOPE=subset20 $CONTRASCF_PY analysis/casf_mutagenesis/scripts/analyze/01_analyze_cofold.py
```

## Latest results (subset20, 2026-05-06)

| model | WT RMSD <2 Å | adversarial <2 Å (rem / pack / inv) |
|---|---|---|
| AF3 + ColabFold MSA | **0.79** (15/19) | 0.37 / 0.37 / 0.26 |
| Boltz-2 (single-seq) | 0.63 (12/19) | 0.26 / 0.37 / 0.21 |
| AF3 (no-MSA, broken baseline) | 0.00 | 0.00 / 0.00 / 0.00 |

Higher WT, lower adversarial = more physics-aware. Both production models
keep the ligand near the WT pose in ~30 % of adversarial cases despite
disrupted pockets — the residual memorisation Masters et al. 2025 quantifies.

See [`docs/casf_mutagenesis.md`](docs/casf_mutagenesis.md) for full
implementation history, gotchas, and the AF3+MSA fix that took six bugs to
land.
