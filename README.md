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

1. **The 16 hand-built adversarial cases** (CDK2 + GDH systems) and the
   analysis pipeline that scores them across four cofolding models plus
   physics-based docking baselines (UniDock2, GNINA, SurfDock, Boltz-2).
   Source under [`analysis/`](analysis/), entry scripts under
   [`analysis/scripts/`](analysis/scripts/).
2. **A new automated CASF-2016 mutagenesis pipeline** that scales the
   paper's binding-site mutagenesis (Fig. 3, n=285) by auto-detecting
   pockets and generating `wt`/`rem`/`pack`/`inv` variants for any
   PDBbind-format complex. Source under
   [`analysis/casf_mutagenesis/`](analysis/casf_mutagenesis/);
   end-to-end docs at [`docs/casf_mutagenesis.md`](docs/casf_mutagenesis.md).

## How to read RMSD numbers from this project

The adversarial variants are **designed to break binding**. A physics-aware
model should place the ligand correctly on **WT** (low RMSD) but somewhere
different on `rem`/`pack`/`inv` (high RMSD), since those pockets no longer
accommodate the original ligand. **High ligand RMSD on adversarial cases is
the desired outcome**, the inverse of typical pose-prediction benchmarks.
"Memorisation rate" = fraction of adversarial cases with RMSD < threshold;
**lower is better**. Always report alongside the WT baseline.

## Project layout

```
analysis/
├── src/                       # 16-case pipeline (config, loaders, native,
│                              #   align, ligand_match, confidence, clashes,
│                              #   pipeline)
├── scripts/                   # 16-case end-to-end (01_fetch → 13_run_surfdock)
├── native/                    # crystal references (1B38, 2VWH, 7XLP)
├── casf_mutagenesis/          # NEW: automated CASF-2016 sweep
│   ├── config.py              # mutation tables (rem→G, pack→F, inv per Miyata)
│   ├── sequence.py / pocket.py
│   ├── mutate.py
│   ├── inputs_{af3,boltz,docking}.py
│   ├── msa_via_boltz.py       # ColabFold piggyback (direct API errors out)
│   ├── analysis.py            # ligand RMSD, memorisation aggregates
│   └── scripts/               # 00_verify → 06_run_af3_msa
└── results/                   # 16-case results (CSV + figures); gitignored
contrasCF/
├── data/                      # paper-shipped predictions (gitignored;
│                              #   download from Zenodo 14749304)
└── Cofolding-Tools-main/      # AF3/Boltz/Chai/RFAA input templates
docking/                       # working dirs for UniDock2 / GNINA / SurfDock
docs/
├── project_notes.md           # 16-case pipeline notes + caveats
└── casf_mutagenesis.md        # CASF sweep details, results, gotchas
```

## Setup

```bash
# Conda envs (already configured on the lab workstation; reproduce on CARC):
#   rdkit_env       — analysis (rdkit, gemmi, biopython, pandas, matplotlib)
#   PyMOL-PoseBench — headless PyMOL for figure rendering
#   unidock2        — UniDock2 + UniDock legacy
#   boltzina_env    — Boltz-2 v2.2.1 binary
#   alphafold3 (CogLigandBench) — AF3 v3.0.1 inference

# Mandatory before any script: prepend env's lib/ to LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/path/to/env/lib:$LD_LIBRARY_PATH
```

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
LIB=/path/to/rdkit_env/lib
PY=/path/to/rdkit_env/bin/python

# 1. gating test (CDK2 must match paper 11/11)
LD_LIBRARY_PATH=$LIB:$LD_LIBRARY_PATH $PY \
    analysis/casf_mutagenesis/scripts/00_verify_reference_systems.py

# 2. build inputs for the 20-PDB subset
LD_LIBRARY_PATH=$LIB:$LD_LIBRARY_PATH $PY \
    analysis/casf_mutagenesis/scripts/01_build_subset20.py

# 3. run Boltz-2 on subset20  (~37 s/job × 76 jobs ≈ 47 min)
LD_LIBRARY_PATH=$LIB:$LD_LIBRARY_PATH $PY \
    analysis/casf_mutagenesis/scripts/03_run_boltz2_subset20.py

# 4. run AF3 with MSA on subset20 (~85 s/job × 76 jobs ≈ 1.8 h)
LD_LIBRARY_PATH=$LIB:$LD_LIBRARY_PATH $PY \
    analysis/casf_mutagenesis/scripts/06_run_af3_msa_subset20.py

# 5. analysis: ligand RMSD + memorisation rates
LD_LIBRARY_PATH=$LIB:$LD_LIBRARY_PATH $PY \
    analysis/casf_mutagenesis/scripts/05_analyze_subset20.py
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
