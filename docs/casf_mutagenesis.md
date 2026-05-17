# CASF-2016 binding-site mutagenesis pipeline

Automated implementation of Masters et al. 2025's binding-site mutagenesis
challenge (Fig. 3, n=285) on the CASF-2016 protein–ligand complex set. Lives
under [`analysis/casf_mutagenesis/`](../analysis/casf_mutagenesis/).

## What the paper does, and what the existing repo had

The paper systematically applies three challenges — **removal**, **packing**,
**inversion** — to every binding pocket in CASF-2016 and tests whether the
co-folding model still places the ligand near the native pose. The shipped
repo contains only the four hand-built CDK2 cases (`bindingsite_{wt,rem,pack,
inv}`) — the systematic CASF-2016 sweep was never reproduced. This module
fills that gap.

## Authoritative rules (Masters et al. 2025, Methods §"Binding site mutations")

**Pocket residue identification:**
> "Residues with side-chain heavy atoms within 3.5 Å of any ligand heavy
> atoms were selected for mutation."

**Mutation tables:**

| variant | rule |
|---|---|
| `wt` | original sequence (control) |
| `rem` | all pocket residues → **Gly** |
| `pack` | all pocket residues → **Phe** |
| `inv` | per-residue Miyata-distance-maximising substitution (Table 2, paper) |

The full Miyata-far table:

| WT 3-letter | mutation |
|---|---|
| TYR, TRP, LYS, ARG | → G |
| VAL, LEU, ILE, MET, PHE | → D |
| ALA, GLY, PRO, HIS, ASP, GLU, CYS, ASN, GLN, THR, SER | → W |

Glycine in WT is excluded from pocket detection (no side chain to clash with
the ligand) and so never appears in `pocket` lists. PHE in WT becomes a no-op
under `pack` (identity substitution); the code silently skips identity
substitutions, so per-system `pack` mutation count = `pocket_size − n_PHE`.

## Data sources

| input | path |
|---|---|
| crystal protein PDB | `data/casf2016/raw/<pdbid>/<pdbid>_protein.pdb` |
| crystal ligand SDF | `data/casf2016/crystal_ligands/<pdbid>_ligand.sdf` |
| HET-code SMILES cache | `data/casf2016/smiles/het_smiles_cache.json` |
| PDBbind metadata (HET → ligand_name) | `data/casf2016/labels/PDBbind_data_dict.json` |
| 20-complex smoke test | `data/casf2016/labels/PDBbind_casf2016_subset20.json` |
| full CASF-2016 (285 PDBs) | `data/casf2016/labels/PDBbind_data_split_cleansplit.json["casf2016"]` |

`data/casf2016` is a symlink to `/home/aoxu/projects/VLS-Benchmark-Dataset/
data/pdbbind_cleansplit` (created in setup).

## Architecture

```
analysis/casf_mutagenesis/
├── config.py             # paths, mutation tables (rem/pack/inv), 3.5 Å cutoff
├── sequence.py           # PDB → per-chain sequence with PDB resnum mapping
├── pocket.py             # 3.5 Å side-chain heavy-atom pocket detection
├── mutate.py             # apply (chain, resnum, new_aa) mutations
├── inputs_af3.py         # AF3 JSON (local "alphafold3" dialect, version 2)
├── inputs_boltz.py       # Boltz YAML (msa: empty for single-sequence mode)
├── inputs_docking.py     # WT-only receptor.pdb + ligand.sdf + box.json
├── build.py              # per-system pipeline orchestrator
├── analysis.py           # per-prediction RMSD + memorization aggregates
├── scripts/
│   ├── 00_verify_reference_systems.py   # gating test (CDK2 + MEK1)
│   ├── 01_build_subset20.py             # build inputs for 20 systems
│   ├── 02_build_full_casf.py            # build inputs for 285 systems (TODO)
│   ├── 03_run_boltz2_subset20.py        # GPU runner — Boltz-2
│   ├── 04_run_af3_subset20.py           # GPU runner — AF3 (no-MSA)
│   └── 05_analyze_subset20.py           # ligand-RMSD + memorization
└── outputs/
    └── <pdbid>/<variant>/{af3.json, boltz.yaml,
                           <prefix>_model_0.cif (Boltz),
                           af3_<prefix>_model_0.cif,
                           confidence_*.json,
                           docking/{receptor.pdb,ligand.sdf,box.json}}
```

## Per-system pipeline (`build.py::build_system`)

For each PDB id:

1. **Resolve files** under `data/casf2016/`.
2. **Extract per-chain sequence** with PDB residue-number mapping
   (`sequence.extract_chain_sequences`, gemmi-based).
3. **Identify pocket residues** (`pocket.detect_pocket`):
   - Load protein (gemmi) and ligand heavy-atom coords (RDKit, from SDF).
   - For each protein residue, take side-chain heavy atoms (skip backbone
     N/CA/C/O/OXT/HA and any `H`).
   - Mark the residue if any side-chain atom is ≤ 3.5 Å from any ligand
     heavy atom.
   - Skip non-standard residues and glycine (no side chain).
4. **Generate the four sequence variants** (`mutate.apply_mutations`):
   `wt` is identity, `rem` and `pack` are blanket substitutions, `inv` looks
   up `INV_TABLE` per-residue. Identity substitutions are silently skipped.
5. **Pick a SMILES.** Prefer crystal-SDF-derived canonical SMILES
   (`Chem.MolToSmiles(mol_from_sdf)`) over the HET cache. Reason: the cached
   SMILES from `het_smiles_cache.json` uses lowercase aromatic notation that
   fails RDKit kekulization for some fused N-heterocycles
   (e.g. 4ih5 pyrazolopyrimidine, 4de1 tetrazole-phthalazine,
   4ivc cyano-fused-heterocycle); the SDF is always parseable. Both
   Boltz-2 and AF3 use RDKit internally so they reject the broken SMILES.
6. **Render per-model inputs:**
   - **AF3:** `dialect: alphafold3`, `version: 2`, with no-MSA fields
     (`unpairedMsa: ">query\n{seq}\n"`, `pairedMsa: ""`, `templates: []`).
   - **Boltz:** YAML with `msa: empty` per protein chain (single-sequence).
   - **Docking (WT only):** `receptor.pdb` (all standard-AA chains preserved
     — differs from the existing 16-case flow that drops everything except
     the largest chain), `ligand.sdf` (copied from crystal SDF, fallback to
     RDKit embedding), `box.json` (centered on crystal-ligand centroid,
     25 × 25 × 25 Å).
7. **Write manifest entry** capturing pocket residues, applied mutations,
   SMILES + source, output paths, warnings.

Mutant-protein docking inputs are deferred — they need an AF3 prediction of
the mutated sequence first, then a follow-on script analogous to
`analysis/scripts/10_prep_docking_inputs.py`.

## Per-model file-naming convention

Both Boltz-2 and AF3 land in the same `outputs/<pdbid>/<variant>/` directory.
To prevent collision, AF3 outputs are prefixed with `af3_`:

| model | structure cif | confidence |
|---|---|---|
| Boltz-2 | `<pdbid>_<variant>_model_0.cif` | `confidence_<pdbid>_<variant>_model_0.json` |
| AF3 | `af3_<pdbid>_<variant>_model_0.cif` | `af3_summary_confidences_<pdbid>_<variant>.json` |

## Implementation gotchas (and workarounds)

| issue | resolution |
|---|---|
| Boltz uses the input YAML's filename stem as the output prefix | runner copies `boltz.yaml` → `<prefix>.yaml` inside a per-job work dir before invoking |
| Cached HET SMILES fails RDKit `MolFromSmiles` for 4ih5 / 4de1 / 4ivc | derive SMILES from crystal SDF via `Chem.MolToSmiles` |
| RDKit warning "molecule is tagged as 2D, but at least one Z coordinate is not zero" | benign — emitted by every PDBbind crystal SDF read; ignore |
| AF3 v3 in single-sequence mode produces structurally bad proteins (~17 Å Cα RMSD) | run with full MSA pipeline OR inject externally-generated MSAs into the JSON; see "Known limitation" below |
| `version: 1` AF3 JSON without `unpairedMsa`/`pairedMsa`/`templates` is rejected | always emit `version: 2` with the no-MSA fields |
| 7XLP's PDB auth_seq numbering is +36 vs the paper's MEK1 numbering | verification reports MEK1 as informational only; CDK2 is the strict gate |
| Crystal PDBs use UniProt-canonical residue numbering (e.g. 4jia chain A is 833-1132); predictions number 1..N | `analysis.superpose_by_index` pairs Cα by chain-position index, not PDB resnum. The repo's `analysis.src.align.superpose_ca` is unsuitable for CASF (only scans ±50 offset) — use the index-based variant for any new PDB-numbering-agnostic system |

## Verification (gating test)

[`scripts/00_verify_reference_systems.py`](../analysis/casf_mutagenesis/scripts/00_verify_reference_systems.py)
runs the auto-detection logic on the two systems Masters et al. enumerate
explicitly:

### CDK2 / 1B38 (strict gate)

Paper-stated 11-residue pocket: `{I10, T14, V18, A31, K33, D86, K129, Q131,
N132, L134, D145}`. Our 3.5 Å side-chain rule on the 1B38 crystal produces
**exactly the same 11 residues**, and the resulting `inv` mutations match
**11/11** the paper's list (`I10D, T14W, V18D, A31W, K33G, D86W, K129G, Q131W,
N132W, L134D, D145W`). This is the authoritative correctness check — the
algorithm is right.

### MEK1 / 7XLP (informational)

The paper's MEK1 numbering is offset by +36 from 7XLP's PDB auth_seq:

> paper A40 → 7XLP A76, paper F173 → 7XLP F209, etc. — all 7 paper residues
> match by AA identity at this offset.

Even after applying the offset, only 2 of the paper's 7 listed residues pass
the strict 3.5 Å side-chain rule. The other 5 (M146/S194/F209 etc.) are at
3.53–3.92 Å side-chain to the FZC inhibitor — within 3.5 Å only when
*backbone-inclusive*. So the paper's MEK1 list is **inconsistent with its
own stated 3.5-Å side-chain rule** for this system. The CDK2 result confirms
the rule itself is correct; MEK1 is reported as informational, not a gate.

### Existing `bindingsite_inv` artifact deviation

The shipped CDK2 `bindingsite_inv` data uses Q131→Q (no change), N132→Q,
D145→Q at three positions — inconsistent with the paper's Table 2, which
maps Q/N/D all to W. This module follows the paper's table, not the
artifact. The artifact predates the systematic Methods description.

## End-to-end: subset20

`01_build_subset20.py` builds inputs for the 20-PDB CASF subset
(`PDBbind_casf2016_subset20.json`).

Per-system invariants checked (all 20 pass):

- Sequence length identical across all 4 variants (only point substitutions).
- `rem` mutation count == `pocket_size` (no GLY in pocket by construction).
- `pack` mutation count == `pocket_size − n_PHE_in_pocket` (PHE→F skipped).
- `inv` mutation count == `pocket_size` (no fixed points in INV_TABLE).

Pocket sizes range 2–12 residues across the subset.

## Boltz-1 vs Boltz-2 — env naming caveat

Both Boltz-1 (v0.4.x) and Boltz-2 (v2.2.x) are pip-installed under the
package name `boltz` and produce a `boltz` binary. The env name doesn't
reliably indicate the version:

| host | env path | Boltz version |
|---|---|---|
| lab | `/home/aoxu/miniconda3/envs/boltzina_env/bin/boltz` | **2.2.1** ✓ |
| lab | `/home/aoxu/miniconda3/envs/rdkit_env/bin/boltz` | 0.4.1 (Boltz-1) |
| CARC | `/project2/katritch_223/aoxu/conda/envs/boltzina_env/bin/boltz` | **2.2.1** ✓ |
| CARC | `/project2/katritch_223/aoxu/.conda/envs/boltz2/bin/boltz` | **2.2.1** ✓ |
| CARC | `/project2/katritch_223/aoxu/.conda/envs/boltzina_env/bin/boltz` | **2.2.1** ✓ |

This module **only uses Boltz-2** (the runners always pass
`--model boltz2`). To prevent accidental Boltz-1 selection, every entry
point calls
[`config.assert_boltz2_binary`](../analysis/casf_mutagenesis/config.py)
which imports `boltz` from the binary's adjacent Python and fails with a
clear error if the major version is < 2. The check is run:
- at import time in `scripts/03_run_boltz2_subset20.py` (the dedicated
  Boltz runner),
- lazily on the first MSA fetch in `msa_via_boltz.fetch_msa_via_boltz`.

dockStrat reuses one `boltz` binary for both `--model boltz1` and
`--model boltz2` (its `dockstrat_config/model/{boltz1,boltz2}_inference
.yaml` both point at the same path); the version selection is by flag,
not by binary. This module **does not use dockStrat's runner** — it
calls `boltz predict` directly with explicit Boltz-2 settings.

## Running on a different host (env-var overrides)

All host-specific paths flow through env vars resolved in
[`config.py`](../analysis/casf_mutagenesis/config.py). The defaults match
the lab workstation; on CARC (or any other host), source the matching
shell file before running anything:

```bash
# Lab workstation
source env/lab.sh

# USC CARC
source env/carc.sh
```

The variables (override individually if needed):

| variable | meaning | lab default |
|---|---|---|
| `CONTRASCF_ROOT` | repo root | `/mnt/katritch_lab2/aoxu/contrasCF` |
| `CONTRASCF_CASF_ROOT` | pdbbind_cleansplit dir | `<root>/data/casf2016` (a symlink) |
| `CONTRASCF_BOLTZ_BIN` | Boltz-2 binary | `/home/aoxu/miniconda3/envs/boltzina_env/bin/boltz` |
| `CONTRASCF_AF3_ENV` | AF3 conda env dir | `/mnt/katritch_lab2/aoxu/CogLigandBench/envs/alphafold3` |
| `CONTRASCF_AF3_DIR` | AF3 source (run_alphafold.py here) | `…/forks/alphafold3/alphafold3` |
| `CONTRASCF_AF3_MODEL_DIR` | AF3 weights (af3.bin) | `…/forks/alphafold3/models` |
| `CONTRASCF_CUDA_DEVICE` | GPU index | `0` |
| `CONTRASCF_PY` | rdkit_env python (set by env/*.sh) | `/home/aoxu/miniconda3/envs/rdkit_env/bin/python` |

The AF3 runner globs `lib/python*/site-packages/nvidia/*/lib` rather than
hard-coding Python 3.12, so it adapts to whatever Python version the
target host's AF3 env uses.

To check what's currently active:

```bash
python -c "from casf_mutagenesis import config as c; print(c.REPO_ROOT, c.BOLTZ_BIN, c.AF3_ENV)"
```

## Running the pipeline

All commands assume `LD_LIBRARY_PATH` points at `rdkit_env`'s `lib/`
(set automatically by `env/lab.sh` / `env/carc.sh`).

```bash
LIB="/home/aoxu/miniconda3/envs/rdkit_env/lib"
PY="/home/aoxu/miniconda3/envs/rdkit_env/bin/python"

# 1. Gating test (CDK2 + MEK1)
LD_LIBRARY_PATH=$LIB:$LD_LIBRARY_PATH $PY \
  analysis/casf_mutagenesis/scripts/00_verify_reference_systems.py

# 2. Build inputs for the 20-system subset
LD_LIBRARY_PATH=$LIB:$LD_LIBRARY_PATH $PY \
  analysis/casf_mutagenesis/scripts/01_build_subset20.py

# 3. Run Boltz-2 (single GPU; ~37 s/job × 76 jobs ≈ 47 min)
LD_LIBRARY_PATH=$LIB:$LD_LIBRARY_PATH $PY \
  analysis/casf_mutagenesis/scripts/03_run_boltz2_subset20.py

# 4. Run AF3 (single GPU; ~85 s/job × 76 jobs ≈ 1 h 47 min)
LD_LIBRARY_PATH=$LIB:$LD_LIBRARY_PATH $PY \
  analysis/casf_mutagenesis/scripts/04_run_af3_subset20.py

# 5. Analyze: ligand RMSD-to-native + memorization rate
LD_LIBRARY_PATH=$LIB:$LD_LIBRARY_PATH $PY \
  analysis/casf_mutagenesis/scripts/05_analyze_subset20.py
```

Outputs:

- `outputs/manifest_subset20.json` — per-system pocket + mutations + warnings.
- `outputs/boltz2_run_log.json`, `outputs/af3_run_log.json` — per-cell
  status, wallclock, error string for failures.
- `outputs/results_subset20.csv` — per (pdbid, variant, model) ligand RMSD,
  Cα RMSD, atoms matched.
- `outputs/memorization_subset20.csv` — per (model, variant) memorization
  rate at 2 Å and 4 Å, median RMSD.

## GPU footprint and skip rules

Single RTX 4090 (24 GB).

| skip | reason |
|---|---|
| 3dx2 (1016 residues, homo-tetramer) | total length > 800-residue budget; would OOM |
| 4ih5 / 4de1 / 4ivc — first run only | RDKit `MolFromSmiles` returned None on cached SMILES; fixed in second run via crystal-SDF SMILES derivation |

Per-job timings on RTX 4090 (subset20 wallclock):

| model | mode | mean | min | max | recycles | diffusion samples |
|---|---|---|---|---|---|---|
| Boltz-2 v2.2.1 | no-MSA | ~37 s | ~36 s | ~58 s | 3 | 5 |
| AF3 v3.0.1 | no-MSA | ~90 s | 84 s | 127 s | 10 | 5 |

## How to read these results (paper framing)

Per Masters et al. 2025, the adversarial variants (`rem`, `pack`, `inv`)
are **designed to break ligand binding** by destroying or inverting the
pocket's chemistry and shape. A **physics-aware** model should:

1. Place the ligand correctly on **WT** (high success rate, low RMSD).
2. Place the ligand **somewhere different** on adversarial variants, since
   the original pocket no longer accommodates it. **High ligand RMSD on
   `rem`/`pack`/`inv` is the desired outcome** — it means the model
   responded to the perturbation rather than recalling the training pose.

Therefore:

- **Memorization rate** = fraction of adversarial cases with ligand RMSD <
  threshold (default 2 Å). **Lower is better.** A high memorization rate
  means the model is recalling the WT-like pose despite the disrupted
  pocket — i.e., it has learned global sequence/structure co-occurrence
  more than the physics of binding.
- **WT placement rate** is the success ceiling: the model must get WT right
  for the adversarial comparison to be meaningful.
- A useful model has **WT rate ≫ adversarial rate**. The gap quantifies
  how much the model lets pocket chemistry determine the pose.

This is the inverse of typical "RMSD-to-native" reporting in standard
docking benchmarks, where lower-is-always-better. Don't reuse those plots
without reframing for this dataset.

## Subset20 results (2026-05-06)

All 19 runnable systems × 4 variants × 3 models = 228 cells produced
predictions. Memorization rates on the **14 adversarial** cases per model:

| model | variant | n | median lig RMSD | median Cα RMSD | <2 Å | <4 Å |
|---|---|---|---|---|---|---|
| **AF3+MSA** | wt   | 19 | **0.58 Å** | **1.21 Å** | **0.89** | 0.89 |
| **AF3+MSA** | rem  | 19 | 6.48 Å | 1.46 Å | 0.37 | 0.37 |
| **AF3+MSA** | pack | 19 | 4.05 Å | 1.23 Å | 0.37 | 0.47 |
| **AF3+MSA** | inv  | 19 | 6.88 Å | 1.46 Å | 0.26 | 0.37 |
| Boltz-2 | wt   | 19 | 1.42 Å | 0.86 Å | 0.63 | 0.79 |
| Boltz-2 | rem  | 19 | 6.47 Å | 1.36 Å | 0.37 | 0.42 |
| Boltz-2 | pack | 19 | 1.96 Å | 1.53 Å | 0.53 | 0.58 |
| Boltz-2 | inv  | 19 | 5.66 Å | 1.54 Å | 0.32 | 0.42 |

(Updated 2026-05-17 after the chain-proximity fix below; see "Pocket-
aligned Cα chain selection" in implementation history. Previous tables
under-counted memorisation because the Cα superposition picked the
largest chain rather than the ligand-binding chain on homo-multimers.)

**WT baseline (higher is better — the model should get unperturbed
systems right):** AF3+MSA places the ligand within 2 Å of the crystal
pose in 79 % of WT systems (median 0.68 Å, Cα 1.21 Å); Boltz-2 in 63 %.
AF3 without MSA is broken at 0 %. These set the success ceiling.

**Adversarial memorization at 2 Å (lower is better — models should
respond to disrupted pockets):**
- **AF3+MSA**: 26–37 % across rem/pack/inv. Cα RMSD stays ~1.2–1.5 Å,
  so the protein folds correctly; the variation is in ligand placement.
- **Boltz-2**: 21–37 %, very similar to AF3+MSA on this subset.
- AF3+MSA's pack and rem rates (37 %) are higher than its inv (26 %)
  — the model is more likely to keep the ligand near the original site
  for blanket-Phe (`pack`) or blanket-Gly (`rem`) substitutions, less so
  for the chemistry-flipping `inv` substitutions. The gap WT (79 %) →
  adversarial (~30 %) is roughly 2× — the model **does** notice the
  perturbation in the majority of cases, but recalls the WT-like pose
  in ~1/3 of adversarial cases anyway.

**Subset20 vs paper:** the paper reports AF3 retains the WT pose on
~50 % of adversarial CASF-285 systems. Our subset20 number (~30 % for
AF3+MSA) is lower; sample size (n=19 vs 285) and ColabFold MSA quality
(vs paper's pipeline-generated MSA) plausibly account for the gap. The
**direction of the result holds**: AF3+MSA shows substantial
memorization — about 1 in 3 adversarial cases keeps the ligand near
the WT pose despite the perturbation. Scaling to full CASF-285 would
tighten the comparison.

**The 3dx2 system (1016 residues, homo-tetramer)** is excluded as too
large for the 24 GB GPU budget — 4 cells skipped per model.

**Boltz-2 WT outliers (>2 Å of crystal):** 7/19 systems; AF3+MSA WT
outliers: 4/19. Both models leave a tail of crystal-pose mismatches even
on unperturbed systems. Worth investigating per-system: alternate binding
modes, our heavy-atom MCS-based ligand selection picking the wrong HET
residue, or genuine model misplacement.

## Known limitation: AF3 single-sequence mode

The AF3 runner uses `--norun_data_pipeline --run_inference=true`, which
means we provide MSAs in the JSON ourselves. We currently provide only the
query sequence as a single-row A3M. **AF3 in single-sequence mode produces
structurally bad protein predictions**: 1pxn-WT (CDK2) gives Cα RMSD vs
crystal of 17 Å, where Boltz-2 in the same single-sequence mode gives 0.96 Å.
Boltz-2 trains for single-sequence; AF3 doesn't. AF3's own confidence
metrics confirm the diagnosis — 1pxn-WT reports `ptm=0.22` (good is ≥0.8),
`iptm=0.39`, `ranking_score=0.36`. The model itself is signaling "I don't
know what this protein looks like."

### Infrastructure built (ready for when ColabFold recovers or DBs land)

- [`msa_fetch.py`](../analysis/casf_mutagenesis/msa_fetch.py) — minimal
  ColabFold MSA-server client (no extra Python deps; uses `urllib`).
  Submits a sequence, polls until COMPLETE, unpacks the tar.gz response,
  and returns merged A3M text. Optional disk cache by sequence hash.
- [`inputs_af3.py::render_af3`](../analysis/casf_mutagenesis/inputs_af3.py)
  now accepts a `chain_msas: dict[chain_id, a3m]` argument; when present
  the JSON's `unpairedMsa` field is set to the supplied A3M.
- [`inputs_af3.py::rewrite_a3m_query`](../analysis/casf_mutagenesis/inputs_af3.py)
  — given a WT-derived A3M and a mutated sequence, swaps the query row's
  letters at the mutation positions while preserving the alignment
  columns. So we fetch one MSA per system (for WT) and reuse it across
  the four variants.

### MSA fetch via Boltz piggyback (2026-05-06 fix)

Direct submission to `https://api.colabfold.com/ticket/msa` returned
persistent `status: ERROR` for every sequence — likely an endpoint /
parameter mismatch with what the server now expects. Boltz's
`--use_msa_server` flag works against the same server but uses the right
internal protocol.

Solution: [`msa_via_boltz.py`](../analysis/casf_mutagenesis/msa_via_boltz.py)
piggy-backs on Boltz's MSA fetch. For each unique WT sequence:

1. Write a minimal protein-only Boltz YAML.
2. Run `boltz predict --use_msa_server` with absolute-minimum settings
   (1 sample, 1 recycle, 5 sampling steps) so structure inference is
   trivial and the dominant cost is MSA download.
3. Read `<out>/boltz_results_<key>/msa/<key>_unpaired_tmp_env/uniref.a3m`
   and the bfd/mgnify counterpart, concatenate, **strip null bytes**
   (one Boltz row in ~10⁴ contains a `\\x00` that AF3 then rejects as
   "Unknown residues"), cache by sequence sha1.

For variants of a given pdbid, we reuse the WT MSA via
[`inputs_af3.rewrite_a3m_query`](../analysis/casf_mutagenesis/inputs_af3.py):
swap only the query row's letters at the mutated positions, leave the
alignment columns unchanged. AF3 sees the variant residues as
non-conserved positions in an otherwise identical MSA.

### A3M cleaning required: per-record shape filter

AF3's `MSA.featurize` raises `Invalid shape — all strings must have the
same number of non-lowercase characters` when a record has more
non-lowercase chars than the query (uppercase + `-` + `X`; lowercase
letters are insertions and don't count). Boltz/MMseqs2 occasionally
returns hits where a longer protein has two regions aligning to the
query, encoded as one A3M record of length ~2*Q. Two further pitfalls:

1. Some hits contain a stray `\\x00` null byte → AF3 rejects as "Unknown
   residues".
2. A3M sequences may wrap to multiple lines (one record's sequence spans
   header + 1..N data lines); a naive line-by-line filter sees a
   correct-shape line and fails to drop the record.

[`inputs_af3.clean_a3m_for_af3`](../analysis/casf_mutagenesis/inputs_af3.py)
parses A3M as records (header + concatenated sequence lines), strips null
bytes, and drops records whose non-lowercase count != Q. Always keeps
the first record (the query). Applied AFTER `rewrite_a3m_query` so the
mutated query is preserved.

### Smoke test: AF3 + MSA on 1pxn-WT (CDK2)

| metric | no-MSA | with MSA | crystal target |
|---|---|---|---|
| Cα RMSD vs crystal | 17.11 Å | **0.96 Å** | — |
| pTM | 0.22 | **0.96** | — |
| iPTM | 0.39 | 0.97 | — |
| ranking_score | 0.36 | 0.97 | — |

The fix is real. AF3 produces paper-quality structures with the MSA
injected; the no-MSA results were genuinely uninterpretable. The
[`scripts/06_run_af3_msa_subset20.py`](../analysis/casf_mutagenesis/scripts/06_run_af3_msa_subset20.py)
runner reproduces this for all 76 (system, variant) cells and writes
`af3msa_<prefix>_*` files alongside the no-MSA `af3_*` baseline so both
can be compared.

Per-system cost: ~2 s MSA fetch (Boltz piggyback) + ~90 s AF3 inference
≈ 95 s. Total subset20 wallclock ~2 h.

## Implementation history (2026-05-05 → 2026-05-06)

The module was built incrementally; six bugs surfaced during integration
and were diagnosed against the paper's stated CDK2/MEK1 ground-truth
residue lists or against Boltz-2 results that *did* work. Recording the
order so future work can reuse the diagnostic patterns.

### Day 1 (2026-05-05) — input generation

1. **Scaffold.** Created `analysis/casf_mutagenesis/` with `config.py`
   (mutation tables from the paper Methods + Table 2), `sequence.py`
   (gemmi-based protein → per-chain sequence), `pocket.py` (3.5 Å
   side-chain heavy-atom rule), `mutate.py` (apply mutations preserving
   length, skip identity substitutions). Symlinked
   `data/casf2016 → /home/aoxu/projects/VLS-Benchmark-Dataset/data/pdbbind_cleansplit`.
2. **CDK2 verification (gating).** Ran auto-detection on 1B38; got
   exactly the paper's 11 residues; `inv` mutations matched 11/11. The
   algorithm was right.
3. **MEK1 / 7XLP discrepancy.** Found a +36 offset between paper's
   numbering and 7XLP's PDB auth_seq (paper A40 = 7XLP A76, etc.) — and
   the paper's MEK1 list is mildly inconsistent with its own 3.5-Å
   side-chain rule (M146/S194/F209 are at 3.5–3.9 Å side-chain, only
   within 3.5 Å backbone-inclusive). Treated MEK1 as informational only;
   CDK2 stays the strict gate.
4. **Renderers.** `inputs_af3.py` (AF3 local dialect; emits `version: 2`
   with no-MSA placeholder fields), `inputs_boltz.py` (`msa: empty` per
   chain), `inputs_docking.py` (WT-only receptor.pdb + ligand.sdf +
   box.json, all-chains-preserved variant of the existing 16-case
   strip-to-protein helper).
5. **Boltz-2 run on subset20.** 64 ok / 4 skip (3dx2 1016 res too large)
   / 12 fail. The 12 failures were 3 systems × 4 variants where Boltz's
   internal `MolFromSmiles` returned None on the cached HET SMILES
   (`het_smiles_cache.json`). Same systems had emitted RDKit
   kekulization warnings during docking-input generation.
6. **AF3 (no-MSA) run on subset20.** 64 ok / 4 skip / 12 fail (same 3
   systems — RDKit fails the same way for both Boltz and AF3).

### Day 2 (2026-05-06) — analysis + AF3 fix

7. **Bug 1 (SMILES kekulization).** Fix: derive SMILES from the crystal
   ligand SDF via `Chem.MolToSmiles(mol_from_sdf)`. The SDF parses fine;
   the round-trip produces a canonical, RDKit-friendly SMILES that both
   Boltz and AF3 accept. Re-ran the 12 failed jobs in each model; 0
   remaining failures.
8. **Analysis pass.** Wrote `analysis.py` (`analyze_prediction`,
   `memorization_stats`) and `scripts/05_analyze_subset20.py`. First run
   produced wildly wrong numbers — Boltz-2 WT median 9 Å, AF3 17 Å.
9. **Bug 2 (Cα superpose).** Diagnosed that `analysis.src.align.
   superpose_ca` only scans ±50 PDB-resnum offset and CASF crystals use
   UniProt-canonical numbering — 4jia chain A is residue **833–1132**
   while predictions are 1..N (offset 832). Five systems silently
   returned 0 paired Cα; others got spurious matches and wrong RMSD.
   Fix: wrote `casf_mutagenesis.analysis.superpose_by_index` that pairs
   by chain-position index (legitimate because the predicted sequence
   *is* the crystal chain's sequence by construction). Re-ran analysis.
10. **AF3 results obviously wrong.** Boltz-2 WT median 1.5 Å (sane);
    AF3 WT median 17 Å, AF3 adversarial 0 % at 2 Å. The user flagged
    this as suspicious before it was clear from the data alone. AF3's
    own confidence reported `ptm = 0.22` (good is ≥ 0.8) on 1pxn-WT —
    the model itself was signalling "I have no idea what this protein
    looks like."
11. **Bug 3 (no-MSA failure).** Diagnosed: AF3 single-sequence mode
    (`--norun_data_pipeline`) doesn't generalise the way Boltz-2 does
    (which trains for single-sequence). The fix is to provide an MSA.
    Genetic databases are not installed locally (~365 GB); started
    looking for a remote MSA source.
12. **Bug 4 (ColabFold direct API → ERROR).** First attempt:
    `casf_mutagenesis.msa_fetch` with direct `POST /ticket/msa`. Every
    submission — including the simplest test sequence (76-AA ubiquitin)
    — returned `status: ERROR`, persistent (the API deduplicates by
    sequence hash). Server-side outage / endpoint mismatch from this
    machine.
13. **Workaround (Boltz piggyback).** Discovered that `boltz predict
    --use_msa_server` works against the same API but via a different
    internal protocol. Wrote `msa_via_boltz.py`: run Boltz with
    minimum-cost settings (1 sample / 1 recycle / 5 sampling steps) and
    pull the resulting `boltz_results_<key>/msa/<key>_unpaired_tmp_env/
    {uniref,bfd.mgnify30...}.a3m`. Cache by sequence sha1.
14. **Bug 5 (`\x00` null byte in MSA).** First AF3 run with MSA crashed
    on row 8823 of the 1pxn uniref MSA — one residue position contained
    a null byte. Fix: `a3m.replace("\x00", "")`. AF3 then accepted the
    cleaned MSA; pTM jumped from 0.22 to **0.96**, Cα RMSD from 17 Å to
    **0.96 Å** on 1pxn-WT. The diagnosis was right.
15. **Bug 6 (multi-line A3M records).** Bulk run failed all 76 jobs
    with "Invalid shape — non-lowercase character count != query
    length". Root cause: some MSA records' sequences wrap to multiple
    lines (one record's full sequence = concatenation of all data lines
    until next `>`); a line-by-line cleaner sees correct-shape lines
    and fails to drop the actual bad record. Fix: rewrite
    `clean_a3m_for_af3` to parse A3M as `(header, full_sequence)`
    records, drop records by full-sequence shape. Re-ran; 76/76 ok.
16. **Re-ran analysis.** AF3+MSA WT 79 % at 2 Å (median 0.68 Å);
    adversarial 26–37 %. Boltz-2 WT 63 %, adversarial 21–37 %. AF3
    no-MSA still 0 % everywhere (kept as the broken-baseline column).

### Day 12 (2026-05-17) — pocket-aligned Cα chain selection

17. **WT outlier audit.** 5 systems (4w9l, 4bkt, 3k5v, 3ao4, 3b5r) had
    WT ligand RMSD ≫ 2 Å in both Boltz-2 and AF3+MSA. The pattern —
    same systems failing in both models AND sometimes-clean Cα RMSD
    (4w9l-WT: AF3+MSA lig 47.7 Å on Cα 0.66 Å) — pointed at an analysis
    bug rather than two independent model failures.
18. **Bug 7 (Cα chain selection on multimers).**
    `analysis.src.loaders.extract_protein_ca` picks the largest polymer
    chain. CASF includes homo-multimers (4w9l is a homo-trimer of
    chains C/K/L); crystal binds chain C, predicted-Boltz binds chain
    L (equivalent pocket, different chain). Picking "largest" gives the
    same chain on both sides (L), the superpose is fine, but the ligand
    sits on the *other* equivalent chain → 47 Å apparent RMSD on a
    correctly-placed pose. Fix:
    `casf_mutagenesis.analysis.extract_protein_ca_near(st, lig_xyz)` —
    pick the chain whose nearest Cα is closest to the ligand centroid.
    Used independently for crystal and predicted; the model may bind a
    different equivalent chain than the crystal, but the question
    "did the ligand land at AN equivalent pocket?" is chain-agnostic
    on symmetric multimers.
19. **Net effect on subset20:** AF3+MSA WT 0.79 → **0.89** (+2 systems
    recovered); Boltz-2 pack 0.37 → 0.53 (+3 systems — pocket-aligned
    reveals more memorisation on `pack` than we previously credited).
    Remaining WT outliers (3k5v, 3ao4, 3b5r) all have either bad Cα
    fold (>5 Å) or ligand placement on the right fold but wrong site —
    real model failures, not pipeline issues.

### Diagnostic patterns to reuse

- **When numbers look wrong, distrust the analysis pipeline before the
  model.** Two of the six bugs (Cα superpose, multi-line A3M cleaner)
  were in our *post-processing*, not in the model run.
- **Always sanity-check at the smallest scale first.** Bug 4 (Cα
  superpose) was caught by checking 4jia's residue ranges directly
  (833–1132 vs 1..N predicted). Bug 6 was caught by inspecting AF3's
  exact error message — it printed the offending row's exact contents,
  which made it possible to count its actual non-lowercase length.
- **When AF3 reports `ptm ≪ 0.5`**, the model is signalling its own
  predictions are unreliable. Don't try to interpret RMSDs on top of
  bad folds.
- **The CDK2/1B38 11-residue paper-stated list is the canonical
  correctness check** for the pocket-detection logic. MEK1's list has
  the +36 offset and is unsuitable as a strict gate.

## Out of scope (deferred)

- **Mutant-protein docking inputs** — need AF3 mutant structures first.
- **Chai / RFAA inputs** — skipped per scope decision.
- **Running the models** is split out from the input-generation step; this
  doc's "Running the pipeline" section covers both, but the build/run
  separation is intentional so input generation is cheap and re-runnable.
- **Full CASF-2016 (285 systems)** — `02_build_full_casf.py` not yet
  written; mirrors `01_` but iterates the 285 IDs in
  `labels/PDBbind_data_split_cleansplit.json["casf2016"]`. Estimated
  wallclock at the current per-job rates: ~3 h Boltz-2 + ~7.5 h AF3.
- **Analysis pipeline integration** — once subset20 results are validated,
  fold the analysis CSVs into `analysis/results/` plotting helpers
  (`04_render_figures.py`, `05_physics_analysis.py`).
