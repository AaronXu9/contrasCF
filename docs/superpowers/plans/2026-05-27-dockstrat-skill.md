# dockstrat Skill Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create a user-global Claude Code skill at `~/.claude/skills/dockstrat/` documenting the dockStrat `dock_engine()` API, per-method conda envs / outputs / gotchas, with a `SKILL.md` + lazy-loaded `references/<method>.md` layout.

**Architecture:** One `SKILL.md` (~250 lines) carrying triggers, decision tree, env conventions, three quickstart snippets, four gotchas, and a method index. Seven `references/*.md` files (~80–150 lines each) loaded only when a specific method is the subject of the task. No code module — the gotchas are documented inline.

**Tech Stack:** Markdown only. Content sourced from `/mnt/katritch_lab2/aoxu/CogLigandBench/` (dockStrat repo CLAUDE.md + `dockstrat/engine.py` + `dockstrat_config/model/*.yaml`) and from working scripts in `/mnt/katritch_lab2/aoxu/contrasCF/analysis/` (especially `scripts/11_run_docking.py`, `scripts/13_run_surfdock.py`, `casf_mutagenesis/scripts/11_run_unidock2_variants.py`).

**Spec:** `docs/superpowers/specs/2026-05-27-dockstrat-skill-design.md`.

**Note on commits:** `~/.claude/` is NOT a git repo on this host, so the skill files themselves cannot be per-file committed there. Only the plan file and the final smoke-test notes (in the contrasCF repo) are committed. Each task ends with a verification step instead of a commit.

---

## File Structure

```
~/.claude/skills/dockstrat/
  SKILL.md                          # Task 2-6 (one section per task)
  references/
    unidock2.md                     # Task 7
    gnina.md                        # Task 8
    surfdock.md                     # Task 9
    boltz.md                        # Task 10 (boltz1+boltz2 share env/CLI)
    af3.md                          # Task 11
    vina.md                         # Task 12
    other-methods.md                # Task 13 (chai/dynamicbind/protenix/icm/icm-rtcnn)
```

Each file has one responsibility (per-method reference, or one SKILL.md section). Files are independent — editing one does not require editing others, except cross-link checks (Task 14).

---

## Task 1: Create directory skeleton

**Files:**
- Create: `~/.claude/skills/dockstrat/`
- Create: `~/.claude/skills/dockstrat/references/`

- [ ] **Step 1: Create directories**

```bash
mkdir -p ~/.claude/skills/dockstrat/references
```

- [ ] **Step 2: Verify**

```bash
ls -la ~/.claude/skills/dockstrat/
```

Expected: `references/` subdir exists. No other files yet.

---

## Task 2: Write SKILL.md — frontmatter + activation triggers

**Files:**
- Create: `~/.claude/skills/dockstrat/SKILL.md`

- [ ] **Step 1: Write the frontmatter and the first section ("When to use this skill")**

Create the file with exactly this content (note: more sections will be appended in subsequent tasks):

```markdown
---
name: dockstrat
description: Use when running protein-ligand docking via the dockStrat framework (/mnt/katritch_lab2/aoxu/CogLigandBench). Covers the dock_engine() Python API for UniDock2, GNINA, SurfDock, AlphaFold3, Boltz-1/2, and Vina, including per-method conda envs, key kwargs, output formats, and known gotchas (GPU desync, SurfDock pocket-center anchoring, rdkit_env import bypass, LD_LIBRARY_PATH). Trigger on mentions of dock_engine, dockStrat, CogLigandBench, or any supported method (UniDock2, GNINA, SurfDock, AlphaFold3, Boltz, Protenix, Chai, DynamicBind), or any task that needs to produce ranked poses for a protein+ligand pair.
---

# dockstrat

A reference for running protein-ligand docking via the **dockStrat** framework
(`/mnt/katritch_lab2/aoxu/CogLigandBench`). The unified entry point is
`dock_engine()` in `dockstrat/engine.py`.

## When to use this skill

Fire on any of:

- A task that needs ranked poses for a protein + ligand pair.
- Any mention of `dock_engine`, `dockStrat`, or `CogLigandBench`.
- Any mention of a supported method: UniDock2, GNINA, SurfDock, AlphaFold3,
  Boltz-1/2, Protenix, Chai, DynamicBind, Vina, ICM.
- Debugging a docking run that failed (timeout, OOM, "pose far from pocket",
  CUDA error).

**Out of scope (point readers elsewhere):**

- Input preparation (PDB cleaning, SMILES → 3D SDF, binding-box detection) —
  see `dockstrat/data/` or `contrasCF/docking/smoketest/smoketest.py`.
- RMSD / scoring / downstream analysis — see `dockstrat/analysis/`.
- Method installation — see `scripts/install_<method>_env.sh` in the dockStrat
  repo.
```

- [ ] **Step 2: Verify the file is well-formed**

```bash
head -5 ~/.claude/skills/dockstrat/SKILL.md
wc -l ~/.claude/skills/dockstrat/SKILL.md
```

Expected: `---` on line 1, frontmatter closes at line 4, body begins; file is ~30 lines.

---

## Task 3: Append sections "Project layout", "Method-selection decision tree", "Env conventions"

**Files:**
- Modify: `~/.claude/skills/dockstrat/SKILL.md` (append)

- [ ] **Step 1: Append these three sections**

Use the Edit tool or `cat >>` to append exactly this content to the end of `SKILL.md`:

````markdown

## Project layout you can assume

```
${COGLIGANDBENCH}=/mnt/katritch_lab2/aoxu/CogLigandBench
${COGLIGANDBENCH}/dockstrat/engine.py            # dock_engine entry point
${COGLIGANDBENCH}/dockstrat/models/              # per-method inference wrappers
${COGLIGANDBENCH}/dockstrat_config/model/        # per-method YAML configs
${COGLIGANDBENCH}/envs/<method>/                 # symlinks to conda envs
${COGLIGANDBENCH}/forks/<method>/                # method source / binaries
```

`dock_engine` is imported as `from dockstrat import dock_engine`. This requires
either `pip install -e ${COGLIGANDBENCH}` into the calling env, or
`${COGLIGANDBENCH}` on `PYTHONPATH`. The default env for this is `dockstrat`
(create with `conda create -n dockstrat python=3.10 && pip install -e .`).

`SUPPORTED_METHODS = ("vina", "gnina", "chai", "dynamicbind", "unidock2", "surfdock", "alphafold3", "boltz1", "boltz2", "protenix")` — note `icm` and `icm-rtcnn` are NOT exposed via `dock_engine`; they ship as standalone scripts in `forks/ICM/`.

## Method-selection decision tree

1. Need physics-only, fast, GPU available → `unidock2`.
2. Need physics + CNN rescoring → `gnina`.
3. Need surface/pocket features (OOD ligand, novel pocket) → `surfdock`.
4. Need a co-folded structure together with the pose:
   - Affinity prediction wanted → `boltz2`.
   - AF3 specifically required → `alphafold3`.
   - Boltz-1 baseline → `boltz1`.
5. CPU-only smoke test, no GPU → `vina`.
6. Anything else (chai, dynamicbind, protenix, icm, icm-rtcnn) → see
   `references/other-methods.md`.

## Env conventions (mandatory on this host)

These are restated here so the skill is self-contained. They mirror
`env_conventions.md` and `tool_locations.md` in the project memory.

### LD_LIBRARY_PATH prepend (every Python invocation)

Before running any Python that imports `rdkit`, `gemmi`, `pymol`, or
`dockstrat`, prepend the env's `lib/` to `LD_LIBRARY_PATH`. Without it you
will hit `libstdc++.so.6: version 'GLIBCXX_3.4.30' not found` or similar.

```bash
export LD_LIBRARY_PATH="/home/aoxu/miniconda3/envs/<env>/lib:$LD_LIBRARY_PATH"
```

### Call Python via absolute path, not `conda run -n`

`conda run -n <env>` has a >1 s activation overhead on this machine. In
scripts and subprocess calls, call the env's Python directly:

```bash
/home/aoxu/miniconda3/envs/<env>/bin/python script.py
```

`conda run -n <env>` is fine for one-off interactive shell commands and is
what `dock_engine` itself uses internally for some methods (e.g., `unidock2`).

### Which env to call `dock_engine` from

`dock_engine` is a thin dispatcher that subprocesses out to the per-method
env. Call it from any env that has `dockstrat` installed (e.g., `dockstrat`).
You do not need to activate the per-method env yourself.

If you want to call `dockstrat` *submodules* from an env that does NOT have
`rootutils` installed (e.g., `rdkit_env`), see gotcha 6b below for the
import bypass.
````

- [ ] **Step 2: Verify**

```bash
grep -c "^## " ~/.claude/skills/dockstrat/SKILL.md
```

Expected: `4` (When to use, Project layout, Method-selection, Env conventions).

---

## Task 4: Append "Quickstart by family" section

**Files:**
- Modify: `~/.claude/skills/dockstrat/SKILL.md` (append)

- [ ] **Step 1: Append the quickstart section**

````markdown

## Quickstart by family

### Physics (UniDock2)

```python
from dockstrat import dock_engine

dock_engine(
    'unidock2',
    protein='receptor.pdb',
    ligand='ligand.sdf',
    output_dir='./out',
    num_poses=20,             # default 9
    box_size=[20, 20, 20],    # default [15, 15, 15] Å
    search_mode='detail',     # detail | balance | fast
)
# Outputs: ./out/unidock2/<ligand-stem>/rank{N}.sdf  (ranked by Vina energy, lowest best)
# GPU required. Wall-clock: ~10-60 s per pose for a typical drug-like ligand.
```

### Hybrid (GNINA)

```python
from dockstrat import dock_engine

dock_engine(
    'gnina',
    protein='receptor.pdb',
    ligand='ligand.sdf',
    output_dir='./out',
    num_modes=9,
    exhaustiveness=8,
    cnn_scoring='rescore',    # rescore | refine | metrorescore | metrorefine | none
    seed=42,
)
# Outputs: ./out/gnina/<ligand-stem>/pose{N}_score{S}.sdf  (ranked by CNNscore, highest best)
# GPU strongly recommended; CPU works but slow. Wall-clock: ~30-120 s.
```

### Co-folding (Boltz-2)

```python
from dockstrat import dock_engine

dock_engine(
    'boltz2',
    protein='receptor.pdb',
    ligand='ligand.sdf',
    output_dir='./out',
    diffusion_samples=5,
    recycling_steps=3,
    sampling_steps=200,
    num_poses_to_keep=5,
    cuda_device_index=0,
)
# Outputs: ./out/boltz2/<system>/rank{N}.sdf  (ranked by confidence_score, highest best)
# Affinity (pred_value + binary probability) attached as SDF properties.
# GPU required (≥24 GB VRAM). Wall-clock: ~3-10 min per pair.
```

For per-method kwargs, output details, and the full set of knobs, **read
the corresponding `references/<method>.md` file before generating code**.
````

- [ ] **Step 2: Verify quickstart blocks parse**

```bash
grep -c '^```python' ~/.claude/skills/dockstrat/SKILL.md
```

Expected: `3` (three code blocks).

---

## Task 5: Append "Gotchas" section with all four gotchas

**Files:**
- Modify: `~/.claude/skills/dockstrat/SKILL.md` (append)

- [ ] **Step 1: Append the gotchas section**

````markdown

## Gotchas

### 6a. LD_LIBRARY_PATH prepend

Already covered under "Env conventions" — repeated here because forgetting it
is the #1 source of `from rdkit import Chem` failures. If you see
`ImportError: libstdc++.so.6: version 'GLIBCXX_*' not found`, you forgot
this.

```bash
export LD_LIBRARY_PATH="/home/aoxu/miniconda3/envs/<env>/lib:$LD_LIBRARY_PATH"
/home/aoxu/miniconda3/envs/<env>/bin/python script.py
```

### 6b. rdkit_env → dockstrat import bypass

`dockstrat/__init__.py` imports `rootutils`. If the caller env doesn't have
`rootutils` installed (e.g., `rdkit_env`), `from dockstrat...` will fail. Use
this importlib trick to load a single dockstrat submodule without touching
`dockstrat/__init__.py`:

```python
import importlib.util
import sys
import types
from pathlib import Path

DOCKSTRAT_ROOT = Path("/mnt/katritch_lab2/aoxu/CogLigandBench")

def _load_surfdock_module():
    # Register stub parent packages so the loaded module's dotted name resolves.
    for name, path in [
        ("dockstrat", DOCKSTRAT_ROOT / "dockstrat"),
        ("dockstrat.models", DOCKSTRAT_ROOT / "dockstrat" / "models"),
    ]:
        if name not in sys.modules:
            m = types.ModuleType(name)
            m.__path__ = [str(path)]
            sys.modules[name] = m
    sf_path = DOCKSTRAT_ROOT / "dockstrat" / "models" / "surfdock_inference.py"
    spec = importlib.util.spec_from_file_location(
        "dockstrat.models.surfdock_inference", str(sf_path),
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules["dockstrat.models.surfdock_inference"] = mod
    spec.loader.exec_module(mod)
    return mod
```

Source: `contrasCF/analysis/scripts/13_run_surfdock.py::_load_surfdock_module`.

### 6c. UniDock2 GPU-desync guard (batch runs only)

When running UniDock2 over dozens of cells in a loop, a single bad cell can
trip the NVIDIA driver into an NVML-mismatch state. Subsequent cells then
fail in <5 s with CUDA-init errors. To avoid burning hundreds of cells:

```python
import subprocess, time

TIMEOUT_S = 120                  # per-cell hard cap
FAST_ERROR_ABORT = 5             # consecutive <5 s errors → abort

n_consec_fast_errors = 0
for cell in cells:
    t0 = time.time()
    try:
        log = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT_S)
        wc = time.time() - t0
        if log.returncode != 0:
            if wc < 5.0:
                n_consec_fast_errors += 1
            else:
                n_consec_fast_errors = 0
            if n_consec_fast_errors >= FAST_ERROR_ABORT:
                raise RuntimeError("GPU desync — recover driver before retrying.")
        else:
            n_consec_fast_errors = 0
    except subprocess.TimeoutExpired:
        pass  # skip this cell; loop continues
```

Source: `contrasCF/analysis/casf_mutagenesis/scripts/11_run_unidock2_variants.py`.
Recovery from the desync state typically requires `sudo nvidia-smi` reset or
a GPU reboot.

### 6d. SurfDock pocket-center anchoring

`dock_engine('surfdock', ...)` runs the 4-step pipeline (surface → CSV →
ESM → diffusion) in a tempdir, which (i) loses SurfDock's score-encoded
filenames, and (ii) does NOT anchor the diffusion's initial position to
the pocket. For OOD ligands (e.g., propyl-ATP, penta-methyl glucose) the
diffusion can place poses thousands of Å away from the pocket.

Fix: drive the 4 steps explicitly with a persistent work dir, inject a
`pocket_center` column into the input CSV, and pass
`--ligand_to_pocket_center` to `inference_accelerate.py`. See
`contrasCF/analysis/scripts/13_run_surfdock.py` for the full worked example;
`references/surfdock.md` has the recipe.

Symptom: rank1 RMSD > 10 Å on a pocket where every other method gives < 5 Å.
````

- [ ] **Step 2: Verify**

```bash
grep -c '^### 6' ~/.claude/skills/dockstrat/SKILL.md
```

Expected: `4`.

---

## Task 6: Append "Per-method reference index" section (closes SKILL.md)

**Files:**
- Modify: `~/.claude/skills/dockstrat/SKILL.md` (append, final section)

- [ ] **Step 1: Append the index**

````markdown

## Per-method reference index

Before generating method-specific code, **read the corresponding
`references/<method>.md` file**.

| Method | Env (under `${COGLIGANDBENCH}/envs/`) | Output filename | "Use when" | Detail |
|--------|------|-----------------|------------|--------|
| `unidock2`   | `unidock2`   | `rank{N}.sdf` | Physics, fast GPU | `references/unidock2.md` |
| `gnina`      | `gnina`      | `pose{N}_score{S}.sdf` | Physics + CNN rescoring | `references/gnina.md` |
| `surfdock`   | `surfdock` (local: `/home/aoxu/miniconda3/envs/SurfDock`) | `rank{N}.sdf` | Surface/pocket features, OOD ligand | `references/surfdock.md` |
| `alphafold3` | `alphafold3` | `rank{N}.sdf` | Cofold (AF3 specifically) | `references/af3.md` |
| `boltz1`     | `boltz`      | `rank{N}.sdf` | Cofold (Boltz-1 baseline) | `references/boltz.md` |
| `boltz2`     | `boltz`      | `rank{N}.sdf` | Cofold + affinity | `references/boltz.md` |
| `vina`       | `vina`       | `pose{N}_score{S:.2f}.sdf` | CPU baseline | `references/vina.md` |
| `chai`, `dynamicbind`, `protenix`, `icm`, `icm-rtcnn` | varies | varies | see full ref | `references/other-methods.md` |

The full source of truth for any method is `${COGLIGANDBENCH}/CLAUDE.md` —
this skill summarizes. If a kwarg or env name looks stale, verify against
that file before suggesting it to the user.
````

- [ ] **Step 2: Verify the full SKILL.md**

```bash
wc -l ~/.claude/skills/dockstrat/SKILL.md
grep -c '^## ' ~/.claude/skills/dockstrat/SKILL.md
```

Expected: 200-300 lines; 7 `## ` headers (When to use, Project layout, Method-selection, Env conventions, Quickstart, Gotchas, Per-method reference index).

---

## Task 7: Write `references/unidock2.md`

**Files:**
- Create: `~/.claude/skills/dockstrat/references/unidock2.md`

- [ ] **Step 1: Write the file**

```markdown
# UniDock2

Modern GPU-accelerated AutoDock-style docking. SDF in, SDF out. No PDBQT
prep needed. Single-pair runs are typically 10-60 s.

## Env

- Conda env: `${COGLIGANDBENCH}/envs/unidock2` → `/home/aoxu/miniconda3/envs/unidock2`
- Binary: `/home/aoxu/miniconda3/envs/unidock2/bin/unidock2`
- GPU required.
- Install: `bash scripts/install_unidock2_env.sh`.

## dock_engine signature

```python
from dockstrat import dock_engine

dock_engine(
    'unidock2',
    protein='receptor.pdb',
    ligand='ligand.sdf',
    output_dir='./out',
    num_poses=9,                 # ranked SDFs to keep
    box_size=[15.0, 15.0, 15.0], # Å; default [15,15,15]
    search_mode='detail',        # detail | balance | fast
    energy_range=15.0,           # kcal/mol window
    timeout=1800,                # per-system hard cap in seconds
    cuda_device_index=0,
)
```

## Output

- `${output_dir}/unidock2/<ligand-stem>/rank{N}.sdf` — one SDF per pose.
- Also writes `{prefix}_poses.sdf` (raw multi-pose).
- Ranked by Vina binding free energy, lowest = best.
- Energy is in the SDF as a property line.

## Knobs that matter

| kwarg | default | typical override |
|-------|---------|-------------------|
| `num_poses` | 9 | 20 for diverse pose set |
| `box_size` | [15,15,15] | [20,20,20] for large/flexible ligands |
| `search_mode` | detail | `fast` for sweeps (loses some accuracy) |
| `energy_range` | 15.0 | rarely changed |
| `cuda_device_index` | 0 | set per-GPU on multi-GPU hosts |
| `timeout` | 1800 | 120 in batch loops (see gotcha 6c) |

## Pitfalls

- **Box size must come from YAML**, not the CLI. `dock_engine` writes a tiny
  per-call YAML to set `Settings.size` because the `unidock2 docking` CLI
  only accepts `-c` (center) on the command line.
- **GPU desync on batch loops** — see SKILL.md gotcha 6c. Always set a
  per-cell timeout and abort on consecutive fast errors.
- **CLI box defaults to [30, 30, 30]** if the YAML override is missing;
  this surprises people writing custom subprocess calls.

## Source of truth

`${COGLIGANDBENCH}/CLAUDE.md` → "UniDock2" section.
`${COGLIGANDBENCH}/dockstrat_config/model/unidock2_inference.yaml`.
```

- [ ] **Step 2: Verify**

```bash
wc -l ~/.claude/skills/dockstrat/references/unidock2.md
```

Expected: 50-100 lines.

---

## Task 8: Write `references/gnina.md`

**Files:**
- Create: `~/.claude/skills/dockstrat/references/gnina.md`

- [ ] **Step 1: Write the file**

```markdown
# GNINA

Hybrid docking: AutoDock-style sampling + CNN rescoring. Reads PDB+SDF
natively (no PDBQT prep). Single-pair runs are ~30-120 s.

## Env

- Conda env: `${COGLIGANDBENCH}/envs/gnina` (thin wrapper around the binary).
- Binary: **hardcoded** at `${COGLIGANDBENCH}/forks/GNINA/gnina`. NOT under
  `envs/gnina/bin/`. On this host an older copy lives at
  `/home/aoxu/projects/PoseBench/forks/GNINA/gnina` — the dockStrat one is
  the canonical path; older scripts may reference the PoseBench path.
- GPU recommended; CPU works but slow.
- Install: `bash scripts/install_gnina_env.sh`.

## dock_engine signature

```python
from dockstrat import dock_engine

dock_engine(
    'gnina',
    protein='receptor.pdb',
    ligand='ligand.sdf',
    output_dir='./out',
    num_modes=9,
    exhaustiveness=8,
    cnn_scoring='rescore',       # rescore | refine | metrorescore | metrorefine | none
    seed=42,
    autobox_ligand='reference_ligand.sdf',  # optional: use reference for box
)
```

## Output

- `${output_dir}/gnina/<ligand-stem>/pose{N}_score{S}.sdf` — one SDF per pose.
- Ranked by CNNscore (highest = best); the score is in the filename.
- Intermediate: `docked.sdf.gz`, `docked.sdf`.

## Knobs that matter

| kwarg | default | typical override |
|-------|---------|-------------------|
| `num_modes` | 9 | 20 for diverse pose set |
| `exhaustiveness` | 8 | 16 for higher-quality sampling |
| `cnn_scoring` | rescore | `refine` for marginal accuracy gain at higher cost |
| `seed` | 42 | for reproducibility |
| `autobox_ligand` | uses input ligand | point at the crystal ligand if known |

## Pitfalls

- **Hardcoded binary path** — `${COGLIGANDBENCH}/forks/GNINA/gnina`. If the
  binary is missing the install script needs to be re-run.
- **Box detection** — `--autobox_ligand` reads the reference ligand's
  bounding box + 4 Å padding. If your input ligand is at the origin (e.g.,
  freshly RDKit-embedded), the box will be at the origin too, not at the
  binding site. Provide an `autobox_ligand=<crystal_ligand.sdf>` or pre-translate
  the ligand into the pocket.
- **CNN models** ship inside the GNINA binary; no separate weights file to
  worry about.

## Source of truth

`${COGLIGANDBENCH}/CLAUDE.md` → "GNINA" section.
`${COGLIGANDBENCH}/dockstrat_config/model/gnina_inference.yaml`.
```

- [ ] **Step 2: Verify**

```bash
wc -l ~/.claude/skills/dockstrat/references/gnina.md
```

Expected: 50-100 lines.

---

## Task 9: Write `references/surfdock.md`

**Files:**
- Create: `~/.claude/skills/dockstrat/references/surfdock.md`

- [ ] **Step 1: Write the file**

```markdown
# SurfDock

Surface-aware diffusion docking. 4-step pipeline (surface → CSV → ESM → diffusion).
Single-pair runs are ~3-10 min. Sensitive to where you put the input ligand —
see Pitfalls.

## Env

- Conda env: `SurfDock` (local: `/home/aoxu/miniconda3/envs/SurfDock`).
  NOT a symlink under `${COGLIGANDBENCH}/envs/`.
- Local install path: `/home/aoxu/projects/SurfDock` — has the editable ESM,
  MSMS/APBS binaries, pdb2pqr. The dockStrat submodule at
  `${COGLIGANDBENCH}/forks/SurfDock/` only contains code + model weights.
  Use the **local install** as `surfdock_dir`, NOT the submodule.
- Required env vars at runtime:
  - `precomputed_arrays=/home/aoxu/projects/precomputed/precomputed_arrays`
  - `PROJECT_ROOT=${COGLIGANDBENCH}` (for OmegaConf interpolation)
- Model weights live at `${COGLIGANDBENCH}/forks/SurfDock/model_weights/{docking,posepredict}`.
- GPU required.

## dock_engine signature

```python
from dockstrat import dock_engine

dock_engine(
    'surfdock',
    protein='receptor.pdb',
    ligand='ligand.sdf',
    output_dir='./out',
    samples_per_complex=40,
    num_poses=10,
    batch_size=40,
    mdn_dist_threshold=3,
    surfdock_dir='/home/aoxu/projects/SurfDock',       # NOT the submodule
    precomputed_arrays='/home/aoxu/projects/precomputed/precomputed_arrays',
    cuda_device_index=0,
)
```

## Output

- `${output_dir}/surfdock/<ligand-stem>/rank{N}.sdf` — one SDF per pose.
- Ranked by confidence (highest = rank1).
- Intermediate files (surface .ply, ESM .pt embeddings, CSV) end up in a
  tempdir by default; pass a persistent `surfdock_work_dir` to keep them.

## Knobs that matter

| kwarg | default | typical override |
|-------|---------|-------------------|
| `samples_per_complex` | 40 | 80 for harder pockets |
| `num_poses` | 10 | 40 to keep everything |
| `batch_size` | 40 | reduce if VRAM-limited |
| `mdn_dist_threshold` | 3 | filter threshold for confidence model |
| `ligand_to_pocket_center` | **false** | set **true** for OOD ligands (see Pitfalls) |

## Pitfalls

### Pose lands ~30 Å from the pocket

`dock_engine`'s default mode does not anchor the diffusion initial position
to the pocket. For OOD ligands (large, charged, halogenated, etc.) the
diffusion can drift thousands of Å. Fix: drive the 4 steps explicitly with a
persistent work dir, inject a `pocket_center` column into the input CSV, and
pass `--ligand_to_pocket_center` to `inference_accelerate.py`.

Worked example: `contrasCF/analysis/scripts/13_run_surfdock.py` — in
particular `_run_inference_with_pocket_center` and the CSV-injection block.
Approximate skeleton:

```python
# After _build_input_csv writes input.csv:
import csv as csv_mod
cx, cy, cz = pocket_center  # from box.json or upstream pocket detector
with open(csv_path) as fh:
    rows = list(csv_mod.DictReader(fh))
if rows and "pocket_center" not in rows[0]:
    for r in rows:
        r["pocket_center"] = f"{cx},{cy},{cz}"
    with open(csv_path, "w", newline="") as fh:
        writer = csv_mod.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
# Then append "--ligand_to_pocket_center" to the accelerate launch cmd.
```

Also translate the input ligand into the pocket BEFORE step 1 (surface
computation), because SurfDock crops an 8 Å pocket around the ligand's
centroid. If the ligand sits at the origin (RDKit default), the pocket crop
is outside the receptor. See `_translate_ligand_to_pocket` in the same file.

### Tempdir loses original filenames

SurfDock encodes scores in its output filenames (e.g.,
`*_rank_1_rmsd_0.5_confidence_3.2.sdf`). `dock_engine`'s tempdir mode
discards these and just keeps `rank{N}.sdf`. For score-aware analysis, drive
the 4 steps yourself with a persistent work dir.

### Local install vs submodule

`surfdock_dir` MUST point at `/home/aoxu/projects/SurfDock`, not
`${COGLIGANDBENCH}/forks/SurfDock`. The submodule lacks the 231 MB
MSMS/APBS/pdb2pqr extract and the editable ESM. Errors look like
"FileNotFoundError: msms" or "No module named esm" if you point at the
submodule.

## Source of truth

`${COGLIGANDBENCH}/CLAUDE.md` → "SurfDock" section.
`${COGLIGANDBENCH}/dockstrat_config/model/surfdock_inference.yaml`.
`contrasCF/analysis/scripts/13_run_surfdock.py` (most complete worked example).
```

- [ ] **Step 2: Verify**

```bash
wc -l ~/.claude/skills/dockstrat/references/surfdock.md
```

Expected: 100-160 lines.

---

## Task 10: Write `references/boltz.md`

**Files:**
- Create: `~/.claude/skills/dockstrat/references/boltz.md`

- [ ] **Step 1: Write the file**

```markdown
# Boltz-1 and Boltz-2

Co-folding model: predicts the protein–ligand complex from sequence + SMILES.
Boltz-2 adds affinity prediction. Both share the `boltz` env and CLI.
Single-pair runs are ~3-10 min.

## Env

- Conda env: `${COGLIGANDBENCH}/envs/boltz` → `/mnt/katritch_lab2/aoxu/envs/boltz`.
- Binary: `${COGLIGANDBENCH}/envs/boltz/bin/boltz` (Boltz-2 v2.2.1+).
- Boltz-1 is also available at `/home/aoxu/miniconda3/envs/rdkit_env/bin/boltz`
  (v0.4.1). Prefer the dedicated env unless you specifically need v0.4.1.
- GPU required (≥24 GB VRAM recommended).
- Install: `bash scripts/install_boltz_env.sh`.

## dock_engine signature

```python
from dockstrat import dock_engine

# Boltz-2 (with affinity)
dock_engine(
    'boltz2',
    protein='receptor.pdb',
    ligand='ligand.sdf',
    output_dir='./out',
    diffusion_samples=5,
    recycling_steps=3,
    sampling_steps=200,
    seed=42,
    num_poses_to_keep=5,
    cuda_device_index=0,
    sampling_steps_affinity=200,
    diffusion_samples_affinity=5,
    affinity_mw_correction=False,
)

# Boltz-1 (pose only)
dock_engine(
    'boltz1',
    protein='receptor.pdb',
    ligand='ligand.sdf',
    output_dir='./out',
    diffusion_samples=5,
    recycling_steps=3,
    sampling_steps=200,
    seed=42,
    num_poses_to_keep=5,
    cuda_device_index=0,
)
```

## Output

- `${output_dir}/boltz{1,2}/<system>/rank{N}.sdf` — top-N poses,
  bond-order-recovered against input SMILES.
- Ranked by `confidence_score` from
  `confidence_{stem}_model_{N}.json` (highest = rank1).
- Boltz-2 attaches affinity properties to each SDF:
  - `affinity_pred_value` (log10 IC50-like)
  - `affinity_probability_binary` (P(binder))
- Intermediate Boltz outputs live in `${output_dir}/boltz_results_<prefix>/`.

## Knobs that matter

| kwarg | default | typical override |
|-------|---------|-------------------|
| `diffusion_samples` | 5 | 10 for higher pose diversity |
| `recycling_steps` | 3 | 6 for more confident folds |
| `sampling_steps` | 200 | 100 for fast sweeps |
| `num_poses_to_keep` | 5 | 1 if you only need the top pose |
| `cuda_device_index` | 0 | per-GPU on multi-GPU hosts |
| `use_msa_server` | False | True to fetch an MSA via Boltz's server |

## Pitfalls

### Single-sequence by default; MSA via Boltz server

Default mode is single-sequence (`msa: "empty"`). For better accuracy on
novel targets, pass `use_msa_server=True` — Boltz piggybacks on its own MSA
endpoint. This requires network access from the calling host.

ColabFold's direct MSA API errored out during 2026-Q1 testing; Boltz's
server is the reliable path. (See contrasCF `casf_mutagenesis_module.md`
memory entry.)

### Affinity is Boltz-2 only

`affinity_mw_correction`, `sampling_steps_affinity`, etc. are silently
ignored for `boltz1`. Don't be surprised if they have no effect there.

### Chain IDs

Chain IDs are assigned A, B, ... for protein chains, next letter for the
ligand. If you're constructing input YAMLs manually, the binder must
reference the ligand's chain letter exactly.

## Source of truth

`${COGLIGANDBENCH}/CLAUDE.md` → "Boltz-1 / Boltz-2" section.
`${COGLIGANDBENCH}/dockstrat_config/model/boltz{1,2}_inference.yaml`.
```

- [ ] **Step 2: Verify**

```bash
wc -l ~/.claude/skills/dockstrat/references/boltz.md
```

Expected: 80-150 lines.

---

## Task 11: Write `references/af3.md`

**Files:**
- Create: `~/.claude/skills/dockstrat/references/af3.md`

- [ ] **Step 1: Write the file**

```markdown
# AlphaFold3

Co-folding model. Single-sequence mode (no MSA databases). Single-pair runs
are ~5-15 min.

## Env

- Conda env: `${COGLIGANDBENCH}/envs/alphafold3` → `/mnt/katritch_lab2/aoxu/envs/alphafold3`.
- AF3 source: `${COGLIGANDBENCH}/forks/alphafold3/alphafold3/`.
- **Decompressed weights**: `${COGLIGANDBENCH}/forks/alphafold3/models/af3.bin`
  (decompressed from `af3.bin.zst` at install time).
- GPU required.
- Install: `bash scripts/install_alphafold3_env.sh`.

## dock_engine signature

```python
from dockstrat import dock_engine

dock_engine(
    'alphafold3',
    protein='receptor.pdb',
    ligand='ligand.sdf',
    output_dir='./out',
    num_samples=5,           # AF3 diffusion samples
    num_seeds=1,             # AF3 random seeds
    num_recycles=10,
    num_poses_to_keep=5,
    cuda_device_index=0,
    timeout_seconds=3600,
)
```

## Output

- `${output_dir}/alphafold3/<prefix_lowercase>/` — full AF3 output (mmCIFs,
  confidence JSONs, `ranking_scores.csv`).
- `${output_dir}/alphafold3/<prefix>/rank{N}.sdf` — top-N ligand poses
  extracted from the mmCIFs and bond-order-recovered against input SMILES.
- Ranked by AF3's `ranking_score` (highest = rank1).

## Knobs that matter

| kwarg | default | typical override |
|-------|---------|-------------------|
| `num_samples` | 5 | 1 for fast smoke tests |
| `num_seeds` | 1 | rarely changed |
| `num_recycles` | 10 | 3 for fast sweeps (less accurate) |
| `num_poses_to_keep` | 5 | 1 for top pose only |
| `cuda_device_index` | 0 | per-GPU |

## Pitfalls

- **Single-sequence mode only.** Uses `--norun_data_pipeline`; no MSA
  databases, no `db_dir`. If you need MSA-driven AF3, use the canonical AF3
  CLI directly (and source MSAs from elsewhere).
- **CLI flag drift.** `num_samples` is mapped internally to AF3's
  `--num_diffusion_samples`. Exact AF3 CLI flag spellings drift across
  releases — verify with `run_alphafold.py --help` after a fresh install.
- **Weights must be decompressed.** If you see `FileNotFoundError: af3.bin`,
  re-run the install script (which decompresses `af3.bin.zst`).

## Source of truth

`${COGLIGANDBENCH}/CLAUDE.md` → "AlphaFold3" section.
`${COGLIGANDBENCH}/dockstrat_config/model/alphafold3_inference.yaml`.
```

- [ ] **Step 2: Verify**

```bash
wc -l ~/.claude/skills/dockstrat/references/af3.md
```

Expected: 50-100 lines.

---

## Task 12: Write `references/vina.md`

**Files:**
- Create: `~/.claude/skills/dockstrat/references/vina.md`

- [ ] **Step 1: Write the file**

```markdown
# AutoDock Vina

CPU baseline. Slower than UniDock2/GNINA but useful for smoke tests and as a
ground-truth physics reference. Wall-clock: 1-5 min per pair on a typical CPU.

## Env

- Conda env: `${COGLIGANDBENCH}/envs/vina` → `/mnt/katritch_lab2/aoxu/envs/vina`.
- Requires `obabel` and `vina` on `$PATH` (both shipped by the env).
- CPU-only.
- Install: `bash scripts/install_vina_env.sh`.

## dock_engine signature

```python
from dockstrat import dock_engine

dock_engine(
    'vina',
    protein='receptor.pdb',
    ligand='ligand.sdf',
    output_dir='./out',
    exhaustiveness=8,
    num_modes=9,
    seed=42,
    box_size=[20, 20, 20],
)
```

## Output

- `${output_dir}/vina/<ligand-stem>/{prefix}_pose{rank}_score{score:.2f}.sdf` —
  one SDF per pose, ranked by Vina score (lowest = best).
- Intermediate files: `protein.pdbqt`, `ligand.pdbqt`, `vina_out.pdbqt`.

## Knobs that matter

| kwarg | default | typical override |
|-------|---------|-------------------|
| `exhaustiveness` | 8 | 16-32 for higher accuracy |
| `num_modes` | 9 | 20 for diverse poses |
| `box_size` | computed from ligand+protein | override for known pockets |
| `seed` | random | set for reproducibility |

## Pitfalls

- **Input preprocessing.** Protein PDB → PDBQT via `obabel -i pdb -o pdbqt -xhk -r`;
  ligand SDF → PDBQT via `obabel -i sdf -o pdbqt -xh`. `dock_engine` does this
  internally; you don't need to prep PDBQT yourself.
- **Box detection.** If `box_size` is not provided, the box is computed from
  the ligand centroid + 50 closest protein atoms (via MDAnalysis). Minimum
  20 Å, 10 Å margin. If the input ligand is at the origin (RDKit-embedded),
  the box will be wrong — translate the ligand into the pocket first, or
  pass an explicit `box_center` + `box_size`.
- **MDAnalysis optional but recommended.** Without it the box detection
  falls back to a cruder algorithm.

## Source of truth

`${COGLIGANDBENCH}/CLAUDE.md` → "Vina" section.
`${COGLIGANDBENCH}/dockstrat_config/model/vina_inference.yaml`.
```

- [ ] **Step 2: Verify**

```bash
wc -l ~/.claude/skills/dockstrat/references/vina.md
```

Expected: 50-90 lines.

---

## Task 13: Write `references/other-methods.md`

**Files:**
- Create: `~/.claude/skills/dockstrat/references/other-methods.md`

- [ ] **Step 1: Write the file**

```markdown
# Other dockStrat methods (quick reference)

Methods covered briefly. For each, the source of truth is the dockStrat
repo's CLAUDE.md section + the per-method YAML config in
`${COGLIGANDBENCH}/dockstrat_config/model/`.

## Chai-1

**Env**: `${COGLIGANDBENCH}/envs/chai`. Requires `chai_lab` package.
**GPU**: yes, ≥24 GB VRAM (proteins > 500 residues may OOM on a 24 GB card).
**Output**: `${output_dir}/chai/<system>/rank{N}.sdf`, ranked by `aggregate_score` (highest = rank1).
**Use when**: co-folding baseline; not your first pick (Boltz-2 typically scores better).
**Pitfall**: requires `python_exec_path` kwarg pointing at the chai-lab env's
Python (e.g., `/path/to/forks/chai-lab/chai-lab/bin/python`), unlike most
other methods.

```python
dock_engine('chai', protein=..., ligand=..., output_dir=...,
            python_exec_path='/path/to/forks/chai-lab/chai-lab/bin/python')
```

Source: `${COGLIGANDBENCH}/CLAUDE.md` → "Chai-1".

## DynamicBind

**Env**: `${COGLIGANDBENCH}/envs/dynamicbind`. Wraps `${COGLIGANDBENCH}/forks/DynamicBind/`.
**GPU**: yes.
**Output**: `${output_dir}/dynamicbind/<system>/rank{N}_ligand_lddt{lDDT}_affinity{aff}.sdf`.
Up to 40 samples (configurable via `samples_per_complex`).
**Use when**: diffusion-based docking with affinity prediction (similar slot to Boltz-2 but different model family).
**Pitfall**: `--no_relax` and `--paper` flags are hardcoded — paper weights only, no relaxation.

Source: `${COGLIGANDBENCH}/CLAUDE.md` → "DynamicBind".

## Protenix

**Env**: `${COGLIGANDBENCH}/envs/protenix`. Requires Python 3.11+. Auto-downloads models to `~/.protenix` on first run.
**GPU**: yes.
**Output**: `${output_dir}/protenix/<prefix>/rank{N}.sdf`, ranked by `ranking_score`.
**Use when**: alternative co-folding model, paper-comparable to AF3/Boltz.
**Pitfall**: input JSON uses `proteinChain` (not `protein`) and `ligand.ligand` for SMILES.
CIF chain IDs are `{letter}0` format (A0, B0, …); ligand discovery uses HETATM scan with bond-order recovery fallback.

Source: `${COGLIGANDBENCH}/CLAUDE.md` → "Protenix".

## ICM

**Env**: none — commercial ICM binary on `$PATH`.
**GPU**: no.
**Use when**: ICM-licensed users want physics-based docking with ICM's MC sampler.
**Pitfall**: NOT exposed via `dock_engine`. Runs via standalone scripts in `${COGLIGANDBENCH}/forks/ICM/`.

Source: `${COGLIGANDBENCH}/CLAUDE.md` → "ICM" section + `forks/ICM/README`.

## ICM-RTCNN

Hybrid: ICM sampling + RTCNN rescoring. Same constraints as ICM (commercial binary, no `dock_engine` entry).

Source: `${COGLIGANDBENCH}/CLAUDE.md` → "ICM-RTCNN".
```

- [ ] **Step 2: Verify**

```bash
wc -l ~/.claude/skills/dockstrat/references/other-methods.md
```

Expected: 50-90 lines.

---

## Task 14: Cross-link audit

**Files:**
- All of `~/.claude/skills/dockstrat/SKILL.md` and `references/*.md`.

- [ ] **Step 1: Check all `references/<name>.md` links in SKILL.md resolve to existing files**

```bash
cd ~/.claude/skills/dockstrat
for ref in $(grep -oE 'references/[a-z0-9_-]+\.md' SKILL.md | sort -u); do
    if [ -f "$ref" ]; then
        echo "OK: $ref"
    else
        echo "MISSING: $ref"
    fi
done
```

Expected: 7 `OK:` lines, 0 `MISSING:`.

- [ ] **Step 2: Check method names referenced in SKILL.md match the index table**

```bash
grep -E "^\| \`(vina|gnina|unidock2|surfdock|alphafold3|boltz1|boltz2|chai|dynamicbind|protenix|icm)" ~/.claude/skills/dockstrat/SKILL.md | wc -l
```

Expected: 8 (one row each for unidock2, gnina, surfdock, alphafold3, boltz1, boltz2, vina, plus the `chai/dynamicbind/protenix/icm/icm-rtcnn` combined row).

- [ ] **Step 3: Total skill size sanity check**

```bash
wc -l ~/.claude/skills/dockstrat/SKILL.md ~/.claude/skills/dockstrat/references/*.md
```

Expected: SKILL.md 200-300 lines; each reference 50-160 lines; total < 1500 lines.

---

## Task 15: Smoke test — single-prompt verification

**Files:**
- No files written; verification only.

- [ ] **Step 1: Run a fresh subagent prompt that should fire the skill**

Dispatch via `Agent` tool with `subagent_type: general-purpose`:

```
Prompt: "Write a Python snippet that uses dock_engine to run UniDock2 on
/tmp/receptor.pdb and /tmp/ligand.sdf, outputting to /tmp/out, with 20 poses
and a 20 Å box. Tell me which conda env to call it from."

Expected behavior:
- Subagent invokes the dockstrat skill.
- Reads references/unidock2.md.
- Produces a snippet using dock_engine('unidock2', ...) with num_poses=20 and
  box_size=[20, 20, 20] (or [20.0, 20.0, 20.0]).
- Mentions the `dockstrat` (or `unidock2`) env.
```

Note for the executor: if the subagent does not auto-fire the skill, the
frontmatter description likely needs more trigger keywords. Inspect the
subagent's output, then revise SKILL.md frontmatter accordingly and re-run.

- [ ] **Step 2: Run a debugging-flavor prompt**

```
Prompt: "My SurfDock run is placing rank1 about 30 Å from the binding pocket.
What's going wrong and how do I fix it?"

Expected behavior:
- Subagent invokes the dockstrat skill.
- Reads references/surfdock.md (or recognizes the pocket-center anchoring
  issue from SKILL.md gotcha 6d).
- Mentions --ligand_to_pocket_center, the pocket_center CSV column, and
  /translate_ligand_to_pocket / 13_run_surfdock.py.
```

- [ ] **Step 3: Record test outcomes**

If both prompts produced the expected behavior, the skill is shippable.
Capture the subagent transcripts (or paraphrased summaries) in a follow-up
note.

If either prompt failed, file the gap (which trigger word was missing? which
section needs to be more prominent?) and iterate on SKILL.md before
declaring done.

---

## Task 16: Commit the plan + smoke-test notes in the contrasCF repo

**Files:**
- Modify: This plan file (`docs/superpowers/plans/2026-05-27-dockstrat-skill.md`) — append a "Smoke test outcomes" section.

- [ ] **Step 1: Append smoke-test outcomes to the plan**

After Task 15, append a section like:

```markdown
## Smoke test outcomes (filled in at end of execution)

- Prompt 1 (UniDock2 quickstart): PASS/FAIL + one-sentence note.
- Prompt 2 (SurfDock debug): PASS/FAIL + one-sentence note.
- Frontmatter adjustments made: [list any].
```

- [ ] **Step 2: Commit the plan + notes in contrasCF**

```bash
cd /mnt/katritch_lab2/aoxu/contrasCF
git add docs/superpowers/plans/2026-05-27-dockstrat-skill.md
git commit -m "$(cat <<'EOF'
docs(plans): dockstrat skill implementation plan + smoke test outcomes

Plan to create ~/.claude/skills/dockstrat/ (SKILL.md + 7 reference files).
Smoke-test outcomes recorded inline. The skill files themselves live outside
this repo (~/.claude is not a git repo).
EOF
)"
```

- [ ] **Step 3: Verify**

```bash
git -C /mnt/katritch_lab2/aoxu/contrasCF log --oneline -1
```

Expected: top commit is the plan-with-outcomes commit.

---

## Success criteria (from spec)

| # | Criterion | Verified by |
|---|-----------|-------------|
| 1 | `~/.claude/skills/dockstrat/SKILL.md` exists, has the frontmatter, is < 300 lines | Task 6 step 2 + Task 14 step 3 |
| 2 | All 7 `references/<method>.md` exist, none > 200 lines | Tasks 7-13 + Task 14 step 3 |
| 3 | Fresh agent given "dock_engine UniDock2 on receptor/ligand" produces a working snippet | Task 15 step 1 |
| 4 | Fresh agent given "SurfDock pose is 30 Å away" surfaces gotcha 6d | Task 15 step 2 |
| 5 | Skill does not duplicate input-prep / RMSD content | Tasks 2-6 by construction (the "Out of scope" block in Task 2 and the absence of those sections enforce this) |
