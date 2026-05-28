# dockstrat skill — design

**Date:** 2026-05-27
**Author:** Aaron Xu (with Claude)
**Status:** design approved, pending implementation plan

## Goal

Create a user-global Claude Code skill that documents the dockStrat
(`/mnt/katritch_lab2/aoxu/CogLigandBench`) docking framework — primarily the
`dock_engine()` Python API and the per-method conda envs / output formats /
known gotchas.

Future Claude sessions should be able to fire this skill on any mention of
docking via dockStrat (or of a supported method like UniDock2, GNINA,
SurfDock, AlphaFold3, Boltz, etc.) and produce correct, working invocations
without rediscovering paths, env names, kwarg names, or pitfalls.

## Non-goals

- Input preparation (PDB cleaning, SMILES → SDF, binding-box detection). The
  skill points readers at `dockstrat/data/` and `contrasCF/docking/smoketest/`
  but does not duplicate that content.
- Downstream analysis (RMSD vs crystal, scoring tables). The skill points at
  `dockstrat/analysis/`.
- Method installation. The skill points at `scripts/install_<method>_env.sh`
  inside the dockStrat repo.
- Documentation of the analysis paper itself, dataset construction, or the
  contrasCF-specific 16-case workflow. That belongs in CLAUDE.md / project
  notes, not a general docking skill.

## Audience and use cases

Audience: Claude Code sessions on Aaron's machine (and any other workstation
with the same dockStrat checkout and conda envs).

Typical triggers:
- "Run UniDock2 on this receptor + ligand."
- "Dock these N ligands with GNINA."
- "Use dock_engine to run Boltz-2."
- "Why is SurfDock placing the pose ~30 Å from the pocket?"
- "Set up a script that loops over cells and runs dock_engine('unidock2', ...)."

## Design choices (recap of brainstorming)

| Decision | Choice | Why |
|----------|--------|-----|
| Primary consumer | Any project that uses dockStrat | Generic, not tied to contrasCF directory layout. |
| Scope | API reference + opinionated quickstart | Skips heavy input-prep / analysis content; focuses on "I have inputs, give me poses." |
| Install location | `~/.claude/skills/dockstrat/` | User-global; available in every session on this machine. |
| Method coverage | 6 deep + 5 brief | Deep: unidock2, gnina, surfdock, alphafold3, boltz1, boltz2, vina (as baseline). Brief: chai, dynamicbind, protenix, icm, icm-rtcnn. |
| Call style | `dock_engine()` Python API first | Matches the unified API in `dockstrat/engine.py`. CLI shown only when the API misses something. |
| Layered structure | `SKILL.md` + lazy-loaded `references/<method>.md` | Keeps initial skill load light (~250 lines); per-method detail loaded only when a specific method is requested. |
| Gotchas included | All 4 selected | UniDock2 GPU desync, SurfDock pocket-center, rdkit_env import bypass, LD_LIBRARY_PATH prepend. |
| Auto-fire on bare method names | Yes | Trigger on any of "dock_engine", "dockStrat", "CogLigandBench", "UniDock2", "GNINA", "SurfDock", "AlphaFold3", "Boltz", "Protenix", "Chai", "DynamicBind". User accepts that this will occasionally over-fire in pure discussion contexts. |
| Restate env conventions inline | Yes | Skill becomes self-contained — readers don't have to chase `env_conventions.md` or `tool_locations.md` memory entries. Slight duplication accepted. |

## Layout

```
~/.claude/skills/dockstrat/
  SKILL.md
  references/
    unidock2.md
    gnina.md
    surfdock.md
    boltz.md           # boltz1 + boltz2 together (same env, same CLI)
    af3.md
    vina.md
    other-methods.md   # chai, dynamicbind, protenix, icm, icm-rtcnn
```

## SKILL.md

### Frontmatter

```yaml
---
name: dockstrat
description: Use when running protein-ligand docking via the dockStrat
  framework (/mnt/katritch_lab2/aoxu/CogLigandBench). Covers the dock_engine()
  Python API for UniDock2, GNINA, SurfDock, AlphaFold3, Boltz-1/2, and Vina,
  including per-method conda envs, key kwargs, output formats, and known
  gotchas (GPU desync, SurfDock pocket-center anchoring, rdkit_env import
  bypass, LD_LIBRARY_PATH). Trigger on mentions of dock_engine, dockStrat,
  CogLigandBench, or any supported method (UniDock2, GNINA, SurfDock,
  AlphaFold3, Boltz, Protenix, Chai, DynamicBind), or any task that needs to
  produce ranked poses for a protein+ligand pair.
---
```

### Sections

1. **When to use this skill**
   - One paragraph listing trigger contexts (running docking, picking a
     method, debugging a method, writing a multi-cell loop).
   - One paragraph on what's out of scope (input prep, RMSD, installs) with
     pointers.

2. **Project layout you can assume**
   - `${COGLIGANDBENCH}=/mnt/katritch_lab2/aoxu/CogLigandBench`
   - Envs symlinked under `${COGLIGANDBENCH}/envs/{method}/`
   - Configs under `${COGLIGANDBENCH}/dockstrat_config/model/`
   - Forks at `${COGLIGANDBENCH}/forks/<method>/`
   - `dock_engine` is `from dockstrat import dock_engine` (requires
     `${COGLIGANDBENCH}` on `PYTHONPATH` or `pip install -e .` in the
     `dockstrat` env).

3. **Method-selection decision tree** (10-line decision block)
   - Physics, fast, GPU available → `unidock2`
   - Physics + CNN rescoring → `gnina`
   - Want surface / pocket features, OOD ligand → `surfdock`
   - Want co-folded structure + ligand pose → `boltz2` (preferred over
     `boltz1`/`af3` for affinity; `af3` if you need AF3 specifically)
   - CPU-only, smoke test, no GPU → `vina`
   - Other (chai/dynamicbind/protenix/icm) → see `references/other-methods.md`

4. **Env conventions** (restated for self-containedness)
   - LD_LIBRARY_PATH prepend pattern (with example).
   - "Call Python via absolute path, not `conda run -n`" rule.
   - Which env to call `dock_engine` from (the `dockstrat` env at
     `/home/aoxu/miniconda3/envs/dockstrat`, if present, or the env the user
     used to `pip install -e .`).
   - Reference to the existing memory entries for users who want to update
     them.

5. **Quickstart by family** — three Python snippets:
   - Physics: `dock_engine('unidock2', protein='r.pdb', ligand='l.sdf', output_dir='./out', num_poses=20)`
   - Hybrid: `dock_engine('gnina', protein='r.pdb', ligand='l.sdf', output_dir='./out')`
   - Co-folding: `dock_engine('boltz2', protein='r.pdb', ligand='l.sdf', output_dir='./out')`
   Each snippet annotated with: where `rank{N}.sdf` lands, how to read the top
   score, whether GPU is required, expected wall-clock.

6. **Gotchas** (the four selected)
   - **6a. LD_LIBRARY_PATH prepend** — `export LD_LIBRARY_PATH=$ENV/lib:$LD_LIBRARY_PATH`
     before Python that imports rdkit/gemmi. With example failure mode.
   - **6b. rdkit_env → dockstrat import bypass** — full snippet:
     `importlib.util.spec_from_file_location` + stub `sys.modules['dockstrat']`,
     lifted from `_load_surfdock_module` in `13_run_surfdock.py`. Use when
     calling dockstrat helpers from an env that doesn't have `rootutils`.
   - **6c. UniDock2 GPU-desync guard** — set per-cell `timeout=120s`; abort
     loop after 5 consecutive errors that returned in `<5s` (CUDA-init
     failure). Symptom: `nvidia-smi` hangs / NVML mismatch. Recovery: usually
     requires a GPU reset.
   - **6d. SurfDock pocket-center anchoring** — `dock_engine('surfdock', ...)`
     uses a tempdir that loses score-encoded filenames AND does not anchor
     the diffusion initial position. For OOD ligands, drive the 4-step
     pipeline manually (see `references/surfdock.md`) and pass
     `--ligand_to_pocket_center` + inject a `pocket_center` column into the
     CSV. Worked example: `contrasCF/analysis/scripts/13_run_surfdock.py`.

7. **Per-method reference index** (table)

   | Method | Env | Output filename | "Use when" | Detail |
   |--------|-----|----------------|------------|--------|
   | unidock2 | `envs/unidock2` | `rank{N}.sdf` | Physics, fast GPU | `references/unidock2.md` |
   | gnina | `envs/gnina` | `rank{N}.sdf` | Physics + CNN | `references/gnina.md` |
   | surfdock | `SurfDock` | `rank{N}.sdf` | Surface/pocket features | `references/surfdock.md` |
   | alphafold3 | `envs/alphafold3` | `rank{N}.sdf` | Cofold (AF3 specifically) | `references/af3.md` |
   | boltz1/2 | `envs/boltz` | `rank{N}.sdf` | Cofold (Boltz, +affinity) | `references/boltz.md` |
   | vina | `envs/vina` | `rank{N}.sdf` | CPU baseline | `references/vina.md` |
   | chai/dynamicbind/protenix/icm/icm-rtcnn | varies | varies | See full ref | `references/other-methods.md` |

   Followed by an explicit instruction to Claude: "Before generating
   method-specific code, Read the relevant `references/<method>.md` file."

## Per-method reference files

Each `references/<method>.md` follows this template (~80–150 lines):

```markdown
# <Method>

## Env
- Conda env: `envs/<method>` (symlink to `/mnt/katritch_lab2/aoxu/envs/<method>`)
- Binary or entrypoint: `...`
- GPU required: yes/no, VRAM minimum
- Install script: `scripts/install_<method>_env.sh`

## dock_engine signature

```python
dock_engine('<method>',
    protein='receptor.pdb',
    ligand='ligand.sdf',
    output_dir='./out',
    # method-specific kwargs:
    ...
)
```

## Output
- Filename pattern (e.g., `rank{N}.sdf`)
- Ranking metric (e.g., "lowest Vina energy", "highest CNNscore")
- SDF tags / properties added (e.g., `ranking_score`, `confidence_score`)
- Intermediate files (kept for post-mortem)

## Knobs that matter
Table of the 4-6 kwargs you actually change, with sensible defaults.

## Pitfalls
Method-specific (e.g., GNINA hardcoded binary path; Boltz needs MSA;
AF3 needs decompressed weights; UniDock2 box size must come from YAML
config, not CLI).
```

### Specific reference content

- **unidock2.md** — box-size YAML quirk; centroid-from-ligand pattern;
  search_mode (`detail`/`balance`/`fast`); GPU desync guard cross-link to
  SKILL.md §6c.
- **gnina.md** — hardcoded binary at
  `${COGLIGANDBENCH}/forks/GNINA/gnina`; `--cnn_scoring rescore` default;
  `--autobox_ligand` reads PDB+SDF natively; GPU optional but slow on CPU.
- **surfdock.md** — 4-step pipeline (surface → CSV → ESM → diffusion); the
  `precomputed_arrays` env var; the local-install-vs-submodule distinction;
  worked example pointer to `13_run_surfdock.py`; the pocket-center
  anchoring fix.
- **boltz.md** — covers both boltz1 and boltz2; affinity-prediction block
  for boltz2; MSA-via-Boltz pattern (single-sequence by default; pass
  `use_msa_server=True` for MSA via Boltz's own server).
- **af3.md** — decompressed `af3.bin` weights at `forks/alphafold3/models/`;
  `num_recycles`, `num_samples`, `num_seeds`; single-sequence mode via
  `--norun_data_pipeline`.
- **vina.md** — obabel preprocessing (PDB→PDBQT, SDF→PDBQT); MDAnalysis
  pocket detection; CPU-only.
- **other-methods.md** — one short paragraph each for chai, dynamicbind,
  protenix, icm, icm-rtcnn. Each paragraph: env, primary use case, one-line
  pitfall, link to dockStrat CLAUDE.md sections for full details.

## Open questions / risks

- **R1. Skill over-fire.** With method names in the trigger, the skill will
  fire even in pure discussion contexts (e.g., "GNINA's CNN is interesting").
  Acceptable per user decision; can tighten later if it becomes noisy.
- **R2. Path coupling.** The skill hardcodes `/mnt/katritch_lab2/aoxu/...`
  paths. If the dockStrat repo moves, the skill needs an edit. Mitigation:
  the SKILL.md opens with a single `COGLIGANDBENCH` "fact" block; everything
  else uses `${COGLIGANDBENCH}`. One edit point.
- **R3. Drift vs. dockStrat CLAUDE.md.** The dockStrat repo's own CLAUDE.md
  is the source of truth for method details; this skill summarizes. Risk
  of drift when new methods are added or kwargs change. Mitigation: each
  `references/<method>.md` ends with a "Source of truth" pointer to the
  dockStrat CLAUDE.md section, so Claude can verify if a kwarg looks
  stale. We do not aim for full faithful mirroring.
- **R4. Skill testing.** Per superpowers:writing-skills, ideally a new skill
  should be validated via TDD-style subagent pressure tests (does an agent
  without the skill misroute UniDock2; does the agent with the skill route
  it correctly). This adds time. The implementation plan will decide
  whether to do this rigorously or rely on smoke-testing the skill on a
  single real prompt ("dock the smoketest case with unidock2") before
  declaring it done.

## Success criteria

1. `~/.claude/skills/dockstrat/SKILL.md` exists, has the frontmatter above,
   and is < 300 lines.
2. All 7 `references/<method>.md` files exist; each follows the template;
   none exceeds 200 lines.
3. A fresh Claude Code session (or subagent), given the prompt "use
   dock_engine to run UniDock2 on receptor.pdb + ligand.sdf into ./out",
   invokes this skill and produces a working snippet without re-reading the
   dockStrat repo.
4. The same session, asked "why is SurfDock placing my pose 30 Å away?",
   surfaces gotcha 6d and points at `13_run_surfdock.py`.
5. The skill does not duplicate input-prep or RMSD content; readers asking
   about those get pointed at the right places.

## Next step

Hand off to `superpowers:writing-plans` to create the implementation plan
(ordering of files, content sourcing per file, smoke-test prompts).
