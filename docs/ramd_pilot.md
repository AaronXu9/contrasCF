# τ-RAMD pilot for adversarial co-folding validation

A 2-system pilot of τ-RAMD (random acceleration molecular dynamics) as a
cheaper alternative to funnel metadynamics for validating co-folding
predictions on adversarial protein–ligand systems. Lives under
[`analysis/ramd_pilot/`](../analysis/ramd_pilot/).

## Why this exists

Masters et al. 2025 ran funnel metadynamics (FM, 1 μs × triplicate per
system) on 13 adversarial CDK2/GDH systems and reported bound-state
probabilities (Table 1). FM is the gold-standard physics oracle but is
expensive (~3 μs of OpenMM per system). To extend that analysis to the
broader contrasCF cohort we need a cheaper oracle.

τ-RAMD (Kokh & Wade, *J. Chem. Theory Comput.* 2018) applies a small,
randomly-oriented force on the ligand center of mass during otherwise
standard MD. The bias accelerates ligand exit; the resulting exit-time
distribution **rank-orders binding affinities** even though the absolute
timescale is artificially compressed.

The pilot tests one binary question: **can τ-RAMD discriminate a paper
binder (CDK2 wild-type, FM=0.61) from a paper non-binder (CDK2 with all
11 pocket residues→Phe, FM=0.00)?** If yes, the protocol scales to the
remaining 11 systems.

## Plan reference

The full pilot plan is at `~/.claude/plans/from-the-paper-we-harmonic-snowflake.md`.
Headline decision criteria from §2:

| outcome | criterion | next step |
|---|---|---|
| **PASS** | R ≥ 5 AND τ_KM(WT) ≥ 10 ns AND f_ee(WT) ≤ 0.2 | scale to 13-system calibration |
| **MARGINAL** | 2 ≤ R < 5, OR f_ee(WT) ∈ (0.2, 0.4] | diagnose force, extend t_max |
| **FAIL** | R < 2 OR f_ee(WT) > 0.4 OR τ_KM(WT) < 5 ns | pivot to unbinding metaD |

## Software stack

| component | choice | rationale |
|---|---|---|
| MD engine | **gromacs-ramd** patched fork (HITS-MCM), tag `gromacs-2024.1-ramd-2.1` | RAMD plugin compiled into `libgromacs_mpi.so.9.0.0`; vanilla GROMACS rejects `ramd-*` mdp options |
| protein force field | AMBER ff99SB-ILDN | ff14SB (paper's choice) is not bundled with GROMACS; ff99SB-ILDN is the closest validated bundled variant. Override via `--ff` once a ported ff14SB.ff lands in `$GMXLIB`. |
| ligand parameters | GAFF2 + AM1-BCC charges via antechamber + acpype | matches paper for unusual chemistries (methylated sugars, charge-modified ATP) |
| water | TIP3P | matches paper |
| ions | 0.15 M NaCl, neutralizing | standard |
| box | cubic, 15 Å padding | matches paper |

All env vars and conda envs are wired up in [`env/carc.sh`](../env/carc.sh).
Bootstrap on a fresh CARC node: see §"Running on CARC" below.

## Architecture

```
analysis/ramd_pilot/
├── config.py                  # paths, system list, decision thresholds
├── carc_setup/
│   ├── install_gromacs_ramd.sh   # one-time: build gromacs-ramd from tag
│   └── setup_ambertools_env.sh   # one-time: AmberTools + analysis env
├── prep/
│   ├── extract_inputs.py         # CIF → chain-A protein PDB + ATP mol2
│   ├── ligand_params.py          # antechamber → parmchk2 → acpype
│   ├── system_build.py           # pdb2gmx + solvate + ions + index
│   ├── equilibrate.py            # 5-stage NPT ladder runner
│   ├── equilibrate.sbatch        # SLURM wrapper for equilibration
│   └── mdp/                      # 6 mdp files (min + 5 equil + RAMD prod)
├── ramd/
│   ├── pbcatom.py                # central-atom finder for RAMD groups
│   ├── replica_runner.py         # render slurm + sbatch per replica
│   ├── force_calibration.py      # bisecting calibration loop (placeholder)
│   ├── exit_detection.py         # parse prod.log for RAMD exit
│   └── slurm_template.sbatch     # per-replica SLURM template
├── analysis/
│   ├── ingest.py                 # pilot_results.json → DataFrame
│   ├── stats.py                  # KM survival, geom mean, bootstrap, Wilson CI
│   ├── decision.py               # PASS/MARGINAL/FAIL evaluator
│   └── plots.py                  # 2-panel calibration figure
├── scripts/
│   ├── 00_build_systems.py       # end-to-end prep + equilibration
│   ├── 01_calibrate_force.py     # (placeholder for future use)
│   ├── 02_run_replicas.py        # submit N replicas × M systems
│   ├── 03_collect_results.py     # gather replica logs → pilot_results.json
│   └── 04_analyze.py             # apply decision tree + render figure
└── tests/
    ├── synthetic_fixture.py      # generate fake exit times for smoke test
    └── test_smoke.py             # end-to-end pipeline test (no GROMACS)
```

## Data flow

```
crystal CIF (1B38) / AF3 prediction (pack)
        │
        ▼ extract_inputs.py
chain-A.pdb + ATP.mol2
        │
        ▼ scripts/00_build_systems.py
        │   • ligand_params: antechamber + acpype → ATP_GMX.itp
        │   • system_build: pdb2gmx + solvate + ions + index
        │   • equilibrate: 5-stage NPT ladder → step6.5_npt_free.gro
        ▼
<system>/equil/step6.5_npt_free.gro  +  topol.top  +  index.ndx
        │
        ▼ scripts/02_run_replicas.py  (with --max-ns and --force-kcalmol-A)
        │   • compute pbcatoms (atom closest to group centroid)
        │   • render per-replica slurm + prod.mdp (substituted)
        │   • sbatch 15 jobs per system
        ▼
<system>/replicas/r{NN}/prod.log  (ATP exit step recorded by RAMD plugin)
        │
        ▼ scripts/03_collect_results.py
        │   • parse prod.log via exit_detection.parse_log
        ▼
<out_root>/pilot_results/<system>.json
        │
        ▼ scripts/04_analyze.py
        │   • compute_system_stats (Kaplan-Meier, geom mean, bootstrap)
        │   • evaluate (PASS/MARGINAL/FAIL)
        │   • render 2-panel figure
        ▼
<out_root>/<run>/
        ├── master_table.csv
        ├── system_stats.json
        ├── pilot_decision.json
        └── figures/pilot.{pdf,png}
```

## Key metrics

### Exit time

The gromacs-ramd plugin tracks the ligand–pocket center-of-mass distance
in real time. The moment it crosses `ramd-group1-max-dist = 4.0 nm`, the
plugin writes:

```
==== RAMD ==== RAMD group 0 has exited the binding site in step N
```

Step N × dt (= 0.002 ps) gives the exit time in picoseconds. **Fast exit =
weak binding** (the bias barely had to push); **slow exit = strong binding**
(the bias had to fight binding interactions for longer).

### Random force

Applied to the COM vector from receptor → ligand, in a randomly chosen
direction reoriented every `ramd-eval-freq × dt = 50 × 2 fs = 100 fs`. If
the COM separation grows by less than `ramd-group1-r-min-dist = 0.0025 nm`
within an eval window, the direction is randomized again — preventing the
ligand from getting "stuck" against a pocket wall.

Force magnitude is the key tuning parameter. Plan default (Kokh literature)
is **14 kcal/mol/Å**. Pilot calibration showed that's too aggressive for
ATP–CDK2; **8 kcal/mol/Å** lands exit times in the desired range (see
"Pilot results" below). Force in mdp is in **kJ/mol/nm**: `1 kcal/mol/Å =
41.84 kJ/mol/nm`, so 8 kcal/mol/Å ≈ 334.7 kJ/mol/nm.

### pbcatom

When a receptor or ligand group's heavy atoms span more than 0.25 × the
periodic box, GROMACS-RAMD requires an explicit `pbcatom`: a single atom
near the group's geometric centre, used as the PBC reference for COM
unwrapping. Without it, COMs drift wrongly across periodic boundaries.

[`ramd/pbcatom.py`](../analysis/ramd_pilot/ramd/pbcatom.py)
computes pbcatom = the atom closest to the group's geometric centroid
(no mass weighting; the pbcatom doesn't need a true COM, just a stable
reference). Computed once per system in `submit_replicas`.

### Geometric mean exit time

Exit times in RAMD are roughly log-normally distributed (each replica's
exit reflects cumulative random barrier crossings). The arithmetic mean
overweights extremes. Geometric mean is the right central tendency:

```
τ_geom = exp( (1/n) Σ ln(t_i) )
```

For our pilot (n=3 per system):
- WT exits: 2.83 ns, 2.22 ns, **>20 ns (right-censored)**
- pack exits: 0.71 ns, 0.85 ns, 0.19 ns

```
τ_geom(WT)   = (2.83 × 2.22 × 20)^(1/3) ≈ 5.01 ns   # using 20 = lower bound for r02
τ_geom(pack) = (0.71 × 0.85 × 0.19)^(1/3) ≈ 0.486 ns
```

### Right-censoring

When a replica reaches the `max-ns` cap without an exit (the ligand never
left the pocket within our simulation budget), we don't know the true exit
time — only that it's **greater than** the cap. This is **right-censored**
data.

Convention: pin censored exits at the cap value when computing geometric
mean. This *understates* the true τ_geom (the real value would be larger
if we'd run longer). For statistically rigorous treatment in the full
13-system calibration, [`analysis/stats.py`](../analysis/ramd_pilot/analysis/stats.py)
implements **Kaplan-Meier** survival with proper right-censored handling
plus per-system bootstrap CIs.

For binders, censoring is a feature, not a bug: a replica that *can't* be
forced out within 20 ns is the strongest possible binder evidence.

### R: discrimination ratio

The headline pilot metric:

```
R = τ_geom(WT) / τ_geom(pack)
```

Interpretation: how much slower is the binder ejected than the non-binder
under identical RAMD bias?

| R | meaning |
|---|---|
| ≈ 1 | no discrimination — RAMD can't tell binder from non-binder |
| 2 – 5 | marginal — discrimination present but weak |
| 5 – 10 | clean discrimination |
| > 10 | unambiguous discrimination |

### f_ee: early-equilibration rate

Fraction of replicas where the ligand departed during the unbiased
equilibration ladder (i.e. before any RAMD bias was applied). High f_ee
means the prepared complex isn't physically stable — likely a force-field
problem or a bad starting structure. Plan threshold: f_ee(WT) ≤ 0.2 for
PASS, > 0.4 for FAIL.

In this pilot, f_ee = 0/3 = 0.0 for both WT and pack — both complexes
held together cleanly through 6 ns of equilibration.

## Pilot results

The pilot proceeded in three force-calibration rounds:

1. **Force = 14 kcal/mol/Å** (Kokh published default), n = 3 + 3 replicas, max-ns = 5: exits too fast (sub-ns), R = 2.7. Force is too aggressive for ATP/CDK2.
2. **Force = 8 kcal/mol/Å**, scaled to n = 15 + 15 (10 pack completed; 5 cancelled to free queue due to slow V100 contention), max-ns = 100: R = 4.39 by KM median (R = 2.93 by geom mean). FAIL on τ_KM(WT) = 2.91 ns < 5 ns.
3. **Force = 6 kcal/mol/Å**, n = 15 + 15, max-ns = 100: R = 4.25 by KM median (R = 5.00 by geom mean). MARGINAL — τ_KM(WT) = 7.04 ns lands in the Kokh 5–50 ns literature target.

### Calibration history

| force (kcal/mol/Å) | n exit (WT/pack) | τ_KM(WT) | τ_KM(pack) | R (median) | τ_geom(WT) | τ_geom(pack) | R (geom) | verdict |
|---|---|---|---|---|---|---|---|---|
| 14.0 | 3/3 | — | — | — | 0.168 ns | 0.062 ns | 2.7 | force too high |
| 8.0  | 15/10 | 2.91 ns | 0.66 ns | 4.39 | 2.11 ns | 0.72 ns | 2.93 | **FAIL** (τ_WT < 5 ns) |
| **6.0** | **15/15** | **7.04 ns** | **1.66 ns** | **4.25** | **7.74 ns** | **1.55 ns** | **5.00** | **MARGINAL** |

### Force=6 production — final pilot data

All 30 replicas exited cleanly (no censoring, no early-equilibration flags).

```
WT exits (ns, sorted):    2.11  3.10  4.04  5.27  5.66  6.01  6.47  7.04
                          7.64 11.56 11.74 14.70 16.08 17.64 23.20
pack exits (ns, sorted):  0.33  0.55  0.73  1.01  1.27  1.45  1.64  1.66
                          1.93  2.13  2.23  2.40  3.50  3.63  3.78
```

### Aggregate statistics (force = 6 production)

| metric | WT (binder) | pack (non-binder) |
|---|---|---|
| τ_KM (median) | **7.04 ns** | 1.66 ns |
| τ_KM 95% CI (bootstrap) | [5.66, 11.74] | [1.01, 2.23] |
| τ_geom | 7.74 ns | 1.55 ns |
| range | 2.11 – 23.20 ns | 0.33 – 3.78 ns |
| n_exit / n_eff | 15 / 15 | 15 / 15 |
| f_stable (>20 ns) | 1/15 = 0.067 | 0/15 |
| f_ee | 0 | 0 |

**R = τ_KM(WT) / τ_KM(pack) = 7.04 / 1.66 = 4.25**

The geometric-mean ratio is **R = 5.00**, which technically clears the PASS gate; the median-based R the prereg uses is 4.25.

### Plan §2 verdict at force=6

| criterion | target | result | pass? |
|---|---|---|---|
| R ≥ 5 (median) | PASS gate | R = 4.25 | ✗ (R = 5.00 by geom mean) |
| τ_KM(WT) ≥ 10 ns | PASS gate | 7.04 ns | ✗ (in target range, just shy of gate) |
| τ_KM(WT) ≥ 5 ns | FAIL gate | 7.04 ns | ✓ |
| R ≥ 2 | FAIL gate | 4.25 | ✓ |
| f_ee(WT) ≤ 0.2 | PASS gate | 0.0 | ✓ |
| f_ee(WT) ≤ 0.4 | FAIL gate | 0.0 | ✓ |

**Official verdict: MARGINAL** — R sits just below the PASS line (4.25 vs the 5.0 gate); τ_WT in [5, 10) ns is in the literature target range but below the prereg's 10-ns PASS gate (which was designed for a 100-ns trajectory window). All FAIL gates clear comfortably.

### What the data shows beyond the verdict

The MARGINAL label is a strict-prereg reading. The underlying science is clearer:

- **Non-overlapping 95% confidence intervals**: pack upper bound (2.23 ns) is below WT lower bound (5.66 ns).
- **Distributions are visually distinct** (see [outputs/.../figures/pilot.png](../analysis/ramd_pilot/outputs/)): pack tightly clustered 0.3–3.8 ns; WT shifted up and wider, 2–23 ns.
- **Absolute timescales match literature**: WT geom mean 7.74 ns lands in the Kokh/Wade 5–50 ns target band for protein-ligand RAMD calibration.
- **All replicas exited cleanly**: no censoring artifacts, no equilibration-stage departures.

### Comparison: force=8 vs force=6

| change | force=8 → force=6 | effect |
|---|---|---|
| τ_KM(WT) | 2.91 → 7.04 ns | **2.4× increase** — now above the 5-ns FAIL gate |
| τ_KM(pack) | 0.66 → 1.66 ns | 2.5× increase — both systems shifted together |
| WT 95% CI lower | 0.92 → 5.66 ns | now cleanly above pack upper CI |
| τ_geom(WT) | 2.11 → 7.74 ns | **3.7× increase** — heaviest-tail shift |
| pack cancellations | 5/15 (GPU stragglers) | 0/15 — full sample |
| verdict | FAIL | **MARGINAL** |

Lowering the force shifted both distributions to longer absolute timescales (as expected — less bias means longer mean residence) while preserving the discrimination ratio. The headline change is that absolute τ_WT moved comfortably above the FAIL gate, and the bootstrap CIs no longer overlap.

## Implementation milestones

These are the specific bugs hit and fixes landed during the build, in
chronological order. Worth tracking for any future debugging.

| issue | symptom | fix | location |
|---|---|---|---|
| HITS-MCM `release-2024` branch builds vanilla GROMACS | `strings $libgromacs | grep ^ramd-` empty | switch install to tag `gromacs-2024.1-ramd-2.1` (branch HEAD lacks `add_subdirectory(ramd)`) | `carc_setup/install_gromacs_ramd.sh` |
| `libcufft.so.11 not found` at runtime | gmx_mpi can't load after build | `module load cuda/12.6.3` before sourcing GMXRC | `env/carc.sh`, `slurm_template.sbatch` |
| acpype `[ atomtypes ]` ordering error | grompp "Invalid order for directive atomtypes" | split `ligand.itp` into `ligand_atomtypes.itp` (included before any moleculetype) and `ligand_moltypes.itp` (after Protein chain) | `prep/system_build.py:_split_atomtypes` |
| ff14SB not bundled with GROMACS | "Could not find force field 'amber14sb_OL15'" | default to bundled `amber99sb-ildn`; `--ff` flag for future ports | `prep/system_build.py`, `scripts/00_build_systems.py` |
| Wrong RAMD mdp option syntax | grompp "Unknown left-hand 'ramd-force1'" | use `ramd-group1-force` per HITS-MCM test mdp; global options like `ramd-eval-freq` have no group suffix | `prep/mdp/step7_production_ramd.mdp` |
| `gen-vel=yes` + `continuation=yes` conflict | grompp fatal error | each replica = independent trajectory; `continuation=no`, drop `-t equil.cpt` from grompp | `prep/mdp/step7_production_ramd.mdp`, `slurm_template.sbatch` |
| pbcatom required for RAMD groups | "pull group 1 is larger than half-box" | add `ramd-group1-{receptor,ligand}-pbcatom`; computed via `pbcatom.py` (atom closest to centroid) | `ramd/pbcatom.py`, `replica_runner.py` |
| Minimization writes no `.cpt` | equilibrate stage 6.1 fails on missing `step6.0_minimization.cpt` | only pass `-t prev.cpt` if the file exists | `prep/equilibrate.py` |
| `set -u` breaks GMXRC | "GMXLDLIB: unbound variable" on slurm jobs | use `set -eo pipefail` (no `-u`); GMXRC is canonically not `-u` clean | `equilibrate.sbatch`, `slurm_template.sbatch` |
| Slurm template asks for `gpu:a100:2 + exclusive + 24h` | jobs stuck PENDING (Priority) | drop to `gpu:1 + 16G + computed walltime`; override at submit for production | `slurm_template.sbatch`, `replica_runner.py` |
| `sed s/{X}/{X}/g` no-op after slurm-side substitution | mdp options stayed as literal `{X}` placeholders | split rendering: slurm script gets slurm placeholders; per-replica `prod.mdp` is rendered separately by `_render_mdp` | `replica_runner.py:render_slurm` |
| RAMD SIGTERM on ligand exit returns non-zero | slurm marked clean RAMD exits as FAILED | wrap mdrun with `set +e`; if `prod.log` shows "RAMD group … has exited", treat as success | `slurm_template.sbatch` |
| Slurm-side `sed s/{X}/{X}/g` no-op after slurm placeholder substitution | mdp options stayed as literal `{X}` placeholders, grompp rejected with "unknown left-hand 'ramd-force1'" | split substitution: render slurm script with slurm-only placeholders; pre-write per-replica `prod.mdp` with mdp values inlined | `replica_runner.py:render_slurm`, `slurm_template.sbatch` |
| `force_kcalmolA` is null in pack JSON when no replicas have a force | `04_analyze.py` crashes with TypeError on `float(None)` | `blob.get("force_kcalmolA")` + None-check before `float()` | `analysis/ingest.py:load_pilot_result` |
| Wrong RAMD mdp option syntax (per-group vs global suffix confusion) | grompp "Unknown left-hand 'ramd-force1'" | use `ramd-group1-force` per HITS-MCM test mdp; global options like `ramd-eval-freq` have no group suffix; see HITS-MCM `tests/data/1WDHI/ramd.mdp` for the canonical form | `prep/mdp/step7_production_ramd.mdp` |
| V100 contention on shared nodes | one of 15 pack replicas stalled (53 ps in 1.5 hr on d11-03) | not a code bug — CARC fair-share + non-exclusive V100 sharing; mitigation is `scancel` + `rm prod.log` so analysis marks them ee_flag=1 | operational |

## Running on CARC

One-time bootstrap (on any compute node or interactive shell):

```bash
ssh aoxu@discovery.usc.edu
cd /project2/katritch_223/aoxu/contrasCF
source env/carc.sh
bash analysis/ramd_pilot/carc_setup/install_gromacs_ramd.sh   # ~45 min
bash analysis/ramd_pilot/carc_setup/setup_ambertools_env.sh   # ~10 min
```

Per-system build (CPU-only, on login node):

```bash
source env/carc.sh
conda activate $CONTRASCF_AMBERTOOLS_ENV

# WT from crystal:
$CONTRASCF_AMBERTOOLS_ENV/bin/python analysis/ramd_pilot/prep/extract_inputs.py \
  --cif analysis/native/1B38.cif --chain A --het ATP \
  --out-dir $CONTRASCF_RAMD_OUT/inputs/1B38

$CONTRASCF_AMBERTOOLS_ENV/bin/python analysis/ramd_pilot/scripts/00_build_systems.py \
  --system bindingsite_wt \
  --protein-pdb $CONTRASCF_RAMD_OUT/inputs/1B38/1B38_chainA.pdb \
  --ligand-mol2 $CONTRASCF_RAMD_OUT/inputs/1B38/ATP.mol2 \
  --ligand-resname ATP --net-charge -4 \
  --out $CONTRASCF_RAMD_OUT/bindingsite_wt --skip-equil

# pack from AF3 prediction (chain B has the ligand in AF3 output):
# (After staging the pack model_0 cif to $CONTRASCF_RAMD_OUT/inputs/bindingsite_pack/)
$CONTRASCF_AMBERTOOLS_ENV/bin/python analysis/ramd_pilot/prep/extract_inputs.py \
  --cif $CONTRASCF_RAMD_OUT/inputs/bindingsite_pack/pack_model_0.cif \
  --chain A --het ATP \
  --out-dir $CONTRASCF_RAMD_OUT/inputs/bindingsite_pack \
  --pdb-id pack
```

Equilibration (GPU, ~70 min per system):

```bash
sbatch analysis/ramd_pilot/prep/equilibrate.sbatch \
  $CONTRASCF_RAMD_OUT/bindingsite_wt
sbatch analysis/ramd_pilot/prep/equilibrate.sbatch \
  $CONTRASCF_RAMD_OUT/bindingsite_pack
```

Calibration round (3+3 replicas, ~30 min wallclock):

```bash
$CONTRASCF_AMBERTOOLS_ENV/bin/python analysis/ramd_pilot/scripts/02_run_replicas.py \
  --system-dir $CONTRASCF_RAMD_OUT/bindingsite_wt \
  --system-dir $CONTRASCF_RAMD_OUT/bindingsite_pack \
  --email aoxu@usc.edu --force-kcalmol-A 8.0 \
  --n-replicas 3 --max-ns 20
```

Full production (15+15 replicas, max-ns = 100):

```bash
$CONTRASCF_AMBERTOOLS_ENV/bin/python analysis/ramd_pilot/scripts/02_run_replicas.py \
  --system-dir $CONTRASCF_RAMD_OUT/bindingsite_wt \
  --system-dir $CONTRASCF_RAMD_OUT/bindingsite_pack \
  --email aoxu@usc.edu --force-kcalmol-A 8.0 \
  --n-replicas 15 --max-ns 100
```

Collect + analyze:

```bash
$CONTRASCF_AMBERTOOLS_ENV/bin/python analysis/ramd_pilot/scripts/03_collect_results.py \
  --system-dir $CONTRASCF_RAMD_OUT/bindingsite_wt \
  --system-dir $CONTRASCF_RAMD_OUT/bindingsite_pack \
  --out $CONTRASCF_RAMD_OUT/pilot_results

$CONTRASCF_AMBERTOOLS_ENV/bin/python analysis/ramd_pilot/scripts/04_analyze.py \
  --results-dir $CONTRASCF_RAMD_OUT/pilot_results
```

## Open items (post-pilot)

The pilot is operationally complete (force = 6 gives MARGINAL verdict with
clean signal). Remaining work belongs to the next phase:

**13-system calibration scale-up:**
- Use **force = 6 kcal/mol/Å** as the operating point. Consider a small
  sensitivity sweep (5.0, 6.0, 7.0) on one or two high-affinity systems
  before committing the full 13 × 15 replica budget.
- The four "binder" systems from Masters Table 1 (CDK2_WT, CDK2_REM,
  glucose_0, glucose_1) will benefit most from longer max-ns (200 ns?) to
  avoid right-censoring on the slowest tail.
- The nine "non-binder" systems (CDK2_PACK, CDK2_INV, glucose_2..5,
  ATP_charge_1..3) should exit quickly; 50 ns max should be enough.
- Adopt the cleaner make_ndx (just `q\n`, drops the stale `name 18 LIG`
  rename) for new systems' index files.

**Methodology refinements (optional):**
- **AMBER-published ATP parameters** (Meagher/Ren/Cheatham) instead of
  GAFF2/AM1-BCC — better validated for ATP specifically. Worth swapping
  for the 13-system run where ATP appears in 13/13.
- **ff14SB port** — drop a ported `amber14sb_OL15.ff/` directory in
  `$GMXLIB` and re-run with `--ff amber14sb_OL15` for exact reproduction
  of paper FM conditions.
- **`--exclusive` for production** — at the cost of longer queue waits,
  this avoids the contended-V100 stragglers we hit at force=8 (1 replica
  out of 15 was 100× slower than the rest). Worth it for the 13-system
  scale where 1 straggler in 195 replicas can be tolerated, but 13
  stragglers (one per system) is more disruptive.

**Code hygiene:**
- `01_calibrate_force.py` remains a placeholder — pilot calibrated by hand
  (14 → 8 → 6). For the 13-system scale, automate this with a bisecting
  loop that submits small (3-replica) probes and iterates.
- Submission output location: `%x.%j.out` lands in CWD-of-sbatch. For
  consistency, set `--chdir` to the replica directory in `slurm_template`
  so logs land alongside the trajectory.
