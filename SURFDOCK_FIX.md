# SurfDock rank-1 positioning — root cause and fix

**Status:** RESOLVED 2026-08-25. Supersedes the 2026-05-28 version of this file,
whose diagnosis was **wrong** (see "The retracted diagnosis" below).

## Symptom

SurfDock rank-1 poses landed 5–8 Å from the crystal ligand on essentially every
system, including wild-type. On the full CASF sweep this read as **0/244 under
2 Å** with a 6.82 Å median — a total failure that was written up as "CASF-2016
sits outside SurfDock's training distribution".

That interpretation was wrong too. The same install scores **rank-1 median
0.98 Å, 80 % under 2 Å** on PoseBusters (428 systems,
`CogLigandBench/notebooks/01_method_analysis/surfdock/surfdock_posebusters_results.csv`).
The model and weights were never the problem.

## Root cause — the missing interface crop

`CogLigandBench/dockstrat/models/_surfdock_surface_helper.py` built the pocket
surface but **skipped SurfDock's ligand-proximity crop**:

```python
faces_to_keep = np.arange(len(faces2))    # keep EVERY face
```

SurfDock's own `comp_surface/prepare_target/computeTargetMesh_test_samples.py`
(lines 100-104, invoked at line 249 with `dist_threshold=8`) instead keeps only
faces whose vertices lie within `dist_threshold - 5` = **3 Å of the ligand**:

```python
kdt = KDTree(atomCoords)
d, r = kdt.query(vertices1)
iface_v = np.where(d <= dist_threshold - 5)[0]
faces_to_keep = [idx for idx, face in enumerate(faces1) if all(v in iface_v for v in face)]
```

So the model was handed the **entire 8 Å pocket surface** instead of the
ligand-proximal interface patch it was trained on — roughly 10× too many
vertices. On the identical 769-atom `1a0q` pocket: **1474 vertices (ours) vs 142
(SurfDock's own reference)**. Across SurfDock's working sample dirs the meshes
run 59–261 vertices; ours ran 1474–1833. A straight train/serve mismatch.

## Secondary defect — `--ligand_to_pocket_center`

`analysis/scripts/13_run_surfdock.py` also appended `--ligand_to_pocket_center`.
SurfDock's own eval scripts never pass it, and `randomize_position()`
(`utils/sampling.py:33-41`) is an **if/else**: the flag *replaces* the trained
`Normal(0, tr_sigma_max=5.0)` translational prior with a deterministic delta at
the pocket centre, which is off-distribution for the reverse diffusion. Removed.

## Measured effect (rank-1 RMSD, SurfDock's own self-reported value)

| system | broken surface + flag | broken, no flag | **fixed surface + flag** | **fixed, no flag** |
|---|---|---|---|---|
| `1a0q` (SurfDock's own test system) | 8.30 Å | 7.74 Å | 0.86 Å | **0.44 Å** |
| `1e66` (CASF WT) | 6.80 Å | 4.87 Å | — | **0.32 Å** |

Both fixes applied, five more CASF WT systems:

| system | before | after |
|---|---|---|
| `1gpk` | 5.15 | **0.42** |
| `1gpn` | 6.65 | **0.46** |
| `1h22` | 7.11 | 3.40 |
| `1k1i` | 7.61 | 3.25 |
| `1h23` | 4.96 | **1.30** |

Median 6.65 → **1.30 Å**; under-2 Å rate 0/5 → **3/5**.

## The retracted diagnosis

The previous version of this file blamed a numpy shape bug at
`inference_accelerate.py:208`:

```python
np.array([(float(x),float(y),float(z))])   # shape (1,3)  -- "buggy"
np.array([float(x),float(y),float(z)])     # shape (3,)   -- "fixed"
```

**This is a no-op.** In the only expression that consumes it
(`utils/sampling.py:36`, `pos - mean(pos, dim=0, keepdim=True) + pocket_center`)
a `(1,3)` broadcasts identically to a `(3,)`: same result shape, same centroids,
verified numerically. The anchoring it was supposed to restore was working the
whole time — the run log shows `Use predict pocket center` firing 40× per
complex against the *unpatched* install.

`apply_surfdock_fix.py`, `fix_and_run_surfdock.sh`, and
`surfdock_pocket_center_shape_fix.patch` implemented that no-op and have been
removed. The shared install at `/home/aoxu/projects/SurfDock` needs **no** patch.

## Consequence for existing results

**All 968 SurfDock cells under `analysis/casf_mutagenesis/outputs/*/*/surfdock/`
were produced with the broken surface step and are invalid.** They must be
regenerated, and the SurfDock row/bar in `docs/casf_overview.md` and
`figures/overview_full.png` must not be read until then.
