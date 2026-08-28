# surfdock-interface-crop

**Type:** postmortem
**Date:** 2026-08-28
**Severity:** experiment-invalidated

## What broke

The published SurfDock CASF-2016 result — **0/244 WT under 2 Å, median 6.82 Å**, written up in `docs/casf_overview.md` as *"CASF-2016 sits outside SurfDock's training distribution"* — was wrong. It was an artifact of **our own** surface preprocessing, not a property of SurfDock. All 968 SurfDock cells were invalid.

## What I saw

Poses were not wild — they were *systematically* displaced, which is what made the "model failure" reading plausible:

- rank1 centroid sat a median **4.48 Å** from the docking-box centre (max 6.1 over 10 WT cells), not the "thousands of Å" the earlier write-up claimed. Exactly one of eleven rank5 poses was a true outlier (119,384 Å).
- My independent symmetry-aware RMSD reproduced the pipeline's number (rank1 median **6.53 Å**, 0/8 under 2 Å), so the analysis code was not at fault.
- SurfDock's own filename-encoded RMSD agreed: `..._rank_1_rmsd_6.7988787227593415_confidence_53.97...sdf` for `1e66/wt`.

The decisive observation: `1a0q`, **SurfDock's own shipped test system**, failed identically through our pipeline — rank1 **8.30 Å** with the anchor flag, **7.74 Å** under SurfDock's own protocol.

## What I expected

SurfDock's published performance is roughly 50–80 % under 2 Å. This same install scores rank-1 median **0.98 Å / 80 % under 2 Å on PoseBusters** (n=428, `CogLigandBench/notebooks/01_method_analysis/surfdock/surfdock_posebusters_results.csv`). So the model, weights and env were demonstrably capable; only the contrasCF path was failing.

## Hypotheses for why

- **numpy shape bug at `inference_accelerate.py:208`** (the previously documented root cause: `np.array([(x,y,z)])` giving `(1,3)` instead of `(3,)`) — **REJECTED, proven a no-op.** In the only expression that consumes it (`utils/sampling.py:36`, `pos - mean(pos, dim=0, keepdim=True) + pocket_center`) a `(1,3)` broadcasts identically to a `(3,)`: same result shape, same centroids, verified numerically. The anchoring it was meant to restore was working all along — the run log shows `Use predict pocket center` firing 40× per complex against the *unpatched* install.
- **`--ligand_to_pocket_center` breaks the diffusion prior** — **PARTIALLY TRUE, not the cause.** `randomize_position` (`utils/sampling.py:33-41`) is an if/else, so the flag *replaces* the trained `Normal(0, tr_sigma_max=5.0)` translational prior with a deterministic delta at the pocket centre. Removing it improved 1e66 6.80→4.87 Å and 1a0q 8.30→7.74 Å — real, but ~2 Å, not the 5–6 Å gap.
- **Broken env / weights** — **REJECTED.** Score model loads 376/377 tensors; the single "missing" key `final_conv.batch_norm.bias` is an empty tensor of shape `(0,)`. Installed versions match SurfDock's pinned `environment.yaml` on every critical package (torch 2.2.2, pyg 2.5.2, torch-scatter/cluster `pt22cu121`, e3nn 0.5.1, numpy 1.24.4, scipy 1.8.1); only rdkit differs (2022.09.5 vs pinned 2023.3.1).
- **Bad inputs** — **REJECTED.** WT box centre equals the crystal ligand centroid to **0.000 Å** across 15 systems; the CSV `pocket_center` is exactly that centroid; the 8 Å pocket is 908 atoms / 58 residues; the surface carries populated charge (±14), hbond and hphob (±4.5 Kyte–Doolittle) channels; the ESM key matches with shape (58, 1280).
- **Missing interface crop** — **CONFIRMED** (below).

## Root cause

`CogLigandBench/dockstrat/models/_surfdock_surface_helper.py:128` kept **every** mesh face:

```python
faces_to_keep = np.arange(len(faces2))
```

SurfDock's own `comp_surface/prepare_target/computeTargetMesh_test_samples.py:100-104` — invoked at its line 249 with `dist_threshold=8` — instead keeps only faces whose vertices lie within `dist_threshold - 5` = **3 Å of the ligand**:

```python
kdt = KDTree(atomCoords)
d, r = kdt.query(vertices1)
iface_v = np.where(d <= dist_threshold - 5)[0]
faces_to_keep = [idx for idx, face in enumerate(faces1) if all(v in iface_v for v in face)]
```

So the model received the **entire 8 Å pocket surface** instead of the ligand-proximal interface patch it was trained on — a straight train/serve mismatch. On the identical 769-atom `1a0q` pocket: **1474 vertices (ours) vs 142 (SurfDock's own reference)**. Across SurfDock's working sample dirs meshes run 59–261 vertices; ours ran 1474–1833.

It degraded every pose silently rather than failing, which is why it survived as a "finding" for three months.

## Fix

- **CogLigandBench `2107de18e9`** — restores the KDTree interface crop. Indices are taken against `vertices2/faces2`, not `vertices1/faces1`, because `remove_isolated_vertices()` re-indexes the mesh; warned fallback to the full mesh if the crop is empty.
- **contrasCF `dfabdaf`** — drops `--ligand_to_pocket_center`, rewrites `SURFDOCK_FIX.md` with the correct diagnosis, retracts the `casf_overview.md` section, deletes the three no-op artifacts (`apply_surfdock_fix.py`, `fix_and_run_surfdock.sh`, `surfdock_pocket_center_shape_fix.patch`).
- **contrasCF `711344e`** — `14_run_surfdock_variants.py` now exports the SurfDock env vars itself (it bypasses `13_`'s `main()`, the only place they were set).

Verification, rank-1 RMSD (SurfDock's own self-reported value):

| system | broken surface | fixed surface |
|---|---|---|
| `1a0q` (SurfDock's own test system) | 7.74 Å | **0.44 Å** |
| `1e66` (CASF WT) | 6.80 Å | **0.32 Å** |
| `1gpk` / `1gpn` / `1h23` (CASF WT) | 5.15 / 6.65 / 4.96 Å | **0.42 / 0.46 / 1.30 Å** |

## Prevention

- [confirmed] Vertex-count health check documented in the `dockstrat` skill and in `SURFDOCK_FIX.md`: `grep "element vertex" <work_dir>/surface/<case>/*_8A.ply` — healthy is 60–260; ~1500+ means the crop regressed.
- [confirmed] `dockstrat` skill gotcha 6d, `references/surfdock.md` and `references/hosts.md` corrected — they previously *recommended* `--ligand_to_pocket_center` as the cure, i.e. the skill would have re-introduced the smaller half of this bug.
- [confirmed] Recorded the PoseBusters reference point (0.98 Å median / 80 % under 2 Å) as the cross-check: a SurfDock sweep far off that on ordinary WT systems is preprocessing, not the model.
- [confirmed] `TMPDIR` redirect + janitor documented — `computeMSMS` writes four scratch files per call (~3.6 MB/cell) into `tempfile.gettempdir()` and never deletes them; on this host `/` is a 49 GB partition with ~1 GB headroom, so a full sweep fills it (see Blast radius).
- [tentative] The surface helper swallows exceptions (`except Exception: print("[WARN]…"); return`), so an MSMS failure exits 0 and surfaces two steps later as an opaque pandas `EmptyDataError`. Tightening it would have saved a reproduction run.
- [open] Cause of the 8 "0 graphs" failures is still unknown — see Blast radius.

## Blast radius

- [confirmed] All **968** SurfDock cells were invalid and have been regenerated. The pre-fix cells are preserved at `analysis/casf_mutagenesis/outputs/_surfdock_invalid_backup/` (159 MB) rather than deleted, so the retracted numbers stay auditable.
- [confirmed] `docs/casf_overview.md`'s SurfDock section and the SurfDock bar in `figures/overview_full.png` were wrong and are retracted; the figure has been regenerated from the corrected sweep.
- [confirmed] `SURFDOCK_FIX.md` documented a fix that does nothing, and the shared install at `/home/aoxu/projects/SurfDock` needs **no** patch. Anyone who applied `apply_surfdock_fix.py` changed nothing.
- [confirmed] No other engine is affected — GNINA and UniDock2 do not use the SurfDock surface path, and their numbers were regenerated independently under `e46d2ea`.
- [open] **8 cells fail with "0 graphs" and I do not know why.** My first explanation — degenerate meshes below SurfDock's ~59-vertex floor — is **refuted** by the mesh census (`analysis/casf_mutagenesis/mesh_census.json`): failed cells span 42–158 vertices while successful ones span 7–166, so mesh size does not predict failure. 7 of 8 are `inv`, so the cause is variant-linked, but that is all that is established.
- [open] The sweep's first attempt died at cell 477/968 with ENOSPC after MSMS scratch filled the 49 GB root partition. Root-fs headroom on this host is a standing hazard for any long job, not just SurfDock; the `TMPDIR` redirect works around it but does not fix the partition.
