# Lab notebook — 2026-06-06 → 06-07 — Pose-swap test: is Boltz-2's affinity head pose-sensitive?

**Project:** contrasCF (do co-folding models learn protein–ligand physics?). **Branch:** `pose-swap`. **Operator:** Aaron Xu (with Claude). **Instrument/env:** lab RTX 4090; `boltzina_env` (Boltz-2 v2.2.1, torch 2.7.0+cu126); analysis via `protenix` env; checkpoints in `/home/aoxu/.boltz` (rsync'd from CARC `~/.boltz`).

## Objective
The CASF-mutagenesis analysis (see `docs/casf_confidence.md`) showed Boltz-2's affinity head barely responds to pocket mutations (median Δlog[IC50] +0.10–0.17 vs a physics expectation of +3 to +6). But a mutation changes **both** the trunk (sequence/MSA → `s`, `z`, `s_inputs`) **and**, downstream, the predicted pose (`x_pred`). **Open question:** is the head *intrinsically* pose-insensitive, or pose-sensitive but the mutations simply don't displace the ligand far enough to register?

## Rationale / design
Source read (spec `docs/superpowers/specs/2026-06-06-pose-swap-test-design.md`): the affinity head reads pose **only** via a distogram of `x_pred` over protein–ligand cross-pairs (`affinity.py:104–108`); the trunk (`s_inputs`, `z`) is pose-independent. So the confound is removable: **hold the trunk fixed, vary only `x_pred`, read affinity.** If a ligand ejected into solvent — imposed directly on the scored coordinates — does not weaken the affinity, the head is functionally pose-insensitive.

## Procedure
- **Hook (script `19_pose_swap_affinity.py`):** monkeypatch `AffinityModule.forward` so that, on the production call, it *also* scores a ladder of decoy `x_pred` with the same modules + captured inputs, then returns the native result. Run via `boltz predict` in-process (`cli(..., standalone_mode=False)`), `recycling_steps=3` to match the casf pipeline. The native call **is** production → identity check is automatic.
- **Readouts per pose:** `affinity_pred_value` (raw ensemble = production value; MW correction is a constant offset across a same-ligand ladder, so it cancels in every Δ), `affinity_probability_binary`, and the **min ligand–protein distance** (added as a verification diagnostic).
- **Decoys:** native; ligand ejected to clearance 5/15/30 Å beyond the protein bounding sphere (radial); random-direction ejection (15 Å); in-pocket rotation; whole-complex translation (sanity).
- **Panel:** 29 diverse CASF WT systems (spread of the 252 affinity-enabled systems). Aggregate (script `20_pose_swap_aggregate.py`): per-system Spearman(separation, aff), slope, native→eject30 gap; panel medians + paired Wilcoxon; figure.

## Problem found & fixed (recorded for the methods)
The **first** decoy scheme translated the ligand along a *local* pocket-exit vector (pocket-centroid → ligand COM). For **buried** pockets that vector is near-arbitrary (e.g. 1bcu: ligand COM is 2.0 Å from its pocket centroid), so a 20 Å translation just **plowed the ligand through the protein** — the min ligand–protein distance stayed **~3 Å** even at "20 Å displacement." That tested *"ligand slid to new contacts,"* **not** *"ligand ejected to solvent."* It was caught by the min-distance diagnostic (added in response to a "which direction?" review question). **Fix:** place the ligand COM **beyond the protein bounding sphere** (`centroid + outward·(R_protein + clearance)`) and **verify** via the recorded min-distance. After the fix: native ~3 Å → eject30 **30–45 Å** (genuinely solvent-exposed). *Lesson: always verify the achieved ligand–protein separation, not the nominal translation.*

## Results (n=29, true ejection)
| metric | value |
|---|---|
| median ligand–protein distance @ `eject30` | **35 Å** (all 30–45 Å; past the 22 Å distogram saturation) |
| median native→eject30 **gap** | **−0.004** log units (IQR −0.036…+0.047) — physics wants +3 to +6 |
| systems weakening by ≥ +1 log unit | **0 / 29** |
| per-system slope | ~0.000–0.002 log/Å (~100× below physics) |
| paired Wilcoxon (native vs ejected, greater) | **p = 0.27** (n.s.) |
| whole-complex-translate sanity (must = 0) | **0.000 on all 29** |
| identity check (native = production affinity) | ✓ exact |
| direction control (random-dir ejection vs radial) | matches (direction-independent) |

Figure: `analysis/casf_mutagenesis/figures/pose_swap_affinity.png` — every per-system curve flat near 0 against the +3 physics floor; `P(binder)` within ±0.03. (The per-system `rho` is scattered +1↔−1 — the *direction* of the sub-0.05 noise is random; the *magnitude* is flat everywhere.)

## Interpretation / conclusion
**Boltz-2's affinity head is functionally, intrinsically pose-insensitive.** A ligand floating 35 Å in solvent (zero protein contacts) is predicted to bind essentially as well as the crystal pose. This is a *construction-based* confirmation (not a correlation): we moved the ligand ourselves, maximally, with the trunk held byte-identical, and the head did not react — so the CASF insensitivity is intrinsic, **not** "the mutations don't move the ligand enough." Architecturally the head *does* consume `x_pred` (distogram); the *learned function* is flat in pose. This closes the diagnosis loop: **structure** half-memorizes, **confidence** registers the break (and "knows it memorized"), **affinity** is blind to ligand geometry.

## Next steps
- **GNINA `--score_only`** reference on the identical decoy poses (a genuine-physics scorer should collapse on the ejected decoys) — *deferred* until the GNINA env is configured; not load-bearing for the core claim (the min-distance + whole-complex controls already establish the decoys are real non-binder geometries).
- Optionally scale to all 252 systems overnight (verdict won't change; tighter CIs only).
- Fold this result + figure into `docs/casf_confidence.md` as the capstone section.

## Provenance
Commits on `pose-swap` (pushed to origin): `5d05e23` (spec), `ab55b0c` (initial hook B + n=11, flawed decoys), `698945a` (buried-pocket fix + n=29). Scripts `19_…` (driver), `20_…` (aggregator+figure). Memory: `pose-swap-result`, `boltz2-affinity-architecture`.

---

## Addendum 2026-06-09 — GNINA reference (the deferred item, now done)

**Setup:** gnina v1.3 (`/mnt/katritch_lab2/aoxu/envs/gnina/bin/gnina`). Per system: take crystal docking inputs (`receptor.pdb` + `ligand.sdf`), eject the ligand radially to the same clearances (5/15/30 Å beyond the protein bounding sphere), rescore each pose with `gnina --score_only` (CNNaffinity pK; Vina kcal/mol). Scripts `21_pose_swap_gnina.py` (driver), `22_pose_swap_contrast.py` (contrast + figure).

**Result (n=29, same panel) — a genuine physics scorer reacts exactly as it should; Boltz-2 does not:**

| metric, ligand ejected ~30+ Å | response |
|---|---|
| **Boltz-2 affinity head** | median Δlog[IC50] = **−0.004** (flat) |
| **GNINA Vina (physics)** | median **−8.7 → 0.0** kcal/mol; **100% (29/29)** collapse to ~0 — all binding energy lost |
| **GNINA CNNaffinity** | median drop **+2.5 pK** (1–4 pK; the CNN affinity has its own non-zero floor) |
| **GNINA CNNscore** (pose plausibility, 0–1) | median **0.95 → 0.54**; collapses (<0.4) for 21% of systems, floors (~0.9) for the rest |

**Hierarchy of "did the scorer notice the ligand left?":** pure physics (Vina) — completely (100% → 0); GNINA's learned CNN — *partially* (CNNscore 0.95→0.54, CNNaff −2.5 pK; floors for most systems); Boltz-2's affinity head — *not at all* (Δ−0.004). All learned heads memorize to some degree; Boltz-2's is the most extreme.

**Command check (Aaron's review):** verified `--score_only` against the lab's standard full-docking command (`gnina --autobox_ligand`, from `CogLigandBench/.../gnina_inference.py`) on 1bcu native — they agree (CNNaffinity 5.34 vs 5.37; Vina −7.19 vs −7.78). `--score_only` is the correct mode here because we score a *fixed imposed* pose; docking would re-search and discard the ejection. gnina's low CNNaffinity (~5.3) is its own calibration, identical from both commands.

Figure: `analysis/casf_mutagenesis/figures/pose_swap_contrast.png` (3 panels: Boltz flat | GNINA Vina collapse | GNINA CNNscore partial). The Vina term hitting exactly 0 for every system confirms the ejected decoys are unambiguous non-binder geometries; Boltz-2's affinity head is alone in not noticing. **Conclusion complete and externally referenced.** Branch merged to `main`.
