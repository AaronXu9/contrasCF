# Pose-swap test — Boltz-2 affinity head — design spec

**Date:** 2026-06-06
**Author:** Aaron Xu (with Claude)
**Status:** design, pending implementation (needs GPU + `boltzina_env`, Boltz-2 v2.2.1)

## Goal

Resolve a single mechanistic question left open by the CASF-mutagenesis affinity
analysis: **is Boltz-2's affinity head functionally pose-insensitive, or is it
pose-sensitive but the binding-site mutations simply don't move the ligand enough
to register?**

The empirical finding (`affinity-head-pose-insensitive.md`, `docs/casf_affinity.md`)
is that under `rem`/`pack`/`inv` pocket mutations the affinity head barely moves:
median Δlog[IC50] is +0.10 to +0.17 (vs a physics expectation of +3 to +6), and
once the WT baseline is partialled out, `partial(ΔAff, |displacement|) ≈ −0.06
(n.s.)` while `partial(ΔAff, ΔRMSD) = +0.12 (p=1e-3)` — a faint, far-sub-physics
residual. The dominant driver of ΔAff is regression-to-mean on WT tightness
(`Spearman(ΔAff, aff_wt) = −0.63`).

That analysis confounds two things, because a mutation changes **both** the trunk
(sequence/MSA → `s`, `z`, `s_inputs`) **and**, downstream, the predicted pose
(`x_pred`). The pose-swap test removes the confound by construction: **hold the
trunk fixed and vary only `x_pred`.** If a 20 Å ligand ejection — imposed directly
on the coordinates the head reads — does not weaken the predicted affinity, the
head is functionally pose-insensitive, full stop. If it does weaken monotonically,
then the head reads pose fine and the mutations were the weak lever.

## Why a clean swap is architecturally feasible

Read of Boltz-2 v2.2.1 source confirms the affinity head **does** consume the
predicted coordinates, and the trunk that feeds it does **not**. The swap is
therefore a surgical intercept, not a re-derivation.

### The affinity head reads pose only via a distogram of `x_pred`

`AffinityModule.forward(s_inputs, z, x_pred, feats)` —
`…/boltz/model/modules/affinity.py:77-139`:

1. `token_to_rep_atom = feats["token_to_rep_atom"]`; `x_pred_repr = torch.bmm(
   token_to_rep_atom.float(), x_pred)` (L95-104) — gathers one representative
   atom per token from the supplied coordinates.
2. `d = torch.cdist(x_pred_repr, x_pred_repr)` (L105) — the **only** place pose
   enters the head: a pairwise distance matrix over token representatives.
3. `distogram = (d.unsqueeze(-1) > self.boundaries).sum(-1).long()` then
   `self.dist_bin_pairwise_embed(distogram)` (L107-108), with
   `boundaries = torch.linspace(2, max_dist=22, num_dist_bins-1=63)` (L49) — 64
   distance bins, max 22 Å. **Any displacement that pushes a protein–ligand
   contact past 22 Å saturates into the top bin** (see "Saturation" caveat).
4. `z = z + self.pairwise_conditioner(z_trunk=z, token_rel_pos_feats=distogram)`
   (L110) — the distogram conditions the pairwise rep.
5. A `PairformerNoSeqModule` runs over the pair rep masked to **protein–ligand
   cross-pairs** (`cross_pair_mask`, L113-125: `mol_type==0` = receptor,
   `affinity_token_mask` = ligand; rec×lig + lig×rec + lig×lig) (L126-130).
6. `AffinityHeadsTransformer` (L142-223) mean-pools `z` over the cross-pair mask
   (`g = sum(z * cross_pair_mask) / sum(cross_pair_mask)`, L208-210) and MLPs to
   `affinity_pred_value` and `affinity_logits_binary`.

So pose enters the head **exclusively** through step 2's `cdist`. Everything else
(`s_inputs`, `z`) is fixed input.

### The trunk that produces `s_inputs` and `z` is pose-independent

In `Boltz2.forward` (`…/boltz/model/models/boltz2.py`):

- `s`, `z` come from the input embedder → MSA module → Pairformer trunk
  (L460-496), computed from sequence/MSA/template/relative-position features.
  None of these is a function of any sampled coordinate.
- The structure module's diffusion `sample()` (L532-544) produces
  `dict_out["sample_atom_coords"]` **from** `s`, `s_inputs`, `z` — i.e. the pose
  is downstream of the trunk, never upstream.
- The affinity block (L608-722) then reads, in order:
  - `best_idx = torch.argsort(dict_out["iptm"], descending=True)[0]` (L621-622) —
    the reported pose is the **highest-ipTM diffusion sample** (not a generic
    "confidence"; it is ipTM specifically).
  - `coords_affinity = dict_out["sample_atom_coords"].detach()[best_idx][None,None]`
    (L623-625) — shape `[1, 1, N, 3]`, matching the 4-D branch in affinity.py
    (L97-99). **This is the `x_pred` we swap.**
  - `s_inputs = self.input_embedder(feats, affinity=True)` (L626) — an
    affinity-specific input embedding, but still pure sequence/feature input.
  - `z_affinity = z * cross_pair_mask[None,:,:,None]` (L614-619) — the trunk pair
    rep, pre-masked to cross-pairs, passed in as `z` (`.detach()`, L631-651).
- Reported value = **ensemble of two modules** `affinity_module1`,
  `affinity_module2` averaged (L629-669), then a **molecular-weight correction**
  (L687-697): `model_coef=1.03525938`, `mw_coef=−0.59992683`, `bias=2.83288489`,
  `mw = feats["affinity_mw"][0] ** 0.3`;
  `pred = model_coef·ensemble_mean + mw_coef·mw + bias`. Applied to the ensemble
  value only; the per-module `affinity_pred_value1/2` are emitted raw.

**Consequence for the swap:** to evaluate any decoy pose for a fixed system, we
need `s_inputs`, `z` (`→ z_affinity`), and `feats` **once**, then call the two
affinity modules with `x_pred = decoy_coords`. The trunk + diffusion are the
expensive part; they are computed once and reused across all poses. The per-pose
cost is two PairformerNoSeq passes over a ≤256-token crop — negligible.

### Where the cropped input lives

`pre_affinity_<id>.npz` (`…/boltz/data/module/inferencev2.py:62-65`) is a
`StructureV2` for the affinity branch. At `__getitem__` (L206-303) it is
tokenized, then cropped by `AffinityCropper` (`…/boltz/data/crop/affinity.py`)
with `max_tokens=256, max_atoms=2048` (inferencev2.py:240-244). The cropper keeps
tokens by **ascending min-distance to the ligand** (`affinity.py:70-81`),
whole-chain-or-contiguous, capped at `max_tokens_protein=200`
(`AffinityCropper.__init__`, neighborhood_size=10). The featurizer then builds
`feats` (`compute_affinity=True`), including `token_to_rep_atom`,
`affinity_token_mask`, `mol_type`, `token_pad_mask`, `affinity_mw`, and a
`coords` tensor (the **input/crystal-frame** coordinates of the crop).

**Important subtlety:** `feats["coords"]` is *not* what the head scores. At the
affinity call site the head is fed `coords_affinity` = the **diffusion-sampled**
pose, not `feats["coords"]`. So editing coordinates inside `pre_affinity_<id>.npz`
changes the crop geometry and the distogram **only on the no-structure path**
(where `x_pred` falls back to `feats["coords"]`); on the normal path the head
ignores npz coords for scoring. This directly informs the implementation choice
below.

## Design

Two orthogonal swaps. Both hold the trunk fixed and keep the **same ligand
molecule**, so molecular weight is constant across each ladder ⇒ the MW correction
(L687-697) contributes an identical additive offset to every pose ⇒ **no MW
confound**, and any change in `affinity_pred_value` is attributable to the
distogram alone.

### Swap 1 — ligand-displacement sweep (primary)

Fix the WT protein + sequence (so trunk, `s_inputs`, `z` are the WT trunk). Take
the model's own best-ipTM WT pose as pose 0 (native-like reference), then generate
a ladder of decoy `x_pred` by rigid-body moving **only the ligand atoms** (the
tokens under `affinity_token_mask`), leaving all protein coordinates untouched:

| pose id        | construction |
|----------------|--------------|
| `native`       | model's best-ipTM WT pose (= `coords_affinity`, unmodified) |
| `exit+2`       | translate ligand COM +2 Å along the pocket-exit vector |
| `exit+5`       | +5 Å along pocket-exit vector |
| `exit+10`      | +10 Å along pocket-exit vector |
| `exit+20`      | +20 Å along pocket-exit vector (fully solvent-exposed) |
| `rand+10`      | +10 Å along a random unit direction (control for exit-axis specificity) |
| `rot_inpocket` | ligand COM held; random rigid rotation about COM (scrambled in-pocket orientation) |

- **Pocket-exit vector:** unit vector from pocket centroid (mean of receptor Cα
  within 8 Å of the native ligand) to the ligand COM, i.e. the natural "out of the
  pocket" direction. Computed from the **fixed protein** coords in `coords_affinity`.
- **Rigid translation only** (COM shift; no internal ligand distortion) keeps the
  ligand's intramolecular distogram block identical — only the protein–ligand
  cross-block of `cdist` changes. This isolates the contact geometry the head is
  supposed to read.
- The `rand+10` arm guards against the head responding to "distance along one
  privileged axis" rather than displacement magnitude.
- `rot_inpocket` probes orientational sensitivity at fixed COM (a healthy scorer
  should disfavor a scrambled binding orientation even without translation).

**Read per pose:** `affinity_pred_value` (MW-corrected ensemble; log[IC50] in µM,
lower = tighter) and `affinity_probability_binary` (= `sigmoid(logits_binary)`,
ensemble-averaged). Also log the raw per-module `affinity_pred_value1/2` to confirm
both ensemble members behave alike.

### Swap 2 — protein-swap at fixed pose (secondary; isolates sequence/contact axis)

Hold the **identical ligand pose** (the native `coords_affinity` ligand block) and
swap the **protein trunk**: run the trunk for WT, `rem`, `pack`, `inv` of the same
system, but feed every one of them the *same* ligand coordinates and the WT protein
coordinates frozen in place. This asks: with geometry held exactly constant, does
the head's affinity change when the sequence/contact identity changes?

- This is the mirror image of the CASF-mutagenesis analysis: there, mutation moved
  *both* trunk and pose; here we move only the trunk and freeze the pose.
- If Swap 1 is flat **and** Swap 2 moves, the head reads sequence/contacts but not
  geometry. If both are flat, the head is reading neither and the WT-baseline RTM
  story is the whole story. If Swap 1 moves, the head reads geometry after all and
  the empirical insensitivity is "mutations don't move the ligand enough."
- Practical note: Swap 2 requires running the trunk per variant (4×), but the
  ligand-pose coordinates are pinned. Use the WT protein frame for all variants so
  the only changing input is the sequence embedding feeding `s_inputs`/`z`.

## Predictions and decision criteria

Let `aff(d)` = MW-corrected `affinity_pred_value` at displacement `d` (Swap 1,
exit axis). Define two summary statistics per system:

- **slope/Spearman:** `ρ_spear(d, aff)` over the exit ladder
  `d ∈ {0, 2, 5, 10, 20}`; and an OLS slope `Δaff/Δd` (log units per Å).
- **native-vs-far gap:** `aff(exit+20) − aff(native)` (positive = weaker when
  ejected, the physically correct sign).

| outcome | `ρ_spear(d, aff)` | gap `aff(+20) − aff(native)` | interpretation |
|---|---|---|---|
| **Functionally pose-insensitive** (our hypothesis) | ≈ 0 (|ρ| < 0.3 median) | ≈ 0 (median < +0.5 log unit) | the head reads `x_pred` architecturally but the learned function is ~flat in pose; the CASF insensitivity is intrinsic, not a mutation-lever artifact |
| **Pose-sensitive, mutation-limited** | ≥ +0.7 (monotone weakening; native tightest) | large (median ≥ +2 log units toward `exit+20`) | the head reads geometry fine; the empirical Δ was small only because rem/pack/inv don't eject the ligand far enough |
| **Weakly pose-sensitive** | +0.3 to +0.7 | +0.5 to +2 | partial reading; reconcile with the faint `+0.12` ΔRMSD residual — this would be the same signal seen directly |

**The test is decisive in the flat case.** A 20 Å ligand ejection imposed directly
on the scored coordinates is an unambiguous non-binder geometry; if `aff` and
`P(binder)` are unmoved, "the mutations didn't move the ligand enough" is excluded
as an explanation, because we moved it ourselves, maximally, with the trunk held
fixed. The `rand+10` and `rot_inpocket` arms back this up: a head that ignores all
three (axis-translation, random translation, rotation) is reading pose in name only.

Quantify across the system panel with: median and IQR of `ρ_spear` and of the
native-vs-far gap; the fraction of systems with gap > +1 log unit; and a paired
comparison of `aff(native)` vs `aff(exit+20)` (Wilcoxon signed-rank across systems).

## Controls and reference

- **GNINA rescore of the identical pose ladder (genuine-physics reference).** A
  bona fide pose-sensitive scorer must penalize displaced decoys hard. Rescore
  every Swap-1 pose with GNINA in **score-only** mode (`gnina --score_only -r
  receptor.pdb -l <decoy_pose>.sdf --cnn_scoring <mode>`), reusing the existing
  invocation pattern in
  `analysis/casf_mutagenesis/scripts/08_run_gnina_variants.py:101-112` (swap
  re-docking for `--score_only`, since here we score a *supplied* pose, not
  search). Expect GNINA's CNNaffinity / Vina score to collapse from `native` to
  `exit+20`; the **GNINA gap minus the Boltz-2 gap** quantifies how much pose
  signal Boltz-2's head is leaving on the table. (Same receptor.pdb + the decoy
  ligand SDF per pose; box centered on the native ligand centroid as in `box.json`.)
- **MW invariance check.** Confirm `feats["affinity_mw"]` is identical across a
  ladder (it must be — same molecule). This makes the MW term a constant offset and
  rules out the "ligand-axis MW memorization" mechanism as a contributor to any
  observed Δ within a ladder.
- **Same-trunk invariant.** Assert that `s_inputs`, `z_affinity`, and `feats`
  (except the swapped `x_pred`) are byte-identical across all poses of a system —
  this is the core experimental control and is trivially enforced by the
  recommended implementation (trunk computed once).
- **Identity/round-trip control.** Feeding the unmodified `native` pose back through
  the intercept must reproduce Boltz-2's normally-reported `affinity_pred_value`
  for that system to within float tolerance. This validates the harness before any
  decoy is scored.

## Implementation

Two viable hooks; recommend (B).

### (A) Edit cropped coordinates in `pre_affinity_<id>.npz`, re-run affinity only

Rewrite the `coords` array of the `StructureV2` in `pre_affinity_<id>.npz` (ligand
atoms shifted), then re-invoke the affinity stage.

- **Rejected as primary.** Per the subtlety above, on the normal (run-structure)
  path the head scores `coords_affinity = sample_atom_coords[best_idx]`, **not**
  `feats["coords"]`. Editing the npz changes the crop/tokenization geometry and the
  *no-structure* fallback `x_pred`, but does not change what the head scores when
  the diffusion path runs. It would only work if we additionally force
  `skip_run_structure`/`--skip_structure`-style behavior so `x_pred` falls back to
  `feats["coords"]` (the `skip_run_structure` branch, boltz2.py:593-594) — fragile,
  version-specific,
  and it also re-runs the (pose-independent) trunk per pose for no benefit. Keep it
  only as a cross-check that the npz path agrees with (B).

### (B) Thin intercept: run the trunk once, call the affinity modules per pose — RECOMMENDED

A small driver that loads the Boltz-2 v2.2.1 model and, per system, executes the
forward pass up to and including the affinity block, but loops the affinity modules
over a supplied pose set. Concretely, reproduce boltz2.py:608-697 with `x_pred`
parameterized:

```
# once per system (WT trunk):
s_inputs = model.input_embedder(feats, affinity=True)          # boltz2.py:626
z_affinity = (z * cross_pair_mask[None,:,:,None]).detach()     # :614-619
native_coords = dict_out["sample_atom_coords"].detach()[best_idx][None, None]  # :623-625

# per pose (decoy x_pred), reusing the two trained ensemble modules:
for x_pred in [native_coords, *decoys]:                        # [1,1,N,3]
    a1 = model.affinity_module1(s_inputs, z_affinity, x_pred, feats, multiplicity=1)
    a2 = model.affinity_module2(s_inputs, z_affinity, x_pred, feats, multiplicity=1)
    val = (a1["affinity_pred_value"] + a2["affinity_pred_value"]) / 2
    val = 1.03525938*val - 0.59992683*(feats["affinity_mw"][0]**0.3) + 2.83288489   # MW corr
    prob = (sigmoid(a1["affinity_logits_binary"]) + sigmoid(a2["affinity_logits_binary"]))/2
```

Tensors/feats that matter (all already present in `feats` after the affinity
featurizer):

- `feats["token_to_rep_atom"]` — `[B, N_tok, N_atom]` gather matrix; `x_pred` must
  be in the **same atom index space** (`bmm(token_to_rep_atom, x_pred)`,
  affinity.py:104). Decoys are built by editing `x_pred` rows, so this stays valid.
- `feats["affinity_token_mask"]` — `[B, N_tok]` ligand-token selector. Use it to map
  ligand **tokens** back to the **atom** rows of `x_pred` to translate (via the
  token→atom association implied by `token_to_rep_atom`, or directly via the
  ligand's atom block in the crop's `mol_type`/`atom_to_token` map).
- `feats["mol_type"]` — `[B, N_tok]`; `==0` = receptor (defines `rec_mask` and the
  pocket centroid for the exit vector).
- `feats["coords"]` — the crop's input coordinates; use as a **fallback** native
  reference and to derive the pocket-exit vector if a structure pose is unavailable,
  but the scored geometry is `x_pred`.
- `feats["affinity_mw"]` — scalar per system; constant across a ladder (assert it).
- `feats["token_pad_mask"]` — masks padding tokens out of all the above.

**Why (B):** it makes the experimental control *exact* (trunk computed once,
provably identical across poses), runs the cheap part per pose, and scores precisely
what production scores (`x_pred`, not npz coords). It needs no monkey-patching of the
diffusion path. A self-contained runner under
`analysis/casf_mutagenesis/scripts/` (e.g. `18_pose_swap_affinity.py`) that imports
the trained checkpoint, builds the WT affinity `feats` via the existing
inference DataModule, and writes one row per (system, pose) is the right shape.
Note the boltz2.py affinity block wraps the modules in
`torch.autocast("cuda", enabled=False)` (:628) — replicate that (fp32) for
numerical agreement with production.

## System selection and scale

- **~10-20 diverse CASF systems × ~6 poses** (the 7-pose ladder above; drop
  `rand+10` or `rot_inpocket` if trimming). Trunk runs once per system; the
  affinity-module loop is ≪1 s/pose. Whole panel is minutes of GPU.
- Draw a deliberate mix from the existing subset20 / full-CASF runs
  (`analysis/casf_mutagenesis/outputs/`), split by the CASF-affinity behavior:
  - **"memorized" systems** — adversarial pose stays near native AND ΔAff ≈ 0
    (`docs/casf_affinity.md`: e.g. the `pack` 70%-still-binder cohort; structure
    RMSD < 2 Å on adversarial). These are where the head was *expected* to flag a
    broken pocket and didn't.
  - **"responded" systems** — adversarial pose ejected (ligand RMSD ≥ 4 Å) — the
    clearest non-binders, where the empirical ΔAff was still only ~+0.07
    (`affinity-head-pose-insensitive.md`). These are the strongest test: the model
    itself moved the ligand far, yet the affinity barely changed.
  If the head is flat on **both** cohorts under a direct 20 Å ejection, the
  pose-insensitivity is intrinsic and cohort-independent.
- Prefer systems whose WT pose is well-placed (best-ipTM WT pose near crystal) so
  pose 0 is a genuine native reference; reuse the WT-quality gate already applied in
  the mutagenesis pipeline.

## Execution steps

Assumes a working `boltzina_env` (Boltz-2 v2.2.1) on a GPU box; env setup is out of
scope for this spec (the user is configuring it).

1. **Pick the panel.** Select 10-20 systems (mix of memorized/responded per above)
   from `analysis/casf_mutagenesis/outputs/paired_affinity_*.csv` +
   `results*.csv`; record each system's `pre_affinity_<id>.npz` location and WT
   best-ipTM pose.
2. **Build the WT affinity `feats` once per system** via the existing Boltz-2
   inference DataModule (`affinity=True`), capturing `s_inputs`, `z_affinity`,
   `native_coords` (= `coords_affinity`), and `feats`.
3. **Generate the pose ladder** per system: compute pocket-exit vector from
   receptor Cα + ligand COM (both from `native_coords`); produce the 6-7 decoy
   `x_pred` tensors by translating/rotating the ligand atom rows only.
4. **Score with Boltz-2** (hook B): loop the two ensemble modules over the ladder,
   apply MW correction, record `affinity_pred_value` + `affinity_probability_binary`
   (+ raw `_value1/2`).
5. **Identity check:** confirm the `native` row reproduces the system's
   production-reported affinity to float tolerance.
6. **GNINA reference:** write each decoy ligand pose to SDF, rescore with
   `gnina --score_only` against the WT receptor (pattern from
   `08_run_gnina_variants.py`), record CNNaffinity / Vina.
7. **Aggregate:** per system compute `ρ_spear(d, aff)`, OLS slope, native-vs-far
   gap, for both Boltz-2 and GNINA; then panel-level medians/IQRs and the paired
   Wilcoxon (native vs exit+20).
8. **Figure + table** (see below).

## Expected outputs

- **Table** `pose_swap_affinity.csv` — one row per (system, pose):
  `pdbid, cohort, pose_id, displacement_A, axis, boltz_aff, boltz_prob,
  boltz_aff1, boltz_aff2, gnina_cnnaff, gnina_vina`.
- **Per-system summary** `pose_swap_summary.csv`:
  `pdbid, cohort, boltz_rho, boltz_slope, boltz_gap_native_far, gnina_rho,
  gnina_gap_native_far, gnina_minus_boltz_gap`.
- **Figure** `figures/pose_swap_affinity.png` — affinity-vs-displacement curves
  (x = displacement Å along exit axis, y = predicted log[IC50]; one faint line per
  system + bold median), Boltz-2 vs GNINA side-by-side. The expected contrast: a
  near-flat Boltz-2 panel against a steeply rising GNINA panel. A second panel for
  `P(binder)` vs displacement.

## Caveats / non-obvious failure modes

- **Distogram saturation at 22 Å (`max_dist`).** Beyond ~22 Å all protein–ligand
  cross-distances land in the top bin, so the head **cannot** distinguish `exit+20`
  from `exit+30`. This is fine for the test (a healthy head should already saturate
  to "non-binder" by 20 Å), but it means the informative dynamic range is
  ~2-20 Å — hence the ladder stops at +20. Do **not** read flatness *above* 22 Å as
  evidence of anything; the architecture is blind there by design.
- **Native ≠ crystal.** Pose 0 is the model's best-ipTM pose, not the crystal pose.
  If the WT pose is itself poor, the "native tightest" prediction is ill-posed —
  hence the WT-quality gate. Optionally add a `crystal` pose (crystal ligand mapped
  into the crop frame) as an extra reference.
- **Token↔atom mapping for the ligand shift.** The translation must move exactly the
  ligand atom rows of `x_pred` (those gathered into ligand tokens via
  `token_to_rep_atom`/`affinity_token_mask`), leaving receptor atoms fixed. Get this
  mapping wrong and you either move nothing or move the whole complex (which leaves
  all relative distances unchanged → spurious flatness). The identity check (step 5)
  plus a "translate everything rigidly → expect exactly zero Δaff" sanity pose catch
  both errors.
- **ensemble fp32 path.** Production runs the affinity modules under
  `autocast(enabled=False)` (boltz2.py:628); mismatching dtype/autocast in the hook
  will shift values and break the identity check.
- **Probability head is weakly calibrated** (WT `P(binder) ≈ 0.54`,
  `docs/casf_affinity.md`); treat `affinity_pred_value` as the primary readout and
  `P(binder)` as corroborating.
- **GNINA score-only requires a valid pose in the receptor frame.** Decoys are built
  in `coords_affinity`'s frame; ensure the SDF written for GNINA is in the same
  frame as `receptor.pdb` (or superpose), else GNINA's gap is an artifact of frame
  mismatch, not displacement.

## Out of scope

- Re-running the CASF-mutagenesis between-case Δ analysis with WT-baseline controls
  (tracked separately in `affinity-head-pose-insensitive.md` "Open").
- The ligand-chemistry axis (`results_ligand.csv`) — a pose-swap there would also
  change MW and reopen that confound; a separate design.
- Retraining / probing internal activations of the affinity head; this spec only
  varies inputs to the frozen v2.2.1 head.
- `boltzina_env` / GPU provisioning.

## Acceptance criteria

1. `pose_swap_affinity.csv` exists with all (system × pose) rows for the panel,
   both Boltz-2 and GNINA columns populated.
2. The `native`-pose identity check passes (Boltz-2 hook reproduces production
   affinity within float tolerance) for every system.
3. The MW-invariance and rigid-translate-everything sanity poses both pass.
4. `pose_swap_summary.csv` + `figures/pose_swap_affinity.png` rendered, with
   panel-level median `ρ_spear`, median native-vs-far gap, and the paired Wilcoxon
   reported for Boltz-2 and GNINA.
5. A one-line verdict recorded against the decision table: **functionally
   pose-insensitive** (flat Boltz-2 under direct ejection, with GNINA steep on the
   same poses) vs **pose-sensitive but mutation-limited** (Boltz-2 also steep).
