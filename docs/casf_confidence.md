# Confidence × Affinity × RMSD on CASF pocket mutagenesis — the three heads dissociate

**Goal.** The memorization story (`casf_mutagenesis.md`, `casf_affinity.md`) is
about **point predictions**: on a broken pocket, the structure head keeps the
ligand near the WT pose and the affinity head keeps ~WT affinity. This doc asks
a different question — about the model's **self-reported uncertainty** and the
**joint structure** of its three outputs:

1. **Q1** Is the affinity invariance because the *structure* memorized (affinity
   correct-given-a-memorized-pose), or because the *affinity module* is
   intrinsically pose-insensitive?
2. **Q2** Does **confidence** drop in step with the structural response
   (mutate → ΔRMSD → Δconfidence)? If it does — even where the structure
   memorized — the model "knows it memorized," and confidence is a usable flag.
3. **Q3** How do confidence, affinity, and RMSD interrelate, *controlling for
   RMSD*?

**Headline — the three heads dissociate.** On a broken pocket Boltz-2's three
outputs give contradictory verdicts:

| head | signal | behavior on rem/pack/inv |
|---|---|---|
| **structure** | ligand RMSD | memorizes ~25–30% of adversarial cells; moves the ligand otherwise |
| **confidence** | interface `iptm` | **registers the damage** — drops, and the drop tracks how far the ligand moved (RTM-robust) |
| **affinity** | log[IC50], P(binder) | **near-blind** — change is dominated by regression-to-mean on WT tightness, not the structure; **pose-swap-confirmed blind to a 35 Å ligand ejection** |

The confidence head and affinity head respond to the same mutation
**independently** (Δ`iptm` ↔ Δaff ≈ 0). The model's uncertainty estimate carries
the perturbation signal that its affinity readout ignores.

---

## Data model and conventions

For Boltz-2, each `(system, variant)` has **5 diffusion-sample poses**. The
three signals live at different levels:

| signal | level | column(s) |
|---|---|---|
| ligand RMSD-to-crystal | **per-pose** | `ligand_rmsd_a` |
| confidence (`iptm` interface; `ptm`/`complex_plddt` global) | **per-pose** | `iptm`, `ptm`, `complex_plddt`, `confidence_score` |
| affinity (log[IC50], P(binder)) | **per-SYSTEM** | `affinity_pred_value`, `affinity_probability_binary` |

Affinity is **one value per case**, ensemble-averaged inside the affinity head and
broadcast over all 5 poses (verified: pose-0 and pose-1 carry identical
`affinity_pred_value`). So any affinity↔RMSD or affinity↔confidence analysis is
necessarily **system-level**; confidence↔RMSD can be done both system-level
(WT vs adv) and pose-level (within a case's 5 samples).

**Pose unit = best-of-5** (min ligand RMSD); top-1 (rank-0) and median are
robustness checks. **Sign conventions** (adversarial vs that system's own WT):

- `ΔRMSD = RMSD(adv) − RMSD(wt)`   →  `+` = ligand moved away from native (good physics)
- `Δconf = iptm(wt) − iptm(adv)`   →  `+` = confidence dropped on the adversarial
- `ΔAff  = logIC50(adv) − logIC50(wt)` → `+` = predicted *weaker* binding (logIC50 lower = tighter)

For "memorization rate" framing and why high adversarial RMSD is the *desired*
outcome, see `casf_mutagenesis.md`. Models with a confidence signal: **Boltz-2**
(full CASF, n≈229/variant) and **AF3+MSA** (confidence at full CASF, n=239, after
rsyncing the CARC `af3msa_summary_confidences` JSONs and ingesting them via
`18_ingest_af3msa_confidence.py`; note WT-RMSD-dependent analyses stay subset20
because the full run predicted only the adversarial structures, not WT). Only
Boltz-2 has an affinity head.

---

## Q1 — affinity invariance: it's the *insensitive-module* answer

Stratifying adversarial cells by their best-of-5 ligand RMSD:

| stratum (adv RMSD) | n | med adv RMSD | median ΔAff | median ΔP(bind) |
|---|---|---|---|---|
| **memorized** (<2 Å) | 189 | 1.07 | **+0.439** [+0.34,+0.66] | −0.191 |
| middle (2–4 Å) | 166 | 2.85 | +0.066 | −0.068 |
| **responded** (≥4 Å) | 332 | 6.33 | **+0.068** [+0.03,+0.12] | −0.057 |

The **decisive test is the responded cells**: when the ligand is ejected ≥4 Å
(the clearest possible non-binder), the affinity head barely moves
(ΔAff +0.07), versus the +3–6 log units biophysics demands. So the module does
**not** register loss of the binding pose.

The naive read of the table — "the affinity head responds *more* when the
structure memorized (+0.44)" — is a **regression-to-the-mean artifact**, not
local-contact sensing:

- Memorized systems are simply much **tighter WT binders** (median WT logIC50
  **−0.26** ≈ 0.6 µM) than responded ones (**+0.80** ≈ 6 µM) — tight, well-defined
  binders are easy to place, so they land in the low-RMSD stratum.
- Tighter WT predictions have more room to drift weaker:
  **Spearman(ΔAff, WT affinity) = −0.63** (p=1e-77) — textbook RTM.
- **Controlling for WT affinity** (partial Spearman, *not* a stratification):
  partial(ΔAff, **absolute** displacement | WT aff) = **−0.06 (n.s.)**.
  partial(ΔAff, **ΔRMSD** | WT aff) = **+0.12 (p=1e-3)** — a faint, ~1.5%-variance
  residual coupling survives.

**Conclusion:** the affinity head is **near pose-insensitive**. Its change under
mutation is governed by regression toward its prior (scaled by WT tightness),
with at most a faint structural signal, and always orders of magnitude below
physics. This resolves "memorized-structure vs inaccurate-module" in favor of
the **inaccurate/insensitive module**. (Whether the faint +0.12 is real
pose-reading is settled by the **pose-swap test** below — it is *not*: the head
ignores a ligand ejected 35 Å into solvent.)

---

## Q2 — confidence tracks the structural response (RTM-robust)

**2a — between-case.** Spearman(ΔRMSD, Δconf): **+0.45** (Boltz-2, p=8e-35),
**+0.52** (AF3+MSA). Stratified Δ`iptm`: memorized +0.012 → middle +0.046 →
responded +0.086 — **monotonic**: the more the ligand moved, the more confidence
dropped. (Treat p as meaningless at n=687; the moderate ρ is the content.)

This drop is **not** regression-to-the-mean — three independent controls:

| control | result | reading |
|---|---|---|
| partial out WT `iptm` | ρ(Δiptm, ΔRMSD): 0.445 → **0.437** | unchanged |
| does the drop persist off the ceiling? | low-WT-`iptm` half still drops **69%** (medΔ +0.046) | RTM would push low-WT *up*; they fall |
| residual of `iptm_adv ~ iptm_wt` (slope 0.80) | corr(residual, displacement) = **−0.41** (p=1e-28) | displacement drives confidence *beyond* RTM |

**AF3+MSA at full CASF (n=239).** After rsyncing the CARC structures + confidence
and re-running `05_analyze`, the AF3 arm is full-CASF and *sharper* than Boltz-2:
`iptm` falls WT 0.96 → adv 0.80–0.85 (median Δ **+0.08 / +0.10 / +0.12** for
rem/pack/inv), 82–90% of systems, **AUROC 0.83–0.88** (vs Boltz-2's 0.67–0.69).
Global `ptm` drops 3–6× less (Δ +0.02–0.03); `inv` hits hardest, matching the
structure side. **Conditional (full CASF):** even on *memorized* cells (ligand
stayed <2 Å), AF3 `iptm` still drops Δ+0.040 (AUROC 0.852, n=168) — 4× the Boltz
effect (Δ+0.011, AUROC 0.65); on *responded* cells Δ+0.160 (AUROC 0.958).
Spearman(ΔRMSD, Δconf) = +0.515 (p=7e-50). The subset20 AF3 result was no fluke.

Plus the controls from `16_confidence_response.py`: the **global** confidence
`ptm` drops 3× less than interface `iptm` (the protein still folds — Cα ~1 Å);
the **broken no-MSA AF3** baseline (ptm~0.28) shows ~zero drop (AUROC 0.54);
and the WT→adv `iptm` drop (0.030) is **3.1× the within-system noise** across
the 5 diffusion samples (SD 0.0096).

**2b — within-case (pose level).** Across a case's 5 samples,
Spearman(`iptm`, RMSD): median **−0.20** (WT), **−0.10** (adv). Confidence is a
genuine *per-pose* quality signal (it ranks lower-RMSD poses higher); on
adversarial cases it still leans toward the WT-like pose (a pose-level
memorization echo, weaker than WT).

**2c — pooled-pose AUROC** (confound-free; no WT subtraction). Confidence
separating near-native (<2 Å) from displaced (≥4 Å) poses: **0.86** (WT),
**0.80** (adv Boltz-2), **0.85** (adv AF3+MSA). Confidence is a usable
per-prediction reliability flag.

**The sharp result — "knows it memorized."** Restricting to cells where the
*structure* memorized (WT correct AND adv RMSD < 2 Å), `iptm` still drops:
Boltz-2 Δ+0.011 (p<1e-4, AUROC 0.65), AF3+MSA Δ+0.030 (AUROC 0.91, every cell).
Where the structure responded (≥4 Å) it drops ~8× more (Δ+0.088, AUROC 0.885).
So confidence registers the broken pocket even when the pose output doesn't move
— attenuated when it memorizes, strong when it responds.

---

## Q3 — joint structure, controlling for RMSD

Spearman correlation matrix (Boltz-2 adversarial cells, best-of-5):

| | RMSD | iptm | plddt | affinity | P(bind) |
|---|---|---|---|---|---|
| RMSD | 1.00 | −0.42 | −0.41 | +0.19 | −0.13 |
| iptm | | 1.00 | +0.59 | −0.15 | +0.02 |
| affinity | | | | 1.00 | −0.22 |

- **RMSD drives confidence** (−0.42) far more than it drives affinity (+0.19).
- **Confidence ↔ affinity is weak** (−0.15) and **half is just shared dependence
  on RMSD**: partial(`iptm`, affinity | RMSD) = **−0.08** (p=0.03).
- **The two heads' *responses* are decoupled:** Δ`iptm` ↔ Δaff = **−0.025
  (n.s.)**; partial|ΔRMSD = −0.067. *Getting less confident does not come with
  predicting weaker binding.*

---

## The pose-swap test — direct confirmation (capstone)

Q1–Q3 are observational and RTM-controlled; they leave one residual ambiguity
(the faint +0.12 ΔAff↔ΔRMSD). The **pose-swap test** removes it by *intervention*:
hold the trunk fixed and hand the affinity head a *supplied* pose with the ligand
rigidly ejected from the pocket — does the predicted affinity react?

**Method** (`19_pose_swap_affinity.py`; design `docs/superpowers/specs/2026-06-06-pose-swap-test-design.md`;
lab notebook `docs/lab_notebook/2026-06-06_pose-swap-test.md`). Boltz-2's affinity
head reads pose *only* via a distogram of `x_pred` over protein–ligand cross-pairs;
the trunk (`s_inputs`, `z`) is pose-independent. So: monkeypatch
`AffinityModule.forward`, capture `(s_inputs, z, x_pred, feats)` on the production
call, then re-invoke it with `x_pred` decoys where the ligand is translated beyond
the protein bounding sphere (clearance 5/15/30 Å, *verified* by the min
ligand–protein distance). The native call **is** production → identity check is
automatic; a whole-complex translation gives Δaff = 0 (translation-invariance
control). [Bug-and-fix in the ELN: the first decoy scheme translated along a
*local* pocket-exit vector that, for buried pockets, plowed the ligand *through*
the protein — caught by the min-distance diagnostic, fixed by ejecting beyond the
bounding sphere.]

**Result (n=29 CASF systems; ligand ejected to a median 35 Å from the protein — zero contacts):**

| scorer | response to a 35 Å ligand ejection | notices? |
|---|---|---|
| **Boltz-2 affinity head** | median native→eject gap = **−0.004** log units; **0/29** weaken by ≥+1; Wilcoxon p=0.27 (n.s.) | ❌ no |
| **GNINA Vina** (pure-physics reference) | **−8.7 → 0.0** kcal/mol; **100%** of systems → ~0 (all binding energy lost) | ✅ completely |
| **GNINA CNNscore** (learned pose-quality) | 0.95 → 0.54 (collapses for 21%, floors for the rest) | ⚠️ partially |

A ligand floating 35 Å in solvent — no protein contacts — is predicted by Boltz-2's
affinity head to bind **essentially as well as the native pose**. So the faint +0.12
was *not* pose-reading: the head is **intrinsically pose-insensitive**, confirmed by
direct intervention rather than correlation. The GNINA contrast sharpens the lesson:
the *physics* term (Vina) collapses completely, while the *learned* heads memorize to
different degrees — gnina's CNN partly, **Boltz-2's affinity head most extremely (flat)**.

Figures: `pose_swap_affinity.png` (flat ejection curves), `pose_swap_contrast.png`
(Boltz vs GNINA, 3-panel).

---

## Interpretation — and the link to CounterFold

The perturbation signal **is present** in the model — the interface/confidence
representation registers the broken pocket cleanly and proportionally. The
structure head acts on it most of the time (moves the ligand) but memorizes
~25–30%; the affinity head essentially ignores it (confirmed directly by the
pose-swap test — capstone above). This is direct evidence for
**CounterFold spec H5** ("the cofolders have the physics features but cannot/do
not propagate them"): the features live in the trunk/interface representation,
so the fix should **reshape the trunk**, not just retrain an output head. It
also rules out the optimistic hope that the affinity head could act as a
physics-aware reranker on top of structure memorization.

Practically: **low / dropped interface confidence is a usable flag** for
memorized or otherwise untrustworthy co-folding predictions (AUROC 0.80–0.86 at
the pose level), whereas the affinity number is not.

---

## Caveats

- **AF3+MSA is now full CASF (n=239), incl. the conditional** — marginal drop,
  conditional ("knows it memorized"), and ΔRMSD-coupling all hold at full scale
  and *stronger* than Boltz-2. (Got there by rsyncing the WT AF3 cifs from CARC
  and re-running `05_analyze` as a CARC CPU job.) *Schema note:* the CARC repo's
  `analysis.py` predates the lab-box `bestfit`/`fullca` RMSD columns, so the
  regenerated `results_full.csv` is missing those (harmless for the confidence
  analyses; original preserved as `results_full.csv.bak`). Re-run the lab box's
  newer `05_analyze` for the canonical full-schema file.
- **Single seed** (seed 42). The 5-sample diffusion spread provides a noise
  floor (used above) but is not a substitute for multi-seed.
- **Affinity is per-system**, so there is no pose-level affinity analysis; the
  faint +0.12 ΔAff↔ΔRMSD residual is the limit of what the system-level data can
  say — now settled by the **pose-swap test** (capstone §): the affinity head is
  intrinsically pose-insensitive (flat under a 35 Å ligand ejection, while a
  pure-physics scorer's Vina term collapses to 0).
- **Pocket axis only.** The ligand axis (`results_ligand.csv`, halo/meth/charge)
  has identical columns and is the next extension.
- **Regression-to-the-mean** affects every WT→adv Δ; claims here are either
  RTM-controlled (partial correlations) or rest on confound-free designs
  (within-case, pooled-pose AUROC, the broken-AF3 null).

---

## Reproduction

```bash
# Pure-stdlib confidence-drop analysis (runs on any python3):
python3 analysis/casf_mutagenesis/scripts/16_confidence_response.py

# Joint conf × affinity × RMSD (needs pandas/scipy/matplotlib/sklearn;
# on this host: /mnt/katritch_lab2/aoxu/envs/protenix/bin/python3):
PY=/mnt/katritch_lab2/aoxu/envs/protenix/bin/python3
$PY analysis/casf_mutagenesis/scripts/17_conf_aff_rmsd.py
```

## Files of interest

| path | what |
|---|---|
| `analysis/casf_mutagenesis/scripts/16_confidence_response.py` | confidence-drop test (Wilcoxon/AUROC/bootstrap, stdlib) |
| `analysis/casf_mutagenesis/scripts/17_conf_aff_rmsd.py` | joint Q1–Q3 analysis + figure |
| `analysis/casf_mutagenesis/outputs/confidence_stats_full.csv` | per (model, variant, metric) drop stats |
| `analysis/casf_mutagenesis/outputs/confidence_conditional_full.csv` | memorized-vs-responded conditional |
| `analysis/casf_mutagenesis/outputs/paired_conf_aff_rmsd.csv` | per-cell WT/adv RMSD + confidence + affinity + deltas |
| `analysis/casf_mutagenesis/outputs/q1_affinity_strata.csv` / `q2_within_case_rho.csv` / `q3_corr_matrix.csv` | per-module tables |
| `analysis/casf_mutagenesis/figures/conf_aff_rmsd_pocket.png` | 4-panel summary figure |
| `analysis/casf_mutagenesis/scripts/19_pose_swap_affinity.py` + `20`–`22` | pose-swap driver, panel aggregator, GNINA reference + contrast |
| `analysis/casf_mutagenesis/scripts/23`/`24_export_*_poses.py` | export native/ejected ligand structures (crystal + predicted frame) for visual inspection |
| `analysis/casf_mutagenesis/figures/pose_swap_affinity.png` / `pose_swap_contrast.png` | pose-swap result + Boltz-vs-GNINA contrast |
| `docs/lab_notebook/2026-06-06_pose-swap-test.md` + `docs/superpowers/specs/2026-06-06-pose-swap-test-design.md` | pose-swap lab-notebook entry + design spec |
