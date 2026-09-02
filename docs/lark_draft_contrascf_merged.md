# contrasCF — Idea, Validation, Experiments

> **Local draft for review.** Merges the main idea doc (`Kg1hdTzdBoHEy9xlxANuIKqxtOh`)
> with the 2026-08 update doc (`TVLddBn1koPfiyxh88UueCU3tNB`).
> **Precedence rule: where the two disagree, the update doc wins** — it is newer and
> every number in it was regenerated after the benchmark corrections.
> Superseded text from the main doc is *deleted*, not kept alongside, so there is one
> current source of truth.
>
> Status markers used below: ✅ verified · ⏳ pending · 🔲 placeholder (not to fill yet)

---

## 1. Idea

*(from main doc §Idea, de-blobbed — its H1 currently contains the entire essay as a title)*

**Core question.** Do deep-learning co-folding models learn the physics of protein–ligand
interaction, or do they memorise sequence/structure co-occurrence? Masters et al. 2025
(*Nat. Commun.* 16:8854) showed AF3, RoseTTAFold-AA, Chai-1 and Boltz-1 keep placing the
ligand in the native pocket even under perturbations that should abolish binding.

**Our extension.** Three moves beyond the paper:

1. **Scale** — full CASF-2016 (251 systems run; 281 receptors available)
   instead of 16 hand-built cases.
2. **Breadth** — add the physics-based/docking arm (GNINA, UniDock2, SurfDock, **ICM**),
   which the paper never tested. If memorisation is a *learned-model* pathology,
   physics-based docking should not show it.
3. **Mechanism** — if the failure is that interaction-specific encoding is not innate to
   a fused protein+ligand representation, then the fix is an explicit interaction-level
   readout. That is the trained-model arm (§4).

**Framing claim.** Interaction-specific encoding is not implicit in co-folding
representations and cannot be learned by fusing ligand and protein alone.

**References**
- Masters, Mahmoud, Lill (2025) *Nat. Commun.* 16:8854 — the reference paper.
- *Do Protein–Ligand Models Learn Binding Sites or Just Binding Likelihood?*
  arXiv 2605.24045 — same concern from the VS side.

---

## 2. Validation

*What had to be true before any result could be trusted. This section absorbs the whole
of the update doc's corrections — they are validation findings, not results.*

### 2.1 Benchmark and denominator ✅

CASF-2016 core = 285 systems. 34 had no receptor under `$CASF/raw/` (HiQBind coverage
gap); **30 recovered from RCSB deposited mmCIF → `raw/` now holds 281.**

Safety argued by measurement, not assertion: all 31 non-peptide ids are absent from
HiQBind's 32,275-row metadata *entirely* (a coverage gap, not a quality veto), and running
the identical recovery on **25 systems HiQBind does have** reproduced the detected 3.5 Å
pocket **exactly on 25/25** (Jaccard 1.00). `build_system()` succeeds 30/30.

Still excluded (4): `1a30`/`3bv9`/`3uri` (peptide ligands — categorical) and `3f3a`.

> ⚠️ **The 281 is receptors, not experiments.** The recovery expanded `raw/` to 281, but
> variants and predictions were never generated for the 30 new systems — verified
> 2026-09-02: they have no output directories on katlab *or* CARC (CARC holds 251 system
> dirs; 0 prediction CIFs newer than 2026-08-20). They do appear in `results_full.csv`,
> but all 360 of their rows are `missing_cif` with no usable RMSD.
>
> **Report it as: 281 receptors available · 251 systems in the current experiment ·
> 239 with mutant docking inputs.** Quoting a bare n=281 overstates what has been run.

### 2.2 Data preparation and the data store ✅

- **Ground truth** = the deposited crystal pose (`crystal_ligands/<id>_ligand.sdf`), *not*
  the re-folded `wt` prediction.
- `crystal_ligands/` is **corpus-wide (14,661 SDFs)**, not CASF-scoped — the `$CASF`
  variable names the PDBbind clean-split root. Scope ladder:
  18,623 → 16,491 → 14,661 → 285 (CASF core) → 281 (with receptor).
- Full lab↔CARC map with the 6 host divergences: `docs/data_store_map.md`,
  machine-verified by `env/verify_data_store_map.sh`.

### 2.3 Integrity corrections ✅

**(a) AF3+MSA emitted a single chain.** The runner rebuilt a one-chain input, so for
multi-chain systems the predicted CIF came back with one chain — sometimes *not the chain
the mutations were on*. Archetype `1bcu`: chains L(26)+H(257), both mutations on H, but the
CIF held only a 26-residue chain, so that cell docked into a stub with **no mutation at all**.

| mutation-presence verdict (690 cells) | before | after |
|---|---|---|
| OK | 601 | **684** |
| PARTIAL | 65 | 0 |
| NOCHAIN / ABSENT | 18 | 0 |
| NOMUT (spec == WT) | 6 | 6 |

Truncated receptors **16 systems → 1**. Effect on AF3+MSA's own rates was material and in
the *conservative* direction — it had been **understating** memorisation.

**(b) SurfDock: retracted, then restored.** The published claim — *0/244 WT under 2 Å,
"CASF sits outside SurfDock's training distribution"* — was an artifact of our own surface
preprocessing skipping SurfDock's interface crop (we kept every mesh face: 1474 vertices
vs SurfDock's own 142 on the same pocket). **SurfDock actually reaches 87.6% WT under 2 Å
at 1.06 Å median.** Falsified three ways, incl. `1a0q` — SurfDock's *own* shipped test
system — failing identically through our pipeline (7.74 Å → 0.44 Å fixed).

**(c) Two earlier explanations were wrong — do not re-derive.** The numpy `(1,3)` vs `(3,)`
shape bug is a **no-op** (broadcasts identically); `--ligand_to_pocket_center` is
**harmful**, not the fix (it replaces the trained translational prior with a delta).

### 2.4 Metrics and failure handling ✅

Full write-up: `docs/rmsd_and_failure_handling.md`. Three things a reader must know:

1. **`ligand_rmsd_a`** = Cα-superpose on the **ligand-nearest chain**, apply R/t to the
   ligand, then symmetry-minimised heavy-atom RMSD. Nearest-chain matters: superposing on
   the largest chain scored `4w9l` WT at 47 Å on a fold that was actually correct.
2. **Failed cells are dropped, not counted as misses.** Every rate is
   *P(RMSD < 2 Å | the prediction succeeded)*, and each method has its own denominator.
3. **WT-conditioning** (§3.4) is the primary aggregation. An unconditioned rate credits a
   method for "responding" on systems it cannot solve at all.

> ⏳ **Two known metric defects, both flattering the headline claim.** (i) The docking arm
> takes a single substructure match rather than minimising over symmetry-equivalent ones —
> measured understatement ≈ 2.4 points. (ii) `ligand_rmsd_bestfit` is 100% NaN for
> `atp_charge` and 83% for `glucose`, because `GetBestRMS` cannot match ligands whose graph
> changed. Tracked as TODO 10/11.

---

## 3. Experiments

### 3.1 Pocket-mutation setup ✅

Binding-site residues (3.5 Å shell) are mutated in three ways, following Masters et al.:

| variant | rule | intent |
|---|---|---|
| `rem` | all selected residues → **GLY** | remove the side-chain chemistry entirely |
| `pack` | → **PHE** | fill the pocket with bulk |
| `inv` | → **Miyata-inverted** residue | invert physicochemical character |

### 3.2 Ligand-mutation setup ✅

| family | rule | note |
|---|---|---|
| methylation | add 1–5 methyls | steric perturbation |
| halogenation | F / Cl / Br swap | electronic + steric |
| charge swap | anionic triphosphate → **neutral alkyl** *or* → **cationic ammonium** | see below |

> **The two charge ladders are not "negative vs positive".** Both amputate the *same*
> anionic triphosphate at the same bond and differ only in the grafted tail:
> `chrg→neutral` removes the charge, `chrg→flipped` reverses it. The WT is the anionic
> one. Rungs 2 and 3 are heavy-atom matched (9 and 13), so the neutral ladder is the
> **isosteric control** for the cationic one. Empirically the two are statistically
> indistinguishable (GNINA p=0.52, UniDock2 p=0.84) — but note neither isolates charge
> alone, since both delete the whole polyphosphate.

### 3.3 Methods, setup and filtering stats

**Experiment set: 251 systems × 4 variants = 1004 cells per method.** Coverage differs by
method, and the reasons are now characterised rather than left as bare n's.

| method | kind | wt | rem | pack | inv | why short of 251 |
|---|---|---|---|---|---|---|
| GNINA | docking | 251 | 239 | 239 | 239 | mutant cells need an AF3 mutant CIF (251→239) |
| UniDock2 | docking | 251 | 239 | 239 | 239 | same |
| ICM | docking | 251 | 232 | 233 | 232 | same, plus 18–19 cells not produced |
| SurfDock | docking | 249 | 238 | 238 | **231** | 12 documented exclusions, variant-linked (7/8 unexplained are `inv`) |
| AF3+MSA | co-folding | 238 | 239 | 239 | 239 | 12–13 systems `missing_cif` |
| Boltz-2 | co-folding | 229 | 229 | 229 | 229 | **the same 22 systems** `missing_cif` throughout |
| AF3 (no MSA) | co-folding | 19 | 19 | 19 | 19 | only ever run on subset20 |

**Two distinct failure modes, and neither is an analysis failure:**

- **Co-folding — 100% `missing_cif`, i.e. the prediction was never generated.** It is
  *system-level*: the same systems are absent from all four variants, so it hits WT and
  mutant equally and does not bias the WT→mutant contrast. These are **re-runnable**.
- **Docking — zero analysis failures.** Every produced cell scores `ok`. WT is essentially
  complete; the deficit is *mutant-only*, because mutant docking inputs depend on the AF3
  mutant CIF existing.

**ICM (new).** Run on CARC 2026-09-01 against `docking/receptor_aligned.pdb` — the
corrected multi-chain receptors (2026-08-25) pre-transformed into the **crystal frame**
(verified: mean Cα distance to crystal 0.80 Å / 1.61 Å, vs 65.95 / 29.75 Å untransformed),
so poses need no superposition. Scored with the *same* matcher as the other engines
(`30_analyze_icm.py`).

- **948/948 cells scored, zero failures.**
- **Parse caveat:** 228 cells (24%) fail strict RDKit sanitisation —
  `AtomValenceException: Explicit valence for atom # 12 P, 7`. Recovered by skipping only
  the valence check; heavy-atom counts still match the crystal. `parse_mode` records the
  path per cell. Without it ICM's WT would read 0.651 on n=189 instead of 0.685 on n=251.
- **Cross-check:** our RMSD tracks ICM's own self-reported value (WT 0.685 vs 0.697).
- ⏳ 251 cells lack a `FINISHED` marker — all WT. They scored fine; clean termination
  unconfirmed.

**Per-system distributions, not just rates.** `29_plot_paired_rmsd.py` renders WT-vs-mutant
scatter per method (`paired_rmsd_{rem,pack,inv}.png`): each point one system, the horizontal
line the 2 Å memorisation threshold, and a grey band marking WT failures — which is exactly
what WT-conditioning discards. ⏳ ICM still to be added as a sixth column.

### 3.4 Cross-method result ✅

#### 3.4a The fair comparison — a common system set (lead with this)

Every method in a per-method table has its own denominator, and co-folding is additionally
credited on systems where it never got the fold right. Docking is *handed* a receptor;
co-folding must *predict* one, so a Cα filter is only definable for co-folding — a
Cα-conditioned table is therefore still not symmetric across arms.

The symmetric construction is a common-system intersection:

- **Set A** — systems all six methods produced for all four variants → **n = 209**
- **Set B** — A, plus **both** co-folding models solved WT (ligand < 2 Å **and** Cα < 2 Å)
  → **n = 96**. In B every method is scored on identical systems and the co-folding models
  are WT-correct by construction, so neither "different denominators" nor "credited on
  systems it cannot solve" applies to anyone.

![matched comparison](../analysis/casf_mutagenesis/figures/matched_comparison.png)

| method | kind | A (n=209) | **B (n=96)** |
|---|---|---|---|
| SurfDock | docking | 0.037 | **0.035** |
| UniDock2 | docking | 0.085 | **0.101** |
| ICM | docking | 0.116 | **0.156** |
| GNINA | docking | 0.131 | **0.163** |
| Boltz-2 | co-folding | 0.219 | **0.333** |
| AF3+MSA | co-folding | 0.325 | **0.420** |

*(adversarial retention = fraction still within 2 Å after the pocket is destroyed; lower =
less memorisation)*

> **Tightening from A to B raises co-folding retention and leaves docking flat.** The
> separation *widens* under the stricter, fairer test — docking **0.035–0.163** vs
> co-folding **0.333–0.420**. That is the opposite of what a critic would predict if the
> effect were an artifact of unequal denominators, which is why this should lead the
> section. Cost: n=96, and that should be stated plainly — it is the price of full symmetry.

#### 3.4b Per-method WT-conditioned rates (larger n, own denominators)

![cross-method conditioned](../analysis/casf_mutagenesis/figures/crossmethod_conditioned.png)

Retention = of the systems a method solved at WT, the fraction still within 2 Å afterwards.

| method | kind | WT ceiling | WT-correct | rem | pack | inv | **retention** |
|---|---|---|---|---|---|---|---|
| SurfDock | docking | 0.876 | 218/249 | 0.024 | 0.043 | 0.030 | **0.032** |
| UniDock2 | docking | 0.578 | 145/251 | 0.141 | 0.141 | 0.104 | **0.129** |
| ICM | docking | 0.685 | 172/251 | 0.133 | 0.170 | 0.114 | **0.139** |
| GNINA | docking | 0.729 | 183/251 | 0.180 | 0.157 | 0.151 | **0.163** |
| Boltz-2 | co-folding | 0.594 | 136/229 | 0.368 | 0.360 | 0.257 | **0.328** |
| AF3+MSA | co-folding | 0.824 | 196/238 | 0.439 | 0.367 | 0.301 | **0.369** |
| AF3 (no MSA) | co-folding | 0.000 | 0/19 | — | — | — | *undefined* |

Adding a protein-structure filter to the co-folding arm alone barely moves it — Boltz-2
0.328 → **0.316**, AF3+MSA 0.369 → **0.368** (WT-correct 136→118 and 196→182). Co-folding
usually gets the fold right when it gets the ligand right (WT Cα < 2 Å on 77% / 87%,
medians 0.75 / 0.60 Å), so the filter removes few cells and they are not preferentially
memorisers.

**The claim.** SurfDock's WT ceiling (0.876) essentially matches AF3+MSA's (0.824) but its
retention is ~10× lower — co-folding-grade accuracy on native structures **does not
require** adversarial retention. **ICM strengthens this specifically:** a mature, purely
physics-based commercial docker, the cleanest physics baseline in the set, behaving like
the other dockers rather than like the co-folding models.

> ⚠️ **Do not use the "gap" column from the update doc.** It is not computed consistently:
> four rows use `1 − retention`, Boltz-2 alone uses `WT_uncond − retention`. Computed
> consistently Boltz-2's gap is **0.672**, not 0.266, which reorders it above AF3+MSA. The
> qualitative conclusion is unaffected. **Report retention, not gap** — gap folds the WT
> ceiling into a memorisation score, which is also why ICM's gap (+0.546) misleadingly
> ranks it below AF3+MSA despite far lower retention.

### 3.5 Deep dives: the three heads dissociate ✅

*(from main doc; Q1 numbers below are the update doc's re-checked, WT-conditioned values)*

- **Q1 — affinity invariance.** Responded-stratum ΔAff is **+0.229** conditioned
  (published: +0.068) — 3.4× larger, because nearly half those cells were systems Boltz-2
  cannot solve. But +0.229 is still **15–25× below** the +3–6 log units biophysics demands.
  *The affinity head remains near pose-insensitive.*
- **Q2 — confidence tracks the structural response** (RTM-robust).
- **Q3 — joint structure, controlling for RMSD.**
- **Pose-swap test.** The strongest evidence, and **immune to the WT-conditioning
  confound**: it holds the trunk byte-identical and ejects the ligand 35 Å on fixed
  systems, so WT-correctness never enters. **Δ = −0.004, p = 0.27.** Lean on this rather
  than the strata.

---

## 4. Trained model 🔲

*Placeholder — do not fill yet.*

Intended content: the interaction-recovery head — an interaction-level readout that is
pose-trapping-robust by construction (removing a side chain zeroes its interaction
regardless of where the ligand sits). Standalone results exist (crystal recover-AUROC
0.955; transfer to co-folded structures Δlost-rate +0.599 pre-registered / +0.942
pose-free with zero false positives), but co-folding integration is not done and results
are single-seed. To be written once that arm is decision-grade.

---

## 5. Open threads

**From the update doc**
- [ ] 8 of SurfDock's 12 exclusions fail with "0 graphs", cause unestablished; 7 of 8 are `inv`.
- [ ] SurfDock is single-seed (seed 42, 40 samples/complex). The WT/adversarial gap is far
      too large to be seed noise, but adversarial rates (0.026 vs 0.046) should not be
      compared to each other without ≥3 seeds.
- [ ] SurfDock has no ligand-mutation arm.
- [ ] Docking WT receptors come from crystal, mutants from AF3-predicted CIFs — the
      WT→mutant drop still conflates mutation with a crystal→predicted change.

**New / carried from the data-prep TODO**
- [ ] Fix or drop the inconsistent "gap" column (§3.4).
- [ ] **Generate variants + predictions for the 30 recovered systems** so the experiment
      actually reaches 281 (currently receptors 281 / experiment 251).
- [ ] Re-run the 22 Boltz-2 and 12–13 AF3+MSA `missing_cif` systems — these are incomplete
      GPU runs, not intrinsic failures, and would lift both co-folding denominators.
- [ ] Add ICM as a sixth column to `29_plot_paired_rmsd.py`.
- [ ] Symmetry-correct the docking matcher (≈ +2.4 points) — TODO 10.
- [ ] `ligand_rmsd_bestfit` NaN for modified ligands — TODO 11.
- [ ] ICM: confirm the 251 WT cells without `FINISHED`.
- [ ] `3mss`/`4eo8` mutant specs identical to WT (impact ≤0.003, but fix the generator).

**Parked ideas from the main doc** (kept, not deleted): τRAMD / funnel-metadynamics layers;
decoy-pocket and pocket-redirection tests ("the paper tests perturbations that destroy
binding but never ones that create or redirect it"); energy minimisation post-prediction;
pocket-based vs blind docking; testing GPCR/kinase targets.
