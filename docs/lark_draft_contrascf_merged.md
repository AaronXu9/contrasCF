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

1. **Scale** — full CASF-2016 (281 systems) instead of 16 hand-built cases.
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

> ⚠️ Anywhere the old doc says **n=251**, the current number is **281**. Results computed
> before 2026-08-23 still use 251.

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

| method | kind | cells | completion | notes |
|---|---|---|---|---|
| Boltz-2 | co-folding | 229 WT-scored | ✅ | affinity + confidence heads available |
| AF3+MSA | co-folding | 238 | ✅ | multi-chain fix + full re-run 2026-08-24; multi-seed confirmed |
| AF3 (no MSA) | co-folding | 19 | ⚠️ | **0/19 WT-correct** — cannot enter a conditioned comparison |
| GNINA | docking | 251 | ✅ | re-run 2026-08-25 on rebuilt receptors |
| UniDock2 | docking | 251 | ✅ | re-run 2026-08-25 on rebuilt receptors |
| SurfDock | docking | 249 | ✅ | 956/968 cells; 12 excluded (1.24%) |
| **ICM** | docking | **251** | ✅ | **new** — scored 2026-09-02, see below |

**ICM (new).** Run on CARC 2026-09-01 against `docking/receptor_aligned.pdb` — the
corrected multi-chain receptors (2026-08-25) pre-transformed into the **crystal frame**
(verified: mean Cα distance to crystal 0.80 Å / 1.61 Å, vs 65.95 / 29.75 Å untransformed),
so poses need no superposition. Scored with the *same* matcher as the other engines
(`30_analyze_icm.py`) for comparability.

- **948/948 cells scored, zero failures.**
- **Parse caveat:** 228 cells (24%) fail strict RDKit sanitisation —
  `AtomValenceException: Explicit valence for atom # 12 P, 7` — ICM writes phosphorus
  valences RDKit rejects. Recovered by skipping only the valence check; heavy-atom counts
  still match the crystal. The `parse_mode` column records which path each cell took.
- **Cross-check:** our RMSD agrees with ICM's own self-reported value
  (WT 0.685 vs 0.697; inv 0.091 vs 0.091).
- ⏳ **251 cells lack a `FINISHED` marker** — all the WT ones. They scored fine, but clean
  termination is unconfirmed; worth asking whether WT ran under a different script.

### 3.4 Cross-method result (WT-conditioned) ✅

![cross-method conditioned](../analysis/casf_mutagenesis/figures/crossmethod_conditioned.png)

**Retention** = of the systems a method solved at WT, the fraction *still* within 2 Å after
the pocket is destroyed. Lower = less memorisation.

| method | kind | WT ceiling | WT-correct | rem | pack | inv | **retention** |
|---|---|---|---|---|---|---|---|
| SurfDock | docking | 0.876 | 218/249 | 0.024 | 0.043 | 0.030 | **0.032** |
| UniDock2 | docking | 0.578 | 145/251 | 0.141 | 0.141 | 0.104 | **0.129** |
| **ICM** | docking | 0.685 | 172/251 | 0.133 | 0.170 | 0.114 | **0.139** |
| GNINA | docking | 0.729 | 183/251 | 0.180 | 0.157 | 0.151 | **0.163** |
| Boltz-2 | co-folding | 0.594 | 136/229 | 0.368 | 0.360 | 0.257 | **0.328** |
| AF3+MSA | co-folding | 0.824 | 196/238 | 0.439 | 0.367 | 0.301 | **0.369** |
| AF3 (no MSA) | co-folding | 0.000 | 0/19 | — | — | — | *undefined* |

**The claim this supports.** SurfDock's WT ceiling (0.876) essentially matches AF3+MSA's
(0.824), but its retention is **~10× lower**. So co-folding-grade accuracy on native
structures **does not require** adversarial retention — the retention is not the price of
the accuracy. All four docking engines cluster at 0.03–0.16; both co-folding models sit at
0.33–0.37.

**ICM strengthens this specifically.** It is a mature, purely physics-based commercial
docker — the cleanest "physics baseline" in the set — and it behaves like the other
dockers (retention 0.139), not like the co-folding models.

> ⚠️ **Do not use the "gap" column from the update doc.** It is not computed consistently:
> four rows use `1 − retention`, Boltz-2 alone uses `WT_uncond − retention`. Computed
> consistently, Boltz-2's gap is **0.672**, not 0.266 — which reorders it *above* AF3+MSA.
> The qualitative conclusion is unaffected, but the number should be fixed or dropped.
> **Recommendation: report retention, not gap** — gap folds the WT ceiling into what is
> supposed to be a memorisation score, which is also why ICM's gap (+0.546) misleadingly
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
- [ ] Symmetry-correct the docking matcher (≈ +2.4 points) — TODO 10.
- [ ] `ligand_rmsd_bestfit` NaN for modified ligands — TODO 11.
- [ ] ICM: confirm the 251 WT cells without `FINISHED`.
- [ ] `3mss`/`4eo8` mutant specs identical to WT (impact ≤0.003, but fix the generator).

**Parked ideas from the main doc** (kept, not deleted): τRAMD / funnel-metadynamics layers;
decoy-pocket and pocket-redirection tests ("the paper tests perturbations that destroy
binding but never ones that create or redirect it"); energy minimisation post-prediction;
pocket-based vs blind docking; testing GPCR/kinase targets.
