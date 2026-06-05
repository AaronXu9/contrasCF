# Project Spec — Contrastive Mutant-Pocket Fine-Tuning of Co-folding Models with Physics-Aware Architectural Priors

**Working title:** *CounterFold: Installing physical sensitivity into co-folding via mutant-counterfactual contrastive supervision*

**Scope:** Single ~3-month focused effort. One paper. Comparative across Boltz-2 and FLOWR.root.

**Status:** Draft v0.1 — internal planning document.

---

## 1. Elevator pitch

Co-folding models (AlphaFold3, Boltz-1/2, Chai-1, FLOWR.root) produce accurate native poses but are insensitive to physically meaningful perturbations of the binding site or ligand (Masters et al., *Nat Commun* 2025; Škrinjar et al., 2025). We hypothesize that this is a **representational** failure (the model has not learned that disrupted pockets should not bind the same ligand the same way), not a **data-volume** failure. The fix is targeted: fine-tune with **counterfactual triplets** — (WT pocket, mutant pocket, same ligand) and (pocket, WT ligand, mutant ligand) — where the WT–mutant pair shares structural scaffold but differs along a single, physics-defined axis. Pair this with an **interaction-typed attention bias** (Axis 1 physics prior) that gives the model an explicit channel for tracking H-bonds, salt bridges, and π-stacks. Compare the same protocol on Boltz-2 (production cofolder, affinity-head hook) and FLOWR.root (generative cofolder, likelihood hook) to test whether contrastive supervision is an architecture-general fix or an architecture-specific patch.

**Headline claim if the project succeeds:** *Mutant-counterfactual contrastive fine-tuning reduces trapped-pose rate on Masters-style adversarial challenges by ≥40 percentage points while preserving native-pose accuracy to within 5 pp, and this result transfers across two structurally different co-folding architectures.*

---

## 2. Background and positioning (one paragraph each)

**The diagnosis.** Masters et al. show that all four current co-folding models retain the WT-ligand pose under three classes of physically destructive perturbation (sidechain removal, packing, charge inversion) and two ligand-side perturbations (methylation, charge flip). Škrinjar et al.'s Runs N' Poses establishes that this is correlated with training-set similarity at the population level. Together: the failure is real, mechanistic, and not fixable by scaling data alone.

**Why existing fixes don't fully address this.** Boltz-1x/2x's inference-time steering potentials reduce *local* geometry failures (clashes, chirality) but do not change the model's underlying response to chemical perturbation. FLOWR.root's affinity head learns from chemically diverse ligands but not from counterfactual mutant pairs of the *same* ligand against perturbed pockets. The Boltz-2 binder/decoy classifier discriminates "is this a binder" but not "given this is a binder for the WT, does it remain a binder for the mutant." None of the released approaches install the response-to-perturbation signal directly.

**What this project adds.** A training-time intervention with two coupled components: (1) a counterfactual contrastive loss using chemistry-defined hard negatives, and (2) an architectural prior (interaction-typed attention bias) that gives the contrastive loss a place to land. Comparative across Boltz-2 and FLOWR.root to establish architecture-generality.

---

## 3. Hypotheses (falsifiable, with predictions)

**H1 — Sensitivity transfer.** Fine-tuning a co-folding model with mutant-counterfactual contrastive supervision reduces trapped-pose rate on Masters-style held-out challenges from the ~50–80% baseline reported to <30%, without dropping native-pose accuracy (CASF / PoseBusters) by more than 5 percentage points.

**H2 — Architecture-generality.** The improvement transfers across Boltz-2 and FLOWR.root when the loss is applied at the structurally-comparable hook (affinity output + structural output). If the improvement appears in only one architecture, the contrastive signal is architecture-specific, which is a publishable but smaller result.

**H3 — Physics-prior additivity.** Interaction-typed attention biases (Axis 1) independently improve perturbation sensitivity. The combination (contrastive + bias) outperforms either alone in ablation. If not additive, the bias is redundant with what the contrastive loss already installs.

**H4 — Generalization across perturbation types.** A model trained on a subset of Masters challenge types (e.g., removal + packing) generalizes to held-out perturbation types (e.g., charge inversion, methylation). Failure to generalize is the signature of shortcut learning ("lots of glycines → displace") and would be a primary abandon signal.

**H5 — The signal is in the cofolder's trunk, not only its head.** Fine-tuning the affinity head alone (frozen trunk) versus also unfreezing the last few trunk layers produces measurably different behavior on the perturbation challenges. Frozen-trunk should perform meaningfully worse if the contrastive signal needs to reshape the structural representation, not just the affinity readout. This is the diagnostic that disambiguates "the model has the physics features but doesn't use them" from "the model lacks the features entirely."

---

## 4. Method

### 4.1 Counterfactual triplet schema

Each training example is a triplet:

```
(anchor: P_WT, L, pose_native)
(positive: P_WT, L, pose ∈ near-native ensemble)
(hard_negative: P', L', δ)
```

where exactly one of {P', L'} differs from {P, L} by a physically-meaningful perturbation, and δ ∈ {0, 1, 2, 3} is the chemistry-oracle expected disruption magnitude.

A **batch** includes WT-only examples for replay (50% of batch) alongside counterfactual examples (50%), to prevent catastrophic forgetting of the native distribution.

### 4.2 Loss formulation (architecture-agnostic)

```
L_total = L_native + λ_aff · L_aff + λ_pose · L_pose + L_reg
```

**Native-pose loss `L_native`** — the model's original training loss on WT examples in the batch (FAPE-style for Boltz-2, flow-matching loss for FLOWR.root). Untouched.

**Counterfactual affinity loss `L_aff`** — margin-based on the affinity-head output `a(·)`:

```
L_aff = Σ over triplets:  max(0,  m_aff(δ) - (a(P_WT, L) - a(P', L')))
```

with margin schedule (log-µM units):

| δ | Interpretation | m_aff(δ) |
|---|---|---|
| 0 | Control mutation (no expected change) | 0.0 (no-op) |
| 1 | Minor disruption (1 H-bond or contact lost) | 0.5 |
| 2 | Moderate disruption (2–3 H-bonds, or 1 salt-bridge, or 1 π-stack lost) | 1.5 |
| 3 | Severe disruption (charge inversion, packing destruction, full pocket loss) | 3.0 |

**Counterfactual pose-divergence loss `L_pose`** — encourages predicted pose distributions to differ when chemistry says they should:

```
L_pose = Σ over triplets with δ ≥ 2:
            max(0, m_pose(δ) - W_2(pose_ensemble(P_WT, L), pose_ensemble(P', L')))
```

where `W_2` is Wasserstein-2 over the predicted pose ensembles (k=5 samples per system). Only δ ≥ 2 triplets contribute to `L_pose` because the pose-divergence signal is too noisy for marginal disruptions. Pose margin schedule: m_pose(2) = 2.0 Å, m_pose(3) = 5.0 Å.

**Regularization `L_reg`** — KL or L2 divergence between the fine-tuned model's WT output distribution and the pre-trained model's WT output distribution. Prevents drift on the native task. Coefficient λ_reg ~ 0.1 initially, tuned on val.

**Hyperparameter starting points:** λ_aff = 1.0, λ_pose = 0.5, λ_reg = 0.1. Sweep over {0.5, 1.0, 2.0} for each in the first training run.

### 4.3 Architectural prior: interaction-typed attention bias (Axis 1)

Six interaction types (PLIP-derived): H-bond (donor → acceptor), H-bond (acceptor → donor), salt-bridge, hydrophobic contact, π-stack, halogen bond. One learnable scalar parameter per type: `w_type ∈ R^6`. Initialized to zero.

For every cross-attention layer between protein and ligand tokens, inject the bias additively into attention logits:

```
attn_logits[i, j] += w_type[ interaction_type(i, j) ]   if cross-attention; else 0
```

Where `interaction_type(i, j)` is computed at **inference time per system** by a fast PLIP-style detector running on the model's current pose estimate (or, during early training, on the input template / a quick docking pose). For mutant systems, the detector returns the interaction types appropriate to the mutant chemistry — so a Lys→Gly mutation removes the "salt-bridge" label for that residue's edges, and the corresponding attention bias is removed.

**This is the architectural channel by which the contrastive loss installs physics:** without the typed bias, the model has to encode the perturbation effect via undifferentiated attention weight changes; with the typed bias, it has a discrete, interpretable handle.

### 4.4 Per-architecture instantiation

**Boltz-2 hook.**
- Updated parameters: affinity module (all weights) + last 3 PairFormer trunk layers + `w_type` (6 params).
- Frozen: MSA module, lower trunk, structure module.
- Pose ensemble: 5 diffusion samples per system at inference.
- Estimated trainable params: ~50–80M out of ~1B total.
- Affinity output: continuous, log-µM scale (per Boltz-2 paper).

**FLOWR.root hook.**
- Updated parameters: pocket encoder (last 2 layers) + ligand decoder (last 2 layers) + affinity head (all) + `w_type` (6 params).
- Frozen: lower pocket encoder, embedding layers.
- Pose ensemble: 10 flow samples per system at inference (FLOWR.root is fast enough to afford this).
- Affinity output: multi-endpoint (pIC50/Ki/Kd/EC50), use pKd head consistently.
- Alternative loss formulation available: replace `L_aff` with a likelihood-margin loss
  `L_lik = max(0, m_lik(δ) - (log p_θ(L | P_WT) - log p_θ(L' | P')))`
  computed along the flow trajectory. Run both `L_aff` and `L_lik` as variants in ablation — `L_lik` is the cleaner contrastive signal mathematically but more expensive to compute (requires likelihood evaluation along the full flow).

**Comparability constraints between Boltz-2 and FLOWR.root runs:**
- Same triplet construction (Section 5).
- Same train / val / test splits (Section 5.5).
- Same margin schedule (m_aff, m_pose).
- Same evaluation protocol (Section 6).
- Same ablation grid.
- Report parameter counts, training compute, and inference compute honestly. *Match the loss formulation and the conceptual hook; do not pretend the architectures are otherwise comparable.*

### 4.5 Out of scope for this 3-month project

- **Axis 2 (decomposed physics-shaped energy head):** would require building a learnable Lennard-Jones / Coulomb / H-bond shaped head with monotonicity constraints. Defer to v2; report as a "natural extension" in discussion.
- **Axis 3 (equivariant force-field local read):** too expensive for this timeline.
- **Multi-target training (different proteins per batch):** keep batches single-target initially to make the contrastive signal cleaner.
- **Generative ligand design** with the fine-tuned model: scope creep; the structural prediction is the deliverable.
- **τRAMD kinetic labels:** these belong to the sibling kinetic-surrogate project (separate spec). The 3-month scope deliberately excludes them.

---

## 5. Data pipeline

### 5.1 Source data

**Primary training source: PLINDER.** Curated, similarity-stratified, ~17k holo protein-ligand systems with quality filters. Use the PLINDER train split.

**Filtering:**
- Pocket has ≥3 sidechain contacts with the ligand (PLIP-detectable).
- Ligand has molecular weight 100–700 Da, ≤7 rotatable bonds (keeps mutant-ligand generation tractable).
- Structure resolution ≤ 2.5 Å.
- No covalent ligands (separate problem).
- Exclude systems where the ligand has no detectable directional interaction (pure hydrophobic binders give no signal for the methylation/charge oracle).

Expected post-filter: ~8k–12k systems. Confirm count in week 2.

**Held-out evaluation sources:**
- The 7 Masters et al. systems (CDK2/ATP, MEK1/inhibitor, glucose DH, FtsE/ATP, CYP109B4/heme, lipid transfer/palmitic acid, plus one in the SI).
- PoseBusters-V2 (for native-pose regression check).
- Runs N' Poses temporal holdout (for similarity-stratified native eval).
- A small in-house ABFE / funnel-MD validation subset (~30 systems) where ground-truth binding free-energy gaps exist.

### 5.2 Mutant generation rules

For each filtered PLINDER system, generate up to **10 mutant variants**, sampled from these categories:

| Category | Specification | Expected δ |
|---|---|---|
| Masters-removal | All sidechains within 3.5 Å of ligand → Gly | 3 |
| Masters-packing | All sidechains within 3.5 Å of ligand → Phe (preserve protein integrity check via predicted pLDDT > 70 on apo) | 3 |
| Masters-inversion | All sidechains within 3.5 Å → physicochemically opposite (charged ↔ opposite charge, polar ↔ hydrophobic, etc.) | 3 |
| Targeted-H-bond | Single H-bond donor/acceptor residue → Ala or Gly | 1–2 |
| Targeted-salt-bridge | Charged residue forming salt-bridge with ligand → Ala | 2–3 |
| Targeted-stack | Aromatic residue stacking with ligand → Ala | 2 |
| Surface-control | Surface residue not in contact with ligand → Ala | 0 |
| Conservative-control | Pocket residue → chemically similar AA (Leu↔Ile, Asp↔Glu, etc.) | 0–1 |

The two **controls** (Surface + Conservative) are crucial — they're how we test whether the model has learned a real perturbation-sensitivity signal vs. a "mutation count" shortcut. A model that drops affinity equally for all categories has learned the wrong invariance.

### 5.3 Ligand modifications

For each system, attempt up to **5 ligand-modification variants** where chemically applicable:

| Modification | Application rule | Expected δ |
|---|---|---|
| Single methylation of H-bond donor OH/NH | If ligand has free O-H or N-H forming detectable H-bond | 1 |
| Multi methylation (2–3 sites) | If ligand has ≥2 OH/NH involved in H-bonding | 2–3 |
| Charge inversion (–COO⁻ → –N⁺(CH3)₃) | If ligand has carboxylate forming salt-bridge | 3 |
| Charge neutralization (–COO⁻ → –C(=O)OCH3) | If ligand has carboxylate forming salt-bridge | 2 |
| Methyl on inert C–H | Control: methyl added to non-interacting position | 0 |

Run modifications through RDKit + Schrödinger LigPrep for protonation/tautomer correctness. Reject modified ligands that violate basic chemistry (e.g., pentavalent carbons).

### 5.4 Chemistry-based oracle for δ (model-independent)

This is the **single most important design decision in the project** (per the risk analysis in the previous turn). The δ label must NOT depend on any co-folding model's output, or we contaminate the training signal with the model we're trying to fix.

**Oracle protocol:**

1. Run PLIP on the original crystal structure of the WT system to inventory interactions: list of (residue_id, atom_id, interaction_type, ligand_atom_id) tuples.
2. For a given mutation (residue R → R'), determine which interactions involving R are lost based on chemistry rules:
   - H-bond lost if R was the donor/acceptor and R' lacks the same chemistry (e.g., Ser → Ala loses the OH).
   - Salt-bridge lost if R was a charged residue and R' is neutral or opposite charge.
   - π-stack lost if R was aromatic and R' is not.
   - Hydrophobic contact: roughly preserved between hydrophobic AAs; lost on hydrophobic → polar mutation.
3. For ligand modifications, same logic from the ligand side: which interactions does the modified ligand atom no longer support?
4. Aggregate to δ via the bin definitions in §4.2.

**Validation of the oracle:** before using oracle labels at scale, validate on a small set (~30 systems) where we run funnel-MD or ABFE for both WT and a representative mutant, and check whether δ-binned mutants give the expected ranking of ΔΔG_bind. If the rank correlation between oracle δ and free-energy ground truth on this subset is below Spearman ρ = 0.6, the oracle is too noisy and we pivot to a quantitative interaction-count oracle.

**Implementation:** PLIP (open-source) for interaction detection; Biopython + RDKit for mutation application; in-house rule table for the chemistry assignments. ~1 week of engineering.

### 5.5 Train / val / test split protocol

Three-axis stratification:

1. **Protein sequence similarity** (mmseqs2 30% sequence identity clusters). No cluster shared between train and test.
2. **Ligand similarity** (Tanimoto on Morgan fingerprints, threshold 0.4). No test ligand within 0.4 Tanimoto of any train ligand.
3. **Pocket similarity** (P2Rank descriptors or PLINDER's `pli_qcov`). Test buckets stratified by quartile of nearest-train similarity.

Splits:
- **Train:** ~80% of PLINDER post-filter (~8k systems × ~5 mutants average = ~40k counterfactual triplets, plus ~8k WT examples).
- **Val:** ~10% (~1k systems).
- **Test held-out (PLINDER-internal):** ~10% with low similarity to train (used for native-pose evaluation under matched distribution).
- **Test held-out (external):** Masters systems, PoseBusters-V2, Runs N' Poses temporal holdout. These are completely external to the training pipeline.

### 5.6 Storage and schema

Single Parquet table for triplets:

```
system_id           str    # original PDB or PLINDER identifier
variant_id          str    # WT / mutant_id
parent_system_id    str    # links mutants to WT
perturbation_kind   enum   # category from §5.2 / §5.3
delta               int    # oracle label, 0–3
protein_seq         str
protein_mutations   list   # [(pos, WT, mut), ...]
ligand_smiles       str
ligand_modifications list  # [(atom, mod_type), ...] or null
crystal_pose_path   str    # null for mutants
holo_path           str    # WT structure
interaction_list    json   # PLIP output on WT, used for oracle
split               enum   # train / val / test_internal / test_external_*
similarity_to_train json   # {seq_id, tanimoto, pli_qcov}
```

Target dataset name: **MutFold-Bench** (or pick something else; this is a placeholder).

---

## 6. Evaluation protocol

### 6.1 Three-tier evaluation

**Tier 1 — Held-out Masters challenges (primary).** Run all four challenge classes (removal, packing, inversion, ligand-modification) on the 7 Masters et al. systems. Report:
- Trapped-pose rate: fraction of systems with ligand RMSD < 2 Å after perturbation (lower = better).
- Affinity-gap recovery: median (a(WT) − a(mut)) on the test mutants, compared to a chemistry-based reference value.
- Pose-divergence: median W_2 distance between WT pose ensemble and mutant pose ensemble.

**Tier 2 — Standard cofolding benchmarks (regression check).**
- PoseBusters-V2: % PB-valid + % within 2 Å of crystal.
- Runs N' Poses similarity-stratified: % within 2 Å in low-low (low protein similarity + low pocket similarity to train) bucket — this is the genuine OOD bucket.
- CASF-2016 scoring power and ranking power for reference (acknowledge contamination, report anyway as baseline).

Critical regression criterion: PoseBusters-V2 native accuracy must not drop by more than 5 percentage points relative to the released checkpoint of each model. If it does, the loss formulation is fighting the original objective.

**Tier 3 — ABFE / funnel-MD subset (high-confidence validation).** ~30 systems where we have free-energy ground truth for WT vs. mutant. Test whether the fine-tuned model's predicted affinity-gap correlates with the ground-truth ΔΔG_bind. Target: Spearman ρ > 0.5 after fine-tuning, vs. roughly zero before.

### 6.2 Ablations (run all on both architectures)

| Ablation | Purpose |
|---|---|
| Baseline (no fine-tuning) | Establishes Masters-replicated failure rate |
| Fine-tune on WT only (no counterfactuals) | Controls for "training longer" effect |
| Counterfactual loss, no attention bias | Tests Axis 1 prior independently (H3) |
| Attention bias only, no counterfactual loss | Tests whether the prior alone is sufficient |
| Counterfactual loss + bias (full method) | Main treatment |
| Hold out perturbation type (train on removal+packing, test on inversion+methylation) | Tests H4, detects shortcut learning |
| Frozen trunk, head only | Tests H5 |
| Random-label oracle control | Sanity check: random δ should give no improvement |

### 6.3 Statistical reporting

- All comparisons report bootstrap 95% CI over systems (n_bootstrap = 1000).
- Per-system significance: paired t-test on ligand RMSD (WT vs. mutant) before and after fine-tuning.
- Architecture comparison: paired across-system test of (Boltz-2 improvement) vs. (FLOWR.root improvement).
- Multi-seed: 3 random seeds per fine-tuning run; report mean ± std. Single-seed results are not acceptable for this paper given the Masters critique itself flagged single-seed as a problem.

---

## 7. Decision criteria (pre-commit before week 9)

### Go signal (proceed to paper)
- Tier-1 trapped-pose rate drops by ≥40 pp on Masters challenges
- Tier-2 PoseBusters native accuracy drops ≤ 5 pp
- Tier-3 ABFE rank correlation Spearman ρ ≥ 0.5
- Held-out-perturbation-type test shows ≥50% of the in-distribution improvement transfers (H4 not violated)
- Result holds for **both** Boltz-2 **and** FLOWR.root, even if magnitudes differ

### Pivot signals (publishable, smaller framing)
- Improvement appears in only one architecture → narrow paper to that architecture; reframe as architecture-specific
- Improvement on Masters but degradation > 10 pp on PoseBusters → publish as "trade-off / Pareto" paper documenting the cost of perturbation sensitivity
- Improvement only in affinity head (frozen trunk wins), not in structural module → publish as diagnostic paper: "the cofolders have the physics features but cannot propagate them" — this is itself a strong finding
- Improvement only on physics-prior ablation (Axis 1 alone matches Axis 1 + contrastive) → publish as "physics priors are sufficient; supervision is the wrong intervention"

### Abandon signals
- No improvement after both architectures attempted with hyperparameter sweep
- Improvement fully explained by shortcut learning (H4 violated badly)
- Catastrophic collapse on standard benchmarks > 15 pp despite hyperparameter tuning
- Oracle validation fails (Spearman ρ < 0.6 against ABFE) AND fixing the oracle would take > 4 weeks

### Pre-commit
Write the decision criteria into the project repo's README and the preregistration document by end of week 4 — before any contrastive training begins. The point is to not let a marginal result rationalize itself into the "go" bucket post-hoc.

---

## 8. Timeline (3 months / 13 weeks)

| Week | Workstream | Deliverable |
|---|---|---|
| 1 | Environment setup. Get Boltz-2 + FLOWR.root running on local GPUs. Reproduce one published number from each (PoseBusters native accuracy). | Reproducible baseline notebooks for both models |
| 2 | Implement Masters challenge replicator on both models. Confirm the published failure rates on CDK2/ATP and one other system. | Baseline failure rates on Masters systems for both models |
| 3 | Build chemistry oracle: PLIP integration, mutation rules, δ assignment logic. | Oracle implementation + 50-system unit-test set |
| 4 | Build mutant generation pipeline + ligand modification pipeline. **Pre-commit decision criteria.** | Pipeline running end-to-end; preregistration doc |
| 5 | Run oracle validation against ABFE/funnel-MD subset (~30 systems). Decide oracle-go/oracle-pivot. | Oracle validation report |
| 6 | Generate full counterfactual dataset (MutFold-Bench v0.1). Apply splits. Compute similarity stratification. | Full Parquet dataset; train/val/test splits with stratification labels |
| 7 | Implement contrastive loss + interaction-typed attention bias for Boltz-2. Dry-run on 100 systems. | Boltz-2 fine-tuning code, dry-run results |
| 8 | Same for FLOWR.root. Dry-run on 100 systems. Full first training runs for both, small data subset. | FLOWR.root fine-tuning code; first comparable training runs |
| 9 | Hyperparameter sweep (λ_aff, λ_pose, λ_reg, m schedule). Full training on full dataset for both. | Best-config trained checkpoints for both architectures |
| 10 | Evaluation: Tier 1 (Masters), Tier 2 (PoseBusters, Runs-N'-Poses), Tier 3 (ABFE subset). | Full evaluation tables |
| 11 | Ablations: all rows in §6.2 table, on both architectures. | Ablation tables |
| 12 | Writing: paper draft + figures. Repo cleanup. Preprint upload. | bioRxiv preprint v1 |
| 13 | Buffer: address gaps surfaced in writing; additional ablations; respond to internal-review feedback. | Preprint v2 |

**Critical-path items:**
- Week 5 oracle validation: if oracle fails, the project pivots before any training compute is spent.
- Week 8 first full training run: if neither architecture trains stably, the loss formulation needs revision.
- Week 10 evaluation: the go/pivot/abandon decision happens here.

**Compute budget estimate:**
- ~4 × A100 GPU-weeks for Boltz-2 fine-tuning (3 seeds × 1 main + 2 hyperparameter variants).
- ~2 × A100 GPU-weeks for FLOWR.root fine-tuning (cheaper architecture, same multipliers).
- ~3 × A100 GPU-weeks for ablations.
- Plus inference at evaluation time: ~1 GPU-week.
- Total: ~10 A100-GPU-weeks. Feasible with one node of 8×A100 over 2 weeks, or one A100 over 10 weeks (won't fit timeline).

---

## 9. Risks and mitigations (ranked by severity)

| # | Risk | Severity | Mitigation |
|---|---|---|---|
| 1 | **Oracle calibration fails** — rule-based δ doesn't correlate with actual ΔΔG_bind | Critical | Validate on ABFE/funnel-MD subset at week 5 *before* any contrastive training. If ρ < 0.6, redesign oracle (consider quantitative interaction-count oracle with learnable weights from the validation subset) or pivot to ABFE labels only on the smaller validation set. |
| 2 | **Shortcut learning** — model learns "Gly-rich pocket → displace" rather than "interaction lost → displace" | High | Mixed-mutation challenges in training (not just Masters-style). Hold out perturbation types between train and test (§6.2). Train control mutations (δ=0) included in batches. |
| 3 | **Catastrophic forgetting** on native distribution | High | Frozen lower layers. WT replay (50% of every batch). KL regularization to original output distribution. PoseBusters as regression check at every checkpoint. |
| 4 | **Architecture-specific hooks behave asymmetrically** — Boltz-2 affinity head and FLOWR.root affinity head don't carry comparable signal | Medium | Report parameter counts and compute honestly. Add the likelihood-margin variant for FLOWR.root as a cleaner alternative. Don't oversell "matched" comparison. |
| 5 | **Compute overrun** | Medium | Pre-budget GPU-weeks. Use LoRA/PEFT adapters on trunk layers if full FT is too expensive. Reduce dataset size first if needed; keep eval set untouched. |
| 6 | **Concurrent publication scoops** | Medium | Post preprint at week 8 with dataset + baselines, before final results. Plant priority. The contrastive-mutant idea is natural enough that someone will publish something similar within 6 months. |
| 7 | **Reviewer pushback on the evaluation** | Medium | Build similarity-stratified evaluation from day 1, not as response-to-reviewer. Include comparison to physics-based docking on mutant pockets. Include the AVE-style bias check on the dataset. |
| 8 | **PoseBusters-V2 native accuracy degrades** by >5 pp | High but recoverable | Tune λ_reg upward. Reduce contrastive loss weight. Frozen-trunk variant from H5 as fallback (smaller intervention, smaller risk). |
| 9 | **FLOWR.root checkpoints aren't released or aren't fine-tuning-ready** | Medium | Check this in week 1, not week 8. If FLOWR.root weights aren't FT-able, pivot to Chai-1 (also open) as the second comparison architecture. |
| 10 | **PLINDER coverage is too narrow** (some protein families overrepresented) | Low | Use Runs N' Poses temporal holdout as primary OOD eval, not PLINDER-internal split. |

**Known landmines from prior planning** (carried over from previous conversation, not in scope but worth noting):
- DUD-E-style decoy bias: not directly relevant since this project uses mutant counterfactuals not DUD-E decoys, but if PLATE-VS is later folded in as auxiliary data, the AVE bias check becomes essential.
- Single-seed predictions: explicitly avoided by 3-seed protocol (§6.3).

---

## 10. Deliverables

**Primary publication:** one paper, target *Nature Machine Intelligence* or *ICML 2026* or *J. Chem. Inf. Model*, framed around the comparative architecture finding.

**Code release:**
- Fine-tuning recipes for both Boltz-2 and FLOWR.root.
- Mutant generation pipeline (reusable for the community).
- Oracle implementation.
- Evaluation suite (Masters challenge replicator + PoseBusters + similarity-stratified eval).

**Data release:** MutFold-Bench (the counterfactual triplet dataset), versioned on Zenodo. Train/val/test splits frozen and documented.

**Models:** Fine-tuned weights for Boltz-2 and FLOWR.root where licenses permit, on HuggingFace.

**Optional:** A small command-line tool that takes (PDB structure, mutation specification) and returns the model's perturbation-sensitivity score — useful for the community independent of fine-tuning.

---

## Appendix A: Chemistry oracle rules (formalized)

PLIP interaction types and their loss-on-mutation rules:

```
H-bond:
    residue side providing donor (Ser-OH, Thr-OH, Tyr-OH, Asn-NH2, Gln-NH2,
                                  Lys-NH3, Arg-guanidinium, His-imidazole)
    lost iff:
        mutation removes the donor atom
        (Ser/Thr → Ala/Val/Leu/Ile/Gly/Pro, Tyr → Phe/Trp/His,
         Asn/Gln → Asp/Glu/Ala, etc.)

    Similar acceptor rules (Asp, Glu, His, Asn, Gln carbonyl O, backbone O).

Salt-bridge:
    requires charged pair: (Asp/Glu) ↔ (Lys/Arg/His+).
    lost iff: charged residue → neutral OR oppositely-charged residue.

π-stack:
    requires aromatic: Phe, Tyr, Trp, His (sometimes).
    lost iff: aromatic → non-aromatic.

Hydrophobic contact:
    requires hydrophobic: Val, Leu, Ile, Met, Phe, Trp, (Ala, Pro).
    lost iff: hydrophobic → polar/charged. Roughly preserved between hydrophobic AAs.

Halogen bond:
    requires halogen on ligand and halogen-acceptor (carbonyl O, His, Met-S).
    typically not relevant for protein-side mutations. Handle on ligand side.
```

δ aggregation (counts apply per mutation event, summed over residues for multi-mutation challenges):

```
n_Hbond_lost, n_salt_lost, n_pi_lost, n_hphob_lost = count_losses(WT_interactions, mutation_spec)

if Masters-removal/packing/inversion:
    δ = 3   # by construction; the entire pocket is destroyed
else:
    severity = 2 * n_salt_lost + 2 * n_pi_lost + n_Hbond_lost + 0.5 * n_hphob_lost
    if severity < 0.5:  δ = 0
    elif severity < 1.5: δ = 1
    elif severity < 3:   δ = 2
    else:                δ = 3
```

This is a starting heuristic. Validate against ABFE in week 5; refine if needed.

---

## Appendix B: Open questions before week 1 kickoff

1. Confirm GPU access (which cluster, queue, # of A100s available end-to-end).
2. Confirm FLOWR.root weights and fine-tuning code are publicly available and FT-compatible. If not, pivot the second architecture to Chai-1 (also open-source).
3. Identify ABFE / funnel-MD subset for oracle validation. ~30 systems where ΔΔG_bind is known for at least one mutant. Candidates: thrombin alanine-scanning, T4 lysozyme mutants, PDE5 mutational series.
4. Decide on dataset versioning host (Zenodo vs. HuggingFace Datasets vs. PLINDER fork).
5. Identify co-authors and approximate contribution split before week 1.
6. Decide whether to attempt anonymized preregistration (OSF) at week 4. Strengthens the credibility of the decision criteria; adds modest overhead.

---

## Appendix C: What this spec deliberately doesn't do

- It does **not** attempt to fix the underlying training objective of co-folding models (e.g., reweighting the diffusion loss). That's a much bigger intervention.
- It does **not** introduce kinetic labels (τRAMD residence times). That's the sibling kinetic-surrogate project.
- It does **not** scale the dataset to 10⁶+ triplets. ~50k is enough to demonstrate the effect; larger scaling is a v2 problem.
- It does **not** include de-novo ligand design as a downstream task. Structural prediction is the deliverable.
- It does **not** propose a new co-folding architecture from scratch. It modifies existing checkpoints.

These exclusions are deliberate — the 3-month scope only works if the project resists the temptation to expand.
