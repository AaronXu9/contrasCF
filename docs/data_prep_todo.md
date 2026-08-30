# Data-prep TODO — `pdbbind_cleansplit` / CASF-2016

Open issues found while auditing the benchmark data store. Items 2-3 are latent
traps or documentation problems that have **not** corrupted a published result;
**item 1 is now resolved** (30 receptors recovered, `raw/` 251 → 281).
**Item 4 was the correctness issue affecting existing GNINA/UniDock2/SurfDock
mutant numbers; 4b is now fixed and the cells rebuilt.** Verified 2026-08-20 on
katlab + CARC; reviewed and extended 2026-08-23.

Companion: [`data_store_map.md`](data_store_map.md) — the corrected data-store
map, machine-verified by `env/verify_data_store_map.sh` (**40/40 PASS**, 2026-08-23).

## Status at a glance

| # | issue | severity | state |
|---|---|---|---|
| 1 | 34 CASF-core systems have no receptor (HiQBind gap) | medium | **DONE** — 30 recovered from RCSB, `raw/` 251 → **281**; 3 peptide-ligand + `3f3a` excluded |
| 2 | `clusters_casf2016.json` O→0 typo (`105b`/`10wh`) | low now | open — 2-char fix not applied |
| 3 | `crystal_ligands/` corpus-wide; `$CASF` misnames the root | docs | **partly done** — map written; notebook+memory stale |
| 4 | mutant docking receptors: AF3 provenance (4a, mild) + lost mutations (4b) | medium | **DONE (4b)** — guard added; cells rebuilt on the AF3 fix: 684/690 OK, truncation 16 → 1 (`2vw5`). 4a provenance still open |
| 4d | `3mss`/`4eo8` mutant spec == WT (generation no-op) | low | **measured** — 2 systems only; rate impact ≤0.003; generator guard still open |
| 4e | AF3+MSA output 1 chain for all 52 multi-chain systems | was HIGH | **DONE** — root cause fixed (`06_run…:77,246`), 4 variants re-run, multi-seed confirmed (inv 0.529±0.020) |
| 4f | 10 unparseable `boltz.yaml` files | low | **new** |
| 5 | 2026-08-14 map claimed lab `crystal_ligands` are symlinks | low | **partly done** — map fixed; notebook+memory stale |
| 6 | figure relabel + stale committed figure | medium | **DONE** — relabelled, re-rendered, committed |
| 7 | uncommitted work; 2026-08-14 notebook entries untracked | medium | **partly done** — committed on `fix/af3-multichain-and-data-audit` (not pushed); notebook entries still untracked + still contain the 2 disproved claims |
| 8 | engine runners: documented idempotency ≠ implemented | low | **new** |
| 9 | SurfDock not re-run — mixed-provenance | medium | **DONE (other session)** — interface-crop bug fixed + re-run 2026-08-26; WT 0.876 (was ~0.001), mutants 0.026–0.046; figure caveat retracted |
| 10 | docking RMSD not symmetry-corrected — rate understated ~2.4 pts | medium | **new** — measured; see `rmsd_and_failure_handling.md` |
| 11 | `ligand_rmsd_bestfit` is 100% NaN for atp_charge, 83% for glucose (GetBestRMS cannot match modified ligands) | medium | **new** — measured |

**Roots**

```bash
LAB_CASF=/home/aoxu/projects/VLS-Benchmark-Dataset/data/pdbbind_cleansplit
CARC_CASF=/project2/katritch_223/aoxu/projects/VLS-Benchmark-Dataset/data/pdbbind_cleansplit
LAB_HIQBIND=/mnt/katritch_lab2/aoxu/data/hiqbind
```

---

## 1. Missing receptors — RESOLVED: 30 of 34 recovered (251 → 281)

**Status:** done (2026-08-23) · **Severity:** was medium (defines the denominator)

34 of the 285 CASF-2016 core ids had no entry under `$CASF/raw/`, because our
receptor source — HiQBind — does not carry them. **30 have now been recovered
from the RCSB deposited mmCIF.** `raw/` holds **281**.

Recovery script: `analysis/casf_mutagenesis/scripts/25_recover_missing_receptors.py`
(`--controls N` to validate, `--recover` to write). Provenance for every
recovered system: `analysis/casf_mutagenesis/recovered_receptors.json`.

### 1a. Why this was safe — the curation question, measured not asserted

The original worry was that HiQBind ships `*_protein_refined.pdb` (PDBBind-Opt:
hydrogens added, missing atoms/residues rebuilt) while an RCSB-derived receptor
is unrefined, so recovering these would mix curation standards mid-benchmark.

Two findings retire that objection:

1. **HiQBind never assessed them.** All 31 non-peptide ids are absent from the
   32,275-row HiQBind metadata entirely. This is a coverage gap, not a quality
   veto being overridden.
2. **The refinement difference does not touch what the pipeline uses.** Running
   the identical recovery on **25 systems HiQBind does have** and diffing
   against the refined receptor: residue counts differ by 0–15 (refinement
   rebuilds missing residues), but the detected 3.5 Å pocket is **EXACT on
   25/25** (Jaccard 1.00, every system). Since pocket detection is the only
   thing `build_system` derives from the receptor geometry, unrefined and
   refined are interchangeable for this benchmark.
   Report: `analysis/casf_mutagenesis/_rcsb_cache/controls_report.json`.

Construction mirrors HiQBind's own protein definition — polymer chains with a
heavy atom within 10 Å of the ligand — so chain selection is identical to the
251. Waters, hydrogens, alternate conformations and non-polymer residues are
dropped; the receptor is anchored to the *same* ligand copy the pipeline will
use (matched on centroid; all 30 matched at **0.0 Å**).

**End-to-end check:** `build_system()` succeeds on **30/30** with zero warnings.
Pocket sizes (median 4, range 1–11) sit inside the pre-existing distribution
(251 systems: median 5, range 0–13, with 8 systems at ≤1 residue), so the
single-residue pocket on `4owm` is normal rather than a defect.

### 1b. The 4 not recovered

| id | reason |
|---|---|
| `1a30` | peptide ligand (PDBbind `3-mer`) |
| `3bv9` | peptide ligand (PDBbind `6-mer`) |
| `3uri` | peptide ligand (PDBbind `8-mer`); **zero** non-polymer entities in the entry |
| `3f3a` | reference ligand SDF matches **no** deposited copy — see below |

The three peptide cases are a **categorical** exclusion, not a coverage gap:
the ligand machinery here (SMILES cache, RDKit MCS matching, the
ligand-mutagenesis rules) assumes small molecules. HiQBind filing two of them
under `raw_data_hiq_poly` was the same categorical call.

- [ ] `[open]` **Scope call for the user — `3f3a`.** Its ligand is a free
      tryptophan and the entry contains **four** TRP copies (A601–A604, 15
      heavy atoms each). The reference SDF has **14** heavy atoms (no OXT) and
      its centroid is **12.71 Å** from the nearest copy, so it corresponds to
      none of them. Recovering it needs a deliberate choice of which copy is
      "the" binding site; burial ranks them A601 (154 protein atoms within
      4.5 Å) > A602 (101) > A603 (79) > A604 (57). Recommendation: recover
      using **A601** and re-fetch its pose via ModelServer
      (`auth_comp_id=TRP&auth_seq_id=601`), or leave excluded at 281.

### 1c. Side finding — `1c5z` had an idealised ligand, and 167 others do too

`1c5z`'s `crystal_ligands` SDF was an **idealised CCD template**, not a crystal
pose: origin-centred coordinates, useless as an RMSD reference. Cause:
`download_crystal_ligands.py` queries RCSB ModelServer by PDBbind's
`ligand_name`, which for this entry is `BAM` while the deposited het code is
`BEN`; the query missed and it silently fell back to `BAM_ideal.sdf`.

Repaired by re-querying with the deposited code (`BEN`); the original is backed
up at `analysis/casf_mutagenesis/_rcsb_cache/ideal_ligand_backup/`. The
recovery script now detects and repairs this class automatically.

**168 entries corpus-wide are `source: "ideal"`. ZERO of them were among the
pre-existing 251** — so no published result is affected. But any future
expansion of `raw/` must screen for it; the script does.

- [ ] `[confirmed]` Fix `download_crystal_ligands.py` upstream in
      VLS-Benchmark-Dataset to fall back to the *deposited* het code before
      falling back to ideal coordinates, so the 168 can be repaired at source.

### 1d. CARC is now out of sync

CARC `raw/` has 285 dirs of which 34 are empty; katlab now has 281 real ones.
The 30 recovered receptors exist **only on katlab**.

- [ ] `[confirmed]` Move the `mkdir` inside the source-exists check at
      `env/setup_casf_data_carc.sh:51` so the sync stops minting empty dirs,
      then push the 30 recovered receptors to CARC and drop the 4 remaining
      empty dirs.
- [ ] `[open]` Stale comment at `setup_casf_data_carc.sh:58-61` cites
      `2vvn, 3arq, 2zy1` as "real-file, not symlink" cases needing the rsync
      second pass. All three were in the missing-34 and are now RCSB-recovered
      real files on katlab — so that pass now genuinely applies to them, but for
      a different reason than the comment states. Re-derive or drop the example.

### 1e. The verifier was too narrow — extended, and a real bug fixed

`env/verify_data_store_map.sh` now asserts the two-provenance invariant
(251 HiQBind symlink pairs + 30 RCSB-recovered real files) instead of a flat
`raw/ = 251`.

While extending it: **`OUT` was never exported**, so the truncated-receptor
check's Python heredoc read `root=""`, globbed nothing, and compared the
expected 16-id list against an empty string. The previous
"38/38 PASS, exit 0" claim was therefore **not true for that assertion**. Fixed
by exporting `OUT`; the check now reproduces the 16 ids. Verifier is
**40/40 PASS, exit 0** (2026-08-23).

---

## 2. `clusters_casf2016.json` has two corrupted PDB ids (letter-O → digit-zero)

**Status:** open · **Severity:** low now, high if anything starts joining on it

`labels/clusters_casf2016.json` — the file that *defines* CASF-2016 core
membership (57 clusters × 5 = 285) — contains `105b` and `10wh` (digit zero),
but the real ids are `1o5b` and `1owh` (letter O):

| id | `raw/` | `crystal_ligands/` | appears in |
|---|---|---|---|
| `1o5b`, `1owh` (letter O) | present | present | `PDBbind_data_dict.json`, `PDBbind_data_split_cleansplit.json`, `PDBbind_data_split_pdbbind.json` |
| `105b`, `10wh` (zero) | absent | absent | `clusters_casf2016.json` **only** |

So `clusters_casf2016.json` is the lone file with the corruption — a classic
O/0 substitution from some text-processing step.

**Why it hasn't bitten:** no Python in this repo reads
`clusters_casf2016.json` (verified by grep), and the pipeline enumerates
`raw/` directly. That is also why `raw/` appears to hold 2 ids "outside the core"
— it doesn't; the label file is misspelling them.

- [ ] `[confirmed]` Fix the two ids in `labels/clusters_casf2016.json`
      (`105b`→`1o5b`, `10wh`→`1owh`). Two-character change; keep a copy of the
      original.
- [ ] `[tentative]` Any future code that selects systems from the cluster file
      (e.g. the full-CASF multi-seed run) would silently drop these 2 systems
      and run 279 instead of 281. Fix before that lands.

---

## 3. `crystal_ligands/` is corpus-wide, not CASF — the `$CASF` name misleads

**Status:** open · **Severity:** documentation

`$CASF/crystal_ligands/` holds **14,661** SDFs, not 285, and that is correct —
the two directories have different scopes:

| directory | scope | count |
|---|---|---|
| `raw/` | CASF-2016 core subset only | **281** (251 HiQBind + 30 RCSB-recovered; 4 excluded — item 1) |
| `crystal_ligands/` | the whole `pdbbind_cleansplit` corpus | 14,661 |

`crystal_ligands/` also contains one non-SDF file, `download_status.json`
(hence 14,662 directory entries), which records **17,272 attempted / 14,661 ok /
2,611 failed** — that is exactly where 14,661 comes from.

The nesting, largest to smallest:

```
18,623  ids in PDBbind_original_train_val_split_f{0..4}
16,491  ids in PDBbind_cleansplit_train_val_split_f{0..4}
14,661  crystal ligands actually on disk   <- crystal_ligands/
   285  CASF-2016 core (labels/clusters_casf2016.json)
   281  ...of which we have a receptor      <- raw/  (251 HiQBind + 30 RCSB)
```

The root is a **general PDBbind clean-split benchmark**; CASF-2016 is one
labeled subset inside it. The misleading part is the variable name — `$CASF`,
`LAB_CASF`, `CONTRASCF_CASF_ROOT` all point at the *PDBbind cleansplit root*,
not a CASF-only tree.

- [x] `[confirmed]` **DONE 2026-08-23.** Scope ladder documented in
      `docs/data_store_map.md`, with an explicit naming warning at the top.
- [ ] `[confirmed]` Still stale: the environment entry
      `journal/2026-08-14-casf-data-store-map.md` and memory
      `remote_locations.md` both imply `crystal_ligands/` is CASF-scoped.
      Correct both (same edit as item 5).
- [ ] `[suggested]` Rename `CONTRASCF_CASF_ROOT` → `CONTRASCF_PDBBIND_ROOT`
      (keeping the old name as an alias) so the name stops implying CASF-only.
      Touches `env/carc.sh`, `env/lab.sh`, `analysis/casf_mutagenesis/config.py`.

---

## 4. Mutant docking receptors: co-folding provenance + lost mutations

**Status:** open · **Severity:** MEDIUM *(downgraded from HIGH on 2026-08-23
after direct measurement — see 4b. The defect is real but its effect on the
published pooled rates is statistically nil, p≈0.99.)*

Gates how a new engine (ICM) should be run.

### Where the docking inputs live

```
$CONTRASCF_ROOT/analysis/casf_mutagenesis/outputs/<pdbid>/<variant>/docking/
    receptor.pdb    # protein only (no ligand, no waters)
    ligand.sdf      # crystal ligand, 3D coords, identical across all 4 variants
    box.json        # {"center":[x,y,z], "size":[25,25,25]}
```

`<variant>` ∈ `wt | rem | pack | inv`. Complete cells: **wt 251**, **rem/pack/inv
239 each** = **968 total**, matching the 968 `gnina/`, `unidock2/`, `surfdock/`
output dirs — this is exactly the input set those engines consumed.

### 4a. WT and mutant receptors have different provenance

| | receptor | box centre | ligand |
|---|---|---|---|
| **wt** (`inputs_docking.py:5-6`) | **crystal protein** | **crystal** ligand centroid | crystal SDF |
| **rem/pack/inv** (`10_build_mutant_docking.py:62,71,74`) | **AF3+MSA predicted CIF**, stripped to protein | **AF3-predicted** ligand centroid | crystal SDF |

So the WT→mutant drop on the docking bars conflates the mutation with a
crystal→predicted structure change, and co-folding influence enters twice
(receptor conformation *and* search-box centre). The ligand side is fine —
these are protein mutations, so the crystal ligand is correctly reused.

**Severity: MILD.** There is *no ground-truth mutant structure to be had.* A
crystal-derived mutant (PyMOL/ICM mutate-residue + repack) is also a model, with
its own rotamer and backbone-relaxation errors, and is not obviously more
trustworthy than AF3's. "Use the crystal instead" trades one modelling
assumption for another rather than removing one.

What *is* cleanly fixable is the **asymmetry**, not the modelling choice: WT uses
a crystal while mutants use a prediction, so the WT bar is not a like-for-like
reference. Cheap fix — a co-folded WT prediction already exists for every
system (**239 `af3msa_<id>_wt_model_0.cif` files on disk**), so a WT docking
receptor can be built from it, putting all four variants on identical
provenance without inventing any new structure.

- [ ] `[suggested]` Build a second WT docking cell from
      `af3msa_<id>_wt_model_0.cif` (e.g. `wt_pred/`) so the docking arm can be
      reported predicted-WT vs predicted-mutant. Keeps the existing crystal-WT
      cell for the co-folding comparison.
- [ ] `[open]` Rebuilding mutants from crystal (4c Arm B) remains worth doing as
      an *independent arm*, not as a replacement — its value is that A−B
      measures the co-folding structural bias, not that B is "correct".

### 4b. Some mutant receptors lose their mutations — measured directly

Residue counts were only a proxy, and a misleading one. The decisive test is
whether the receptor actually **carries the mutations the spec asked for**:
recover them by diffing wt-vs-variant sequences in `boltz.yaml`, align the
receptor onto the spec chain (difflib opcodes, tolerant of AF3's routine 2-3
residue terminal trimming), and check each mutated position.

Script `analysis/casf_mutagenesis/scripts/audit_mutation_presence.py`
→ `analysis/casf_mutagenesis/mutation_presence_audit.json`.

> **Method warning.** The first version of this audit required the *full* spec
> chain to sit inside the receptor sequence and so misread normal terminal
> trimming as "chain missing" — it reported 90 dead cells and 30 broken systems.
> Those were **false positives**. The numbers below are from the corrected
> alignment. Do not quote the v1 figures.

**Result — 690 mutant cells, 3,882 specified mutations:**

| verdict | cells | share |
|---|---|---|
| OK (all mutations present) | 601 | 87.1 % |
| PARTIAL (some lost) | 65 | 9.4 % |
| ABSENT / NOCHAIN (none survived) | 18 | 2.6 % |
| NOMUT (spec had no mutation at all — see 4d) | 6 | 0.9 % |

**3,658 / 3,882 mutations (94.2 %) survive into the receptor.** Of the 224 lost:
`no_chain` 45, `not_modeled` 179, **`wrong_aa` 0** — AF3 never writes the wrong
residue; loss is always a position it did not model. 202 systems are fully clean
across all three variants.

Only **6 systems** have a cell where *no* mutation survived:
`1bcu 1lpg 1oyt 3ag9 3utu 5c2h` (`1bcu/1lpg/1oyt/3utu` lose the whole mutated
chain; `3ag9/5c2h` have every mutated position unmodelled). 22 systems show
partial loss, mostly mild (23 cells keep ≥75 %, 36 keep 50-75 %, 6 keep <50 %).

#### Impact on the published numbers: negligible, and self-cancelling

Recomputed from `outputs/docking_results.csv`:

| engine | pooled rate, all cells | excl. invalid | Δ | p |
|---|---|---|---|---|
| GNINA | 0.1241 (n=717) | 0.1241 (n=693) | −0.0000 | 0.999 |
| UniDock2 | 0.0837 (n=717) | 0.0851 (n=693) | +0.0015 | 0.922 |
| SurfDock | 0.0014 (n=704) | 0.0015 (n=680) | +0.0001 | 0.980 |

The reason the exclusion moves nothing is that **the two defect classes bias in
opposite directions and cancel**:

| group | n | rate <2 Å | median RMSD |
|---|---|---|---|
| all mutant cells | 2138 | 0.0702 | 7.50 Å |
| tiny-stub receptors (`1bcu 1lpg 1oyt 3utu`) | 36 | **0.0000** | 19.24 Å |
| WT-duplicate spec (`3mss 4eo8`, see 4d) | 18 | **0.2222** | 5.28 Å |
| valid cells | 2066 | 0.0707 | 7.44 Å |

WT-duplicates inflate apparent memorization exactly as theory predicts (3× the
overall rate — docking into a genuinely wild-type pocket recovers the native
pose more often), while truncated stubs drive it to zero. **The published
"docking does not memorize" conclusion is not inflated by this defect.**

- [x] `[confirmed]` **DONE 2026-08-23.** Measured directly; impact quantified.
- [ ] `[confirmed]` Add a build-time guard to `10_build_mutant_docking.py`:
      after writing `receptor.pdb`, verify every specified mutation is present;
      refuse (or mark) the cell otherwise. This is the durable fix — a residue-
      count tolerance would *not* have caught `3ag9`/`5c2h`.
- [ ] `[suggested]` Exclude the 18 dead cells from mutant docking analyses and
      report the exclusion; the pooled numbers do not change, so this is for
      correctness of the record, not to alter a result.

### 4d. `3mss` and `4eo8` have mutant specs identical to WT (generation bug)

For both systems, `rem`/`pack`/`inv` `boltz.yaml` protein sequences are
**byte-identical to `wt`** — the mutation generator produced no mutation, so
these 6 cells are WT duplicates masquerading as adversarial variants. They are
also the highest-scoring "mutant" cells in the whole set (0.222 vs 0.070), which
is exactly the contamination direction that matters.

### Checked 2026-08-25 — scope and impact measured

**Scope: exactly 2 systems, 6 cells, no wider spread.** Swept all 756 mutant
specs via `af3.json` (parses even where `boltz.yaml` is broken, see 4f):
only `3mss` and `4eo8` are byte-identical to WT.
List: `analysis/casf_mutagenesis/wt_duplicate_specs.json`.

**Yes, it reaches the co-folding arm**, and `4eo8` scores as memorized for both
models — but the effect on the headline rates is negligible:

| model | cell | rmsd (rem / pack / inv) | counts as memorized? |
|---|---|---|---|
| Boltz-2 | `3mss` | 2.15 / 2.15 / 2.15 | no (just misses) |
| Boltz-2 | `4eo8` | 0.67 / 0.67 / 0.67 | **yes** |
| AF3+MSA | `3mss` | 2.33 / 2.33 / 2.34 | no |
| AF3+MSA | `4eo8` | 0.93 / 0.93 / 0.93 | **yes** |

Excluding both systems moves every rate by **≤ 0.003** (Boltz-2 rem
0.231→0.229, AF3+MSA inv 0.259→0.257). No conclusion depends on it.

**Useful tell:** the RMSD is *identical across rem/pack/inv* — the signature of
one input producing one output. That is a cheap detector for this whole class of
silent no-op, needing no sequence comparison.

- [ ] `[confirmed]` Fix the generator so it cannot silently emit a mutation-free
      variant: raise (or mark the cell) when the rendered mutant sequence equals
      WT. `3mss` (274 aa) and `4eo8` (562 aa) presumably had no pocket residue
      pass the selection filter.
- [ ] `[suggested]` Add the identical-RMSD-across-variants check to the analysis
      as a standing guard — it catches no-op mutants without re-reading specs.
- [ ] `[suggested]` Until fixed, exclude `3mss`/`4eo8` from memorization rates
      and say so; the correction is ≤0.003, so this is for correctness of the
      record rather than to change a result.

### 4c. Recommended protocol for a new engine (ICM)

Run **both** arms; the difference is itself the measurement.

- **Arm A (comparability):** reuse `docking/` as-is → directly comparable to the
  existing GNINA/UniDock2 bars, inherits 4a and 4b caveats.
- **Arm B (clean):** build the mutant from `$CASF/raw/<id>/<id>_protein.pdb` by
  applying the same point mutations (ICM mutate-residue or PyMOL), keeping the
  full construct and backbone; box on the **crystal** ligand centroid. Removes
  co-folding bias and fixes 4a + 4b in one move.
- **A − B measures the co-folding structural bias** on docking outcomes.

- [ ] `[confirmed]` Write a mutation-extraction script: diff WT vs mutant
      protein sequences in `outputs/<id>/<variant>/boltz.yaml` per chain to
      recover the mutation list (there is **no** `applied_mutations` file — the
      mutation is implicit in the sequence). Verified on `1bcu/rem` → chain H
      `D199G, C231G`.
- [ ] `[open]` **Numbering caveat:** those indices are 1-based *within the chain
      sequence*, not crystal PDB residue numbers. A sequence→PDB-numbering
      alignment (gaps, insertion codes) is required before applying them in ICM.
      Validate on 2-3 systems before scaling.

---

## 4e. ROOT CAUSE — AF3+MSA emits **single-chain** predictions for every multi-chain system

**Status:** open · **Severity:** HIGH — this is upstream of 4b *and* it
contaminates AF3+MSA's own co-folding numbers, which 4b does not.

Tracing where the chains are lost, for `2qnq` / `3ag9` / `1bcu`:

| stage | chains |
|---|---|
| `boltz.yaml` (input) | 2 ✅ |
| `af3.json` (input) | 2 ✅ |
| **`af3msa_<id>_<var>_model_0.cif` (output)** | **1 ❌** |
| `docking/receptor.pdb` | 1 (faithfully mirrors the CIF) |

**The inputs are correct and our stripping code is correct** — the loss is in the
AF3+MSA prediction output itself. Sweep over the 52 multi-chain systems:

| model | output chain counts |
|---|---|
| Boltz-2 | 2 chains ×44, 3 chains ×6, missing ×2 — **preserved** |
| AF3+MSA | **1 chain × 51 of 52** — collapsed |

100 % collapse is systematic, not stochastic.

### ROOT CAUSE CONFIRMED (2026-08-23): our AF3+MSA runner rebuilds a 1-chain input

**It is not AF3 and not the CIF ingest — it is input construction in
`analysis/casf_mutagenesis/scripts/06_run_af3_msa_subset20.py`.** The correct
multi-chain `af3.json` is *never passed to AF3*; it is only mined for one
sequence and the ligand SMILES, and a fresh single-chain JSON is rendered:

```python
# line 77-83 — returns the FIRST protein and silently drops the rest
def _read_af3_protein_sequence(af3_json_path: Path) -> str:
    """Extract the (single) protein sequence from a no-MSA af3.json."""
    for s in d["sequences"]:
        if "protein" in s:
            return s["protein"]["sequence"]      # <-- first chain only

# line 234, then 246-250 — hardcodes a single chain "A"
variant_seq = _read_af3_protein_sequence(af3_json)
render_af3(name=prefix, chain_seqs=[("A", variant_seq)],
           ligand_smiles=smiles, out_path=json_msa,
           chain_msas={"A": variant_a3m})
```

Proof on disk: **80 `af3_msa.json` files (the file actually handed to AF3),
0 with more than one protein chain.** And it explains `1bcu` exactly — its
`af3.json` lists `L`(26) *before* `H`(257), so the first-protein rule kept the
26-residue L chain and discarded the 257-residue H chain that carried both
mutations.

**Fix and its complication.** `_read_af3_protein_sequence` must return all
`(id, seq)` pairs and `render_af3` must receive all of them. The non-trivial
part is the MSA: `chain_msas` is currently one a3m for one chain, so a
multi-chain rerun needs an MSA **per distinct sequence** (a homodimer can reuse
one; a hetero-complex like `1bcu` L/H needs two). The MSA-fetch path upstream
assumes a single chain too, so this is a real change, not a one-liner.

### It depresses AF3+MSA's own accuracy

Top-1 ligand RMSD < 2 Å from `results_full.csv`, split by system type:

| model | variant | single-chain | multi-chain | p |
|---|---|---|---|---|
| **AF3+MSA** | wt | **0.855** (n=179) | **0.627** (n=51) | **0.000** |
| AF3+MSA | rem | 0.330 (n=179) | 0.216 (n=51) | 0.119 |
| Boltz-2 (control) | wt | 0.609 (n=179) | 0.540 (n=50) | 0.380 |
| Boltz-2 | rem | 0.196 (n=179) | 0.360 (n=50) | 0.015 |

AF3+MSA loses **23 points** of WT accuracy on multi-chain systems while Boltz-2
— which keeps its chains — shows no significant gap. The natural reading is that
the missing chain is depressing AF3+MSA's score, i.e. **the reported AF3+MSA WT
rate is an under-estimate**. (Correlational: multi-chain systems could be
intrinsically harder, but the Boltz-2 control argues against that.)

- [x] `[confirmed]` **DONE 2026-08-23.** Traced to
      `06_run_af3_msa_subset20.py:77-83` (first-protein-only) +
      `:246-250` (`chain_seqs=[("A", variant_seq)]` hardcoded). Verified: 80/80
      `af3_msa.json` are single-chain.
- [x] `[confirmed]` **FIXED 2026-08-23** in `06_run_af3_msa_subset20.py`:
      - `_read_af3_protein_sequence` → **`_read_af3_protein_chains`**, returning
        every `(chain_id, sequence)` (scalar or list `id` both expanded).
      - Phase 1 now caches an MSA **per protein chain**, fetching once per
        *distinct* sequence (a homodimer reuses one; `1bcu` L/H gets two) and
        length-gating on the **total** across chains.
      - Phase 2 rewrites each chain's own a3m onto that chain's (possibly
        mutated) sequence and calls
        `render_af3(chain_seqs=variant_chains, chain_msas=chain_msas)`.
      `render_af3` already supported multi-chain + per-chain MSAs, so no change
      was needed there. Run log now records `n_chains` per cell.

      **Validation gate (no GPU):**
      `analysis/casf_mutagenesis/scripts/26_verify_af3_multichain_input.py`
      re-renders every job and asserts chain count, per-chain sequence identity,
      and that each MSA's query row equals its chain's sequence.
      **756 cells / 165 multi-chain: 0 failures, exit 0.**

      Demonstrated on the previously broken cells:
      | system | rendered chains | mutations now inside the job |
      |---|---|---|
      | `1bcu/rem` | `L`(26) + `H`(257) | `H:D199G, H:C231G` (were both lost) |
      | `2qnq/rem` | `A`(99) + `B`(99) | 9 across both chains (were 4 lost) |
      | `3ag9/rem` | `A`(325) + `B`(319) | 9 on chain B (were all 9 lost) |
- [x] `[confirmed]` **RE-RUN DONE 2026-08-24.** All 4 variants re-run for the 53
      runnable multi-chain systems: **202/213 cells ok, 6.0 h AF3 wall**
      (~108 s/cell; runtime is near-flat in sequence length, not quadratic).
      Driver `27_pilot_af3_multichain_wt.py --variants wt,rem,pack,inv`;
      single-chain baselines preserved under `outputs/_af3_mc_pilot_backup/`.
      Failures: `1o5b` (all 4, AF3 exit 1 — also the O/0-typo id and one of the
      10 unparseable YAMLs), `2c3i` (all 4), `2xbv/wt` (transient MSA server).

### RESULT — the defect was **masking memorization**, not inflating it

On the affected systems (paired, McNemar exact):

| variant | n | old <2 Å | new <2 Å | Δ | old median | new median | p |
|---|---|---|---|---|---|---|---|
| wt | 50 | 0.600 | 0.700 | +0.100 | 1.48 | 1.04 | 0.302 |
| **rem** | 51 | 0.255 | **0.588** | **+0.333** | 6.20 | 1.63 | **<0.001** |
| **pack** | 50 | 0.300 | **0.560** | **+0.260** | 4.19 | 1.65 | **0.002** |
| **inv** | 51 | 0.078 | **0.549** | **+0.471** | 6.24 | 1.93 | **<0.001** |

Effect on the published AF3+MSA bars (all 239 systems):

| variant | published | corrected | Δ |
|---|---|---|---|
| wt | 0.803 | 0.824 | +0.021 |
| rem | 0.305 | 0.377 | +0.071 |
| pack | 0.261 | 0.315 | +0.055 |
| inv | **0.159** | **0.259** | **+0.100** |

**Interpretation.** A single-chain prediction is structurally broken, so the
ligand landed far from native and scored as "not memorized". That was an
artifact of a mutilated input, not the model recognising the mutation. Given
the complete complex, AF3+MSA places the ligand near-native *despite* the
pocket mutations — i.e. it memorises **more** than published, and the
`inv` bar is understated by a factor of ~1.6.

This **strengthens** the paper's core claim (co-folding models memorise) while
showing the published AF3+MSA numbers were wrong in the conservative direction.
It also contradicts my own pre-run prediction that memorization would drop.

- [ ] `[confirmed]` Regenerate `results_full.csv` / `memorization_full.csv` and
      the overview figure from the corrected CIFs — the AF3+MSA bars in
      `figures/overview_full.png` are now stale in a way that matters.
- [x] `[confirmed]` **MULTI-SEED CONFIRMED 2026-08-24.** `inv` re-run under
      seeds 2 and 3 (102 jobs, `28_af3_multiseed_check.py`; seed runs isolated
      in `outputs/_af3_seeds/`, RMSD via the pipeline's own `analyze_prediction`
      through a try/finally CIF swap).

      | seed | n | rate <2 Å | median |
      |---|---|---|---|
      | 1 | 51 | 0.549 | 1.93 |
      | 2 | 51 | 0.510 | 1.97 |
      | 3 | 51 | 0.529 | 1.84 |

      **mean 0.529, sd 0.020, range [0.510, 0.549].** Against the single-chain
      baseline of 0.078 the correction is **+0.45, about 22× the across-seed
      sd** — not a seed artifact.

- [ ] `[confirmed]` **Aggregate rates are stable; individual systems are NOT.**
      Same data, per-system: the <2 Å verdict is consistent across all 3 seeds
      for 43/51 (84 %) but **flips for 8/51 (16 %)**. Per-system RMSD spread has
      median sd 0.19 Å but a 90th percentile of 2.86 Å and a max of 13.54 Å
      (`2wn9/inv`: 7.68 / 7.96 / 31.28 Å; `4bkt/inv`: 7.20 / 0.37 / 0.33 Å).
      Consequence for the whole benchmark: **aggregate rates may be quoted from
      a single seed (sd ≈ 0.02), but any per-system or paired claim must not
      be** — including the "5 WT systems regressed" observation above, which is
      most likely seed noise rather than a real effect.
- [ ] `[open]` `2vw5` (856 aa) and `2j78` (888 aa) still exceed
      `MAX_TOTAL_LENGTH=800` once chains are summed honestly — decide whether to
      raise the gate or exclude them explicitly.
- [ ] `[open]` Once fixed, re-run AF3+MSA for the 52 multi-chain systems and
      re-check the panel (a) WT bar — currently 0.803 (n=239) and likely low.
- [ ] `[open]` Blast radius beyond docking: any AF3+MSA-derived number
      (memorization, confidence, affinity) for those 52 systems is computed on
      a half-structure.

### 4f. 10 `boltz.yaml` files are unparseable

`1o3f 1o5b 1uto 2y5h 2yge 3uuo 4e5w 4m0z 4x6p 5dwr` — `yaml.safe_load` raises
(e.g. `1o3f/rem/boltz.yaml` line 8: "expected \<block end>, but found scalar").

- [ ] `[open]` Determine whether Boltz-2 ever consumed these (a malformed input
      may have silently failed) or whether the corruption post-dates the run.
      Note `1o5b` is also the O/0-typo id from item 2.

---

## 5. Correction to the 2026-08-14 data-store map

**Status:** open · **Severity:** low (documentation)

The environment entry `journal/2026-08-14-casf-data-store-map.md` states lab
`crystal_ligands/<id>_ligand.sdf` entries are symlinks into HiQBind and that
CARC holds real files. **The lab side is wrong** — a full scan of all 14,661
entries gives **14,661 real files, 0 symlinks**. Only `raw/` uses symlinks on
lab (502/502 symlinks, all resolve). Real files on *both* hosts.

- [x] `[confirmed]` **DONE 2026-08-23.** Corrected in `docs/data_store_map.md`
      (marked `[CHANGED]`) and now asserted by the verifier.
- [ ] `[confirmed]` Still to fix: the same wrong claim in the environment entry
      `journal/2026-08-14-casf-data-store-map.md` and in memory
      `remote_locations.md`. Both still say lab = symlink.

---

## 6. Overview figure: charge relabel done but uncommitted; committed PNG is stale

**Status:** partly done · **Severity:** medium (a published figure is wrong)

### 6a. Charge-group labels were misleading — fixed, not committed

Panel (b) legend used the raw dict keys `chrg-` / `chrg+`, which read as a
symmetric negative-vs-positive dichotomy. They are not. Both ladders amputate the
**same anionic triphosphate** at the same bond and differ only in the grafted
tail (`ligand_mutagenesis/rules/charge_swap.py:36-45`):

| key | variants | replacement tail | charge effect |
|---|---|---|---|
| `chrg-` | `chrg_neu_{methyl,ethyl,propyl}` | branched alkane | anion → **neutral** (charge *removed*) |
| `chrg+` | `chrg_pos_{1,2,3}` | 1-3 quaternary ammonium | anion → **cation** (charge *flipped*) |

The WT is the anionic one; these are two rungs going the *same* direction. Rungs
2 and 3 are heavy-atom matched (9/9, 13/13), so the neutral ladder is the
**isosteric control** for the cationic one — that, not "two severity levels", is
why both exist. (Rung 1 is unmatched: 7 vs 5.)

Empirically the two rungs are **statistically indistinguishable** — GNINA
0.182 vs 0.214 (z=−0.64, p=0.52); UniDock2 0.065 vs 0.071 (z=−0.20, p=0.84).
Caveat on interpretation: neither arm isolates charge alone, because *both*
delete the whole polyphosphate. The honest reading is "given the anchor is
destroyed, sign does not further matter" — **not** "charge is irrelevant".

- [x] `[confirmed]` **DONE 2026-08-23.** Relabelled to `chrg→neutral` /
      `chrg→flipped (+)` with a subtitle stating both start from the anionic WT
      (`13_plot_overview.py`, new `LIG_LABELS` map; matches the pattern already
      in `ligand_mutagenesis/scripts/06_plot_affinity.py:52-56`). Re-rendered.
- [ ] `[confirmed]` Not committed. Also the counterfold worktree still holds the
      **old** copy of `13_plot_overview.py` — the two have diverged.

### 6b. The committed `overview_full.png` was ~2 weeks stale

`figures/overview_full.png` was last committed **2026-05-23** (`e702fa5`), but
`outputs/results_full.csv` was modified **2026-06-06**. Re-rendering moved
AF3+MSA WT in panel (a) from **n=19, rate 0.897** to **n=239, rate 0.803** —
caused by data added in June, *not* by the label edit.

- [ ] `[confirmed]` Any slide/report quoting AF3+MSA WT ≈ 0.90 (n=19) is stale;
      the current value is ≈ 0.80 (n=239). Re-check downstream text.
- [ ] `[open]` Why does panel (a) show AF3+MSA WT n=239 when
      `results_full.csv` holds **285** AF3+MSA WT rows (and AF3 shows 285)?
      `load_wt_rates_from_results` is filtering 46 systems — confirm that is
      intended before the number is quoted anywhere.

---

## 7. Nothing is committed; the 2026-08-14 notebook entries were never tracked

**Status:** open · **Severity:** medium (work is unprotected and partly wrong on disk)

Uncommitted on `main`:

```
 M analysis/casf_mutagenesis/scripts/13_plot_overview.py   # item 6a
 M analysis/casf_mutagenesis/figures/overview_full.png     # item 6b
 ?? analysis/casf_mutagenesis/receptor_size_audit.json     # item 4b
 ?? docs/data_prep_todo.md                                 # this file
 ?? docs/data_store_map.md                                 # verified map
 ?? env/verify_data_store_map.sh                           # 38-check verifier
```

Untracked in the counterfold worktree — the 2026-08-14 lab-notebook entries were
**completed but never committed**, so they exist only on this disk *and* now
contain claims items 3 and 5 prove wrong:

```
 ?? journal/2026-08-14-casf-data-store-map.md   (+ .meta.yaml, .results.yaml)
 ?? experiments/2026-08-14-casf-data-store-map-verify/
```

- [ ] `[confirmed]` `main` is the default branch — branch before committing the
      six paths above rather than committing straight to `main`.
- [ ] `[confirmed]` Fix the two wrong claims in the environment entry
      (items 3, 5) **before** committing it, so the corrected version is what
      lands.
- [ ] `[suggested]` Sync `13_plot_overview.py` from `main` into the counterfold
      worktree so the relabel is not lost on the next worktree run.

---

## Verification commands

```bash
R=$LAB_CASF
ls -1 $R/raw | wc -l                 # 281  (251 HiQBind + 30 RCSB-recovered)
ls -1 $R/crystal_ligands | wc -l     # 14662 (14661 sdf + download_status.json)
# cluster ids (57 clusters x 5 = 285), and the O/0 typo:
python -c "
import json;d=json.load(open('$R/labels/clusters_casf2016.json'))
ids=[r[0].lower() for c in d.values() for r in c]
print(len(d),'clusters',len(ids),'ids')
print([i for i in ids if i in ('105b','10wh')])"
```

---

## 8. Engine runners: documented idempotency ≠ implemented idempotency

**Status:** open · **Severity:** low (but it costs a silent no-op run)

`08_run_gnina_variants.py` and `11_run_unidock2_variants.py` both document
"skips cells whose `poses.sdf` already exists". They actually gate on their own
run log:

```python
seen = {(r["system"], r["variant"]) for r in runs if r.get("status") == "ok"}
...
if (system, variant) in seen: continue        # the poses.sdf check never runs
```

So clearing `poses.sdf` to force a re-run does nothing — the job reports
`new ok=0 skip=0 fail=0` and exits. Re-running a subset requires removing those
cells from `{gnina,unidock2}_variants_casf_mutagenesis_run_log.json` as well.

- [ ] `[confirmed]` Either honour `poses.sdf` (drop the log gate) or fix the
      docstrings, and add a `--force` / `--ids` flag like the one added to
      `10_build_mutant_docking.py`.

---

## 9. SurfDock cannot be re-run on this host — panel (b) is mixed-provenance

**Status:** open · **Severity:** medium (a figure implies a comparison it does not have)

The conda env `surfdock` is absent on katlab, so SurfDock was **not** re-run
after the receptors were rebuilt. Its bars come from the pre-fix single-chain
receptors while GNINA and UniDock2 come from corrected ones.

Its results were briefly moved aside with the other engines during the re-run
and had to be **restored from `_docking_prefix_backup/`** (153 dirs) — without
that, its bar would have been computed on 183 cells instead of 234, silently.

A provenance note is now printed on `overview_full.png` itself
(`13_plot_overview.py`), so the figure cannot be read as like-for-like.

- [ ] `[confirmed]` Re-run SurfDock once its env is available (or on CARC) and
      drop the caveat.
- [ ] `[open]` SurfDock's rate is ~0.000 across all variants, so it contributes
      nothing to the memorization conclusion either way — decide whether it
      earns a place in the figure at all.
