# Data-store map — CASF-2016 / PDBbind-cleansplit (katlab ↔ CARC)

Supersedes the 2026-08-14 version. Every claim below is mechanically verified by
`env/verify_data_store_map.sh` (**40/40 PASS, exit 0**, re-run 2026-08-23 against
both hosts). Lines marked **[CHANGED]** were wrong in the previous map. Open defects
are tracked in [`data_prep_todo.md`](data_prep_todo.md).

> **Naming warning.** `$CASF` / `LAB_CASF` / `CONTRASCF_CASF_ROOT` point at the
> **PDBbind clean-split root**, *not* a CASF-only tree. CASF-2016 is one labeled
> subset (285 ids) inside a ~14.6k-complex corpus. See the scope ladder below.

## Roots

**KatLab (this workstation)**

```bash
export LAB_REPO=/mnt/katritch_lab2/aoxu/contrasCF
export LAB_CASF=/home/aoxu/projects/VLS-Benchmark-Dataset/data/pdbbind_cleansplit  # == $LAB_REPO/data/casf2016 (symlink)
export LAB_HIQBIND=/mnt/katritch_lab2/aoxu/data/hiqbind                            # ultimate source of protein+ligand
export LAB_CF=$LAB_REPO/.claude/worktrees/counterfold                              # worktree-counterfold branch
export LAB_CACHE=$LAB_CF/counterfold/_cache                                        # CounterFold-derived data (lab only)
```

**CARC (discovery.usc.edu)** — `source $CARC_REPO/env/carc.sh` sets these:

```bash
export CARC_REPO=/project2/katritch_223/aoxu/contrasCF                                              # == $CONTRASCF_ROOT
export CARC_CASF=/project2/katritch_223/aoxu/projects/VLS-Benchmark-Dataset/data/pdbbind_cleansplit # == $CONTRASCF_CASF_ROOT
export CARC_HIQBIND=/project2/katritch_223/aoxu/data/hiqbind
export CARC_SCRATCH=/scratch1/aoxu/contrasCF
# NOTE: CounterFold code lives at $CARC_REPO/counterfold (branch checked out at root; no .claude/worktrees/)
```

## Scope ladder — what each id count means **[CHANGED]**

```
18,623  ids in PDBbind_original_train_val_split_f{0..4}
16,491  ids in PDBbind_cleansplit_train_val_split_f{0..4}
14,661  crystal ligands actually on disk        <- $CASF/crystal_ligands/
   285  CASF-2016 core (57 clusters x 5)        <- $CASF/labels/clusters_casf2016.json
   281  ...of which a receptor exists           <- $CASF/raw/
        = 251 HiQBind-backed + 30 RCSB-recovered (2026-08-23; see divergence 4)
```

`crystal_ligands/` is **corpus-wide, not CASF-scoped** — 14,661 SDFs (one file
per PDB entry, exactly one molecule each) plus `download_status.json`
(17,272 attempted / 14,661 ok / 2,611 failed) = 14,662 directory entries.
Consumers index it by `<id>`; nothing globs it, so only the ~281 CASF entries
are ever read.

## Data — paths relative to `$CASF` and `$REPO`

Substitute `$CASF`→`$LAB_CASF`/`$CARC_CASF`, `$REPO`→`$LAB_REPO`/`$CARC_REPO`.
`<id>` = 4-char PDB code (e.g. `1bcu`).

```
# ── 1. ORIGINAL / GROUND-TRUTH crystal structures (the benchmark reference) ──
$CASF/crystal_ligands/<id>_ligand.sdf     # GT crystal ligand pose  [BOTH hosts: REAL FILE]   **[CHANGED]**
                                          #   lab 14,661/14,661 real, 0 symlinks (full scan)
$CASF/raw/<id>/<id>_protein.pdb           # GT receptor. TWO provenances:
                                          #   251 symlink→$HIQBIND (PDBBind-Opt refined), all resolve
                                          #    30 REAL file, RCSB deposited coords, UNREFINED [lab only]
$CASF/raw/<id>/<id>_ligand.sdf            # GT ligand     [symlink→$HIQBIND; only the 251 HiQBind ids]
$CASF/labels/clusters_casf2016.json       # CASF-2016 membership -- CONTAINS 2 CORRUPT IDS, see divergence 4
$CASF/labels/PDBbind_casf2016_subset20.json
$CASF/labels/PDBbind_cleansplit_train_val_split_f{0..4}.json

# ── 2. MUTATED structures + ALL-METHOD predictions (the benchmark itself) ──
$REPO/analysis/casf_mutagenesis/outputs/<id>/{wt,rem,pack,inv}/
    ├─ boltz.yaml, af3.json                 # mutation spec (input). NO applied_mutations file --
    │                                       #   the mutation is implicit; recover by sequence diff
    ├─ <id>_<variant>_model_{0..4}.cif      # Boltz-2 co-folded poses
    ├─ af3msa_<id>_<variant>_model_0.cif    # AF3+MSA co-folded pose
    ├─ confidence_*.json, affinity_*.json   # co-folding confidence/affinity
    └─ docking/{receptor.pdb, ligand.sdf, box.json}   # engine inputs -- see §2b
       gnina/  unidock2/  surfdock/         # engine outputs
$REPO/analysis/casf_mutagenesis/{gnina,unidock2}_variants_casf_mutagenesis_run_log.json
$REPO/analysis/ligand_mutagenesis/outputs/<id>/...   # ligand-side variants (methyl/charge/halogen)

# ── 3. COUNTERFOLD-derived data (labels, models, generated poses) — LAB ONLY ──
$LAB_CACHE/plip_labels/<id>_complex.pdb   # crystal complexes for PLIP interaction labels
$LAB_CACHE/plip/*.xml                     # per-system PLIP interaction analysis
$LAB_CACHE/flowr/{dry12,dataset100}/      # FLOWR datasets (.smol)
$LAB_CACHE/flowr/carc_heldout/{gen_finetuned,gen_frozen}/   # S_ft vs S_frozen poses (L_pose test)
$LAB_CACHE/affinity/*.ckpt                # trained heads (interaction_head*, typed_contrastive*)
```

### §2b. The `docking/` contract — what GNINA / UniDock2 / SurfDock consumed **[NEW]**

`receptor.pdb` (protein only) + `ligand.sdf` (crystal ligand, identical across
all 4 variants) + `box.json` (`{"center":[x,y,z],"size":[25,25,25]}`).
Directly reusable by a new engine (e.g. ICM).

**Receptor provenance differs between WT and mutants** — this is a confound,
not a convention:

| variant | receptor | box centre | source |
|---|---|---|---|
| `wt` | **crystal protein** | **crystal** ligand centroid | `inputs_docking.py:5-6` |
| `rem`/`pack`/`inv` | **AF3+MSA predicted CIF**, stripped to protein | **AF3-predicted** ligand centroid | `10_build_mutant_docking.py:62,71,74` |

Cell counts — **katlab**: wt 251, rem/pack/inv 239 each = **968**, matching the
968 `gnina`/`unidock2`/`surfdock` output dirs. **CARC**: wt 251, mutants **0**
(see divergence 5).

**16 systems have truncated mutant receptors** (`< 50 %` of WT residues; 56 have
any mismatch). Audit: `analysis/casf_mutagenesis/receptor_size_audit.json`.

```
1bcu 1lpg 1oyt 2vw5 2wn9 3n7a 3n86 3utu 4bkt 4f2w 4u4s 4w9c 4w9h 4w9i 4w9l 5c2h
```

Worst case `1bcu`: input had chains L(26)+H(257) with **both mutations on chain
H**; the predicted CIF has only chain A = 26 residues, so the receptor is a
26-residue stub containing **no mutation at all**.

## Where the two hosts diverge (6 caveats) **[EXPANDED from 3]**

1. `$LAB_REPO/data/casf2016` is a symlink → `$LAB_CASF`. On CARC,
   `$CARC_REPO/data/casf2016` is a real dir with only `labels/` — reach
   crystal/raw via **`$CARC_CASF`**, never `$CARC_REPO/data/casf2016/...`.
2. Mutagenesis outputs: **293** entries on lab vs **273** on CARC; gap = loose
   logs/CSVs + one system dir `7xlp_mek1` (absent on CARC).
3. `$LAB_CACHE` (§3) is **lab-only** — on CARC bootstrapped via
   `$LAB_CACHE/carc_runkit/` bundles, not a byte-for-byte mirror.
4. **[UPDATED 2026-08-23]** `raw/` holds **281** on lab, not 285. HiQBind lacked
   34 CASF-core systems; **30 were recovered from the RCSB deposited mmCIF** by
   `analysis/casf_mutagenesis/scripts/25_recover_missing_receptors.py` (chain
   selection uses HiQBind's own 10 Å rule; validated on 25 HiQBind-backed
   controls, which reproduce the 3.5 Å pocket **exactly 25/25**). The 4 still
   absent are `1a30`/`3bv9`/`3uri` (peptide ligands — categorical exclusion) and
   `3f3a` (reference SDF matches none of its 4 TRP copies). The recovered 30 are
   **lab-only**; on **CARC all 285 dirs exist and 34 are still empty**, because
   `setup_casf_data_carc.sh:51` runs `mkdir -p` before the source-exists check.
   Provenance per system: `analysis/casf_mutagenesis/recovered_receptors.json`.
5. **[NEW]** Docking inputs: CARC has **wt only** (251). Mutant `docking/` cells
   are **absent** (0), and the 968 `gnina`/`unidock2` dirs there are **empty
   scaffolding, not results** — every docking result to date was produced on
   katlab. The 239 AF3+MSA mutant CIFs *are* on CARC, so mutant inputs can be
   regenerated in place with `10_build_mutant_docking.py` rather than rsynced
   (but that reproduces the same 16 truncated receptors — guard first).
6. **[NEW]** `labels/clusters_casf2016.json` contains **`105b` and `10wh`**
   (digit zero) where the real ids are **`1o5b` and `1owh`** (letter O). It is
   the only file with the corruption; the other three labels files and both data
   dirs use the correct spelling. No repo code reads it today, so nothing is
   broken yet — but a naive set-difference against `raw/` reports 36 missing
   systems instead of the true 34.

## Quick verification

```bash
ls -1 $LAB_CASF/raw | wc -l                       # 281 (251 HiQBind + 30 RCSB)
ls -1 $LAB_CASF/crystal_ligands | wc -l           # 14662 (14661 sdf + download_status.json)
find $LAB_REPO/analysis/casf_mutagenesis/outputs -path "*/docking/receptor.pdb" | wc -l   # 968
```
