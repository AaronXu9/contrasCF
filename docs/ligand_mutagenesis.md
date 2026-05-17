# CASF-2016 ligand-mutagenesis pipeline

Automated implementation of Masters et al. 2025's *ligand-side* adversarial
challenges (Fig. 4 glucose methylation, Fig. 5 ATP charge swap) on the
CASF-2016 protein–ligand complex set. Lives under
[`analysis/ligand_mutagenesis/`](../analysis/ligand_mutagenesis/).

This module is the **ligand-side sibling** of
[`analysis/casf_mutagenesis/`](../analysis/casf_mutagenesis/). The two
modules share data discovery, input renderers, and docking helpers; they
differ in what they perturb (protein pocket vs ligand chemistry) and in
the manifest shape.

## What the paper does, and what the existing repo had

The paper hand-built **12 ligand variants on 2 systems**:

- **Fig. 4 (n=6):** progressive O-methylation of glucose hydroxyls (one
  more methylated OH per rung), tested on GDH (2VWH).
- **Fig. 5 (n=6):** ATP triphosphate replaced by either neutral alkyl
  chains (methyl / ethyl / propyl gem-dimethyl ladder) or quaternary
  ammonium chains (+1 / +2 / +3 formal-charge ladder), tested on CDK2
  (1B38).

The paper did **not** run a CASF-scale ligand sweep analogous to Fig. 3's
pocket sweep. The shipped repo contains only the 12 hand-built SMILES
(see [`analysis/src/config.py::ATP_CHARGE_SMILES`](../analysis/src/config.py)
and `GLUCOSE_SMILES`). This module **encodes the paper's transformation
rules** so they apply systematically to every CASF-2016 ligand carrying
the relevant functional group — extending the paper's two systems × 12
variants to the full dataset.

## Authoritative rules

Paper § "Glucose methylation" and § "ATP charge variants" only describe
the transformations narratively; the algorithmic specification below is
inferred from the literal SMILES they shipped, decoded carefully so the
generator reproduces all 12 hand-built SMILES exactly (canonical-SMILES
equality).

### Hydroxyl methylation (paper Fig. 4)

| variant | rule |
|---|---|
| `wt` | original ligand SMILES (control) |
| `meth_k` (k=1..N) | first **k** aliphatic free hydroxyls methylated (-OH → -OCH₃) |

Eligibility SMARTS: `[OX2H][CX4]` — free hydroxyl on sp3 carbon. Ring
oxygens (no implicit H) and phenolic hydroxyls (sp2 C) are excluded by
construction.

**Order of methylation across the ladder** is chemistry-aware so the
paper's glucose ladder reproduces exactly:

1. **Class priority:** hemiacetal OH → secondary OH → primary OH →
   methanol-like.
2. **Within-class tiebreak:** shortest bond-path distance from the
   hemiacetal C (so secondary OHs are walked around the ring from the
   anomeric C — for glucose: C2, C3, C4 in pyranose order).
3. **Final tiebreak:** RDKit canonical atom rank.

Pure canonical-rank ordering picks the wrong OH first on glucose (it
picks the 6'-CH₂OH before the anomeric OH). The chemistry-aware ordering
+ path-distance tiebreak is required for paper reproduction.

Ladder length is capped at `MAX_METHYLATIONS = 6` to bound combinatorial
blow-up on polyhydroxyl ligands.

### Phosphate / carboxylate charge swap (paper Fig. 5)

Two ladders, both attached at a single anchor:

| variant | replacement attached to bridging-O |
|---|---|
| `chrg_neu_methyl` | `-CH₂-C(CH₃)₂-CH₂-CH₂-CH₃` (7 atoms) |
| `chrg_neu_ethyl`  | `-CH₂-C(CH₃)₂-CH₂-C(CH₃)₂-CH₃` (9 atoms) |
| `chrg_neu_propyl` | `-CH₂-C(CH₃)₂-CH₂-C(CH₃)₂-CH₂-C(CH₃)₂-CH₃` (13 atoms) |
| `chrg_pos_1` | `-CH₂-N⁺(CH₃)₃` (5 atoms, +1 charge) |
| `chrg_pos_2` | `-CH₂-N⁺(CH₃)₂-CH₂-N⁺(CH₃)₃` (9 atoms, +2 charge) |
| `chrg_pos_3` | `-CH₂-N⁺(CH₃)₂-CH₂-N⁺(CH₃)₂-CH₂-N⁺(CH₃)₃` (13 atoms, +3 charge) |

Increment from `chrg_neu_methyl` → `chrg_neu_ethyl` is uneven (+2 atoms,
not +3) because the methyl variant ends in a -CH₂-CH₂-CH₃ propyl tail
while the ethyl/propyl variants end in a -C(CH₃)₂-CH₃ tert-butyl tail.
The replacement SMILES are decoded directly from the paper's
`ATP_CHARGE_SMILES` strings — see code comments in
[`rules/charge_swap.py`](../analysis/ligand_mutagenesis/rules/charge_swap.py).

**Eligibility (priority: phosphate > carboxylate):**

| path | SMARTS | paper-faithful? |
|---|---|---|
| phosphate | `[#6][OX2;!R][P]` (-C-O-P linkage; bridging O kept) | yes (paper Fig. 5 target) |
| carboxylate | `[CX4][CX3](=O)[OX2H,OX1-]` (alpha-C – COOH/COO⁻) | no (soft extension) |

**Algorithm.** Match the SMARTS; cut the bond at the anchor → leaving-group
boundary (O-P bond for phosphate; alpha-C – carboxyl-C bond for
carboxylate). Drop the disconnected leaving-group fragment via
`Chem.GetMolFrags`. Attach the replacement fragment (atom 0 = attachment
point) at the kept anchor via `Chem.CombineMols` + `AddBond`. Atom-index
tracking through `RemoveAtom` calls uses `Atom.SetAtomMapNum(99)` on the
kept anchor before deletion.

**Stereochemistry preservation:** the kept anchor C and the rest of the
sugar are untouched by the rule, so RDKit preserves their `[C@H]`/
`[C@@H]` tags through the round-trip.

### Aromatic halogenation (informational extension — NOT in the paper)

| variant | rule |
|---|---|
| `halo_F_1`, `halo_Cl_1`, `halo_Br_1` | single aromatic-CH → C-X at the lowest-canonical-rank aromatic CH |

Eligibility SMARTS: `[c;H1]`. RDKit handles aromatic-CH → aromatic-CX
natively; no kekulization concern. Single substitution only (no
combinatorial ladder — that's deliberate to bound the manifest size).
Variants are flagged `paper_faithful=False`.

## Variant naming scheme

Names are filename-safe (alphanumeric + underscores). Unlike the
binding-site sibling — which has a fixed `VARIANTS = ("wt","rem","pack",
"inv")` tuple — the variant set here is **dynamic per system** because
the eligible rules depend on the ligand's functional groups. Downstream
model-runner scripts should read variant names from the manifest, not
from a hard-coded tuple.

| pattern | rule | meaning |
|---|---|---|
| `wt` | – | unmodified ligand |
| `meth_<k>` | methylation | k aliphatic OHs methylated, k ∈ {1..min(N_OH, 6)} |
| `chrg_neu_{methyl,ethyl,propyl}` | charge_swap | neutral-alkyl ladder |
| `chrg_pos_{1,2,3}` | charge_swap | cationic-ammonium ladder |
| `halo_{F,Cl,Br}_1` | halogenation | single aromatic halogen substitution |

## Data sources

Same as the binding-site sibling — all CASF paths flow through
`casf_mutagenesis.config` (re-exported by `ligand_mutagenesis.config`):

| input | path |
|---|---|
| crystal protein PDB | `data/casf2016/raw/<pdbid>/<pdbid>_protein.pdb` |
| crystal ligand SDF | `data/casf2016/crystal_ligands/<pdbid>_ligand.sdf` |
| HET-code SMILES cache | `data/casf2016/smiles/het_smiles_cache.json` |
| PDBbind metadata | `data/casf2016/labels/PDBbind_data_dict.json` |
| 20-complex smoke test | `data/casf2016/labels/PDBbind_casf2016_subset20.json` |
| full CASF-2016 (285 PDBs) | `data/casf2016/labels/PDBbind_data_split_cleansplit.json["casf2016"]` |

## Architecture

```
analysis/ligand_mutagenesis/
├── config.py             # re-exports CASF paths; adds OUTPUT_ROOT, MAX_METHYLATIONS, HALOGENS
├── ligand_io.py          # SMILES from crystal SDF (cross-imports casf_mutagenesis)
├── rules/
│   ├── base.py           # LigandRule protocol + LigandVariant dataclass
│   ├── methylation.py    # paper Fig. 4 — chemistry-aware OH ladder
│   ├── charge_swap.py    # paper Fig. 5 — phosphate / carboxylate swap
│   └── halogenation.py   # informational extension
├── build.py              # per-system pipeline orchestrator
├── scripts/
│   ├── 00_verify_reference_systems.py   # gating test (paper 11/11 SMILES)
│   ├── 01_build_subset20.py             # 20-system smoke build
│   └── 02_build_full_casf.py            # 285-system full build (--start/--limit)
└── outputs/
    ├── <pdbid>/<variant>/{af3.json, boltz.yaml,
    │                      docking/{receptor.pdb, ligand.sdf, box.json}}
    ├── manifest_subset20.json
    └── manifest_full.json
```

## Per-system pipeline (`build.py::build_system`)

For each PDB id:

1. **Resolve files** — `protein.pdb` and `_ligand.sdf` under `data/casf2016/`.
2. **Get WT SMILES.** Prefer crystal-SDF-derived canonical SMILES
   (`Chem.MolToSmiles(mol_from_sdf)`) over the HET cache — same reasoning
   as the binding-site sibling (kekulization safety).
3. **Extract protein chain sequences** via
   `casf_mutagenesis.sequence.extract_chain_sequences` (gemmi-based).
   Sequences are **constant across all variants** in this module since
   we only perturb the ligand.
4. **Run rules.** For each rule (methylation, charge_swap, halogenation):
   record per-rule eligibility in the manifest. If eligible, generate
   the ladder variants.
5. **Render per-variant inputs:**
   - **AF3:** `casf_mutagenesis.inputs_af3.render_af3` (local "alphafold3"
     dialect, version 2, no-MSA placeholder).
   - **Boltz:** `casf_mutagenesis.inputs_boltz.render_boltz` (YAML,
     `msa: empty`).
   - **Docking:** `receptor.pdb` + `box.json` are constant per system —
     written once for `wt` and **hardlinked or copied** to every other
     variant's `docking/` directory. `ligand.sdf` is freshly
     RDKit-embedded per variant via
     `casf_mutagenesis.inputs_docking._embed_ligand`. Box is centred on
     the WT crystal ligand centroid (25 × 25 × 25 Å, same as
     binding-site sibling).
6. **Write manifest entry** capturing wt SMILES + source, per-rule
   eligibility report, per-variant {rule, smiles, applied changes,
   paper_faithful flag, out_dir}, warnings.

## Verification (gating test)

[`scripts/00_verify_reference_systems.py`](../analysis/ligand_mutagenesis/scripts/00_verify_reference_systems.py)
asserts that the rules reproduce all 11 paper SMILES on the two paper
systems. Without this gate the systematic CASF sweep is meaningless.

| ladder | n | parent SMILES | source in repo |
|---|---|---|---|
| glucose methylation | 5 | `GLUCOSE_SMILES["glucose_0"]` | `analysis/src/config.py:72` |
| ATP charge swap (neutral) | 3 | `ATP_SMILES` | `analysis/src/config.py:52-55` |
| ATP charge swap (cationic) | 3 | `ATP_SMILES` | `analysis/src/config.py:52-55` |

The test canonicalises both the generated and the expected SMILES via
`Chem.MolToSmiles(canonical=True, isomericSmiles=True)` before comparing,
so equivalent SMILES writings (different starting atom, different
canonicalisation seed) compare equal. Stereochemistry is verified by
direct atom-by-atom comparison of the canonical strings.

**Negative-eligibility check** included: glucose has no phosphate or
carboxylate, so `charge_swap.is_eligible` must refuse it.

**Last known good run (2026-05-07):**

```
=== glucose hydroxyl methylation (paper Fig. 4) ===
  ✓ meth_1..5 all match glucose_1..5

=== ATP triphosphate charge swap (paper Fig. 5) ===
  ✓ chrg_neu_{methyl,ethyl,propyl} all match atp_charge_{methyl,ethyl,propyl}
  ✓ chrg_pos_{1,2,3} all match atp_charge_{1,2,3}

=== negative eligibility checks ===
  ✓ charge_swap refuses glucose

=== summary ===
  glucose methylation : PASS
  ATP charge swap     : PASS
  negative eligibility: PASS
```

## Pilot run — subset20 (2026-05-07)

`01_build_subset20.py` on the 20-PDB CASF subset
(`PDBbind_casf2016_subset20.json`).

**Outcome:** 20/20 systems built, 0 errors, 0 warnings, **93 total variants**.

**Per-rule eligibility distribution:**

| rule | n eligible / 20 | systems |
|---|---|---|
| methylation | 5 | 4bkt, 4w9l, 3dx2, 3b5r, 4ivc |
| charge_swap | 2 | 1vso (phosphate), 1z9g (carboxylate) |
| halogenation | 19 | all except 1vso |

Most CASF ligands have at least one aromatic CH (halogenation eligible).
Methylation eligibility is much rarer (5 of 20 have a free aliphatic
hydroxyl). Charge-swap eligibility is rarest (2 of 20: one phosphate,
one carboxylate).

**Cross-module sanity:** `outputs/1pxn/wt/boltz.yaml` ligand SMILES
matches `casf_mutagenesis/outputs/1pxn/wt/boltz.yaml` exactly — both come
from the same `_smiles_from_crystal_sdf` path.

**Per-variant chemistry spot-check (4bkt):** WT ligand SMILES
`CNC(=O)[C@@H]1C[C@@H](O)CN1C(=O)Cc1cc(C)no1` → `meth_1` SMILES
`CNC(=O)[C@@H]1C[C@@H](OC)CN1C(=O)Cc1cc(C)no1` — exactly one extra `-OC`
group, stereochemistry preserved.

**Results path:**

```
analysis/ligand_mutagenesis/outputs/manifest_subset20.json
analysis/ligand_mutagenesis/outputs/<pdbid>/<variant>/{af3.json, boltz.yaml,
                                                       docking/{receptor.pdb,
                                                                ligand.sdf,
                                                                box.json}}
```

Manifest entries record per-variant `{rule, smiles, applied, paper_faithful,
out_dir}` plus a top-level `rules_eligible` dict; downstream model-runner
scripts read the variant set from this manifest (not from a fixed tuple).

## Running on a different host (env-var overrides)

This module re-exports all CASF / Boltz / AF3 paths from
`casf_mutagenesis.config`, so the same env-var overrides apply
([`docs/casf_mutagenesis.md` § "Running on a different host"](casf_mutagenesis.md#running-on-a-different-host-env-var-overrides)).

The one module-specific override:

| variable | meaning | lab default |
|---|---|---|
| `CONTRASCF_LIGAND_OUT` | output root for this module | `analysis/ligand_mutagenesis/outputs` |

## Running the pipeline

All commands assume `LD_LIBRARY_PATH` points at `rdkit_env`'s `lib/`
(set automatically by `env/lab.sh` / `env/carc.sh`).

```bash
LIB="/home/aoxu/miniconda3/envs/rdkit_env/lib"
PY="/home/aoxu/miniconda3/envs/rdkit_env/bin/python"

# 1. Gating test — must pass before doing anything else
LD_LIBRARY_PATH=$LIB:$LD_LIBRARY_PATH $PY \
  analysis/ligand_mutagenesis/scripts/00_verify_reference_systems.py

# 2. Build inputs for the 20-system subset
LD_LIBRARY_PATH=$LIB:$LD_LIBRARY_PATH $PY \
  analysis/ligand_mutagenesis/scripts/01_build_subset20.py

# 3. Build inputs for the full CASF-2016 (285 systems).
#    Supports --start / --limit / --manifest for chunked re-runs.
LD_LIBRARY_PATH=$LIB:$LD_LIBRARY_PATH $PY \
  analysis/ligand_mutagenesis/scripts/02_build_full_casf.py
```

## Implementation gotchas (and workarounds)

| issue | resolution |
|---|---|
| Paper's `glucose_0` SMILES is a different canonical writing of the same α-D-glucopyranose as `GLC_SMILES` | verify always compares canonicalised SMILES; never compare raw strings |
| Pure canonical-rank ordering picks 6'-CH₂OH before anomeric OH on glucose | chemistry-aware priority class (hemiacetal → secondary → primary) + path-distance tiebreak required for paper reproduction |
| `Chem.GetShortestPath(mol, a, a)` errors with "Invariant Violation aid1 != aid2" | guard `c_idx != hemi_c_idx` before calling |
| Atom indices shift after `RWMol.RemoveAtom` | tag the anchor with `SetAtomMapNum(99)` before deletion, then look up by map number |
| Charge-swap chain-length off-by-one easy to make | replacement SMILES decoded by hand from `analysis/src/config.py::ATP_CHARGE_SMILES`; covered by 6/6 gating-test assertions |
| RDKit "molecule is tagged as 2D" warning when reading crystal SDFs | benign — emitted by most PDBbind crystal SDFs; ignore |

## Out of scope (deferred to follow-ons)

- **Running the models** (AF3 / Boltz-2 / docking) on the generated
  inputs. Mirrors `casf_mutagenesis/scripts/03_*` and `04_*`. Variant
  names in the manifest are filename-safe so the wrappers can be written
  by reading variant names from the manifest (not from a fixed tuple).
- **MSA injection** for the AF3 inputs. Protein sequences are constant
  per system in this module, so one MSA per system can be reused across
  all variants of that system. Reuse `casf_mutagenesis.msa_via_boltz`
  (Boltz piggyback) and `inputs_af3.clean_a3m_for_af3`.
- **Per-variant ligand-RMSD scoring.** The binding-site sibling's
  `casf_mutagenesis.analysis.superpose_by_index` will work directly
  (protein numbering is identical across variants since protein is
  constant), but per-rule common-substructure SMARTS — analogous to
  `analysis/src/config.py::COMMON_SUBSETS["ADENOSINE"]` and
  `["GLC_CORE"]` — are needed to pair "the same atoms" across a WT
  ligand and its variants for RMSD purposes. Design that in a follow-on.
- **More chemistry rules.** Deisostere swap (C=O → C=S, etc.) and
  combinatorial halogenation are candidate extensions; current scope is
  the paper's two rules + single-position halogenation only.
