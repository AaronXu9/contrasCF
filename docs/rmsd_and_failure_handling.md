# RMSD computation and failure handling

How every RMSD in this benchmark is computed, and what happens to cells that
fail. Written 2026-08-28 by reading the implementation, not from memory; every
claim cites `file:line`.

There are **two independent RMSD code paths** — one for co-folding predictions,
one for docking poses. They differ in ways that matter for cross-method
comparison (see §4).

---

## 1. Co-folding path — `analysis/casf_mutagenesis/analysis.py`

Three RMSD variants are recorded per prediction, all heavy-atom only
(`analysis.py:121-138`).

### 1a. `ligand_rmsd_a` — the canonical metric

**This is the number behind every memorization rate and every `<2 Å` bar.**

1. Load the crystal ligand and the predicted complex.
2. `extract_protein_ca_near` on **both** sides — Cα of the chain *closest to the
   ligand*, not the largest chain (`analysis.py:600-607`).
3. `superpose_by_index(pred_ca, native_ca)` → rotation+translation mapping
   predicted → crystal frame; record `ca_rmsd_a`, `n_ca_paired`
   (`analysis.py:610-612`).
4. Apply that same R,t to the predicted ligand's heavy coordinates
   (`analysis.py:615`).
5. Symmetry-corrected heavy-atom RMSD against the crystal ligand
   (`_matched_rmsd`, `analysis.py:631`).

**Why the ligand-near chain and not the largest** (`analysis.py:600-604`): CASF
contains homo-multimers where the crystal and the prediction place the ligand on
*different homologous chains*. Superposing on the wrong copy makes a correct
prediction look catastrophic — `4w9l` WT scored **47 Å** ligand RMSD on a fold
that was actually right (Cα 0.66 Å). This choice is load-bearing, not cosmetic.

### 1b. `ligand_rmsd_fullca_a` — all-chain alignment

Same procedure but superposing on **all** chains' Cα
(`analysis.py:637-658`). Larger than `ligand_rmsd_a` on multi-chain systems when
the model gets inter-chain orientation wrong. Wrapped in a bare `except:
pass` — a failure here silently leaves the column `None` (`analysis.py:659-660`).

### 1c. `bestfit_rmsd_a` — pocket-blind, the paper's metric

Symmetry-aware **Kabsch** fit of the predicted ligand onto the crystal ligand,
ignoring the protein entirely (`_bestfit_rmsd`, `analysis.py:440-486`). Answers
"is the internal geometry right?", not "is it in the right pocket". Also
wrapped in `except: pass`.

### 1d. Symmetry correction — `_matched_rmsd` (`analysis.py:361-437`)

Enumerates substructure matches between predicted and crystal heavy-atom
graphs with `uniquify=False, maxMatches=200`, and returns the **minimum RMSD
over all correspondences** — so a symmetric ligand (a phenyl flip, equivalent
carboxylate oxygens) is not penalised for atom relabelling.

Fallback chain, in order — each silently degrades to a **same-order** mapping:

| condition | code | behaviour |
|---|---|---|
| crystal/pred heavy counts differ | `:375-379` | truncate **both to `min(n)`** and compare index-by-index |
| `pred_mol is None` | `:381-384` | assume identical atom order |
| no substructure match either direction | `:397-408` | assume identical atom order |
| heavy count ≠ match-atom count | `:420-422` | assume identical atom order |

**These fallbacks are silent.** They return a finite RMSD with no status flag,
so a nonsense correspondence is indistinguishable in the output from a real
symmetry-minimised match. The truncation case (`min(n)`) is the most dangerous:
it compares the first *n* atoms in file order, which is arbitrary.

### 1e. Failure statuses

Only three (`analysis.py:570, 591, 671`):

| status | meaning |
|---|---|
| `missing_cif` | no prediction file for that cell |
| `no_ligand` | the CIF had no recoverable ligand |
| `error` | any exception; message stored in `rec.error` |

Confidence and affinity sidecars are merged **regardless** of RMSD success
(`analysis.py:674-683`) — a failed-RMSD cell still reports its model's
confidence.

---

## 2. THE KEY POINT — failures are dropped, not penalised

`memorization_stats` (`analysis.py:750`) and the pose selectors
(`:794`, `:905`) all filter:

```python
if r.status != "ok" or r.ligand_rmsd_a is None:
    continue
```

and then set `n_total = len(surviving)` (`analysis.py:758`).

**A cell that fails leaves the denominator entirely.** It is not counted as a
miss. So every rate in this benchmark is *conditional on the model having
produced a parsable prediction*:

> rate = P(RMSD < 2 Å **given** the prediction succeeded)

Consequences to keep in mind:

- **`n` varies by model and variant** — this is why the bars carry per-method
  `n=` labels (Boltz-2 229, AF3+MSA 239, GNINA/UniDock2 251, SurfDock 244).
  They are *different denominators*, not a shared one.
- A model that crashes often is **not** punished; it is scored on the subset it
  managed. A method failing on its hardest cases would look artificially good.
- Comparisons across methods with different `n` are not strictly like-for-like.
  For a strict comparison, restrict to the intersection of successful cells.

This is a defensible convention (it separates "wrong pose" from "no output"),
but it must be stated whenever the rates are published.

---

## 3. Docking path — `analysis/casf_mutagenesis/gnina_analysis.py`

Different code, different assumptions.

### 3a. No superposition by default

`_aligned_rmsd` (`gnina_analysis.py:109-121`) is a **direct atom-by-atom heavy
RMSD with no superposition** — docked poses are already in the receptor frame.
Returns `NaN` on shape mismatch.

### 3b. …except for mutant cells, which need a frame change

For `rem/pack/inv`, docking ran against the **AF3-predicted mutant receptor**,
so poses are in the AF3 frame and must be mapped to the crystal frame first
(`_ca_transform`, `gnina_analysis.py:122-138`). For WT and all of
`ligand_mutagenesis` the receptor *is* the crystal protein, so the transform is
the identity and is skipped.

### 3c. Atom matching — `_mcs_match_indices` (`gnina_analysis.py:75-107`)

1. `pose.GetSubstructMatch(crystal)` — fast path; pose is a superstructure
   (halogenation/methylation **add** atoms).
2. reverse direction.
3. `rdFMCS.FindMCS` with `CompareElements` / `CompareAny` bonds, `timeout=10` —
   needed for `charge_swap`, where neither molecule contains the other.

Returns `[], [], 0` if MCS finds nothing → status `no_match`.

### 3d. Failure statuses (`gnina_analysis.py:150-185`)

`missing_pose`, `missing_crystal`, `parse_error`, `no_match`, `align_failed`,
`error` — a richer set than the co-folding path's three, and the same
drop-from-denominator rule applies.

---

## 4. ⚠ The two paths are not symmetry-equivalent

| | co-folding (`_matched_rmsd`) | docking (`_mcs_match_indices`) |
|---|---|---|
| matches enumerated | up to **200**, `uniquify=False` | **one** (`GetSubstructMatch`, singular) |
| RMSD reported | **minimum over all mappings** | the single match found |
| symmetric ligands | correctly handled | **can be inflated** |

The docking path takes the *first* substructure match rather than minimising
over equivalent ones. For a ligand with symmetry (para-substituted rings,
equivalent carboxylate oxygens, freely rotating terminal groups) the docking
RMSD is therefore **larger than it should be**, while the co-folding RMSD for
the same molecule is symmetry-minimised.

### Measured, not assumed (2026-08-28)

**71 %** of a 250-ligand sample of CASF crystal ligands have more than one graph
automorphism, so most ligands are exposed to this. Re-scoring **244 WT GNINA
cells** both ways (WT chosen because the pose is already in the crystal frame,
so no superposition step confounds the comparison):

| | single match (current) | symmetry-minimised |
|---|---|---|
| mean RMSD | 1.943 Å | **1.760 Å** |
| median RMSD | 1.086 Å | **0.815 Å** |
| rate < 2 Å | 0.738 | **0.762** |

Inflation (single − minimised): mean **0.182 Å**, median 0.000 Å, max
**8.543 Å**. **29 %** of cells are inflated by >0.1 Å and **11 %** by >0.5 Å.
Worst cases: `1h23` 11.37 → 2.82 Å, `2qnq` 7.27 → 5.29 Å (32 matches),
`4mgd` 2.03 → 0.63 Å (48 matches).

**Net effect: the docking `<2 Å` rate is understated by ≈ 2.4 points**
(+6 of 244 cells) purely from the matcher.

That is the same direction as the paper's headline claim — "co-folding
memorises, docking does not" — so the metric is mildly flattering to the
conclusion. The correction is far smaller than the reported co-folding/docking
gap (≈ 0.38 vs ≈ 0.13), so the claim survives comfortably; but the number
should be fixed before publication rather than defended.

- [x] `[confirmed]` Quantified 2026-08-28 (above).
- [ ] `[confirmed]` Port `_matched_rmsd`'s match enumeration into
      `gnina_analysis.py::_mcs_match_indices` (swap `GetSubstructMatch` for
      `GetSubstructMatches(..., uniquify=False)` and minimise) so both arms use
      one matcher, then re-run `12_analyze_docking_engines.py`. Expect docking
      rates to rise ~2 points.

---

## 5. Practical notes

- All RMSDs are **heavy-atom only**; hydrogens are stripped
  (`Chem.RemoveHs`) on both sides.
- Everything is rounded to 3 decimals on write.
- Bond orders are recovered from a **crystal-derived SMILES**, the same string
  fed to the model when the inputs were generated (`analysis.py:583-588`).
- `pose_idx = 0` is the model's own top-ranked pose; all headline rates use it.
- Thresholds: `< 2 Å` is the headline; `< 4 Å` is also recorded
  (`MemorizationStats.n_below_4A`).
