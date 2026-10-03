# surfdock-protonation-sensitivity

**Type:** validation
**Kind:** validation
**Date:** 2026-09-04
**Status:** refuted

## What this validates
**Plan step:** n/a
**Validates iteration:** surfdock-receptor-provenance-artifact
**Validates upcoming experiment:** re-run of the SurfDock mutant arm on protonated receptors

## Validation kind
sensitivity

## Assertion

SurfDock consumes an MSMS molecular surface plus MaSIF-style chemical channels. If that
representation is **insensitive** to receptor protonation, then the protonation asymmetry in
our inputs (wt = protonated HiQBind crystal, mutants = AF3 predictions with zero hydrogens)
is cosmetic and can be ignored. Concretely: building the surface for the *same pocket* from a
protonated vs an H-stripped receptor should give the same vertex count (± a few %), vertices
in the same places (median nearest-neighbour displacement ≪ 0.5 Å, i.e. below MSMS's own
discretisation), and unchanged hbond / charge / hydrophobicity channels.

Falsified if vertex counts, vertex positions or channel statistics move materially — in which
case the SurfDock mutant cells were scored off-distribution and their numbers are unsafe.

## Setup
- **Code:** `f2529b8cea19d89c36524d76efa8edf879a28403` on `fix/af3-multichain-and-data-audit` (dirty)
- **Env:** `/home/aoxu/miniconda3/envs/SurfDock` (python 3.8, pymesh 0.3.1, plyfile, numpy 1.24.4, scipy 1.8.1)
- **Host:** katlab
- `TMPDIR=/home/aoxu/tmp/msms_htest` — MSMS leaves ~3.6 MB scratch per call and fills `/` otherwise.

## Model & data provenance

| artifact | exact id / name | size | source | version / sha256 | role |
|---|---|---|---|---|---|
| surface builder | `_surfdock_surface_helper.py` | 216 lines | `CogLigandBench` @ `2107de18e9` (post interface-crop fix) | dist_threshold=8, crop at 3 Å | builds the .ply SurfDock consumes |
| SurfDock install | SurfDock (MaSIF-derived surface stack) | n/a | `/home/aoxu/projects/SurfDock` | MSMS 2.6.1, APBS 3.4.1, pdb2pqr 2.1.1 | `computeMSMS`, `computeCharges` |
| wt receptors | `outputs/<sys>/wt/docking/receptor.pdb` | 3 systems | HiQBind / PDBBind-Opt curated crystal | 1e66 8356 atoms (4106 H), 3u5j 2121 (1059 H), 2brb 4463 (2233 H) | protonated arm |
| mutant-like receptors | same PDBs, hydrogens deleted | 3 systems | derived in-place by element filter | 1e66 4250, 3u5j 1062, 2brb 2230 atoms | stands in for AF3 predictions (0 H) |
| SurfDock reference input | `1a0q_protein_processed.pdb` | 6284 atoms | `SurfDock/data/eval_sample_dirs/test_samples/` | **3101 H (49 %)** | the distribution SurfDock expects |

## Validation code
`validation-snapshots/compare.py` — reads both .ply sets, reports vertex/face counts, per-channel
mean ± sd, and a KD-tree nearest-neighbour displacement between the two vertex clouds.

## Commands & run log
```
python _surfdock_surface_helper.py --data_dir <withH|noH> --out_dir out_<mode>
python compare.py
```

## Expected outcome
Vertex counts within a few percent, median NN displacement ≪ 0.5 Å, channels unchanged.

## Actual outcome

```
system  mode     verts   faces   charge            hbond            hphob           si
1e66    withH      154     258   -1.892+-0.664    -0.050+-0.198    -0.099+-1.891   +0.426
1e66    noH        148     236   -1.424+-1.114    -0.080+-0.213    -0.012+-2.048   +0.402
3u5j    withH      125     181   -0.307+-0.738    -0.030+-0.124    +1.151+-2.808   +0.293
3u5j    noH         80     118   -0.280+-0.799    -0.026+-0.112    +0.570+-2.726   +0.344
2brb    withH      137     204   -0.697+-0.590    -0.039+-0.238    +1.101+-2.901   +0.323
2brb    noH         97     137   -0.503+-0.957    -0.078+-0.194    +1.754+-2.387   -0.226

system    d_verts       %   median NN   p95 NN   hbond mean
1e66           -6   -3.9%      0.552A   1.029A   -0.050 -> -0.080
3u5j          -45  -36.0%      0.574A   0.900A   -0.030 -> -0.026
2brb          -40  -29.2%      0.589A   0.947A   -0.039 -> -0.078
```

Vertex loss up to **36 %**; median surface displacement **0.55–0.59 Å** (p95 ≈ 0.9–1.0 Å);
hbond mean changes by up to **2×**; `2brb` shape index **flips sign** (+0.323 → −0.226).

Mechanism, traced in the source rather than inferred:
- `computeMSMS.py:20` — the `protonate=True` argument **adds no hydrogens**. It only selects
  `output_pdb_as_xyzrn` over the deprecated `pdb2xyzrn`. There is **no protonation step anywhere
  in the pipeline**; the input PDB's H content passes straight through to MSMS.
- `triangulation/xyzrn.py:37` — hydrogens are first-class: given a radius, written to the xyzrn,
  and classified against `polarHydrogens[resname]`.
- `computeCharges.py:127` — `if atom_name in polarHydrogens[res.get_resname()]` is how the
  **H-bond channel** is built. With no hydrogens in the file that branch never fires.

## Pass / fail
**fail** — the assertion is falsified. SurfDock's input representation is materially
protonation-dependent, and its own reference input is 49 % hydrogen, so unprotonated
AF3 mutant receptors are off-distribution.

## Code ↔ theory alignment

| claim / property checked | code (file:line / function) | justified by (ref ↓) | match? |
|---|---|---|---|
| hydrogens enter the solvent-excluded surface via per-atom radii | `SurfDock/comp_surface/prepare_target/triangulation/xyzrn.py:19-39` | Sanner 1996 (MSMS: SES from atomic radii) | exact |
| H-bond channel is defined by explicit polar hydrogens | `computeCharges.py:127` | Gainza 2020 (MaSIF chemical features) | exact |
| `protonate=True` does not protonate | `computeMSMS.py:20-29` | n/a — misleading legacy flag; verified by reading both branches | **divergent from its own name** |
| expected input distribution is protonated | `SurfDock/data/eval_sample_dirs/test_samples/1a0q/` (49 % H) | MaSIF prep protonates with Reduce upstream | exact |
| Vina-family scoring is protonation-invariant (contrast) | GNINA `--score_only`, 3 systems | Trott & Olson 2010 (united-atom Vina scoring) | exact — affinity bit-identical |

## References
- Sanner, M. F., Olson, A. J., Spehner, J.-C. (1996). Reduced surface: an efficient way to compute molecular surfaces. *Biopolymers* 38, 305. DOI 10.1002/(SICI)1097-0282(199603)38:3<305::AID-BIP4>3.0.CO;2-Y
- Gainza, P. et al. (2020). Deciphering interaction fingerprints from protein molecular surfaces using geometric deep learning (MaSIF). *Nat. Methods* 17, 184. DOI 10.1038/s41592-019-0666-6
- Word, J. M. et al. (1999). Asparagine and glutamine: using hydrogen atom contacts in the choice of side-chain amide orientation (Reduce). *J. Mol. Biol.* 285, 1735. DOI 10.1006/jmbi.1998.2401
- Olsson, M. H. M. et al. (2011). PROPKA3: consistent treatment of internal and surface residues. *J. Chem. Theory Comput.* 7, 525. DOI 10.1021/ct100578z
- Trott, O., Olson, A. J. (2010). AutoDock Vina. *J. Comput. Chem.* 31, 455. DOI 10.1002/jcc.21334
- Code: `CogLigandBench@2107de18e9`; SurfDock install `/home/aoxu/projects/SurfDock`.

## Optional: If failed — what to fix

- [confirmed] The SurfDock mutant arm ran off-distribution. Its adversarial rates (0.024–0.046) and its WT→adversarial gap (+0.842) are not safe to publish until re-run.
- [confirmed] Fix is to **protonate the AF3 mutant receptors** (Reduce, or PDB2PQR+PROPKA at pH 7.4) before the surface step — not to strip hydrogens from the crystal, since SurfDock expects them.
- [open] Re-run the ~717 SurfDock mutant cells on protonated receptors and re-measure the gap.
- [open] Quantify how much of the 5.2–5.5 Å zero-mutation drift this accounts for: rebuild those two control cells protonated and see whether the drift closes.
- [tentative] GNINA and UniDock2 need no re-run — Vina affinity is bit-identical with and without H, and the GNINA CNN moves ≤ 0.006.
- [suggested] Add a receptor-protonation assertion to the SurfDock runner so an unprotonated receptor fails loudly instead of silently degrading.
