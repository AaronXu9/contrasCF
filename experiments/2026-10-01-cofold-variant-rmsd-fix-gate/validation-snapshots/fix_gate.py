"""Gate for the _atom_correspondences fix in casf_mutagenesis/analysis.py.
F1 known-answer: crystal ligand vs (crystal + extra F) / (crystal - 1 atom) /
   (crystal with one atom's element changed) placed at the crystal coords.
   True in-place and bestfit RMSD are 0 on the shared atoms -> must be < 1e-6.
   The OLD code must be shown to FAIL F1 on the +F case (proves the test bites).
F2 protein-arm regression: re-score a random sample of protein-arm cells with
   the NEW code; ligand_rmsd_a must equal results_full.csv (old code) exactly.
F3 ligand arm: re-scored halo_F_1 median bestfit RMSD must fall to within
   0.3 A of WT's (a one-atom swap cannot change scaffold geometry by Å).
"""
import sys, importlib.util, random, math
from pathlib import Path
import numpy as np, pandas as pd
from rdkit import Chem
REPO = Path("/mnt/katritch_lab2/aoxu/contrasCF"); sys.path.insert(0, str(REPO/"analysis"))
import casf_mutagenesis.analysis as NEW
spec = importlib.util.spec_from_file_location("casf_mutagenesis.analysis_old", sys.argv[1])
OLD = importlib.util.module_from_spec(spec); sys.modules[spec.name] = OLD; spec.loader.exec_module(OLD)
ok = True
def check(n, c, d):
    global ok; ok &= bool(c); print(f"[{'PASS' if c else 'FAIL'}] {n}: {d}")

# ---------- F1 synthetic known-answer ----------
CASF = Path("/home/aoxu/projects/VLS-Benchmark-Dataset/data/pdbbind_cleansplit/crystal_ligands")
xtal = Chem.MolFromMolFile(str(CASF/"1bcu_ligand.sdf"))
xh = Chem.RemoveHs(xtal); xyz = xh.GetConformer().GetPositions()
def with_coords(rw, coords):
    m = rw.GetMol(); c = Chem.Conformer(m.GetNumAtoms())
    for i, p in enumerate(coords): c.SetAtomPosition(i, p.tolist())
    m.RemoveAllConformers(); m.AddConformer(c, assignId=True); return m
# +F on an aromatic carbon that has an H
ar = next(a.GetIdx() for a in xh.GetAtoms() if a.GetIsAromatic() and a.GetSymbol()=="C" and a.GetTotalNumHs()>0)
rw = Chem.RWMol(xh); f = rw.AddAtom(Chem.Atom(9)); rw.AddBond(ar, f, Chem.BondType.SINGLE)
plusF = with_coords(rw, np.vstack([xyz, xyz[ar] + [1.35, 0, 0]])); Chem.SanitizeMol(plusF)
# -1 atom: drop a terminal heavy atom
term = next(a.GetIdx() for a in xh.GetAtoms() if a.GetDegree()==1)
rw2 = Chem.RWMol(xh); rw2.RemoveAtom(term)
minus1 = with_coords(rw2, np.delete(xyz, term, axis=0))
# element change on a terminal atom -> neither is a substructure -> MCS path
rw3 = Chem.RWMol(xh); a = rw3.GetAtomWithIdx(term); a.SetAtomicNum(17 if a.GetAtomicNum()!=17 else 35)
swap = with_coords(rw3, xyz)
for name, pm in [("+F", plusF), ("-1 atom", minus1), ("element swap", swap)]:
    pxyz = Chem.RemoveHs(pm).GetConformer().GetPositions()
    r, n = NEW._matched_rmsd(xtal, pxyz, pm); b = NEW._bestfit_rmsd(xtal, pxyz, pm)
    mode = NEW._atom_correspondences(xh, Chem.RemoveHs(pm))[1]
    check(f"F1 new {name}", r < 1e-6 and b < 1e-6, f"in-place {r:.2e}  bestfit {b:.2e}  matched {n}  mode {mode}")
pxyz = plusF.GetConformer().GetPositions()
r_old, _ = OLD._matched_rmsd(xtal, pxyz, plusF)
# shuffle pred atom order the way an independently written SDF would be ordered
perm = list(range(plusF.GetNumAtoms())); random.Random(0).shuffle(perm)
shuf = Chem.RenumberAtoms(plusF, perm); sx = shuf.GetConformer().GetPositions()
r_old_s, _ = OLD._matched_rmsd(xtal, sx, shuf); r_new_s, _ = NEW._matched_rmsd(xtal, sx, shuf)
check("F1 test bites: OLD code wrong on reordered +F", r_old_s > 0.5, f"old {r_old_s:.2f} A vs new {r_new_s:.2e} A (old, file-order +F: {r_old:.2e})")

# ---------- F2 protein-arm regression ----------
rf = pd.read_csv(REPO/"analysis/casf_mutagenesis/outputs/results_full.csv")
rf = rf[rf.status=="ok"]
samp = rf.groupby(["model","variant"], group_keys=False).apply(lambda g: g.sample(min(len(g), 25), random_state=0))
OUT = REPO/"analysis/casf_mutagenesis/outputs"; diffs = []; n=0
for _, r in samp.iterrows():
    sp = NEW.MODEL_FILES[r.model]; vdir = OUT/r.pdbid/r.variant; prefix=f"{r.pdbid}_{r.variant}"
    cifs = sorted(vdir.glob(sp["cif_glob"].format(prefix=prefix)))
    if r.pose_idx >= len(cifs): continue
    rec = NEW._analyze_single_pose(r.pdbid, r.variant, r.model, int(r.pose_idx), cifs[int(r.pose_idx)], sp, vdir)
    n += 1; diffs.append(abs((rec.ligand_rmsd_a if rec.ligand_rmsd_a is not None else math.nan) - r.ligand_rmsd_a))
diffs = np.array(diffs)
check("F2 protein arm unchanged", n >= 200 and np.nanmax(diffs) < 1e-3 and not np.isnan(diffs).any(),
      f"{n} cells re-scored, max |Δ| {np.nanmax(diffs):.4f} A, nan {int(np.isnan(diffs).sum())}")
print("GATE-PARTIAL", "PASS" if ok else "FAIL")
