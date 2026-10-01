"""Gate for S1: is the new ligand-arm scoring trustworthy?
G1 regression : Boltz-2 rows in results_ligand.csv unchanged by the AF3+MSA extension.
G2 reference  : ligand-arm AF3+MSA WT top-1 reproduces the protein-arm AF3+MSA WT
                (same chains, same native SMILES, same MSA cache) -> rate within 0.08
                and per-system pass/fail concordance >= 0.85 on the overlap.
G3 regression : docking_results.csv rows for module=casf and ligand gnina/unidock2
                unchanged; ligand surfdock rows newly present.
"""
import sys, pandas as pd, numpy as np
SP = sys.argv[1]
ok = True
def check(name, cond, detail):
    global ok; ok &= bool(cond); print(f"[{'PASS' if cond else 'FAIL'}] {name}: {detail}")

new = pd.read_csv("analysis/ligand_mutagenesis/outputs/results_ligand.csv")
old = pd.read_csv(f"{SP}/results_ligand.pre_af3.csv")
key = ["pdbid","variant","model","pose_idx"]
nb = new[new.model=="Boltz2"].sort_values(key).reset_index(drop=True)
ob = old.sort_values(key).reset_index(drop=True)
cols = ["status","ligand_rmsd_a","confidence_score","affinity_pred_value"]
same = nb.shape==ob.shape and all(np.allclose(nb[c].astype(float), ob[c].astype(float), equal_nan=True) if c!="status" else (nb[c]==ob[c]).all() for c in cols)
check("G1 Boltz-2 rows unchanged", same, f"new {nb.shape} vs old {ob.shape}")

prot = pd.read_csv("analysis/casf_mutagenesis/outputs/results_full.csv")
p = prot[(prot.model=="AF3+MSA")&(prot.variant=="wt")&(prot.pose_idx==0)&(prot.status=="ok")].set_index("pdbid").ligand_rmsd_a
l = new[(new.model=="AF3+MSA")&(new.variant=="wt")&(new.pose_idx==0)&(new.status=="ok")].set_index("pdbid").ligand_rmsd_a
ov = p.index.intersection(l.index)
rp, rl = (p[ov]<2).mean(), (l[ov]<2).mean()
conc = ((p[ov]<2)==(l[ov]<2)).mean()
check("G2a WT rate reproduces", abs(rp-rl)<=0.08, f"overlap n={len(ov)}  protein-arm {rp:.3f}  ligand-arm {rl:.3f}  (full ligand-arm n={len(l)} rate {(l<2).mean():.3f})")
check("G2b per-system concordance", conc>=0.85, f"{conc:.3f}  median |dRMSD|={np.median(np.abs(p[ov]-l[ov])):.2f} A")

d_new = pd.read_csv("analysis/casf_mutagenesis/outputs/docking_results.csv")
d_old = pd.read_csv(f"{SP}/frozen_pre_s1/docking_results.csv")
k2 = ["module","engine","system","variant"]
def sub(df): return df[~((df.module=="ligand")&(df.engine=="surfdock"))].sort_values(k2).reset_index(drop=True)
a, b = sub(d_new), sub(d_old)
same2 = a.shape==b.shape and (a[k2+["status"]]==b[k2+["status"]]).all().all() and np.allclose(a.rmsd_a, b.rmsd_a, equal_nan=True)
check("G3a pre-existing docking rows unchanged", same2, f"{a.shape} vs {b.shape}")
ls = d_new[(d_new.module=="ligand")&(d_new.engine=="surfdock")]
check("G3b ligand SurfDock ingested", len(ls)>=1260, f"{len(ls)} rows, {ls.system.nunique()} systems, ok={ (ls.status=='ok').sum() }")
print("GATE", "PASS" if ok else "FAIL")
