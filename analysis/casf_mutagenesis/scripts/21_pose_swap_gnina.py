#!/usr/bin/env python3
"""GNINA reference for the pose-swap test: rescore the SAME ligand-ejection
ladder with `gnina --score_only`. A genuine physics scorer should collapse as
the ligand leaves the pocket — the steep (GNINA) vs flat (Boltz-2, scripts
19/20) contrast lands the conclusion that Boltz-2's affinity head ignores pose.

Per system: crystal docking inputs (receptor.pdb + ligand.sdf); eject the ligand
radially beyond the protein bounding sphere (clearance 5/15/30 A, matching
19_pose_swap_affinity.py); gnina --score_only each pose; record CNNaffinity (pKd,
higher = tighter) + Vina (kcal/mol, more negative = tighter) + min lig-protein
distance. Writes gnina_<system>.csv.

Run (rdkit_env): .../rdkit_env/bin/python 21_pose_swap_gnina.py --system 1bcu --out /tmp/gnina_ps/1bcu
"""
from __future__ import annotations
import argparse
import csv
import re
import subprocess
from pathlib import Path

import numpy as np
from rdkit import Chem

GNINA = "/mnt/katritch_lab2/aoxu/envs/gnina/bin/gnina"
CLEAR = [0.0, 5.0, 15.0, 30.0]
REPO = Path("/mnt/katritch_lab2/aoxu/contrasCF")


def protein_geom(pdb):
    xyz = []
    for line in open(pdb):
        if line.startswith(("ATOM", "HETATM")):
            try:
                xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                pass
    P = np.array(xyz)
    return P, P.mean(0), float(np.linalg.norm(P - P.mean(0), axis=1).max())


def gnina_score(receptor, sdf):
    out = subprocess.run([GNINA, "-r", receptor, "-l", sdf, "--score_only"],
                         capture_output=True, text=True, timeout=300).stdout

    def g(k):
        m = re.search(rf"{k}:\s*(-?\d+\.?\d*)", out)
        return float(m.group(1)) if m else float("nan")
    return g("Affinity"), g("CNNaffinity"), g("CNNscore")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    D = REPO / "analysis/casf_mutagenesis/outputs" / args.system / "wt" / "docking"
    rec = str(D / "receptor.pdb")
    mol = Chem.MolFromMolFile(str(D / "ligand.sdf"), sanitize=False)
    if mol is None or mol.GetNumConformers() == 0:
        print(f"[gnina] {args.system}: unreadable ligand.sdf — skip"); return
    conf = mol.GetConformer()
    lig = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    P, c, R = protein_geom(rec)
    lcom = lig.mean(0)
    outward = lcom - c
    outward = outward / np.linalg.norm(outward) if np.linalg.norm(outward) > 2 else np.array([0., 0., 1.])

    work = Path(args.out); work.mkdir(parents=True, exist_ok=True)
    rows = []
    for cl in CLEAR:
        new = lig + ((c + outward * (R + cl)) - lcom) if cl > 0 else lig.copy()
        md = float(np.linalg.norm(new[:, None, :] - P[None, :, :], axis=2).min())
        m2 = Chem.Mol(mol); cf = m2.GetConformer()
        for i in range(m2.GetNumAtoms()):
            cf.SetAtomPosition(i, (float(new[i, 0]), float(new[i, 1]), float(new[i, 2])))
        psdf = work / f"{args.system}_eject{int(cl)}.sdf"
        Chem.MolToMolFile(m2, str(psdf), kekulize=False)
        vina, cnnaff, cnnsc = gnina_score(rec, str(psdf))
        rows.append(dict(system=args.system, pose="native" if cl == 0 else f"eject{int(cl)}",
                         clearance=cl, min_lig_prot=round(md, 1), vina=vina, cnnaffinity=cnnaff, cnnscore=cnnsc))
    with open(work / f"gnina_{args.system}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    nat = rows[0]
    print(f"=== GNINA {args.system} (CNNaffinity pKd higher=tighter; Vina kcal/mol) ===")
    for r in rows:
        print(f"  {r['pose']:8s} minD={r['min_lig_prot']:5.1f}A  CNNaff={r['cnnaffinity']:.2f}  Vina={r['vina']:+.2f}")
    print(f"  native->eject30: CNNaff drop = {nat['cnnaffinity']-rows[-1]['cnnaffinity']:+.2f} pK ; "
          f"Vina change = {rows[-1]['vina']-nat['vina']:+.2f} kcal/mol")


if __name__ == "__main__":
    main()
