#!/usr/bin/env python3
"""Export the pose-swap ligand poses for visual inspection: target protein +
native ligand + ejected ligands at increasing distance, all in the SAME frame.
This is the geometry both the Boltz and GNINA pose-swap tests apply (radial
ejection beyond the protein bounding sphere); use it to eyeball that the ligand
is genuinely pulled out of the pocket into solvent.

Writes (per system, under --out/<system>/):
  <sys>_target.pdb         the receptor (crystal protein)
  <sys>_ligand_poses.sdf   multi-model SDF: native + eject*, each titled with its
                           min ligand-protein distance
  <sys>_view.pml           PyMOL script (surface + all ligand poses as sticks)

Run (rdkit_env): .../rdkit_env/bin/python 23_export_poses.py --system 1bcu --out /tmp/poseswap_viz
"""
from __future__ import annotations
import argparse
import shutil
from pathlib import Path

import numpy as np
from rdkit import Chem

REPO = Path("/mnt/katritch_lab2/aoxu/contrasCF")
CLEAR = [0.0, 5.0, 15.0, 30.0]   # clearance beyond protein bounding sphere


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", required=True)
    ap.add_argument("--out", default="/tmp/poseswap_viz")
    args = ap.parse_args()
    D = REPO / "analysis/casf_mutagenesis/outputs" / args.system / "wt" / "docking"
    work = Path(args.out) / args.system
    work.mkdir(parents=True, exist_ok=True)
    shutil.copy(D / "receptor.pdb", work / f"{args.system}_target.pdb")

    mol = Chem.MolFromMolFile(str(D / "ligand.sdf"), sanitize=False)
    conf = mol.GetConformer()
    lig = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    P, c, R = protein_geom(str(D / "receptor.pdb"))
    lcom = lig.mean(0)
    outward = lcom - c
    outward = outward / np.linalg.norm(outward) if np.linalg.norm(outward) > 2 else np.array([0., 0., 1.])

    writer = Chem.SDWriter(str(work / f"{args.system}_ligand_poses.sdf"))
    print(f"=== {args.system} ligand poses ===")
    for cl in CLEAR:
        new = lig + ((c + outward * (R + cl)) - lcom) if cl > 0 else lig.copy()
        md = float(np.linalg.norm(new[:, None, :] - P[None, :, :], axis=2).min())
        m2 = Chem.Mol(mol); cf = m2.GetConformer()
        for i in range(m2.GetNumAtoms()):
            cf.SetAtomPosition(i, (float(new[i, 0]), float(new[i, 1]), float(new[i, 2])))
        name = "native" if cl == 0 else f"eject{int(cl)}"
        m2.SetProp("_Name", f"{name}_minD{md:.0f}A")
        m2.SetProp("min_lig_prot_A", f"{md:.1f}")
        writer.write(m2)
        print(f"  {name:8s} clearance={cl:4.0f}A  ->  min ligand-protein distance = {md:.1f} A")
    writer.close()

    pml = work / f"{args.system}_view.pml"
    pml.write_text(
        f"load {args.system}_target.pdb, target\n"
        f"load {args.system}_ligand_poses.sdf, poses\n"
        "hide everything\n"
        "show surface, target\nset transparency, 0.45\ncolor grey80, target\n"
        "set all_states, on\nshow sticks, poses\nutil.cbao poses\nset stick_radius, 0.25, poses\n"
        "bg_color white\norient\n"
        "# states (step with the state slider): 1=native(in pocket) -> 4=ejected 30 A clearance\n"
    )
    print(f"[wrote] {work}/  ({args.system}_target.pdb, _ligand_poses.sdf, _view.pml)")


if __name__ == "__main__":
    main()
