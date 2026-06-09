#!/usr/bin/env python3
"""Export Boltz-PREDICTED-frame pose-swap structures for visual inspection: the
predicted protein + the predicted native ligand pose + the predicted ligand
ejected radially (same geometry the pose-swap test applies), all in the
predicted frame. Faithful to what the affinity pose-swap actually scored (the
Boltz prediction), vs 23_export_poses.py which uses the crystal docking inputs.

Reads the predicted complex CIF from a pose-swap panel run:
  <panel>/<sys>/boltz_results_<sys>/predictions/<sys>/<sys>_model_0.cif
Writes (under --out/<sys>/):
  <sys>_pred_target.pdb         predicted protein
  <sys>_pred_ligand_poses.pdb   multi-MODEL: native + ejected ligand (HETATM)
  <sys>_pred_view.pml

Run (rdkit_env, has gemmi): .../rdkit_env/bin/python 24_export_predicted_poses.py --system 1bcu
"""
from __future__ import annotations
import argparse
from pathlib import Path

import numpy as np
import gemmi

AA = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL "
         "MSE SEC PYL".split())
CLEAR = [0.0, 5.0, 15.0, 30.0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", required=True)
    ap.add_argument("--panel", default="/tmp/poseswap_panel2")
    ap.add_argument("--out", default="/tmp/poseswap_viz_pred")
    args = ap.parse_args()
    cif = next(Path(args.panel, args.system).glob(
        f"boltz_results_{args.system}/predictions/{args.system}/{args.system}_model_0.cif"))
    st = gemmi.read_structure(str(cif))
    st.setup_entities()

    lig_atoms, lig_orig, prot_xyz = [], [], []
    for ch in st[0]:
        for res in ch:
            if res.name == "HOH":
                continue
            is_lig = res.name not in AA
            for at in res:
                if is_lig:
                    lig_atoms.append(at); lig_orig.append((at.pos.x, at.pos.y, at.pos.z))
                else:
                    prot_xyz.append((at.pos.x, at.pos.y, at.pos.z))
    P = np.array(prot_xyz); c = P.mean(0); R = float(np.linalg.norm(P - c, axis=1).max())
    L0 = np.array(lig_orig); lcom = L0.mean(0)
    outward = lcom - c
    outward = outward / np.linalg.norm(outward) if np.linalg.norm(outward) > 2 else np.array([0., 0., 1.])

    work = Path(args.out, args.system); work.mkdir(parents=True, exist_ok=True)
    prot = st.clone(); prot.remove_ligands_and_waters(); prot.remove_empty_chains()
    prot.write_pdb(str(work / f"{args.system}_pred_target.pdb"))

    print(f"=== {args.system} (predicted frame; {len(lig_atoms)} ligand atoms) ===")
    models = []
    for cl in CLEAR:
        shift = ((c + outward * (R + cl)) - lcom) if cl > 0 else np.zeros(3)
        for at, (ox, oy, oz) in zip(lig_atoms, lig_orig):
            at.pos = gemmi.Position(ox + shift[0], oy + shift[1], oz + shift[2])
        new = L0 + shift
        md = float(np.linalg.norm(new[:, None, :] - P[None, :, :], axis=2).min())
        tmp = work / "_tmp.pdb"; st.write_pdb(str(tmp))
        het = [ln for ln in tmp.read_text().splitlines() if ln.startswith("HETATM")]
        models.append(het); tmp.unlink()
        print(f"  {'native' if cl == 0 else f'eject{int(cl)}':8s} clearance={cl:4.0f} A -> min ligand-protein distance = {md:.1f} A")

    with open(work / f"{args.system}_pred_ligand_poses.pdb", "w") as f:
        for mi, het in enumerate(models, 1):
            f.write(f"MODEL     {mi:>4d}\n"); f.write("\n".join(het) + "\n"); f.write("ENDMDL\n")
        f.write("END\n")

    (work / f"{args.system}_pred_view.pml").write_text(
        f"load {args.system}_pred_target.pdb, target\n"
        f"load {args.system}_pred_ligand_poses.pdb, poses\n"
        "hide everything\nshow surface, target\nset transparency, 0.45\ncolor grey80, target\n"
        "set all_states, on\nshow sticks, poses\nutil.cbag poses\nset stick_radius, 0.25, poses\n"
        "bg_color white\norient\n"
        "# states 1..4: native (in pocket) -> ejected ~30 A clearance (predicted frame)\n")
    print(f"[wrote] {work}/  (predicted-frame: target.pdb + multi-model ligand poses + pml)")


if __name__ == "__main__":
    main()
