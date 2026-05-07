"""Generate docking inputs (receptor PDB + ligand SDF + box) for all 16 cases.

For each case, we take the AF3 prediction (model_0) as the starting geometry:
  - strip everything except the protein chain                 -> receptor.pdb
  - build a 3D SDF from the case's target SMILES              -> ligand.sdf
  - centre the docking box on AF3's predicted ligand centroid -> box.json

This mirrors the experiment design approved in the plan: docking engines are
handed AF3's predicted receptor and the case-specific ligand, and we ask whether
they place the ligand where AF3 did (memorisation) or elsewhere (physics).

Run:
    LD_LIBRARY_PATH=/home/aoxu/miniconda3/envs/rdkit_env/lib:$LD_LIBRARY_PATH \
        /home/aoxu/miniconda3/envs/rdkit_env/bin/python \
        analysis/scripts/10_prep_docking_inputs.py
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

REPO_ROOT = Path("/mnt/katritch_lab2/aoxu/contrasCF")
sys.path.insert(0, str(REPO_ROOT / "analysis" / "src"))

from config import CASES, MODELS, case_dir, find_first  # noqa: E402
from loaders import read_structure, extract_protein_ca, select_target_ligand, STANDARD_AA  # noqa: E402
from native import load_native  # noqa: E402
from align import superpose_ca  # noqa: E402


DOCKING_ROOT = REPO_ROOT / "docking" / "inputs"
BOX_SIZE = (25.0, 25.0, 25.0)

# Pocket residue numbers used to centre the box when the AF3 ligand has been
# ejected. CDK2 (1B38) is hard-coded from the paper; GDH (2VWH) is derived
# dynamically from the native structure.
CDK2_POCKET_RESNUMS = [10, 14, 18, 31, 33, 86, 129, 131, 132, 134, 145]


def _gdh_pocket_resnums() -> list[int]:
    """Residue numbers within 6 Å of the native BGC ligand in 2VWH."""
    nat = load_native("gdh_glc")
    lig_xyz = nat.ligand.heavy_coords
    out = []
    for ch, num, name, xyz in nat.ca.residues:
        if np.min(np.linalg.norm(lig_xyz - xyz, axis=1)) <= 6.0:
            out.append(num)
    return out


def strip_to_protein_pdb(st: gemmi.Structure) -> str:
    """Return a PDB string with ONLY standard-AA residues of the largest chain."""
    model = st[0]
    best, best_n = None, 0
    for ch in model:
        n = sum(1 for r in ch if r.name in STANDARD_AA)
        if n > best_n:
            best, best_n = ch, n
    if best is None:
        raise RuntimeError("no protein chain found")

    lines = ["REMARK   1 receptor stripped for docking"]
    serial = 1
    chain_ch = best.name[0]
    for r in best:
        if r.name not in STANDARD_AA:
            continue
        for a in r:
            el = (a.element.name or "C").upper()
            name = a.name
            name_col = (" " + name[:3].ljust(3)) if len(el) == 1 else name[:4].ljust(4)
            lines.append(
                f"ATOM  {serial:5d} {name_col}"
                f" "
                f"{r.name[:3]:>3s}"
                f" {chain_ch}"
                f"{r.seqid.num:4d}"
                f"    "
                f"{a.pos.x:8.3f}{a.pos.y:8.3f}{a.pos.z:8.3f}"
                f"  1.00{a.b_iso:6.2f}"
                f"          "
                f"{el:>2s}"
            )
            serial += 1
    lines.append("TER")
    lines.append("END")
    return "\n".join(lines) + "\n"


def make_ligand_sdf(smiles: str, path: Path, seed: int = 42) -> None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RuntimeError(f"bad SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        # Retry with a non-random strategy.
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise RuntimeError("ETKDG embed failed")
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        pass
    mol.SetProp("_Name", path.stem)
    w = Chem.SDWriter(str(path))
    w.write(mol)
    w.close()


def _pocket_centroid(st, pocket_resnums) -> np.ndarray:
    """Cα centroid of the given residue numbers in `st` (AF3 frame)."""
    ca = extract_protein_ca(st)
    hits = [r[3] for r in ca.residues if r[1] in pocket_resnums]
    if not hits:
        raise RuntimeError(f"no pocket residues ({pocket_resnums}) found in receptor")
    return np.mean(np.stack(hits), axis=0)


def prep_case(case: str, gdh_pocket: list[int]) -> dict:
    """Prepare receptor + ligand + box for one case.

    Box centre strategy: we use the native crystal ligand's centroid mapped into
    the AF3 receptor frame. This gives a uniform "dock into the real binding
    pocket of the AF3 receptor" target across all 16 cases, independent of
    whether AF3 put its own ligand prediction in that pocket or ejected it into
    solvent. The mapping is done by Cα-superposing the AF3 protein onto the
    native crystal protein (`analysis.src.align.superpose_ca`) and then
    inverse-transforming the native ligand centroid.
    """
    spec = CASES[case]
    out = DOCKING_ROOT / case
    out.mkdir(parents=True, exist_ok=True)

    d = case_dir("AF3", case)
    cif_path = find_first(d, MODELS["AF3"].structure_globs)
    if cif_path is None:
        raise RuntimeError(f"AF3 CIF not found for {case}")
    st = read_structure(cif_path)

    (out / "receptor.pdb").write_text(strip_to_protein_pdb(st))

    # Align AF3 receptor Cα onto native Cα; use inverse to map native ligand
    # centroid into the AF3 frame.
    nat = load_native(spec.native_key)
    pred_ca = extract_protein_ca(st)
    sup = superpose_ca(pred_ca, nat.ca)
    nat_lig_cen = nat.ligand.heavy_coords.mean(axis=0)
    # superpose_ca.apply(pred) = pred @ R.T + t  (puts pred in native frame)
    # inverse: native -> AF3 frame is  (x - t) @ R
    cen = (nat_lig_cen - sup.t) @ sup.R

    # Also record diagnostics: where AF3 put its own ligand, and how far that is
    # from the native-based target centre (tells us if AF3 memorised the pocket).
    pocket_resnums = CDK2_POCKET_RESNUMS if spec.family != "glucose" else gdh_pocket
    pocket_cen = _pocket_centroid(st, pocket_resnums)
    af3_lig, dbg = select_target_ligand(st, spec.target_smiles)
    if af3_lig is None:
        raise RuntimeError(f"AF3 predicted ligand not found for {case}: {dbg}")
    af3_cen = af3_lig.heavy_coords.mean(axis=0)
    af3_dist = float(np.linalg.norm(af3_cen - cen))

    box = {
        "center": [float(x) for x in cen],
        "size": list(BOX_SIZE),
        "centre_source": "native_ligand_centroid_in_af3_frame",
        "af3_ligand_centroid": [float(x) for x in af3_cen],
        "pocket_ca_centroid_af3_frame": [float(x) for x in pocket_cen],
        "af3_to_target_dist": af3_dist,
        "ca_alignment_rmsd": float(sup.rmsd),
        "n_ca_aligned": int(sup.n_paired),
        "pocket_resnums": pocket_resnums,
        "source_cif": str(cif_path),
        "af3_ligand": {
            "chain": af3_lig.chain,
            "resname": af3_lig.resname,
            "n_heavy": int(len(af3_lig.heavy_coords)),
        },
    }
    (out / "box.json").write_text(json.dumps(box, indent=2))

    make_ligand_sdf(spec.target_smiles, out / "ligand.sdf")

    return {
        "case": case,
        "family": spec.family,
        "center": box["center"],
        "af3_to_target_dist": af3_dist,
        "ca_alignment_rmsd": float(sup.rmsd),
        "lig_chain": af3_lig.chain,
        "lig_resname": af3_lig.resname,
        "n_heavy_af3": int(len(af3_lig.heavy_coords)),
    }


def main():
    gdh_pocket = _gdh_pocket_resnums()
    print(f"[prep] GDH pocket residues (from 2VWH, within 6 Å of BGC): {gdh_pocket}")

    results = []
    for case in CASES:
        print(f"[prep] {case} ...", flush=True)
        try:
            info = prep_case(case, gdh_pocket)
            results.append(info)
            print(f"[prep]   ca_align_rmsd={info['ca_alignment_rmsd']:.2f} Å  "
                  f"af3_lig→target={info['af3_to_target_dist']:.1f} Å  "
                  f"n_heavy_af3={info['n_heavy_af3']}")
        except Exception as e:
            print(f"[prep]   FAIL: {type(e).__name__}: {e}")
            results.append({"case": case, "error": f"{type(e).__name__}: {e}"})

    # Summary manifest at the root of docking/inputs.
    (DOCKING_ROOT / "manifest.json").write_text(json.dumps(results, indent=2))
    fails = [r for r in results if "error" in r]
    print(f"[prep] done. {len(results)-len(fails)}/{len(results)} prepared, {len(fails)} failures.")
    if fails:
        for r in fails:
            print(f"         - {r['case']}: {r['error']}")


if __name__ == "__main__":
    main()
