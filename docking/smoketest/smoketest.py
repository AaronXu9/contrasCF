"""One-case smoke test for the docking pipeline.

Reads AF3 `bindingsite_wt`, strips everything except the protein chain into a
receptor PDB, computes the AF3-predicted ligand centroid as the box centre,
builds an ATP SDF from SMILES via RDKit, and exits ready for docking runs.

Run with:
    LD_LIBRARY_PATH=/home/aoxu/miniconda3/envs/rdkit_env/lib:$LD_LIBRARY_PATH \
        /home/aoxu/miniconda3/envs/rdkit_env/bin/python docking/smoketest/smoketest.py
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# Make analysis.src importable.
REPO_ROOT = Path("/mnt/katritch_lab2/aoxu/contrasCF")
sys.path.insert(0, str(REPO_ROOT / "analysis" / "src"))

from config import CASES, MODELS, case_dir, find_first   # noqa: E402
from loaders import read_structure, extract_protein_ca, select_target_ligand, STANDARD_AA  # noqa: E402

CASE = "bindingsite_wt"
OUT = REPO_ROOT / "docking" / "inputs" / CASE
OUT.mkdir(parents=True, exist_ok=True)


def strip_to_protein_pdb(st: gemmi.Structure, chain_hint: str | None = None) -> str:
    """Return a PDB-format string containing ONLY standard-AA residues of one chain."""
    model = st[0]
    # Find the largest protein chain.
    best, best_n = None, 0
    for ch in model:
        n = sum(1 for r in ch if r.name in STANDARD_AA)
        if n > best_n:
            best, best_n = ch, n
    assert best is not None, "no protein chain found"

    lines = ["REMARK   1 receptor stripped for docking"]
    serial = 1
    for r in best:
        if r.name not in STANDARD_AA:
            continue
        for a in r:
            el = a.element.name.upper() or "C"
            # Proper PDB formatting
            name = a.name
            name_col = (" " + name[:3].ljust(3)) if len(el) == 1 else name[:4].ljust(4)
            lines.append(
                f"ATOM  {serial:5d} {name_col}"
                f" "
                f"{r.name[:3]:>3s}"
                f" {best.name[0]}"
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
    assert mol is not None, f"bad SMILES: {smiles}"
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    assert AllChem.EmbedMolecule(mol, params) == 0, "embed failed"
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        pass
    mol.SetProp("_Name", path.stem)
    w = Chem.SDWriter(str(path))
    w.write(mol)
    w.close()


def main():
    spec = CASES[CASE]
    d = case_dir("AF3", CASE)
    cif_path = find_first(d, MODELS["AF3"].structure_globs)
    print(f"[smoke] AF3 CIF: {cif_path}")
    assert cif_path is not None

    st = read_structure(cif_path)

    # Receptor PDB: protein only.
    receptor_pdb = strip_to_protein_pdb(st)
    (OUT / "receptor.pdb").write_text(receptor_pdb)
    n_atoms = receptor_pdb.count("\nATOM  ") + receptor_pdb.startswith("ATOM  ")
    print(f"[smoke] receptor.pdb written, {n_atoms} ATOM lines")

    # Box centre: AF3 predicted ligand centroid.
    lig, dbg = select_target_ligand(st, spec.target_smiles)
    assert lig is not None, f"failed to find predicted ligand: {dbg}"
    cen = lig.heavy_coords.mean(axis=0)
    print(f"[smoke] AF3 predicted ligand: chain={lig.chain} resname={lig.resname} n_heavy={len(lig.heavy_coords)}")
    print(f"[smoke] box centre: {cen.tolist()}")

    box = {
        "center": [float(x) for x in cen],
        "size": [25.0, 25.0, 25.0],
        "source_ligand": {"chain": lig.chain, "resname": lig.resname, "n_heavy": len(lig.heavy_coords)},
    }
    (OUT / "box.json").write_text(json.dumps(box, indent=2))
    print(f"[smoke] box.json written")

    # Ligand SDF from SMILES.
    make_ligand_sdf(spec.target_smiles, OUT / "ligand.sdf")
    print(f"[smoke] ligand.sdf written")

    print("[smoke] prep done. Ready to run docking.")
    print(f"[smoke] OUT = {OUT}")


if __name__ == "__main__":
    main()
