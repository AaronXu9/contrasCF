"""Generate docking inputs (receptor.pdb + ligand.sdf + box.json) for the WT
variant of a CASF system. Mutant variants are deferred — they require an
AF3 prediction of the mutated sequence, which is downstream of this pipeline.

The receptor is the crystal protein (all standard-AA residues across all
chains preserved). The box is centred on the crystal ligand centroid with the
contrasCF default size (25 × 25 × 25 Å).
"""
from __future__ import annotations
import json
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from .config import AA3_STANDARD, DOCKING_BOX_A
from .pocket import load_ligand_heavy_coords


def _strip_to_protein_pdb(st: gemmi.Structure) -> str:
    """PDB string with only standard-AA residues from all protein chains.

    Differs from analysis/scripts/10_prep_docking_inputs.py:strip_to_protein_pdb
    which picks the single largest chain — here we keep all chains, since CASF
    ligand-binding sites can sit at chain interfaces.
    """
    lines = ["REMARK   1 receptor stripped for docking (all chains)"]
    serial = 1
    for ch in st[0]:
        any_aa = False
        for r in ch:
            if r.name not in AA3_STANDARD:
                continue
            any_aa = True
            for a in r:
                el = (a.element.name or "C").upper()
                name = a.name
                name_col = (" " + name[:3].ljust(3)) if len(el) == 1 else name[:4].ljust(4)
                lines.append(
                    f"ATOM  {serial:5d} {name_col}"
                    f" "
                    f"{r.name[:3]:>3s}"
                    f" {ch.name[0]}"
                    f"{r.seqid.num:4d}"
                    f"    "
                    f"{a.pos.x:8.3f}{a.pos.y:8.3f}{a.pos.z:8.3f}"
                    f"  1.00{a.b_iso:6.2f}"
                    f"          "
                    f"{el:>2s}"
                )
                serial += 1
        if any_aa:
            lines.append("TER")
    lines.append("END")
    return "\n".join(lines) + "\n"


def _embed_ligand(smiles: str, out_sdf: Path, seed: int = 42) -> None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RuntimeError(f"bad SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        # fallback: random coords
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol)
    mol = Chem.RemoveHs(mol)
    w = Chem.SDWriter(str(out_sdf))
    w.write(mol)
    w.close()


def render_docking_wt(
    pdb_path: Path | str,
    crystal_ligand_sdf: Path | str,
    target_smiles: str,
    out_dir: Path | str,
    seed: int = 42,
) -> None:
    """Write receptor.pdb, ligand.sdf, box.json for the WT variant.

    For WT, the ligand chemistry is unchanged from the crystal — use the
    crystal SDF directly (it has 3D coordinates already). Fall back to RDKit
    embedding from SMILES only if the crystal SDF cannot be parsed.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    st = gemmi.read_structure(str(pdb_path))
    (out_dir / "receptor.pdb").write_text(_strip_to_protein_pdb(st))

    crystal_ligand_sdf = Path(crystal_ligand_sdf)
    out_sdf = out_dir / "ligand.sdf"
    try:
        # crystal SDF has 3D coords already; just sanitize and copy
        suppl = Chem.SDMolSupplier(str(crystal_ligand_sdf), removeHs=False, sanitize=False)
        mol = next((m for m in suppl if m is not None), None)
        if mol is None:
            raise RuntimeError("could not parse crystal SDF")
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            pass  # tolerate sanitize failures (some PDBbind ligands are exotic)
        w = Chem.SDWriter(str(out_sdf))
        w.write(mol)
        w.close()
    except Exception:
        # last-ditch fallback: embed from SMILES
        _embed_ligand(target_smiles, out_sdf, seed=seed)

    lig_xyz = load_ligand_heavy_coords(crystal_ligand_sdf)
    centre = lig_xyz.mean(axis=0)
    box = {
        "center": [float(centre[0]), float(centre[1]), float(centre[2])],
        "size": list(DOCKING_BOX_A),
    }
    (out_dir / "box.json").write_text(json.dumps(box, indent=2))
