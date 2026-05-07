"""Identify pocket residues by side-chain heavy-atom proximity to ligand.

Implements the rule from Masters et al. 2025 §"Removing interacting residues":
    "Residues with side-chain heavy atoms within 3.5 Å of any ligand heavy
     atoms were selected for mutation."
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem

from .config import AA3_STANDARD, BACKBONE_ATOMS, POCKET_CUTOFF_A


@dataclass(frozen=True)
class PocketResidue:
    chain: str
    resnum: int
    ins_code: str
    aa3: str            # WT 3-letter; only standard AAs are returned


def load_ligand_heavy_coords(sdf_path: Path | str) -> np.ndarray:
    """Read ligand SDF (first molecule, take heavy atoms only)."""
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=True, sanitize=False)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        raise RuntimeError(f"could not parse ligand: {sdf_path}")
    conf = mol.GetConformer()
    coords = []
    for i, atom in enumerate(mol.GetAtoms()):
        if atom.GetAtomicNum() == 1:
            continue
        p = conf.GetAtomPosition(i)
        coords.append((p.x, p.y, p.z))
    if not coords:
        raise RuntimeError(f"no heavy atoms in ligand: {sdf_path}")
    return np.asarray(coords, dtype=float)


def detect_pocket(
    pdb_path: Path | str,
    ligand_sdf: Path | str,
    cutoff: float = POCKET_CUTOFF_A,
) -> list[PocketResidue]:
    """Pocket residues with any side-chain heavy atom within `cutoff` Å of any
    ligand heavy atom. Excludes non-standard residues and glycine."""
    lig_xyz = load_ligand_heavy_coords(ligand_sdf)
    return detect_pocket_from_coords(pdb_path, lig_xyz, cutoff=cutoff)


def detect_pocket_from_coords(
    pdb_path: Path | str,
    lig_xyz: np.ndarray,
    cutoff: float = POCKET_CUTOFF_A,
) -> list[PocketResidue]:
    """Coords-based variant of `detect_pocket` — takes pre-loaded ligand heavy
    coords instead of an SDF path. Useful when the ligand lives inside a CIF/
    PDB rather than a separate SDF (e.g. paper-reference systems).
    """
    cutoff2 = cutoff * cutoff
    st = gemmi.read_structure(str(pdb_path))
    model = st[0]
    out: list[PocketResidue] = []
    seen_keys = set()
    for chain in model:
        for res in chain:
            aa3 = res.name.strip().upper()
            if aa3 not in AA3_STANDARD or aa3 == "GLY":
                continue
            # collect side-chain heavy-atom coords
            sc_coords = []
            for a in res:
                if a.name in BACKBONE_ATOMS:
                    continue
                if a.element.name == "H":
                    continue
                sc_coords.append((a.pos.x, a.pos.y, a.pos.z))
            if not sc_coords:
                continue
            sc = np.asarray(sc_coords, dtype=float)
            # min squared distance to any ligand heavy atom
            d2 = np.sum((sc[:, None, :] - lig_xyz[None, :, :]) ** 2, axis=-1)
            if d2.min() <= cutoff2:
                key = (chain.name, res.seqid.num, res.seqid.icode.strip())
                if key in seen_keys:
                    continue
                seen_keys.add(key)
                out.append(PocketResidue(
                    chain=chain.name,
                    resnum=res.seqid.num,
                    ins_code=res.seqid.icode.strip(),
                    aa3=aa3,
                ))
    out.sort(key=lambda r: (r.chain, r.resnum, r.ins_code))
    return out
