"""RMSD utilities."""
from __future__ import annotations
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from rdkit.Geometry import Point3D


def paired_rmsd(a: np.ndarray, b: np.ndarray) -> float:
    """Heavy-atom RMSD between two equal-length coord arrays (no fitting)."""
    if a.shape != b.shape or a.size == 0:
        return float("nan")
    diff = a - b
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def _mol_with_coords(template: Chem.Mol, heavy_xyz: np.ndarray) -> Chem.Mol:
    """Return a copy of `template` whose conformer places heavy atoms at heavy_xyz."""
    m = Chem.Mol(template)
    if m.GetNumConformers() == 0:
        conf = Chem.Conformer(m.GetNumAtoms())
        m.AddConformer(conf, assignId=True)
    conf = m.GetConformer()
    heavy_idx = [a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum() > 1]
    if len(heavy_idx) != len(heavy_xyz):
        raise ValueError(
            f"heavy-atom count mismatch: template has {len(heavy_idx)}, coords has {len(heavy_xyz)}"
        )
    for i, hidx in enumerate(heavy_idx):
        conf.SetAtomPosition(hidx, Point3D(*heavy_xyz[i]))
    return m


def bestfit_rmsd(
    pred_mol: Chem.Mol, nat_mol: Chem.Mol,
    pred_heavy_xyz: np.ndarray, nat_heavy_xyz: np.ndarray,
) -> float:
    """Paper-style ligand RMSD: symmetry-corrected optimal rigid-body fit.

    Uses rdkit's `rdMolAlign.GetBestRMS`, which:
      1. Enumerates all symmetric heavy-atom correspondences between mols.
      2. Applies an optimal Kabsch (rigid-body) superposition of the ligands.
      3. Returns the minimum heavy-atom RMSD over all those correspondences.

    This is *conformation-only*: it ignores where the ligand sits relative to
    the protein, answering "is the ligand's internal geometry correct?" rather
    than "is the ligand placed correctly in the pocket?". The Masters et al.
    paper uses this convention for its quoted WT values (e.g. AF3 0.2 Å).

    Returns NaN on any rdkit exception (e.g. mismatched atom counts).
    """
    if pred_mol is None or nat_mol is None:
        return float("nan")
    try:
        p = _mol_with_coords(pred_mol, pred_heavy_xyz)
        n = _mol_with_coords(nat_mol, nat_heavy_xyz)
        return float(rdMolAlign.GetBestRMS(p, n))
    except Exception:
        return float("nan")
