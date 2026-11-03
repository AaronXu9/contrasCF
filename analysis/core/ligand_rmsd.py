"""Ligand atom correspondence and RMSD — shared by co-folding and docking scoring.

Moved verbatim from casf_mutagenesis/analysis.py (co-folding path) and
casf_mutagenesis/gnina_analysis.py (docking path) in S3 (2026-10). The two matchers
are deliberately NOT unified yet: the docking matcher takes a single match without
symmetry enumeration, so unifying it changes docking RMSDs (data_prep_todo item 10).
Co-location makes the divergence visible; a silent fallback in only one of them is
what corrupted every ligand-variant RMSD until 2026-10-01 (item 16).
"""
from __future__ import annotations

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFMCS


# ---- co-folding path (from casf_mutagenesis/analysis.py) -------------------

def _heavy_coords(mol: Chem.Mol) -> np.ndarray:
    conf = mol.GetConformer()
    out = []
    for i, a in enumerate(mol.GetAtoms()):
        if a.GetAtomicNum() == 1:
            continue
        p = conf.GetAtomPosition(i)
        out.append((p.x, p.y, p.z))
    return np.asarray(out, dtype=float)


def _heavy_indices(mol: Chem.Mol) -> list[int]:
    return [i for i, a in enumerate(mol.GetAtoms()) if a.GetAtomicNum() != 1]


_CORR_CACHE: dict = {}


def _topology_key(m: Chem.Mol) -> tuple:
    """Exact, index-preserving topology fingerprint (coords excluded).

    Two mols with the same key have identical atom indexing, so a cached index
    correspondence is valid for both. A cell's poses share topology, so the
    expensive MCS runs once per cell instead of 15 times."""
    atoms = tuple((a.GetAtomicNum(), a.GetIsAromatic(), a.GetFormalCharge())
                  for a in m.GetAtoms())
    bonds = tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx(),
                          b.GetBondTypeAsDouble()) for b in m.GetBonds()))
    return atoms, bonds


def _atom_correspondences(
    crystal_heavy: Chem.Mol, pred_heavy: Chem.Mol, max_matches: int = 200,
) -> tuple[list[tuple[list[int], list[int]]], str]:
    key = (_topology_key(crystal_heavy), _topology_key(pred_heavy), max_matches)
    hit = _CORR_CACHE.get(key)
    if hit is None:
        hit = _atom_correspondences_uncached(crystal_heavy, pred_heavy, max_matches)
        _CORR_CACHE[key] = hit
    return hit


def _atom_correspondences_uncached(
    crystal_heavy: Chem.Mol, pred_heavy: Chem.Mol, max_matches: int = 200,
) -> tuple[list[tuple[list[int], list[int]]], str]:
    """All symmetry-equivalent heavy-atom correspondences crystal ↔ pred.

    Returns ([(crystal_idx, pred_idx), ...], mode). Each pair of equal-length
    index lists says crystal_heavy atom crystal_idx[k] corresponds to
    pred_heavy atom pred_idx[k]. mode is one of
      "pred_superset"    crystal is a substructure of pred (identical ligand,
                         or a variant that ADDS atoms: halogenation, methylation)
      "crystal_superset" pred is a substructure of crystal (variant REMOVES atoms)
      "mcs"              neither; maximum common substructure (charge swaps)
      "none"             no usable correspondence.

    Replaces logic that, once the two ligands differed in heavy-atom count,
    discarded every valid substructure match and paired atoms by FILE ORDER
    while reporting a full match. That silently corrupted every ligand-
    mutagenesis variant (top-1 RMSD ~5-6 Å on a one-atom fluorine swap whose
    true RMSD is ~0.5 Å). Postmortem: journal/2026-10-01-cofold-ligand-variant-
    rmsd-mapping.md. Identical-ligand cells (the whole protein arm) take the
    "pred_superset" path with equal sizes and get exactly the matches they got
    before.
    """
    nc, npd = crystal_heavy.GetNumAtoms(), pred_heavy.GetNumAtoms()
    kw = dict(useChirality=False, uniquify=False, maxMatches=max_matches)
    m = pred_heavy.GetSubstructMatches(crystal_heavy, **kw)
    if m:
        return [(list(range(nc)), list(t)) for t in m], "pred_superset"
    m = crystal_heavy.GetSubstructMatches(pred_heavy, **kw)
    if m:
        return [(list(t), list(range(npd))) for t in m], "crystal_superset"
    res = rdFMCS.FindMCS(
        [crystal_heavy, pred_heavy],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        timeout=10, completeRingsOnly=False,
    )
    if res.numAtoms < 3 or not res.smartsString:
        return [], "none"
    patt = Chem.MolFromSmarts(res.smartsString)
    if patt is None:
        return [], "none"
    cm = crystal_heavy.GetSubstructMatches(patt, **kw)
    pm = pred_heavy.GetSubstructMatch(patt)
    if not cm or not pm:
        return [], "none"
    # Fix the pred side, enumerate crystal symmetry: covers the equivalent
    # mappings without a combinatorial product.
    return [(list(c), list(pm)) for c in cm], "mcs"


def _matched_rmsd(
    crystal_mol: Chem.Mol, pred_xyz: np.ndarray, pred_mol: Chem.Mol,
) -> tuple[float, int]:
    """Return (RMSD, n_matched_heavy) between the crystal ligand and the
    predicted ligand whose heavy-atom coords have been transformed into the
    crystal frame. In-place: no superposition of the ligand itself.

    Minimum over all symmetry-equivalent correspondences from
    `_atom_correspondences`, computed on the matched atoms only, so variant
    ligands with extra or missing atoms are scored on their shared scaffold.
    Raises ValueError when no correspondence exists, so the cell is marked
    `error` rather than given a meaningless number.
    """
    if pred_mol is None:
        # Legacy path: no RDKit mol for the prediction. Same-order assumption
        # is only meaningful for identical ligands; refuse otherwise.
        crystal_xyz = _heavy_coords(crystal_mol)
        if len(crystal_xyz) != len(pred_xyz):
            raise ValueError("no RDKit mol for prediction and heavy-atom counts differ")
        d = pred_xyz - crystal_xyz
        return float(np.sqrt((d * d).sum() / len(d))), len(d)

    crystal_heavy = Chem.RemoveHs(crystal_mol)
    pred_heavy = Chem.RemoveHs(pred_mol)
    crystal_pts = _heavy_coords(crystal_heavy)
    pred_pts = _heavy_coords(pred_heavy)
    pairs, mode = _atom_correspondences(crystal_heavy, pred_heavy)
    if not pairs:
        raise ValueError("no atom correspondence between crystal and predicted ligand")
    best, n_best = float("inf"), 0
    for ci, pi in pairs:
        d = pred_pts[pi] - crystal_pts[ci]
        r = float(np.sqrt((d * d).sum() / len(d)))
        if r < best:
            best, n_best = r, len(ci)
    return best, n_best


def _bestfit_rmsd(
    crystal_mol: Chem.Mol, pred_heavy_xyz: np.ndarray, pred_mol: Chem.Mol,
) -> float:
    """Symmetry-corrected Kabsch RMSD between predicted and crystal ligand.

    Same correspondences as `_matched_rmsd`, but each candidate mapping gets an
    independent rigid-body Kabsch superposition of the matched predicted atoms
    onto the matched crystal atoms; returns the minimum. Paper-style
    "conformation-only" RMSD: pocket-blind. NaN when no correspondence exists.
    """
    crystal_heavy = Chem.RemoveHs(crystal_mol)
    crystal_xyz = _heavy_coords(crystal_heavy)
    if pred_mol is None:
        if len(crystal_xyz) != len(pred_heavy_xyz):
            return float("nan")
        return _kabsch_rmsd(pred_heavy_xyz, crystal_xyz)
    pred_heavy = Chem.RemoveHs(pred_mol)
    pred_xyz = _heavy_coords(pred_heavy)
    pairs, _ = _atom_correspondences(crystal_heavy, pred_heavy)
    if not pairs:
        return float("nan")
    return min(_kabsch_rmsd(pred_xyz[pi], crystal_xyz[ci]) for ci, pi in pairs)


def _kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """Minimum RMSD of P onto Q under a rigid-body (rotation + translation)."""
    if P.shape != Q.shape or len(P) < 3:
        return float("nan")
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    H = Pc.T @ Qc
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    fitted = Pc @ R.T
    diff = fitted - Qc
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


# ---- docking path (from casf_mutagenesis/gnina_analysis.py) ----------------

def _mcs_match_indices(
    crystal: Chem.Mol, pose: Chem.Mol
) -> tuple[list[int], list[int], int]:
    """Find the maximum common substructure and return (crystal_idx, pose_idx, n)."""
    # Quick path: if pose is a superstructure of crystal (halogenation,
    # methylation add atoms), GetSubstructMatch is fast.
    direct = pose.GetSubstructMatch(crystal)
    if direct:
        return list(range(crystal.GetNumAtoms())), list(direct), crystal.GetNumAtoms()
    # Reverse direction (crystal contains pose — unusual but try).
    direct = crystal.GetSubstructMatch(pose)
    if direct:
        return list(direct), list(range(pose.GetNumAtoms())), pose.GetNumAtoms()
    # Fall back to MCS (handles charge_swap etc. where neither is a
    # substructure of the other).
    res = rdFMCS.FindMCS(
        [crystal, pose],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        timeout=10,
        completeRingsOnly=False,
    )
    if res.numAtoms == 0 or res.smartsString == "":
        return [], [], 0
    pattern = Chem.MolFromSmarts(res.smartsString)
    if pattern is None:
        return [], [], 0
    crystal_match = crystal.GetSubstructMatch(pattern)
    pose_match = pose.GetSubstructMatch(pattern)
    if not crystal_match or not pose_match:
        return [], [], 0
    return list(crystal_match), list(pose_match), len(crystal_match)
