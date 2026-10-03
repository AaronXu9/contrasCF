"""Ligand atom pairing: native <-> predicted.

Two pairings are produced per (case, prediction):
  * `all`    : MCS over full molecules (the predicted ligand may be bigger or
               smaller than native; we match whatever atoms are in common).
  * `common` : restricted to a case-family SMARTS core so values are comparable
               across different variants within a family.
"""
from __future__ import annotations
import itertools

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFMCS


def _mcs_pairs(pred_mol: Chem.Mol, nat_mol: Chem.Mol, timeout: int = 10) -> tuple[list[int], list[int], object]:
    """Return heavy-atom index lists (pred, native) paired by MCS, and the query patt."""
    res = rdFMCS.FindMCS(
        [pred_mol, nat_mol],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        ringMatchesRingOnly=True,
        completeRingsOnly=True,
        timeout=timeout,
        matchValences=False,
    )
    if res.numAtoms == 0:
        return [], [], None
    patt = Chem.MolFromSmarts(res.smartsString)
    pred_match = pred_mol.GetSubstructMatch(patt)
    nat_match = nat_mol.GetSubstructMatch(patt)
    return list(pred_match), list(nat_match), patt


def _best_symmetric_pairing(
    pred_mol: Chem.Mol, nat_mol: Chem.Mol,
    pred_xyz: np.ndarray, nat_xyz: np.ndarray,
    query: Chem.Mol,
) -> tuple[list[int], list[int], float]:
    """Given SMARTS `query`, enumerate all symmetric matches on both sides and
    pick the (pred_hit, nat_hit) pair that minimizes aligned heavy RMSD."""
    pred_hits = list(pred_mol.GetSubstructMatches(query, useChirality=False))
    nat_hits = list(nat_mol.GetSubstructMatches(query, useChirality=False))
    if not pred_hits or not nat_hits:
        return [], [], float("nan")

    best = (None, None, float("inf"))
    for ph, nh in itertools.product(pred_hits, nat_hits):
        if len(ph) != len(nh):
            continue
        p = pred_xyz[list(ph)]
        n = nat_xyz[list(nh)]
        diff = p - n
        rmsd = float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))
        if rmsd < best[2]:
            best = (list(ph), list(nh), rmsd)
    return best[0] or [], best[1] or [], best[2]


def match_all(
    pred_mol: Chem.Mol, nat_mol: Chem.Mol,
    pred_xyz: np.ndarray, nat_xyz: np.ndarray,
) -> tuple[list[int], list[int], float]:
    """Full-molecule MCS pairing with symmetry-minimizing selection.

    Enumerates all symmetric matches of the MCS pattern on both predicted and
    native molecules and returns the (pred_hit, nat_hit) pair with the lowest
    pre-fitted RMSD on the already-aligned coordinates. This is important for
    atoms with permutational symmetry (e.g. phosphate oxygens) that a single
    GetSubstructMatch() would pick arbitrarily.
    """
    _, _, patt = _mcs_pairs(pred_mol, nat_mol)
    if patt is None:
        return [], [], float("nan")
    return _best_symmetric_pairing(pred_mol, nat_mol, pred_xyz, nat_xyz, patt)


def match_common(
    pred_mol: Chem.Mol, nat_mol: Chem.Mol,
    pred_xyz: np.ndarray, nat_xyz: np.ndarray,
    common_smarts: str,
) -> tuple[list[int], list[int], float]:
    """SMARTS-core pairing with symmetry enumeration; returns (pred_idx, nat_idx, best_rmsd_under_that_pairing)."""
    q = Chem.MolFromSmarts(common_smarts)
    if q is None:
        # Fall back to MCS if SMARTS is malformed.
        p, n = _mcs_pairs(pred_mol, nat_mol)
        if not p:
            return [], [], float("nan")
        rmsd = float(np.sqrt(np.mean(np.sum((pred_xyz[p] - nat_xyz[n]) ** 2, axis=1))))
        return p, n, rmsd
    return _best_symmetric_pairing(pred_mol, nat_mol, pred_xyz, nat_xyz, q)
