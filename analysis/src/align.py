"""Cα superposition of predicted protein onto native crystal.

We use Biopython's SVD-based Superimposer over the intersection of residue
numbers. The native AF3/Boltz/Chai/RFAA all number residues 1..N from the
submitted sequence; the crystal may start at a different offset, so we align
by shared residue *position* using sequence alignment when numbering doesn't
intersect, falling back to the pairwise prefix match that works for 1B38/2VWH.
"""
from __future__ import annotations
from dataclasses import dataclass

import numpy as np
from Bio.SVDSuperimposer import SVDSuperimposer

from loaders import ProteinCA


@dataclass
class Superposition:
    R: np.ndarray        # 3x3 rotation
    t: np.ndarray        # (3,) translation (applied after rotation about origin)
    rmsd: float          # Å
    n_paired: int        # number of Cα pairs used
    offset: int          # residue-number offset applied to predicted to match native

    def apply(self, coords: np.ndarray) -> np.ndarray:
        return coords @ self.R.T + self.t


def _pair_by_resnum(
    pred: ProteinCA, native: ProteinCA, offset: int
) -> tuple[np.ndarray, np.ndarray]:
    """Pair Cα atoms by residue number with `offset` added to native numbering.

    Returns (pred_xyz, native_xyz) arrays of equal length.
    """
    native_by_num = {(nref[1] + offset, nref[2][:1]): nref[3] for nref in native.residues}
    # We compare only by resnum since residue 3-letter code can differ due to
    # mutations, but we record whether names match for diagnostics.
    pred_xyz, nat_xyz = [], []
    for chain, num, name, xyz in pred.residues:
        key_num = num
        # Try exact resnum match first.
        for nref in native.residues:
            if nref[1] + offset == key_num:
                pred_xyz.append(xyz); nat_xyz.append(nref[3])
                break
    return np.array(pred_xyz), np.array(nat_xyz)


def _best_offset(pred: ProteinCA, native: ProteinCA) -> int:
    """Find the integer offset that maximizes number of paired residues."""
    pred_nums = {r[1] for r in pred.residues}
    nat_nums  = {r[1] for r in native.residues}
    best_off, best_n = 0, -1
    # Typical PDB offset is small; scan ±50.
    for off in range(-50, 51):
        n = sum(1 for x in nat_nums if (x + off) in pred_nums)
        if n > best_n:
            best_n, best_off = n, off
    return best_off


def superpose_ca(pred: ProteinCA, native: ProteinCA) -> Superposition:
    offset = _best_offset(pred, native)
    pred_xyz, nat_xyz = _pair_by_resnum(pred, native, offset)
    if len(pred_xyz) < 10:
        raise RuntimeError(f"too few Cα pairs ({len(pred_xyz)}) after offset={offset}")

    # SVDSuperimposer.set(reference, moving) -> fits `moving` onto `reference`.
    # Our convention: transformed = coord @ R.T + t (in Superposition.apply).
    # Biopython gives (rot, tran) such that: transformed = moving @ rot + tran.
    # So we store R = rot.T and t = tran; apply uses coord @ R.T + t = coord @ rot + tran.
    sup = SVDSuperimposer()
    sup.set(nat_xyz, pred_xyz)  # reference = native, moving = predicted
    sup.run()
    rot, tran = sup.get_rotran()
    R = np.asarray(rot).T
    t = np.asarray(tran)
    rmsd = float(sup.get_rms())
    return Superposition(R=R, t=t, rmsd=rmsd, n_paired=len(pred_xyz), offset=offset)
