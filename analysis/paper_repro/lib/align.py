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


_OFFSET_RANGE = range(-150, 151)
_MIN_PAIRS = 10


def _pair_by_resnum(
    pred: ProteinCA, native: ProteinCA, offset: int
) -> tuple[np.ndarray, np.ndarray]:
    """Pair Cα atoms by residue number with `offset` added to native numbering."""
    native_by_num = {nref[1] + offset: nref[3] for nref in native.residues}
    pred_xyz, nat_xyz = [], []
    for _, num, _, xyz in pred.residues:
        nxyz = native_by_num.get(num)
        if nxyz is not None:
            pred_xyz.append(xyz)
            nat_xyz.append(nxyz)
    return np.array(pred_xyz), np.array(nat_xyz)


def _aa_identity_fraction(
    pred: ProteinCA, native: ProteinCA, offset: int,
) -> float:
    """Fraction of paired residues whose AA3 names match at this offset.

    At the correct offset, predicted and native sequences should match
    almost everywhere (except at deliberately mutated pocket residues).
    At a residue-shifted offset, most pairs are different residue types
    so this fraction drops sharply — even when the pair COUNT is identical.
    """
    native_by_num = {nref[1] + offset: nref[2] for nref in native.residues}
    matches = total = 0
    for _, num, name, _ in pred.residues:
        nname = native_by_num.get(num)
        if nname is None:
            continue
        total += 1
        if nname == name:
            matches += 1
    return matches / total if total else 0.0


def _svd_fit(
    pred_xyz: np.ndarray, nat_xyz: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Run SVD fit of pred onto nat. Returns (R, t, rmsd).

    SVDSuperimposer.set(reference, moving) fits `moving` onto `reference`.
    Our convention `transformed = coord @ R.T + t` matches Biopython's
    `transformed = moving @ rot + tran` with R = rot.T and t = tran.
    """
    sup = SVDSuperimposer()
    sup.set(nat_xyz, pred_xyz)
    sup.run()
    rot, tran = sup.get_rotran()
    return np.asarray(rot).T, np.asarray(tran), float(sup.get_rms())


def superpose_ca(pred: ProteinCA, native: ProteinCA) -> Superposition:
    """Find the best resnum offset + Cα superposition in one pass.

    Picking by RMSD alone fails on conformationally-divergent predictions
    (e.g. Boltz-1 on MEK1: a small mis-paired subset happens to align with
    lower RMSD than the correct ~278-pair alignment). Picking by count
    alone fails on systems with internal gaps where a boundary offset ties
    at the same count but mis-shifts the pairings (the 7XLP case).
    Strategy: only accept offsets whose overlap is within 90 % of the
    maximum-overlap offset, then rank that band by SVD RMSD.
    """
    pred_nums = {r[1] for r in pred.residues}
    nat_nums = {r[1] for r in native.residues}
    min_pairs_required = max(_MIN_PAIRS, len(pred_nums) // 4)

    counts = {off: sum(1 for x in nat_nums if (x + off) in pred_nums)
              for off in _OFFSET_RANGE}
    max_count = max(counts.values()) if counts else 0
    overlap_floor = max(min_pairs_required, max_count - 2)
    candidates = [off for off, n in counts.items() if n >= overlap_floor]

    # 7XLP's internal gaps (221-224 + 277-306) make multiple offsets tie at
    # max-overlap; RMSD alone then picks shifted-but-coincidentally-tight
    # alignments on conformationally-divergent predictions (Boltz-1 on MEK1).
    # Disambiguate by AA3 identity — at the correct offset paired residues
    # share their type almost everywhere; at a shifted offset they don't.
    if len(candidates) > 1:
        id_fracs = {off: _aa_identity_fraction(pred, native, off) for off in candidates}
        max_id = max(id_fracs.values())
        if max_id > 0.5:
            candidates = [off for off, f in id_fracs.items() if f >= max_id - 0.05]

    best: tuple[float, int, np.ndarray, np.ndarray, float] | None = None
    for off in candidates:
        pred_xyz, nat_xyz = _pair_by_resnum(pred, native, off)
        if len(pred_xyz) < _MIN_PAIRS:
            continue
        try:
            R, t, rmsd = _svd_fit(pred_xyz, nat_xyz)
        except (np.linalg.LinAlgError, ValueError):
            continue
        key = (rmsd, -counts[off])
        if best is None or key[0] < best[0] or (key[0] == best[0] and key[1] < best[1]):
            best = (rmsd, off, R, t, float(len(pred_xyz)))

    if best is None:
        # No offset gave a usable fit. Fall back to the offset with the
        # highest overlap (uses counts already computed, no re-scan).
        off = max(counts, key=counts.get)
        pred_xyz, nat_xyz = _pair_by_resnum(pred, native, off)
        if len(pred_xyz) < _MIN_PAIRS:
            raise RuntimeError(
                f"too few Cα pairs ({len(pred_xyz)}) at any offset")
        R, t, rmsd = _svd_fit(pred_xyz, nat_xyz)
        return Superposition(R=R, t=t, rmsd=rmsd, n_paired=len(pred_xyz),
                              offset=off)

    rmsd, offset, R, t, n_paired = best
    return Superposition(R=R, t=t, rmsd=rmsd, n_paired=int(n_paired),
                          offset=offset)
