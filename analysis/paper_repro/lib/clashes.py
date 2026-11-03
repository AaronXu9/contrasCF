"""Protein-ligand steric clash detection: heavy-atom pairs closer than cutoff."""
from __future__ import annotations
import numpy as np


def count_clashes(
    protein_heavy: np.ndarray, ligand_heavy: np.ndarray, cutoff: float = 2.0
) -> tuple[int, float]:
    """Return (n_clash_pairs, min_distance)."""
    if protein_heavy.size == 0 or ligand_heavy.size == 0:
        return 0, float("nan")
    diff = protein_heavy[None, :, :] - ligand_heavy[:, None, :]  # (L, P, 3)
    d = np.sqrt(np.sum(diff * diff, axis=-1))
    n = int(np.sum(d < cutoff))
    return n, float(d.min())
