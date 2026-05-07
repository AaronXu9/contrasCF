"""Ligand SMILES discovery for the CASF dataset.

Thin wrapper around `casf_mutagenesis.build._smiles_from_crystal_sdf` so the
two modules share one canonicalization path. Cross-import keeps the SMILES
that lands in this module's outputs identical to the SMILES that the
binding-site sibling module emits for the same system.
"""
from __future__ import annotations
import json
from pathlib import Path

from casf_mutagenesis.build import _smiles_from_crystal_sdf

from .config import DATA_DICT_JSON, SMILES_CACHE


_DATA_DICT_CACHE: dict | None = None
_SMILES_CACHE: dict | None = None


def data_dict() -> dict:
    global _DATA_DICT_CACHE
    if _DATA_DICT_CACHE is None:
        _DATA_DICT_CACHE = json.loads(DATA_DICT_JSON.read_text())
    return _DATA_DICT_CACHE


def het_smiles_cache() -> dict[str, str]:
    global _SMILES_CACHE
    if _SMILES_CACHE is None:
        _SMILES_CACHE = json.loads(SMILES_CACHE.read_text())
    return _SMILES_CACHE


def resolve_wt_smiles(ligand_sdf: Path, het_code: str | None) -> tuple[str | None, str | None]:
    """Return (smiles, source) where source ∈ {'crystal_sdf', 'het_cache', None}."""
    smi = _smiles_from_crystal_sdf(ligand_sdf)
    if smi is not None:
        return smi, "crystal_sdf"
    if het_code:
        smi = het_smiles_cache().get(het_code)
        if smi is not None:
            return smi, "het_cache"
    return None, None
