"""Load cached native references (1B38 CDK2-ATP, 2VWH GDH-BGC)."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from rdkit import Chem

from config import NATIVES, NATIVE_DIR, CASES
from loaders import (
    read_structure, extract_protein_ca, select_target_ligand,
    extract_all_protein_heavy, ProteinCA, LigandBlock,
)


@dataclass
class NativeLoaded:
    pdb_id: str
    ca: ProteinCA
    protein_chain: str
    protein_heavy: np.ndarray
    ligand: LigandBlock


_CACHE: dict[str, NativeLoaded] = {}


def load_native(native_key: str) -> NativeLoaded:
    if native_key in _CACHE:
        return _CACHE[native_key]
    nref = NATIVES[native_key]
    st = read_structure(NATIVE_DIR / f"{nref.pdb_id}.cif")
    ca = extract_protein_ca(st)
    protein_heavy = extract_all_protein_heavy(st, ca.residues[0][0])
    # Use the target SMILES of any case in this family as the template for bond-order
    # assignment on the native ligand. For natives we want a clean parent form:
    # ATP for cdk2_atp, BGC/glucose for gdh_glc.
    case_for_template = next(c for c, sp in CASES.items()
                             if sp.native_key == native_key and sp.family != "atp_charge" and sp.family != "glucose" or
                                (sp.native_key == native_key and sp.name in ("bindingsite_wt","glucose_0")))
    template_smiles = CASES[case_for_template].target_smiles
    lig, _ = select_target_ligand(
        st, template_smiles,
        target_is_native=True, native_resname=nref.ligand_resname,
    )
    if lig is None:
        raise RuntimeError(f"native ligand {nref.ligand_resname} not found in {nref.pdb_id}")
    loaded = NativeLoaded(
        pdb_id=nref.pdb_id, ca=ca,
        protein_chain=ca.residues[0][0],
        protein_heavy=protein_heavy,
        ligand=lig,
    )
    _CACHE[native_key] = loaded
    return loaded
