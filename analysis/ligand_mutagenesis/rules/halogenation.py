"""Aromatic halogenation rule (informational extension — NOT in the paper).

Adds a single F / Cl / Br substitution at one aromatic-CH position per
variant, picking the lowest-canonical-rank aromatic CH so the choice is
deterministic. Single substitution only — no combinatorial ladder.

Variants emitted: halo_F_1, halo_Cl_1, halo_Br_1.

This rule is `paper_faithful=False` — a chemistry-aware extension we
include because the user asked for it. Aromatic CH→CX preserves
aromaticity (RDKit handles it natively); no kekulization concern.
"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem.rdchem import Atom, BondType, Mol, RWMol

from ..config import HALOGENS
from .base import LigandVariant


_HALOGEN_Z = {"F": 9, "Cl": 17, "Br": 35, "I": 53}


class HalogenationRule:
    name = "halogenation"

    SMARTS = Chem.MolFromSmarts("[c;H1]")

    def is_eligible(self, mol: Mol) -> tuple[bool, str, int]:
        matches = mol.GetSubstructMatches(self.SMARTS)
        if not matches:
            return (False, "no aromatic CH", 0)
        return (True, "", len(matches))

    def generate(self, mol: Mol) -> list[LigandVariant]:
        matches = mol.GetSubstructMatches(self.SMARTS)
        if not matches:
            return []
        ar_atoms = sorted({m[0] for m in matches})
        ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=True))
        ar_atoms.sort(key=lambda i: ranks[i])
        target_idx = ar_atoms[0]

        variants: list[LigandVariant] = []
        for x in HALOGENS:
            new_smi = _halogenate(mol, target_idx, x)
            if new_smi is None:
                continue
            variants.append(LigandVariant(
                name=f"halo_{x}_1",
                rule=self.name,
                smiles=new_smi,
                applied=[f"aromatic-C@{target_idx}H→{x}"],
                paper_faithful=False,
            ))
        return variants


def _halogenate(mol: Mol, c_idx: int, halogen: str) -> str | None:
    rwmol = RWMol(mol)
    c_atom = rwmol.GetAtomWithIdx(c_idx)
    c_atom.SetNumExplicitHs(0)
    c_atom.SetNoImplicit(True)
    x = rwmol.AddAtom(Atom(_HALOGEN_Z[halogen]))
    rwmol.AddBond(c_idx, x, BondType.SINGLE)
    try:
        out = rwmol.GetMol()
        Chem.SanitizeMol(out)
        return Chem.MolToSmiles(out, canonical=True, isomericSmiles=True)
    except Exception:
        return None
