"""Hydroxyl methylation rule (paper Fig. 4).

Progressively methylates aliphatic free hydroxyls (-OH on sp3 C). Ring
oxygens and phenolic oxygens are excluded by SMARTS construction. Order of
methylation across the ladder is chemistry-aware so the paper's glucose
ladder reproduces exactly:

    1. hemiacetal OH first (anomeric in sugars — C also bonded to ring O)
    2. secondary alcohols next, in canonical atom-rank order
    3. primary alcohols last (e.g. the 6'-CH2OH of glucose)
"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem.rdchem import Atom, BondType, Mol, RWMol

from ..config import MAX_METHYLATIONS
from .base import LigandVariant


class MethylationRule:
    name = "methylation"

    # Aliphatic-attached free hydroxyl: -O-H bonded to an sp3 C.
    # Excludes ring oxygens (no implicit H) and phenolic OHs (sp2 C).
    SMARTS = Chem.MolFromSmarts("[OX2H][CX4]")
    # Hemiacetal: -O-H on a C that is also bonded to another oxygen.
    HEMIACETAL_SMARTS = Chem.MolFromSmarts("[OX2H][CX4][OX2,OX1-]")

    def is_eligible(self, mol: Mol) -> tuple[bool, str, int]:
        matches = mol.GetSubstructMatches(self.SMARTS)
        if not matches:
            return (False, "no aliphatic free hydroxyl", 0)
        return (True, "", len(matches))

    def generate(self, mol: Mol) -> list[LigandVariant]:
        matches = mol.GetSubstructMatches(self.SMARTS)
        if not matches:
            return []
        oh_atoms = sorted({m[0] for m in matches})
        ordered = _order_hydroxyls(mol, oh_atoms)

        n_max = min(len(ordered), MAX_METHYLATIONS)
        variants: list[LigandVariant] = []
        for k in range(1, n_max + 1):
            new_smi, applied = _methylate(mol, ordered[:k])
            if new_smi is None:
                continue
            variants.append(LigandVariant(
                name=f"meth_{k}",
                rule=self.name,
                smiles=new_smi,
                applied=applied,
                paper_faithful=True,
            ))
        return variants


def _order_hydroxyls(mol: Mol, o_atoms: list[int]) -> list[int]:
    """Return o_atoms sorted by chemistry-aware priority.

    Priority class (lowest first methylated):
      0 = hemiacetal OH (the C bonded to this O also has another O neighbour)
      1 = secondary alcohol (C bonded to ≥2 other carbons)
      2 = primary alcohol (C bonded to 1 other carbon)
      3 = methanol-like (C with no other C neighbours; rare)
    Tiebreaks within class:
      a. shortest-bond-path distance from the hemiacetal C (so secondary
         OHs are walked around the ring from the anomeric C — C2, C3, C4
         in pyranose order, matching paper Fig. 4)
      b. RDKit canonical atom rank
    """
    hemi_matches = mol.GetSubstructMatches(MethylationRule.HEMIACETAL_SMARTS)
    hemi_o_set = {m[0] for m in hemi_matches}
    hemi_c_idx: int | None = None
    if hemi_matches:
        hemi_c_idx = hemi_matches[0][1]
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=True))

    def key(o_idx: int) -> tuple[int, int, int]:
        atom = mol.GetAtomWithIdx(o_idx)
        c_neighbours = [n for n in atom.GetNeighbors() if n.GetSymbol() == "C"]
        if not c_neighbours:
            return (3, 0, ranks[o_idx])
        c_atom = c_neighbours[0]
        c_idx = c_atom.GetIdx()
        if hemi_c_idx is not None and c_idx != hemi_c_idx:
            path = rdmolops.GetShortestPath(mol, hemi_c_idx, c_idx)
            dist = len(path) - 1 if path else 99
        else:
            dist = 0
        if o_idx in hemi_o_set:
            return (0, dist, ranks[o_idx])
        n_c = sum(1 for n in c_atom.GetNeighbors() if n.GetSymbol() == "C")
        if n_c >= 2:
            return (1, dist, ranks[o_idx])
        return (2, dist, ranks[o_idx])

    return sorted(o_atoms, key=key)


def _methylate(mol: Mol, o_atoms: list[int]) -> tuple[str | None, list[str]]:
    rwmol = RWMol(mol)
    applied: list[str] = []
    for o_idx in o_atoms:
        o_atom = rwmol.GetAtomWithIdx(o_idx)
        o_atom.SetNumExplicitHs(0)
        o_atom.SetNoImplicit(True)
        c_new = rwmol.AddAtom(Atom(6))
        rwmol.AddBond(o_idx, c_new, BondType.SINGLE)
        applied.append(f"O@{o_idx}→OCH3")
    try:
        new_mol = rwmol.GetMol()
        Chem.SanitizeMol(new_mol)
        smi = Chem.MolToSmiles(new_mol, canonical=True, isomericSmiles=True)
        return smi, applied
    except Exception as exc:
        return None, [f"failed: {exc}"]
