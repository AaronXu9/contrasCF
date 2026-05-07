"""Phosphate / carboxylate charge-swap rule (paper Fig. 5).

Two ladders, mirroring the paper's hand-built ATP variants:

  Neutral-alkyl ladder   (anionic triphosphate → neutral alkyl tail):
    chrg_neu_methyl, chrg_neu_ethyl, chrg_neu_propyl

  Cationic-ammonium ladder (anionic triphosphate → +1 / +2 / +3 amine tail):
    chrg_pos_1, chrg_pos_2, chrg_pos_3

Algorithm (phosphate path — paper-faithful):
  1. Match `[#6][OX2;!R][P]` to identify the anchor C (kept), the bridging O
     (kept), and the first P (start of the leaving group).
  2. Cut the O-P bond. The molecule splits into the "kept" fragment (sugar +
     base) and the "leaving" fragment (entire polyphosphate chain).
  3. Delete the leaving fragment.
  4. Parse the replacement SMILES; atom 0 is the attachment point.
  5. Combine kept + replacement and bond bridging-O to replacement atom 0.

Algorithm (carboxylate path — informational extension, not in paper):
  Cut the alpha-C – carboxyl-C bond, drop the carboxyl group entirely, and
  attach the replacement fragment directly to the alpha C. Flagged as
  `paper_faithful=False`.
"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem.rdchem import BondType, Mol, RWMol

from .base import LigandVariant


# Replacement fragments. Atom 0 of each SMILES is bonded to the kept
# attachment atom (bridging O for phosphate path, alpha C for carboxylate).
# Verified by hand-decoding the paper's atp_charge_* SMILES.
_NEUTRAL_ALKYL = {
    "chrg_neu_methyl": "CC(C)(C)CCC",
    "chrg_neu_ethyl":  "CC(C)(C)CC(C)(C)C",
    "chrg_neu_propyl": "CC(C)(C)CC(C)(C)CC(C)(C)C",
}
_CATIONIC = {
    "chrg_pos_1": "C[N+](C)(C)C",
    "chrg_pos_2": "C[N+](C)(C)C[N+](C)(C)C",
    "chrg_pos_3": "C[N+](C)(C)C[N+](C)(C)C[N+](C)(C)C",
}


class ChargeSwapRule:
    name = "charge_swap"

    PHOSPHATE_SMARTS = Chem.MolFromSmarts("[#6][OX2;!R][P]")
    CARBOXYLATE_SMARTS = Chem.MolFromSmarts("[CX4][CX3](=O)[OX2H,OX1-]")

    def is_eligible(self, mol: Mol) -> tuple[bool, str, int]:
        m_phos = mol.GetSubstructMatches(self.PHOSPHATE_SMARTS)
        if m_phos:
            return (True, "phosphate", len(m_phos))
        m_carb = mol.GetSubstructMatches(self.CARBOXYLATE_SMARTS)
        if m_carb:
            return (True, "carboxylate", len(m_carb))
        return (False, "no phosphate or carboxylate", 0)

    def generate(self, mol: Mol) -> list[LigandVariant]:
        m_phos = mol.GetSubstructMatches(self.PHOSPHATE_SMARTS)
        if m_phos:
            anchor_c, attach_o, root = m_phos[0]
            attach_idx = attach_o
            paper_faithful = True
            cut_pair = (attach_o, root)
            applied_prefix = "triphosphate"
        else:
            m_carb = mol.GetSubstructMatches(self.CARBOXYLATE_SMARTS)
            if not m_carb:
                return []
            alpha_c, carboxyl_c, _, _ = m_carb[0]
            attach_idx = alpha_c
            paper_faithful = False
            cut_pair = (alpha_c, carboxyl_c)
            applied_prefix = "carboxylate"

        variants: list[LigandVariant] = []
        all_fragments = list(_NEUTRAL_ALKYL.items()) + list(_CATIONIC.items())
        for name, frag_smi in all_fragments:
            new_smi = _swap(mol, attach_idx, cut_pair, frag_smi)
            if new_smi is None:
                continue
            variants.append(LigandVariant(
                name=name,
                rule=self.name,
                smiles=new_smi,
                applied=[f"{applied_prefix}→{frag_smi}"],
                paper_faithful=paper_faithful,
            ))
        return variants


def _swap(
    mol: Mol,
    attach_idx: int,
    cut_pair: tuple[int, int],
    frag_smi: str,
) -> str | None:
    """Cut the bond at cut_pair, drop the leaving fragment, attach `frag_smi`.

    `attach_idx` is the original-mol index of the atom that stays in the kept
    fragment AND is the attachment point for the replacement. It must be one
    of the two atoms in `cut_pair`.
    """
    rwmol = RWMol(mol)
    a, b = cut_pair
    rwmol.RemoveBond(a, b)

    # Tag attach atom with a unique map number to track it through deletions.
    rwmol.GetAtomWithIdx(attach_idx).SetAtomMapNum(99)

    frags = Chem.GetMolFrags(rwmol, asMols=False, sanitizeFrags=False)
    if len(frags) != 2:
        return None
    drop_frag = next((f for f in frags if attach_idx not in f), None)
    if drop_frag is None:
        return None

    for idx in sorted(drop_frag, reverse=True):
        rwmol.RemoveAtom(idx)

    new_attach = None
    for a_obj in rwmol.GetAtoms():
        if a_obj.GetAtomMapNum() == 99:
            new_attach = a_obj.GetIdx()
            a_obj.SetAtomMapNum(0)
            break
    if new_attach is None:
        return None

    frag_mol = Chem.MolFromSmiles(frag_smi)
    if frag_mol is None:
        return None
    n_kept = rwmol.GetNumAtoms()
    combined = RWMol(Chem.CombineMols(rwmol.GetMol(), frag_mol))
    combined.AddBond(new_attach, n_kept, BondType.SINGLE)

    try:
        out = combined.GetMol()
        Chem.SanitizeMol(out)
        return Chem.MolToSmiles(out, canonical=True, isomericSmiles=True)
    except Exception:
        return None
