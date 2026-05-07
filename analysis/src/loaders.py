"""Loading predicted / native structures with gemmi; extracting ligands to RDKit.

Strategy for target-ligand selection in predictions
--------------------------------------------------
Predicted outputs may contain multiple non-protein residues (target ligand,
Mg/Zn ions, NADP cofactor, waters). We do NOT assume "chain B = target".
Instead, for each non-protein residue we:

1. Skip obvious non-targets: metals (MG, ZN, NA, CA, etc.), waters, AA-only
   residues, and NADP (CCD NAP, NDP, APR).
2. Count heavy atoms; keep residues whose heavy-atom count is close (±10%
   or ±3 atoms) to the target SMILES heavy-atom count.
3. Among remaining candidates, the one with the largest MCS against the
   target SMILES wins.

This works regardless of chain-ID differences across AF3 / Boltz / Chai / RFAA.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS

# Residue 3-letter codes to exclude outright when searching for a ligand.
METAL_RESNAMES = {
    "MG", "ZN", "CA", "NA", "K", "MN", "FE", "CU", "NI", "CO", "CD", "SR",
}
COFACTOR_RESNAMES = {
    "NAP", "NDP", "NAD", "NAI", "APR",   # NAD/NADP and analogs
    "FAD", "FMN",
    "SAH", "SAM",
    "COA",
    "ACE",   # acetyl cap (1B38)
}
WATER_RESNAMES = {"HOH", "H2O", "WAT", "DOD"}

# Standard amino-acid 3-letter codes (used to decide "is polymer").
STANDARD_AA = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS",
    "MET","PHE","PRO","SER","THR","TRP","TYR","VAL","MSE","SEC","PYL",
}


# ---------------------------------------------------------------------------

@dataclass
class ProteinCA:
    residues: list[tuple[str, int, str, np.ndarray]]   # (chain, resnum, resname, xyz)

    @property
    def coords(self) -> np.ndarray:
        return np.array([r[3] for r in self.residues], dtype=np.float64)


@dataclass
class LigandBlock:
    chain: str
    resname: str
    resnum: int
    atoms: list[tuple[str, str, np.ndarray, float]]   # (name, element, xyz, bfactor)
    mol: Chem.Mol | None                               # RDKit Mol with bond orders assigned, or None

    @property
    def heavy_coords(self) -> np.ndarray:
        return np.array(
            [a[2] for a in self.atoms if a[1].upper() != "H"],
            dtype=np.float64,
        )

    @property
    def heavy_elements(self) -> list[str]:
        return [a[1].upper() for a in self.atoms if a[1].upper() != "H"]

    @property
    def mean_bfactor(self) -> float:
        bs = [a[3] for a in self.atoms if a[1].upper() != "H"]
        return float(np.mean(bs)) if bs else float("nan")


def read_structure(path: Path) -> gemmi.Structure:
    """Read a CIF or PDB into a gemmi Structure."""
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    return st


def extract_protein_ca(st: gemmi.Structure, chain_hint: str | None = None) -> ProteinCA:
    """Return Cα coordinates of the (first) protein chain.

    If chain_hint is provided and that chain exists, use it; otherwise pick the
    chain with the largest number of polymer residues.
    """
    model = st[0]
    # Candidate chains ranked by polymer residue count.
    candidates = []
    for ch in model:
        n_poly = sum(1 for r in ch if r.name in STANDARD_AA)
        if n_poly > 20:
            candidates.append((n_poly, ch))
    if not candidates:
        raise RuntimeError("no protein chain found")
    candidates.sort(key=lambda x: -x[0])
    chain = None
    if chain_hint:
        for n, ch in candidates:
            if ch.name == chain_hint:
                chain = ch; break
    if chain is None:
        chain = candidates[0][1]

    residues = []
    for r in chain:
        if r.name not in STANDARD_AA:
            continue
        ca = None
        for atom in r:
            if atom.name == "CA":
                ca = atom; break
        if ca is None:
            continue
        pos = np.array([ca.pos.x, ca.pos.y, ca.pos.z], dtype=np.float64)
        residues.append((chain.name, r.seqid.num, r.name, pos))
    return ProteinCA(residues=residues)


def _collect_candidate_ligands(
    st: gemmi.Structure,
    exclude_resnames: frozenset[str] = frozenset(),
) -> list[LigandBlock]:
    """Enumerate HETATM-style residues that could be the target ligand."""
    out: list[LigandBlock] = []
    model = st[0]
    for chain in model:
        for r in chain:
            if r.name in STANDARD_AA:
                continue
            if r.name in WATER_RESNAMES:
                continue
            if r.name in METAL_RESNAMES:
                continue
            if r.name in COFACTOR_RESNAMES:
                continue
            if r.name in exclude_resnames:
                continue
            atoms = []
            for a in r:
                xyz = np.array([a.pos.x, a.pos.y, a.pos.z], dtype=np.float64)
                atoms.append((a.name, a.element.name, xyz, float(a.b_iso)))
            if len(atoms) < 3:
                continue
            out.append(LigandBlock(
                chain=chain.name,
                resname=r.name,
                resnum=r.seqid.num,
                atoms=atoms,
                mol=None,
            ))
    return out


def _write_ligand_pdb_block(lig: LigandBlock) -> str:
    """Serialize a LigandBlock to a PDB HETATM string for RDKit.

    Fixed-column PDB record:
      cols 1-6   record name (HETATM)
      cols 7-11  serial
      col  12    space
      cols 13-16 atom name (left-justified if starts with letter)
      col  17    altLoc (space)
      cols 18-20 resname
      col  21    space
      col  22    chain
      cols 23-26 resSeq
      cols 27-30 iCode + padding
      cols 31-38 x (8.3f)
      cols 39-46 y (8.3f)
      cols 47-54 z (8.3f)
      cols 55-60 occupancy (6.2f)
      cols 61-66 temp/B (6.2f)
      cols 77-78 element (right-justified)
    """
    lines = []
    # PDB atom names: 4 chars. If element is 1 char, name starts at col 14 (space in col 13).
    for i, (name, elem, xyz, b) in enumerate(lig.atoms, start=1):
        elem_u = (elem or "").upper()
        if not elem_u:
            # Derive element from atom name by stripping trailing digits.
            elem_u = "".join(ch for ch in name if ch.isalpha())[:2].upper()
        if len(elem_u) == 1:
            name_col = " " + name[:3].ljust(3)          # space + 3-char name
        else:
            name_col = name[:4].ljust(4)                # 4-char name
        chain_ch = (lig.chain or "A")[0]
        resn = lig.resname[:3].rjust(3)
        line = (
            f"HETATM{i:5d} {name_col}"                  # 1-12, 13-16
            f" "                                         # 17 altLoc
            f"{resn}"                                    # 18-20
            f" "                                         # 21
            f"{chain_ch}"                                # 22
            f"{lig.resnum:4d}"                           # 23-26
            f"    "                                      # 27-30
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"  # 31-54
            f"  1.00{b:6.2f}"                            # 55-66
            f"          "                                # 67-76
            f"{elem_u:>2s}"                              # 77-78
        )
        lines.append(line)
    lines.append("END")
    return "\n".join(lines)


def _ligand_to_mol(lig: LigandBlock, template_smiles: str | None) -> Chem.Mol | None:
    """Build an RDKit Mol from a LigandBlock, assigning bond orders from SMILES template."""
    pdb_block = _write_ligand_pdb_block(lig)
    mol = Chem.MolFromPDBBlock(pdb_block, removeHs=True, proximityBonding=True, sanitize=False)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE ^ Chem.SANITIZE_SETAROMATICITY)
    except Exception:
        pass
    # Force RingInfo initialization even when sanitize fails — downstream
    # GetSubstructMatch / MCS calls require it.
    try:
        mol.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(mol)
    except Exception:
        pass
    if template_smiles:
        try:
            tmpl = Chem.MolFromSmiles(template_smiles)
            if tmpl is not None:
                tmpl = Chem.RemoveHs(tmpl)
                mol_h = Chem.RemoveHs(mol)
                mol2 = AllChem.AssignBondOrdersFromTemplate(tmpl, mol_h)
                if mol2 is not None:
                    return mol2
        except Exception:
            pass
    return mol


def select_target_ligand(
    st: gemmi.Structure,
    target_smiles: str,
    target_is_native: bool = False,
    native_resname: str | None = None,
) -> tuple[LigandBlock | None, dict]:
    """Pick the HETATM residue most similar to `target_smiles`.

    If `target_is_native`, we first try to locate a residue matching
    `native_resname` directly (crystal CCD code like ATP / BGC).

    Returns: (best_ligand_block_or_None, debug_info).
    """
    candidates = _collect_candidate_ligands(st)
    debug = {"n_candidates": len(candidates), "cands": [(c.chain, c.resname, len(c.atoms)) for c in candidates]}

    if target_is_native and native_resname:
        for c in candidates:
            if c.resname == native_resname:
                c.mol = _ligand_to_mol(c, target_smiles)
                debug["picked_by"] = "resname"
                return c, debug

    # Heavy-atom count of target.
    tmpl = Chem.MolFromSmiles(target_smiles)
    if tmpl is None:
        return None, {**debug, "err": "bad template SMILES"}
    target_heavy = tmpl.GetNumHeavyAtoms()

    scored: list[tuple[float, LigandBlock, Chem.Mol | None]] = []
    for c in candidates:
        heavy = sum(1 for a in c.atoms if a[1].upper() != "H")
        size_diff = abs(heavy - target_heavy)
        # Soft size window; over-allow and score by MCS below.
        if heavy < 5:
            continue
        mol = _ligand_to_mol(c, target_smiles)
        mcs_n = 0
        if mol is not None:
            try:
                res = rdFMCS.FindMCS(
                    [mol, tmpl],
                    atomCompare=rdFMCS.AtomCompare.CompareElements,
                    bondCompare=rdFMCS.BondCompare.CompareAny,
                    ringMatchesRingOnly=False,
                    completeRingsOnly=False,
                    timeout=5,
                )
                mcs_n = res.numAtoms
            except Exception:
                mcs_n = 0
        # Score: big MCS is good; big size mismatch is bad.
        score = mcs_n - 0.5 * size_diff
        scored.append((score, c, mol))

    if not scored:
        return None, {**debug, "err": "no viable candidates after filtering"}

    scored.sort(key=lambda x: -x[0])
    top_score, best, best_mol = scored[0]
    best.mol = best_mol
    debug["picked_by"] = "mcs_score"
    debug["top_score"] = top_score
    debug["ranking"] = [(c.resname, c.chain, sc) for sc, c, _ in scored[:5]]
    return best, debug


def extract_all_protein_heavy(st: gemmi.Structure, chain: str) -> np.ndarray:
    """All heavy atom coords of protein chain (for clash detection)."""
    model = st[0]
    coords = []
    for ch in model:
        if ch.name != chain:
            continue
        for r in ch:
            if r.name not in STANDARD_AA:
                continue
            for a in r:
                if a.element.name.upper() == "H":
                    continue
                coords.append([a.pos.x, a.pos.y, a.pos.z])
    return np.array(coords, dtype=np.float64) if coords else np.zeros((0, 3))
