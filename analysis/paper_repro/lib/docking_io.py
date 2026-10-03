"""Shared I/O helpers for docking-style methods (UniDock2, GNINA, SurfDock).

Produces `combined_top1.pdb` = receptor.pdb (AF3-prepped) + top-1 docked pose as
a HETATM block (resname LIG, chain Z). This is the canonical structure file
consumed by the analysis pipeline for all docking-style methods.
"""
from __future__ import annotations
from pathlib import Path

from rdkit import Chem


def read_top_pose(sdf_path: Path) -> Chem.Mol:
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=False)
    for mol in suppl:
        if mol is not None:
            return mol
    raise RuntimeError(f"no parseable pose in {sdf_path}")


def mol_to_hetatm_lines(mol: Chem.Mol, chain: str = "Z", resname: str = "LIG",
                        resnum: int = 1, start_serial: int = 1) -> list[str]:
    """Serialize an RDKit mol's conformer as PDB HETATM lines."""
    lines = []
    conf = mol.GetConformer()
    serial = start_serial
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        pos = conf.GetAtomPosition(atom.GetIdx())
        el = atom.GetSymbol().upper()
        name = f"{el}{serial}"[:4]
        name_col = (" " + name[:3].ljust(3)) if len(el) == 1 else name[:4].ljust(4)
        lines.append(
            f"HETATM{serial:5d} {name_col}"
            f" "
            f"{resname[:3]:>3s}"
            f" {chain[0]}"
            f"{resnum:4d}"
            f"    "
            f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}"
            f"  1.00  0.00"
            f"          "
            f"{el:>2s}"
        )
        serial += 1
    return lines


def write_combined_top1(receptor_pdb: Path, poses_sdf: Path, out_dir: Path) -> Path:
    """Produce receptor + top-1 pose as a single PDB (HETATM ligand)."""
    receptor_text = receptor_pdb.read_text()
    lines = receptor_text.rstrip().splitlines()
    while lines and lines[-1].strip() in ("END", ""):
        lines.pop()
    if not lines or not lines[-1].startswith("TER"):
        lines.append("TER")
    last_serial = 0
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            try:
                last_serial = max(last_serial, int(line[6:11]))
            except ValueError:
                pass
    top_mol = read_top_pose(poses_sdf)
    lines.extend(mol_to_hetatm_lines(top_mol, start_serial=last_serial + 1))
    lines.append("END")
    out_pdb = out_dir / "combined_top1.pdb"
    out_pdb.write_text("\n".join(lines) + "\n")
    return out_pdb
