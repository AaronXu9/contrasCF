"""Extract protein chain + ATP mol2 from a crystal CIF.

Used to seed `00_build_systems.py`:
    extract_inputs.py --cif analysis/native/1B38.cif --chain A \\
        --het ATP --out-dir <pilot_inputs>/

Outputs:
    <out_dir>/<pdb_id>_chain<C>.pdb   chain-A protein only (no waters/het)
    <out_dir>/<het>.pdb               raw het record (gemmi)
    <out_dir>/<het>.mol2              ready for ligand_params.parameterize

Requires gemmi (in boltzina_env / rdkit_env) and `obabel` (in AmberTools env).
For the mol2 conversion specifically, run this from the AmberTools env so
obabel is on PATH; gemmi-only PDBs are written first so the protein side
works even if obabel is unavailable.
"""
from __future__ import annotations
import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import gemmi


def extract_chain_protein(cif_path: Path, chain_id: str, out_pdb: Path) -> Path:
    s = gemmi.read_structure(str(cif_path))
    s.setup_entities()
    out = gemmi.Structure()
    out.spacegroup_hm = s.spacegroup_hm
    out.cell = s.cell
    model = gemmi.Model(s[0].name)
    src = s[0][chain_id]
    new_chain = gemmi.Chain(chain_id)
    for residue in src:
        # keep only canonical amino acids
        if residue.het_flag == "H":
            continue
        if residue.name == "HOH":
            continue
        new_chain.add_residue(residue)
    model.add_chain(new_chain)
    out.add_model(model)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    out.write_pdb(str(out_pdb))
    return out_pdb


def extract_het_pdb(cif_path: Path, het_name: str, chain_id: str | None,
                    out_pdb: Path) -> Path:
    """Pull all atoms of the specified HET residue (e.g. 'ATP') into a PDB.

    If `chain_id` is None, searches ALL chains (HET ligands are commonly in
    a separate chain from the protein — e.g. AF3 puts ligands in chain B
    while protein is in A). Pass an explicit chain_id only when there are
    multiple HETs of the same name in different chains and you need to pick
    a specific one.
    """
    s = gemmi.read_structure(str(cif_path))
    out = gemmi.Structure()
    out.cell = s.cell
    model = gemmi.Model(s[0].name)
    new_chain = gemmi.Chain("X")
    found = 0
    for chain in s[0]:
        if chain_id is not None and chain.name != chain_id:
            continue
        for residue in chain:
            if residue.name.strip().upper() == het_name.upper():
                new_chain.add_residue(residue)
                found += 1
    if found == 0:
        scope = f"chain {chain_id}" if chain_id else "any chain"
        raise RuntimeError(f"HET {het_name!r} not found in {scope} of {cif_path}")
    model.add_chain(new_chain)
    out.add_model(model)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    out.write_pdb(str(out_pdb))
    return out_pdb


def pdb_to_mol2(pdb_path: Path, mol2_path: Path) -> Path:
    """Convert PDB ligand to mol2 via openbabel. Adds hydrogens at pH 7.4."""
    if shutil.which("obabel") is None:
        raise FileNotFoundError(
            "`obabel` not on PATH. Activate $CONTRASCF_AMBERTOOLS_ENV first."
        )
    subprocess.run([
        "obabel", str(pdb_path), "-O", str(mol2_path),
        "-p", "7.4",  # add hydrogens at physiological pH
    ], check=True)
    return mol2_path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--cif", type=Path, required=True)
    p.add_argument("--chain", default="A",
                   help="protein chain to extract (default A)")
    p.add_argument("--het", default="ATP", help="HET residue name to extract")
    p.add_argument("--het-chain", default=None,
                   help="restrict HET search to this chain; default is all chains")
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--pdb-id", default=None,
                   help="defaults to CIF stem (e.g. 1B38)")
    p.add_argument("--skip-mol2", action="store_true",
                   help="only emit PDBs; convert with obabel later")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    pdb_id = (args.pdb_id or args.cif.stem).upper()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    protein_pdb = args.out_dir / f"{pdb_id}_chain{args.chain}.pdb"
    extract_chain_protein(args.cif, args.chain, protein_pdb)
    print(f"[extract] protein  -> {protein_pdb}")

    het_pdb = args.out_dir / f"{args.het}.pdb"
    extract_het_pdb(args.cif, args.het, args.het_chain, het_pdb)
    print(f"[extract] het PDB  -> {het_pdb}")

    if not args.skip_mol2:
        het_mol2 = args.out_dir / f"{args.het}.mol2"
        pdb_to_mol2(het_pdb, het_mol2)
        print(f"[extract] het mol2 -> {het_mol2}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
