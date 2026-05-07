"""Single source of truth for CASF-2016 mutagenesis pipeline."""
from __future__ import annotations
from pathlib import Path

REPO_ROOT = Path("/mnt/katritch_lab2/aoxu/contrasCF")
MODULE_ROOT = REPO_ROOT / "analysis" / "casf_mutagenesis"
OUTPUT_ROOT = MODULE_ROOT / "outputs"
NATIVE_DIR = REPO_ROOT / "analysis" / "native"

CASF_ROOT = REPO_ROOT / "data" / "casf2016"  # symlink → pdbbind_cleansplit
CASF_RAW = CASF_ROOT / "raw"
CASF_LIGANDS = CASF_ROOT / "crystal_ligands"
CASF_LABELS = CASF_ROOT / "labels"
SUBSET20_JSON = CASF_LABELS / "PDBbind_casf2016_subset20.json"
DATA_DICT_JSON = CASF_LABELS / "PDBbind_data_dict.json"
SPLIT_JSON = CASF_LABELS / "PDBbind_data_split_cleansplit.json"
SMILES_CACHE = CASF_ROOT / "smiles" / "het_smiles_cache.json"

AF3_TEMPLATE = REPO_ROOT / "contrasCF/Cofolding-Tools-main/examples/af_input/fold_input.json"
BOLTZ_TEMPLATE = REPO_ROOT / "contrasCF/Cofolding-Tools-main/examples/boltz_input/example.yaml"

POCKET_CUTOFF_A = 3.5  # Å, side-chain heavy-atom → ligand heavy-atom
DOCKING_BOX_A = (25.0, 25.0, 25.0)

VARIANTS = ("wt", "rem", "pack", "inv")

# Backbone atoms excluded when measuring "side-chain heavy atom" distances.
BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O", "OXT", "H", "HA"})

# Standard 20 amino acids (3-letter); used to filter pocket residues.
AA3_STANDARD = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
})

AA3_TO_AA1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # Common modified-residue mappings (for sequence reconstruction; will not
    # be auto-mutated — pocket detection skips non-standard residues).
    "MSE": "M", "HSD": "H", "HSE": "H", "HSP": "H", "CSO": "C", "SEC": "U",
    "PYL": "O", "HYP": "P",
}

# Inversion table — Miyata-distance-maximising substitution per WT residue.
# Source: Masters et al. 2025, Table 2.
INV_TABLE = {
    "TYR": "G", "TRP": "G", "LYS": "G", "ARG": "G",
    "VAL": "D", "LEU": "D", "ILE": "D", "MET": "D", "PHE": "D",
    "ALA": "W", "GLY": "W", "PRO": "W", "HIS": "W", "ASP": "W",
    "GLU": "W", "CYS": "W", "ASN": "W", "GLN": "W", "THR": "W",
    "SER": "W",
}
assert AA3_STANDARD <= set(INV_TABLE), "INV_TABLE must cover all 20 standard AAs"

REM_AA1 = "G"   # all-glycine
PACK_AA1 = "F"  # all-phenylalanine
