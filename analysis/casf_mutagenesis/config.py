"""Single source of truth for CASF-2016 mutagenesis pipeline.

All host-specific paths are read from environment variables with lab-default
fallbacks. To run on a different host (e.g. CARC), set the env vars listed in
`env/carc.sh` (or `env/lab.sh`) before invoking any script. See
[`docs/casf_mutagenesis.md`](../../docs/casf_mutagenesis.md#running-on-carc)
for the full list and the recommended workflow.
"""
from __future__ import annotations
import os
from pathlib import Path


def _path_env(name: str, default: str | Path) -> Path:
    """Return Path from $name if set, else `default`. Always Path-typed."""
    val = os.environ.get(name)
    return Path(val) if val else Path(default)


# --- project root ---------------------------------------------------------

REPO_ROOT = _path_env("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF")
MODULE_ROOT = REPO_ROOT / "analysis" / "casf_mutagenesis"
OUTPUT_ROOT = MODULE_ROOT / "outputs"
NATIVE_DIR = REPO_ROOT / "analysis" / "native"
ANALYSIS_SRC = REPO_ROOT / "analysis" / "src"

# --- CASF dataset ---------------------------------------------------------
# Default points at the in-repo symlink `<repo>/data/casf2016`. To override
# (e.g. on a host where the dataset lives elsewhere on disk), set
# CONTRASCF_CASF_ROOT to the pdbbind_cleansplit directory directly.

CASF_ROOT = _path_env("CONTRASCF_CASF_ROOT", REPO_ROOT / "data" / "casf2016")
CASF_RAW = CASF_ROOT / "raw"
CASF_LIGANDS = CASF_ROOT / "crystal_ligands"
CASF_LABELS = CASF_ROOT / "labels"
SUBSET20_JSON = CASF_LABELS / "PDBbind_casf2016_subset20.json"
DATA_DICT_JSON = CASF_LABELS / "PDBbind_data_dict.json"
SPLIT_JSON = CASF_LABELS / "PDBbind_data_split_cleansplit.json"
SMILES_CACHE = CASF_ROOT / "smiles" / "het_smiles_cache.json"

# --- input templates (live in the repo; no override needed) ---------------

AF3_TEMPLATE = REPO_ROOT / "contrasCF/Cofolding-Tools-main/examples/af_input/fold_input.json"
BOLTZ_TEMPLATE = REPO_ROOT / "contrasCF/Cofolding-Tools-main/examples/boltz_input/example.yaml"

# --- tool locations -------------------------------------------------------

BOLTZ_BIN = _path_env(
    "CONTRASCF_BOLTZ_BIN",
    "/home/aoxu/miniconda3/envs/boltzina_env/bin/boltz",
)

# AF3 v3.0.1 + dockstrat-managed env layout.
AF3_ENV = _path_env(
    "CONTRASCF_AF3_ENV",
    "/mnt/katritch_lab2/aoxu/CogLigandBench/envs/alphafold3",
)
AF3_DIR = _path_env(
    "CONTRASCF_AF3_DIR",
    "/mnt/katritch_lab2/aoxu/CogLigandBench/forks/alphafold3/alphafold3",
)
AF3_MODEL_DIR = _path_env(
    "CONTRASCF_AF3_MODEL_DIR",
    "/mnt/katritch_lab2/aoxu/CogLigandBench/forks/alphafold3/models",
)

# Override which CUDA device the runners pin to. Default 0 = first GPU.
CUDA_DEVICE = os.environ.get("CONTRASCF_CUDA_DEVICE", "0")

# --- algorithm parameters --------------------------------------------------

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
