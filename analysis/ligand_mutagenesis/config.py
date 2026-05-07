"""Single source of truth for the ligand-mutation pipeline.

Path constants are reused from `casf_mutagenesis.config` (same CASF-2016
dataset, same Boltz-2 binary, same docking-box parameters) so that any host
override (env var, CARC bootstrap) propagates uniformly across both modules.
"""
from __future__ import annotations
import os
from pathlib import Path

from casf_mutagenesis.config import (  # noqa: F401  (re-exports)
    AF3_DIR,
    AF3_ENV,
    AF3_MODEL_DIR,
    AF3_TEMPLATE,
    BOLTZ_BIN,
    BOLTZ_TEMPLATE,
    CASF_LABELS,
    CASF_LIGANDS,
    CASF_RAW,
    CASF_ROOT,
    CUDA_DEVICE,
    DATA_DICT_JSON,
    DOCKING_BOX_A,
    REPO_ROOT,
    SMILES_CACHE,
    SPLIT_JSON,
    SUBSET20_JSON,
    assert_boltz2_binary,
)


def _path_env(name: str, default: str | Path) -> Path:
    val = os.environ.get(name)
    return Path(val) if val else Path(default)


# --- module-specific paths ------------------------------------------------

MODULE_ROOT = REPO_ROOT / "analysis" / "ligand_mutagenesis"
OUTPUT_ROOT = _path_env("CONTRASCF_LIGAND_OUT", MODULE_ROOT / "outputs")

# --- algorithm parameters --------------------------------------------------

# Cap the methylation ladder length so multi-OH ligands don't explode the
# manifest. Paper's glucose ladder caps at 5; 6 is a safe upper bound.
MAX_METHYLATIONS = 6

# Halogens emitted by the halogenation rule (single-substitution per halogen).
HALOGENS = ("F", "Cl", "Br")
