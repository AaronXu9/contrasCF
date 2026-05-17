#!/bin/bash
# Create the AmberTools conda env on USC CARC.
#
# AmberTools provides antechamber, parmchk2, tleap — needed to derive GAFF2
# parameters and AM1-BCC charges for non-standard ligands (ATP charge
# variants, methylated glucoses). acpype converts the AMBER topology to
# GROMACS format.
#
# Run on a CARC interactive shell:
#     bash $CONTRASCF_ROOT/analysis/ramd_pilot/carc_setup/setup_ambertools_env.sh
#
# Verify:
#     conda activate $CONTRASCF_AMBERTOOLS_ENV
#     antechamber -h | head -3
#     acpype --help | head -3
set -euo pipefail

ENV_PATH=${CONTRASCF_AMBERTOOLS_ENV:-/project2/katritch_223/aoxu/conda/envs/ambertools}

if [[ -d "$ENV_PATH" ]]; then
    echo "[setup_ambertools] env already exists at $ENV_PATH"
    exit 0
fi

# Resolve the conda base. ENV_PATH lives at <user_envs_root>/envs/ambertools,
# but the actual conda install is the CARC system miniforge — not under
# /project2/.../conda/. Ask conda directly.
if ! command -v conda >/dev/null 2>&1; then
    echo "[setup_ambertools] conda not on PATH. module load conda first?"
    exit 1
fi
CONDA_ROOT=$(conda info --base 2>/dev/null)
if [[ -z "$CONDA_ROOT" || ! -f "$CONDA_ROOT/etc/profile.d/conda.sh" ]]; then
    echo "[setup_ambertools] could not resolve conda base via 'conda info --base'."
    exit 1
fi
echo "[setup_ambertools] conda base: $CONDA_ROOT"
source "$CONDA_ROOT/etc/profile.d/conda.sh"

# Prefer mamba over conda — ambertools pulls a large dep tree and the
# classic conda solver routinely takes 15-25 min on CARC. mamba cuts
# that to 2-3 min. Falls back to conda if mamba isn't on PATH.
SOLVER=conda
if command -v mamba >/dev/null 2>&1; then
    SOLVER=mamba
fi
echo "[setup_ambertools] solver: $SOLVER"

conda create -y -p "$ENV_PATH" python=3.11
conda activate "$ENV_PATH"

# AmberTools 23+ from conda-forge bundles antechamber, parmchk2, tleap.
# Also add the analysis stack so this env can run scripts/04_analyze.py
# and tests/test_smoke.py without needing a separate Python env.
"$SOLVER" install -y -p "$ENV_PATH" -c conda-forge ambertools=23 acpype openbabel rdkit gemmi \
    numpy scipy pandas matplotlib biopython

echo "[setup_ambertools] DONE. To use: conda activate $ENV_PATH"
