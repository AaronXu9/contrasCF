# USC CARC (discovery.usc.edu) environment for contrasCF.
# Source before running any analysis/casf_mutagenesis script:
#     source env/carc.sh
#
# The CARC clone of this repo lives at $CONTRASCF_ROOT and the analysis
# pipeline expects all the tool envs / dataset paths below.

export CONTRASCF_ROOT=/project2/katritch_223/aoxu/contrasCF
export CONTRASCF_CASF_ROOT=/project2/katritch_223/aoxu/projects/VLS-Benchmark-Dataset/data/pdbbind_cleansplit

# --- tool binaries / envs --------------------------------------------------
# Boltz-2 v2.2.1 lives in the CARC-shared boltzina_env conda env.
export CONTRASCF_BOLTZ_BIN=/project2/katritch_223/aoxu/conda/envs/boltzina_env/bin/boltz
# AF3 v3.0.1 + dockstrat-managed install.
export CONTRASCF_AF3_ENV=/project2/katritch_223/aoxu/envs/alphafold3
export CONTRASCF_AF3_DIR=/project2/katritch_223/aoxu/dockStrat/forks/alphafold3/alphafold3
export CONTRASCF_AF3_MODEL_DIR=/project2/katritch_223/aoxu/dockStrat/forks/alphafold3/models

export CONTRASCF_CUDA_DEVICE=0

# rdkit_env on CARC — verify the path before relying on this.
# If the rdkit_env isn't installed yet on this machine, install with:
#     conda env create -f env/rdkit_env.yml -p /project2/katritch_223/aoxu/conda/envs/rdkit_env
export _RDKIT_ENV=/project2/katritch_223/aoxu/conda/envs/rdkit_env
if [ -d "$_RDKIT_ENV" ]; then
    export LD_LIBRARY_PATH="$_RDKIT_ENV/lib:${LD_LIBRARY_PATH:-}"
    export CONTRASCF_PY="$_RDKIT_ENV/bin/python"
else
    echo "[contrasCF] WARNING: rdkit_env not found at $_RDKIT_ENV"
    echo "             Set CONTRASCF_PY to your rdkit/gemmi-equipped python before running."
fi
unset _RDKIT_ENV

echo "[contrasCF] env: USC CARC. CONTRASCF_ROOT=$CONTRASCF_ROOT"
