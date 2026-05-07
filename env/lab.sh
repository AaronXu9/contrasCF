# Lab workstation environment for contrasCF.
# Source before running any analysis/casf_mutagenesis script:
#     source env/lab.sh
#
# These match the defaults baked into config.py but exporting them
# explicitly is cheap and makes the active host visible in `env | grep
# CONTRASCF`.

export CONTRASCF_ROOT=/mnt/katritch_lab2/aoxu/contrasCF
export CONTRASCF_CASF_ROOT=/home/aoxu/projects/VLS-Benchmark-Dataset/data/pdbbind_cleansplit

# --- tool binaries / envs --------------------------------------------------
export CONTRASCF_BOLTZ_BIN=/home/aoxu/miniconda3/envs/boltzina_env/bin/boltz
export CONTRASCF_AF3_ENV=/mnt/katritch_lab2/aoxu/CogLigandBench/envs/alphafold3
export CONTRASCF_AF3_DIR=/mnt/katritch_lab2/aoxu/CogLigandBench/forks/alphafold3/alphafold3
export CONTRASCF_AF3_MODEL_DIR=/mnt/katritch_lab2/aoxu/CogLigandBench/forks/alphafold3/models

export CONTRASCF_CUDA_DEVICE=0

# Python: needs rdkit_env's lib on LD_LIBRARY_PATH (gemmi/PyMOL libstdc++ quirk).
export LD_LIBRARY_PATH=/home/aoxu/miniconda3/envs/rdkit_env/lib:${LD_LIBRARY_PATH:-}
export CONTRASCF_PY=/home/aoxu/miniconda3/envs/rdkit_env/bin/python

echo "[contrasCF] env: lab workstation. Use \$CONTRASCF_PY for analysis scripts."
