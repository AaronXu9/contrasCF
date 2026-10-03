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

# --- SurfDock on the lab box (these are also the code defaults) -----------
export CONTRASCF_SURFDOCK_ENV=/home/aoxu/miniconda3/envs/SurfDock
export CONTRASCF_SURFDOCK_DIR=/home/aoxu/projects/SurfDock
export CONTRASCF_SURFDOCK_WEIGHTS=/mnt/katritch_lab2/aoxu/CogLigandBench/forks/SurfDock/model_weights
export CONTRASCF_SURFDOCK_PRECOMPUTED=/home/aoxu/projects/precomputed/precomputed_arrays
export CONTRASCF_DOCKSTRAT_ROOT=/mnt/katritch_lab2/aoxu/CogLigandBench

# Python: needs rdkit_env's lib on LD_LIBRARY_PATH (gemmi/PyMOL libstdc++ quirk).
export LD_LIBRARY_PATH=/home/aoxu/miniconda3/envs/rdkit_env/lib:${LD_LIBRARY_PATH:-}
export CONTRASCF_PY=/home/aoxu/miniconda3/envs/rdkit_env/bin/python

# --- τ-RAMD pilot ---------------------------------------------------------
# Output root for the ramd_pilot module (per-system sim dirs + analysis).
export CONTRASCF_RAMD_OUT=$CONTRASCF_ROOT/analysis/ramd_pilot/outputs
# Path to the gromacs-ramd build (HITS-MCM patched fork). Empty on the lab
# workstation since the pilot runs on CARC.
export CONTRASCF_GMX_RAMD=
# AmberTools env (antechamber, parmchk2, acpype) for ligand parameterization.
export CONTRASCF_AMBERTOOLS_ENV=

echo "[contrasCF] env: lab workstation. Use \$CONTRASCF_PY for analysis scripts."
