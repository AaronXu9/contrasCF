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

# --- SurfDock on CARC -----------------------------------------------------
# SurfDock IS runnable here (an earlier note claiming otherwise was wrong: it
# only ever checked /project2/.../conda/envs, /project2/.../envs and
# ~/miniconda3/envs, and the env actually lives under /home1).
#   env      : /home1/aoxu/.conda/envs/SurfDock_CARC (torch 2.2.2/cu121, pyg 2.6.1)
#   source   : /project2/katritch_223/aoxu/projects/SurfDock
#   weights  : that tree's model_weights/{docking,posepredict,screen} (140 MB)
#   MSMS/APBS: bundled at comp_surface/tools/transfer/APBS-3.4.1.Linux/bin
#   arrays   : 6 HIDDEN .npy dotfiles, rsynced from lab 2026-09-04 (413 MB).
#              Copy them with `rsync -a src/ dst/`; a `cp src/*` glob silently
#              copies NOTHING because every file starts with a dot.
export CONTRASCF_SURFDOCK_ENV=/home1/aoxu/.conda/envs/SurfDock_CARC
export CONTRASCF_SURFDOCK_DIR=/project2/katritch_223/aoxu/projects/SurfDock
export CONTRASCF_SURFDOCK_WEIGHTS=/project2/katritch_223/aoxu/projects/SurfDock/model_weights
export CONTRASCF_SURFDOCK_PRECOMPUTED=/project2/katritch_223/aoxu/projects/precomputed/precomputed_arrays
export CONTRASCF_DOCKSTRAT_ROOT=/project2/katritch_223/aoxu/dockStrat

# Analysis Python: CARC's `boltzina_env` already has rdkit + gemmi +
# biopython + numpy + pandas + yaml AND the Boltz-2 binary, so we use it
# for both analysis and Boltz inference. The dedicated rdkit_env at
# /project2/katritch_223/aoxu/conda/envs/rdkit_env is missing gemmi as
# of 2026-05-06; if it gets gemmi later, switch CONTRASCF_PY back to it.
export _ANALYSIS_ENV=/project2/katritch_223/aoxu/conda/envs/boltzina_env
if [ -d "$_ANALYSIS_ENV" ]; then
    export LD_LIBRARY_PATH="$_ANALYSIS_ENV/lib:${LD_LIBRARY_PATH:-}"
    export CONTRASCF_PY="$_ANALYSIS_ENV/bin/python"
else
    echo "[contrasCF] WARNING: analysis env not found at $_ANALYSIS_ENV"
    echo "             Set CONTRASCF_PY manually."
fi
unset _ANALYSIS_ENV

# --- τ-RAMD pilot ---------------------------------------------------------
# Output root for the ramd_pilot module. Sim working dirs land under
# $CONTRASCF_RAMD_OUT/<system>/replicas/r{NN}; analysis under
# $CONTRASCF_RAMD_OUT/<run>/.
export CONTRASCF_RAMD_OUT=/scratch1/aoxu/contrasCF/ramd_pilot
# Path to the gromacs-ramd build (HITS-MCM patched fork). Built per
# analysis/ramd_pilot/carc_setup/install_gromacs_ramd.sh — that script
# leaves the binaries under $CONTRASCF_GMX_RAMD/{bin,share}.
export CONTRASCF_GMX_RAMD=/project2/katritch_223/aoxu/gromacs-ramd-2024.3
# Auto-source GMXRC so gmx_mpi is on PATH after `source env/carc.sh`.
# Skipped silently if the build hasn't run yet (first-time setup).
# CUDA must be loaded before GMXRC for libcufft.so.11 to resolve at runtime.
if [ -f "$CONTRASCF_GMX_RAMD/bin/GMXRC" ]; then
    if command -v module >/dev/null 2>&1; then
        module load cuda/12.6.3 2>/dev/null || true
    fi
    source "$CONTRASCF_GMX_RAMD/bin/GMXRC"
fi
# AmberTools env (antechamber, parmchk2, acpype). Built per
# analysis/ramd_pilot/carc_setup/setup_ambertools_env.sh.
export CONTRASCF_AMBERTOOLS_ENV=/project2/katritch_223/aoxu/conda/envs/ambertools

echo "[contrasCF] env: USC CARC. CONTRASCF_ROOT=$CONTRASCF_ROOT"
