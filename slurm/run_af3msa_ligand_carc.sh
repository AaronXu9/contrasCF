#!/usr/bin/env bash
# AF3+MSA over the ligand-mutagenesis arm (251 systems / 1300 cells) on CARC.
#
# Sizing: AF3+MSA measures a 167 s median per cell (n=76, lab 4090), p90 318 s,
# so the whole arm is ~60 GPU-hours. That does not fit one job, hence a job
# array of 10-system chunks (~50 cells, ~2.3 h median each). Short tasks also
# backfill far better than long ones on a partition sitting at ~270 pending.
#
# GPU type is pinned deliberately. AF3 v3 needs compute 7.0+; CARC's P100 is
# Pascal (6.0) and fails, so NEVER use a bare `--gres=gpu:1` here — SLURM may
# hand you a P100. a40 (48 GB) has the most non-drained nodes; l40s and a100
# are equivalent fallbacks.
#
# Submit from the worktree:
#   cd /project2/katritch_223/aoxu/contrasCF-af3 && mkdir -p slurm/logs
#   sbatch --array=0        slurm/run_af3msa_ligand_carc.sh   # canary first
#   sbatch --array=1-25%8   slurm/run_af3msa_ligand_carc.sh   # then the rest
#
#SBATCH --job-name=af3msa_lig
#SBATCH --account=katritch_223
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --array=0-25%8
#SBATCH --output=slurm/logs/%x_%A_%a.out
#SBATCH --error=slurm/logs/%x_%A_%a.err

set -o pipefail

WORKTREE=/project2/katritch_223/aoxu/contrasCF-af3
MAIN=/project2/katritch_223/aoxu/contrasCF

cd "$WORKTREE"
# `set -u` must be OFF across this source. env/carc.sh sources GROMACS's
# GMXRC for the tau-RAMD pilot, and GMXRC is not -u clean: it reads `shell`
# and `GMXLDLIB` while unset, so the job dies in 1 s before reaching Python.
set +u
source env/carc.sh
set -u

# Code from the worktree, data from the main clone. env/carc.sh points
# CONTRASCF_ROOT at the main clone, but that clone sits on the CounterFold
# branch and does NOT carry the AF3 multi-chain fix (0169b42) or the
# module-agnostic runner — so override it, or every import silently comes
# from the wrong checkout. outputs/ is gitignored and exists only in MAIN.
export CONTRASCF_ROOT="$WORKTREE"
export CONTRASCF_OUTPUTS_ROOT="$MAIN/analysis/ligand_mutagenesis/outputs"
export CONTRASCF_MSA_CACHE="$MAIN/analysis/casf_mutagenesis/outputs/_msa_cache"
export CONTRASCF_SCOPE=disk

module load cuda/12.6.3 2>/dev/null || true

CHUNK=${CONTRASCF_CHUNK:-10}
export CONTRASCF_START=$(( SLURM_ARRAY_TASK_ID * CHUNK ))
export CONTRASCF_END=$(( CONTRASCF_START + CHUNK ))
export CONTRASCF_RUN_LOG="$CONTRASCF_OUTPUTS_ROOT/_af3_logs/af3_msa_run_log_${SLURM_ARRAY_TASK_ID}.json"

echo "=== task ${SLURM_ARRAY_TASK_ID} systems [${CONTRASCF_START}:${CONTRASCF_END}] on $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
echo "code:  $CONTRASCF_ROOT"
echo "data:  $CONTRASCF_OUTPUTS_ROOT"
echo "cache: $CONTRASCF_MSA_CACHE ($(ls "$CONTRASCF_MSA_CACHE" 2>/dev/null | wc -l) a3m)"

$CONTRASCF_PY analysis/casf_mutagenesis/scripts/06_run_af3_msa_subset20.py
echo "=== task ${SLURM_ARRAY_TASK_ID} exit=$? at $(date -Is) ==="
