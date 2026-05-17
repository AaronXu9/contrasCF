#!/usr/bin/env bash
#SBATCH --job-name=contrasCF_full
#SBATCH --partition=gpu
# AF3 needs Volta+ for tokamax; bare gpu:1 may schedule a P100.
#SBATCH --gres=gpu:v100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=23:00:00
#SBATCH --array=0-9%5
#SBATCH --output=slurm/logs/full_%A_%a.out
#SBATCH --error=slurm/logs/full_%A_%a.err

# Full CASF-2016 (n=285) split into 10 chunks of ~29 systems each.
# Each task processes a contiguous slice via CONTRASCF_START / _END.
#
# Outputs land in the shared outputs/<pdbid>/<variant>/ tree; tasks
# write disjoint subtrees so no concurrency issues. skip_existing in
# both runners means resubmitting a failed array task only redoes the
# missing cells.
#
# %5 caps to 5 simultaneously-running tasks (CARC fair-share); raise if
# the queue is light.
#
# Submit from the repo root on CARC:
#     sbatch slurm/run_full_casf_carc.sh

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-/project2/katritch_223/aoxu/contrasCF}"
source env/carc.sh

# Configure scope + slice from SLURM array index
TOTAL=285
N_CHUNKS=10
CHUNK_SIZE=$(( (TOTAL + N_CHUNKS - 1) / N_CHUNKS ))   # 29
START=$(( SLURM_ARRAY_TASK_ID * CHUNK_SIZE ))
END=$(( START + CHUNK_SIZE ))
[ "$END" -gt "$TOTAL" ] && END=$TOTAL

export CONTRASCF_SCOPE=full
export CONTRASCF_START=$START
export CONTRASCF_END=$END

echo "=== contrasCF FULL CASF-2016 chunk $SLURM_ARRAY_TASK_ID ==="
echo "Host: $(hostname)"
echo "GPU:  $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo "PDB slice: $START..$END (of $TOTAL)"
echo "Time: $(date -Iseconds)"
echo

echo "=== 1/4 build inputs ==="
$CONTRASCF_PY analysis/casf_mutagenesis/scripts/02_build_full_casf.py

echo
echo "=== 2/4 Boltz-2 ==="
$CONTRASCF_PY analysis/casf_mutagenesis/scripts/03_run_boltz2_subset20.py

echo
echo "=== 3/4 AF3 + MSA ==="
$CONTRASCF_PY analysis/casf_mutagenesis/scripts/06_run_af3_msa_subset20.py

# Step 4 (analysis) is intentionally OUT of the array — run once after
# all chunks finish, not per-chunk (analysis is fast + needs the full
# results).
echo
echo "=== chunk $SLURM_ARRAY_TASK_ID done at $(date -Iseconds) ==="
