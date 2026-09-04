#!/usr/bin/env bash
# SurfDock over the ligand-mutagenesis arm on CARC.
#
# SurfDock IS runnable on CARC. An earlier note in this repo and in the
# dockstrat skill said it was not; that check only looked in
# /project2/.../conda/envs, /project2/.../envs and ~/miniconda3/envs, and
# missed the env under /home1. Everything needed is present:
#   env     : /home1/aoxu/.conda/envs/SurfDock_CARC (torch 2.2.2/cu121)
#   source  : /project2/katritch_223/aoxu/projects/SurfDock
#   weights : that tree's model_weights/{docking,posepredict} (140 MB)
#   MSMS    : bundled, comp_surface/tools/transfer/APBS-3.4.1.Linux/bin/msms
#   arrays  : /project2/katritch_223/aoxu/projects/precomputed/precomputed_arrays
#
# CRITICAL: CARC's dockStrat checkout shipped the PRE-FIX surface helper
# (`faces_to_keep = np.arange(len(faces2))`, keep every face). Run with that
# and every pose is silently degraded — it is exactly the bug behind the
# retracted "SurfDock fails on CASF" result. The fixed
# dockstrat/models/_surfdock_surface_helper.py was copied over on 2026-09-04
# (original kept as *.prefix_backup). Verify before trusting a sweep:
#   grep -n "faces_to_keep" $CONTRASCF_DOCKSTRAT_ROOT/dockstrat/models/_surfdock_surface_helper.py
# and healthy meshes are 60-260 vertices, not ~1500.
#
# Smoke test one system:
#   sbatch --array=0 --export=ALL,CONTRASCF_PDBID_FILTER=1bcu slurm/run_surfdock_ligand_carc.sh
# Back half of the arm (lab takes the front half):
#   sbatch --array=0-4 slurm/run_surfdock_ligand_carc.sh
#
#SBATCH --job-name=sd_lig
#SBATCH --account=katritch_223
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=slurm/logs/%x_%A_%a.out
#SBATCH --error=slurm/logs/%x_%A_%a.err

set -o pipefail

WORKTREE=/project2/katritch_223/aoxu/contrasCF-af3
MAIN=/project2/katritch_223/aoxu/contrasCF

cd "$WORKTREE"
# -u must be off across this source: env/carc.sh pulls in GROMACS's GMXRC,
# which reads `shell` and `GMXLDLIB` while unset and would kill the job in 1 s.
set +u
source env/carc.sh
set -u

# Code from the worktree, data from the main clone (outputs/ is gitignored).
export CONTRASCF_ROOT="$WORKTREE"
export CONTRASCF_OUTPUTS_ROOT="$MAIN/analysis/ligand_mutagenesis/outputs"

# MSMS leaks ~1 MB of scratch per cell and never cleans up. /tmp on a compute
# node is small; use the per-user NVMe scratch.
export TMPDIR=/scratch1/aoxu/tmp/msms
mkdir -p "$TMPDIR"

# Split the arm: lab covers systems [0:126], CARC covers [126:251].
BASE=${CONTRASCF_SPLIT_BASE:-126}
CHUNK=${CONTRASCF_CHUNK:-25}
export CONTRASCF_SYSTEM_START=$(( BASE + SLURM_ARRAY_TASK_ID * CHUNK ))
export CONTRASCF_SYSTEM_LIMIT=$CHUNK

echo "=== sd task ${SLURM_ARRAY_TASK_ID} systems [${CONTRASCF_SYSTEM_START}:+${CHUNK}] on $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
echo "surfdock env : $CONTRASCF_SURFDOCK_ENV"
echo "weights      : $CONTRASCF_SURFDOCK_WEIGHTS"
echo "arrays       : $CONTRASCF_SURFDOCK_PRECOMPUTED ($(ls -a "$CONTRASCF_SURFDOCK_PRECOMPUTED" 2>/dev/null | grep -c npy) npy)"
echo -n "crop fix     : "
grep -q "iface_v" "$CONTRASCF_DOCKSTRAT_ROOT/dockstrat/models/_surfdock_surface_helper.py" \
  && echo "PRESENT" || { echo "MISSING - refusing to run"; exit 2; }

$CONTRASCF_PY analysis/casf_mutagenesis/scripts/14_run_surfdock_variants.py
echo "=== sd task ${SLURM_ARRAY_TASK_ID} exit=$? at $(date -Is) ==="
