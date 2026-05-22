#!/usr/bin/env bash
#SBATCH --job-name=contrasCF_boltz_lig
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=07:00:00
#SBATCH --array=0-9%5
#SBATCH --output=slurm/logs/boltz_lig_%A_%a.out
#SBATCH --error=slurm/logs/boltz_lig_%A_%a.err

# Boltz-2 (with affinity head) on the ligand_mutagenesis module — sibling
# of run_boltz2_affinity_full_casf_carc.sh but for the ligand-side variants
# (halogenation, charge-swap, methylation; no rem/pack/inv here).
#
# 255 systems × ~5 variants/system = ~1300 cells. Split into 10 chunks of
# ~26 systems (~130 cells) each. V100 @ ~150s/cell ⇒ ~5.5h per chunk; 7h
# time limit gives headroom.
#
# Pre-flight (lab box):
#   1. Regenerate ligand_mutagenesis YAMLs so they have the
#      properties.affinity block:
#        $CONTRASCF_PY analysis/ligand_mutagenesis/scripts/02_build_full_casf.py
#   2. Rsync the YAMLs to CARC (~50 KB):
#        rsync -avz --include='*/' --include='boltz.yaml' --exclude='*' \
#          analysis/ligand_mutagenesis/outputs/                                \
#          discovery.usc.edu:/project2/katritch_223/aoxu/contrasCF/analysis/ligand_mutagenesis/outputs/
#
# On CARC:
#   cd /project2/katritch_223/aoxu/contrasCF
#   mkdir -p slurm/logs
#   sbatch slurm/run_boltz2_affinity_ligand_carc.sh
#
# Build (02) NOT run on CARC — same rationale as casf_mutagenesis: the 42 GB
# CASF data tree is gitignored, ship YAMLs pre-rendered.
# Analysis NOT run on CARC — crystal structures + ligands needed; run on lab
# box after rsync of model_*.cif + affinity_*.json back.

set -eo pipefail
cd "${SLURM_SUBMIT_DIR:-/project2/katritch_223/aoxu/contrasCF}"
source env/carc.sh

# Slice the SORTED list of ligand_mutagenesis system dirs by array index.
# The 03_run_boltz2.py runner honors CONTRASCF_START / CONTRASCF_END to
# select a pdbid slice.
N_CHUNKS=10
LIG_OUT="$PWD/analysis/ligand_mutagenesis/outputs"
TOTAL=$($CONTRASCF_PY - <<PY
from pathlib import Path
p = Path("$LIG_OUT")
n = sum(1 for d in p.iterdir() if d.is_dir() and not d.name.startswith("_"))
print(n)
PY
)
CHUNK_SIZE=$(( (TOTAL + N_CHUNKS - 1) / N_CHUNKS ))
START=$(( SLURM_ARRAY_TASK_ID * CHUNK_SIZE ))
END=$(( START + CHUNK_SIZE ))
[ "$END" -gt "$TOTAL" ] && END=$TOTAL
export CONTRASCF_START=$START
export CONTRASCF_END=$END

echo "=== contrasCF ligand_mutagenesis Boltz-2 chunk $SLURM_ARRAY_TASK_ID ==="
echo "Host: $(hostname)"
echo "GPU:  $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo "system slice: $START..$END (of $TOTAL)"
echo "Time: $(date -Iseconds)"
echo

$CONTRASCF_PY analysis/ligand_mutagenesis/scripts/03_run_boltz2.py

echo "=== Done chunk $SLURM_ARRAY_TASK_ID at $(date -Iseconds) ==="
