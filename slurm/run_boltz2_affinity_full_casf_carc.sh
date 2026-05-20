#!/usr/bin/env bash
#SBATCH --job-name=contrasCF_boltz_aff_full
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=07:00:00
#SBATCH --array=0-9%5
#SBATCH --output=slurm/logs/boltz_aff_full_%A_%a.out
#SBATCH --error=slurm/logs/boltz_aff_full_%A_%a.err

# Full CASF-2016 (n=285; 251 with usable inputs) × 4 variants = ~1004 cells.
# Split into 10 chunks of ~29 systems each (~116 cells per chunk). V100
# affinity-enabled inference runs ~150s/cell, so each chunk takes ~5h.
# Time limit 7h gives headroom.
#
# Pre-flight (run from the lab box BEFORE submitting):
#   1. Regenerate YAMLs locally so they have the `properties: - affinity:`
#      block:  $CONTRASCF_PY analysis/casf_mutagenesis/scripts/02_build_full_casf.py
#   2. Rsync labels + the 1004 YAMLs to CARC:
#      rsync -avz data/casf2016/labels/                      \
#          discovery.usc.edu:/project2/katritch_223/aoxu/contrasCF/data/casf2016/labels/
#      rsync -avz --include='*/' --include='boltz.yaml' --exclude='*'  \
#          analysis/casf_mutagenesis/outputs/                \
#          discovery.usc.edu:/project2/katritch_223/aoxu/contrasCF/analysis/casf_mutagenesis/outputs/
#
# Then on CARC:
#   cd /project2/katritch_223/aoxu/contrasCF
#   mkdir -p slurm/logs
#   sbatch slurm/run_boltz2_affinity_full_casf_carc.sh
#
# Per-chunk steps:
#   - Backup existing Boltz-2 CIFs+confidence JSONs to _backup_pre_affinity/
#     for cells in this chunk's pdbid slice (idempotent).
#   - Run Boltz-2 on this chunk; skip_existing on model_0.cif still works
#     because the backup moved them out of the way.
#
# Build (01) is NOT run on CARC — the data tree is gitignored and 42 GB,
# we ship YAMLs pre-rendered via rsync.
# Analyze (05) is NOT run on CARC — it needs the crystal data; run on the
# lab box after rsyncing affinity_*.json and model_*.cif back.

set -eo pipefail
cd "${SLURM_SUBMIT_DIR:-/project2/katritch_223/aoxu/contrasCF}"
source env/carc.sh

TOTAL=285
N_CHUNKS=10
CHUNK_SIZE=$(( (TOTAL + N_CHUNKS - 1) / N_CHUNKS ))   # 29
START=$(( SLURM_ARRAY_TASK_ID * CHUNK_SIZE ))
END=$(( START + CHUNK_SIZE ))
[ "$END" -gt "$TOTAL" ] && END=$TOTAL
export CONTRASCF_SCOPE=full
export CONTRASCF_START=$START
export CONTRASCF_END=$END

echo "=== contrasCF Boltz-2 affinity full-CASF chunk $SLURM_ARRAY_TASK_ID ==="
echo "Host: $(hostname)"
echo "GPU:  $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo "PDB slice: $START..$END (of $TOTAL)"
echo "Time: $(date -Iseconds)"
echo

# Step 1: back up existing Boltz-2 outputs for THIS chunk's slice.
# Idempotent — already-moved files are simply absent on the next run.
echo "=== 1/2 backup existing Boltz-2 outputs for chunk slice ==="
$CONTRASCF_PY - <<'PY'
import json, shutil, os
from pathlib import Path
ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/project2/katritch_223/aoxu/contrasCF"))
out = ROOT / "analysis/casf_mutagenesis/outputs"
backup = out / "_backup_pre_affinity"
split = ROOT / "data/casf2016/labels/PDBbind_data_split_cleansplit.json"
all_ids = json.loads(split.read_text())["casf2016"]
start = int(os.environ.get("CONTRASCF_START", "0"))
end = int(os.environ.get("CONTRASCF_END", str(len(all_ids))))
ids = all_ids[start:end]
n_cif = n_conf = 0
for pid in ids:
    for v in ("wt", "rem", "pack", "inv"):
        sys_dir = out / pid / v
        if not sys_dir.is_dir():
            continue
        dst = backup / pid / v
        for f in sys_dir.glob(f"{pid}_{v}_model_*.cif"):
            dst.mkdir(parents=True, exist_ok=True)
            shutil.move(str(f), str(dst / f.name))
            n_cif += 1
        for f in sys_dir.glob(f"confidence_{pid}_{v}_model_*.json"):
            dst.mkdir(parents=True, exist_ok=True)
            shutil.move(str(f), str(dst / f.name))
            n_conf += 1
print(f"Chunk slice [{start}:{end}] ({len(ids)} systems): moved {n_cif} CIFs, {n_conf} confidence JSONs")
PY
echo

echo "=== 2/2 Boltz-2 (with affinity) ==="
$CONTRASCF_PY analysis/casf_mutagenesis/scripts/03_run_boltz2_subset20.py
echo

echo "=== Done chunk $SLURM_ARRAY_TASK_ID at $(date -Iseconds) ==="
