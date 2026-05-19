#!/usr/bin/env bash
#SBATCH --job-name=contrasCF_boltz_aff
#SBATCH --partition=gpu
# Boltz-2 affinity head is small; V100 (or any Volta+) is plenty. Don't ask
# for bare gpu:1 — see notes in run_subset20_carc.sh.
#SBATCH --gres=gpu:v100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:30:00
#SBATCH --output=slurm/logs/boltz_aff_%j.out
#SBATCH --error=slurm/logs/boltz_aff_%j.err

# Re-runs Boltz-2 on the CASF-mutagenesis subset20 with the binding-affinity
# head enabled. The prior Boltz YAMLs didn't request
#   properties:
#     - affinity:
#         binder: <lig>
# so no affinity_<prefix>.json files were ever produced. After commit
# 4b697ef the renderer always adds that block, but existing on-disk
# predictions don't have the sidecar — they need to be re-run.
#
# This script:
#   1. Moves the existing Boltz-2 CIFs + confidence JSONs into
#      outputs/_backup_pre_affinity/ (NON-destructive; structures will
#      be deterministic with seed=42 in the runner, so the new CIFs
#      should match these — keep them around to verify).
#   2. Re-renders boltz.yaml for every subset20 (pdbid, variant) with
#      the affinity property included.
#   3. Re-runs Boltz-2 — produces new CIFs, confidence JSONs, AND
#      the affinity_<prefix>.json sidecar.
#   4. Re-runs the analyze driver (subset20 scope) — paired_affinity_
#      subset20.csv is the new output of interest.
#
# AF3 / AF3+MSA are intentionally NOT re-run: they don't predict affinity,
# and their CIFs are unchanged by this work.
#
# Submit from the repo root on CARC:
#     cd /project2/katritch_223/aoxu/contrasCF
#     mkdir -p slurm/logs
#     sbatch slurm/run_boltz2_affinity_subset20_carc.sh

set -eo pipefail
cd "${SLURM_SUBMIT_DIR:-/project2/katritch_223/aoxu/contrasCF}"
source env/carc.sh

echo "=== contrasCF Boltz-2 affinity re-run on CARC ==="
echo "Host: $(hostname)"
echo "GPU:  $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo "Time: $(date -Iseconds)"
echo

# Step 1: back up existing Boltz-2 outputs (idempotent — no-op on a second
# submission once the backup is in place).
echo "=== 1/4 backup existing Boltz-2 outputs ==="
$CONTRASCF_PY - <<'PY'
import json, shutil
from pathlib import Path
import os
ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/project2/katritch_223/aoxu/contrasCF"))
out = ROOT / "analysis/casf_mutagenesis/outputs"
backup = out / "_backup_pre_affinity"
ids = json.loads((ROOT / "data/casf2016/labels/PDBbind_casf2016_subset20.json").read_text())["casf2016"]
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
print(f"Moved to _backup_pre_affinity: {n_cif} CIFs, {n_conf} confidence JSONs")
PY
echo

echo "=== 2/2 Boltz-2 (with affinity) ==="
# Note: input regen (01_build_subset20.py) intentionally NOT run here —
# CARC doesn't have the 42 GB CASF data tree (raw PDBs + ligand SDFs).
# The affinity-enabled boltz.yaml files are rsync'd from the lab box
# before submission, and Boltz-2 only needs the YAMLs as input.
#
# Analysis (05_analyze_subset20.py) also intentionally NOT run here —
# it needs the crystal data for RMSD. Run analysis back on the lab box
# after rsyncing the affinity JSONs back.
$CONTRASCF_PY analysis/casf_mutagenesis/scripts/03_run_boltz2_subset20.py
echo

echo "=== Done at $(date -Iseconds) ==="
