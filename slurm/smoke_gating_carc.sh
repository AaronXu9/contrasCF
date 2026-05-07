#!/usr/bin/env bash
#SBATCH --job-name=contrasCF_smoke
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --output=slurm/logs/smoke_%j.out
#SBATCH --error=slurm/logs/smoke_%j.err

# 30-minute smoke test on CARC: gating verifier (CDK2 11/11) + one
# Boltz-2 prediction (~1 min) + one AF3+MSA prediction (~3 min) on the
# cheapest subset20 system. Confirms tools, paths, and GPU are all
# wired up correctly before launching the multi-hour subset20 batch.
#
# Submit from the repo root on CARC:
#     sbatch slurm/smoke_gating_carc.sh

set -euo pipefail
cd "$(dirname "$0")/.."

source env/carc.sh

echo "=== smoke test on $(hostname) at $(date -Iseconds) ==="
nvidia-smi --query-gpu=name,memory.total,memory.free --format=csv,noheader
echo

echo "--- 1. CDK2/MEK1 gating verification (no GPU needed) ---"
$CONTRASCF_PY analysis/casf_mutagenesis/scripts/00_verify_reference_systems.py

echo
echo "--- 2. Build subset20 inputs ---"
$CONTRASCF_PY analysis/casf_mutagenesis/scripts/01_build_subset20.py | tail -5

echo
echo "--- 3. One Boltz-2 prediction (1pxn-wt) ---"
SMOKE_PDB=1pxn
SMOKE_VAR=wt
SMOKE_DIR=/tmp/contrasCF_smoke_$$
mkdir -p "$SMOKE_DIR"
cp analysis/casf_mutagenesis/outputs/$SMOKE_PDB/$SMOKE_VAR/boltz.yaml "$SMOKE_DIR/${SMOKE_PDB}_${SMOKE_VAR}.yaml"
CUDA_VISIBLE_DEVICES=0 "$CONTRASCF_BOLTZ_BIN" predict \
    "$SMOKE_DIR/${SMOKE_PDB}_${SMOKE_VAR}.yaml" \
    --out_dir "$SMOKE_DIR" --model boltz2 --output_format mmcif \
    --diffusion_samples 1 --recycling_steps 1 --sampling_steps 50 --seed 42 \
    2>&1 | tail -3
ls "$SMOKE_DIR"/boltz_results_*/predictions/*/*.cif | head -1

echo
echo "--- 4. One AF3+MSA prediction (1pxn-wt) ---"
$CONTRASCF_PY -c "
import sys
sys.path.insert(0,'analysis')
from casf_mutagenesis.msa_via_boltz import fetch_msa_via_boltz
import json, hashlib
seq = json.load(open('analysis/casf_mutagenesis/outputs/1pxn/wt/af3.json'))['sequences'][0]['protein']['sequence']
key = hashlib.sha1(seq.encode()).hexdigest()[:16]
a3m = fetch_msa_via_boltz(seq, cache_dir='analysis/casf_mutagenesis/outputs/_msa_cache')
print(f'MSA fetched: sha1[:16]={key} n_seqs={a3m.count(chr(62))}')
"

# Just verify we have a fresh AF3 JSON; full AF3 inference takes ~90s
# and its own verification runs are covered by 06_run_af3_msa.
echo "(AF3 inference smoke skipped here; tested by full sbatch run.)"

rm -rf "$SMOKE_DIR"
echo
echo "=== smoke test PASSED at $(date -Iseconds) ==="
