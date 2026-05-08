#!/usr/bin/env bash
#SBATCH --job-name=contrasCF_smoke
#SBATCH --partition=gpu
# AF3 v3 needs Volta+ for tokamax kernels — pin to V100. See
# slurm/run_subset20_carc.sh for the longer note.
#SBATCH --gres=gpu:v100:1
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

# SLURM copies the script to /var/spool/.../slurm_script, so $0 isn't the
# repo path. Use SLURM_SUBMIT_DIR (sbatch's CWD) when available, else fall
# back to the canonical CARC repo path.
cd "${SLURM_SUBMIT_DIR:-/project2/katritch_223/aoxu/contrasCF}"
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
echo "--- 4. One AF3+MSA prediction (1pxn-wt; ~90s) ---"
# Run the actual AF3 inference on a single system so we surface
# host-specific issues (missing XLA_FLAGS, missing
# --flash_attention_implementation flag, GPU-arch mismatch, etc.) BEFORE
# the multi-hour full sbatch run. This is more expensive but catches
# bugs early — the alternative is wasting 4+ hours per missed flag.
$CONTRASCF_PY -c "
import sys, json, hashlib, shutil, subprocess
sys.path.insert(0,'analysis')
from casf_mutagenesis.msa_via_boltz import fetch_msa_via_boltz
from casf_mutagenesis.inputs_af3 import (
    clean_a3m_for_af3, render_af3, rewrite_a3m_query,
)
import importlib.util
spec = importlib.util.spec_from_file_location(
    'runner06', 'analysis/casf_mutagenesis/scripts/06_run_af3_msa_subset20.py')
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)

# Fetch + clean the MSA, build a fresh JSON, run AF3, check the cif lands.
seq = json.load(open('analysis/casf_mutagenesis/outputs/1pxn/wt/af3.json'))['sequences'][0]['protein']['sequence']
smiles = json.load(open('analysis/casf_mutagenesis/outputs/1pxn/wt/af3.json'))['sequences'][1]['ligand']['smiles']
a3m = fetch_msa_via_boltz(seq, cache_dir='analysis/casf_mutagenesis/outputs/_msa_cache')
a3m = clean_a3m_for_af3(rewrite_a3m_query(a3m, seq), len(seq))

import tempfile
with tempfile.TemporaryDirectory(prefix='af3_smoke_') as tmp:
    from pathlib import Path
    json_path = Path(tmp)/'af3.json'
    render_af3('1pxn_wt_smoke', [('A', seq)], smiles, json_path, chain_msas={'A': a3m})
    sys_dir = mod._run_af3(json_path, Path(tmp)/'work')
    # Find produced cif
    cifs = list(sys_dir.glob('seed-*/*model.cif'))
    print(f'AF3 produced {len(cifs)} cif files; first: {cifs[0] if cifs else \"(none)\"}')
    assert cifs, 'AF3 produced no cif'
print('AF3 smoke OK')
"

rm -rf "$SMOKE_DIR"
echo
echo "=== smoke test PASSED at $(date -Iseconds) ==="
