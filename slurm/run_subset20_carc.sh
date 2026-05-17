#!/usr/bin/env bash
#SBATCH --job-name=contrasCF_subset20
#SBATCH --partition=gpu
# AF3 v3 requires Volta+ (compute 7.0+) for tokamax attention kernels.
# CARC's P100 is Pascal (compute 6.0) and raises
# `NotImplementedError: Not supported on Tesla P100-PCIE-16GB.`
# v100 is the most plentiful Volta+ on CARC (29 nodes); a100/a40/l40s
# also work. Don't use bare `gpu:1` here — SLURM may schedule a P100.
#SBATCH --gres=gpu:v100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=slurm/logs/subset20_%j.out
#SBATCH --error=slurm/logs/subset20_%j.err

# Submit from the repo root on CARC:
#     cd /project2/katritch_223/aoxu/contrasCF
#     mkdir -p slurm/logs
#     sbatch slurm/run_subset20_carc.sh
#
# Runs the full CASF-mutagenesis subset20 pipeline:
#   1. Build inputs (cheap, ~10 s).
#   2. Boltz-2 on all 80 cells (~47 min on RTX 4090; CARC nodes vary).
#   3. AF3 + ColabFold MSA on all 76 runnable cells (~1.8 h).
#   4. Analysis: ligand RMSD + memorisation rates.
#
# Total wallclock: ~3 h on a single GPU. Adjust --time if your CARC GPU
# is slower than an RTX 4090.

# Note: -u dropped because env/carc.sh now sources GMXRC (for τ-RAMD),
# which uses unset $shell / $GMXLDLIB internally. Keep -e + pipefail.
set -eo pipefail

# SLURM copies the script to /var/spool/.../slurm_script, so $0 isn't the
# repo path. Use SLURM_SUBMIT_DIR (sbatch's CWD) when available, else fall
# back to the canonical CARC repo path.
cd "${SLURM_SUBMIT_DIR:-/project2/katritch_223/aoxu/contrasCF}"
source env/carc.sh

echo "=== contrasCF subset20 run on CARC ==="
echo "Host: $(hostname)"
echo "GPU:  $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo "Time: $(date -Iseconds)"
echo

echo "=== 1/4 build inputs ==="
$CONTRASCF_PY analysis/casf_mutagenesis/scripts/01_build_subset20.py
echo

echo "=== 2/4 Boltz-2 ==="
$CONTRASCF_PY analysis/casf_mutagenesis/scripts/03_run_boltz2_subset20.py
echo

echo "=== 3/4 AF3 + MSA ==="
$CONTRASCF_PY analysis/casf_mutagenesis/scripts/06_run_af3_msa_subset20.py
echo

echo "=== 4/4 Analysis ==="
$CONTRASCF_PY analysis/casf_mutagenesis/scripts/05_analyze_subset20.py
echo

echo "=== Done at $(date -Iseconds) ==="
