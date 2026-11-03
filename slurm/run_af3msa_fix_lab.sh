#!/usr/bin/env bash
# Re-fold 1o5b + 2c3i (AF3+MSA) in both arms after the MSA merge fix (lab GPU).
set -uo pipefail
cd /mnt/katritch_lab2/aoxu/contrasCF; source env/lab.sh >/dev/null 2>&1
export CONTRASCF_SCOPE=disk CONTRASCF_PDBID_FILTER=1o5b,2c3i
for ARM in casf_mutagenesis ligand_mutagenesis; do
  export CONTRASCF_OUTPUTS_ROOT=$CONTRASCF_ROOT/analysis/$ARM/outputs
  export CONTRASCF_RUN_LOG=$CONTRASCF_OUTPUTS_ROOT/_af3_logs/af3_msa_run_log_msafix_20261102.json
  echo "=== $ARM $(date -Is)"; $CONTRASCF_PY analysis/casf_mutagenesis/scripts/run/04_run_af3_msa.py
done; echo "=== finished $(date -Is)"
