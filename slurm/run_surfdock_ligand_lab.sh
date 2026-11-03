#!/usr/bin/env bash
# Ligand-arm SurfDock re-run (data_prep_todo item 15), LAB BOX ONLY (CARC has no PyMesh).
# Input ligand is re-embedded in the crystal frame (CONTRASCF_SURFDOCK_LIGAND=crystal_frame)
# so SurfDock's surface crop sees the crystal footprint, as on the protein arm.
# Usage: slurm/run_surfdock_ligand_lab.sh wt        # phase 1: WT gate
#        slurm/run_surfdock_ligand_lab.sh variants  # phase 2: everything else
set -uo pipefail
cd /mnt/katritch_lab2/aoxu/contrasCF
source env/lab.sh >/dev/null 2>&1
export TMPDIR=/home/aoxu/tmp/msms
export CONTRASCF_OUTPUTS_ROOT=$CONTRASCF_ROOT/analysis/ligand_mutagenesis/outputs
export CONTRASCF_SURFDOCK_LIGAND=crystal_frame
if [ "${1:-}" = wt ]; then export CONTRASCF_VARIANT_FILTER=wt; else unset CONTRASCF_VARIANT_FILTER; fi
echo "=== SurfDock ligand phase=${1:-all} started $(date -Is) ==="
$CONTRASCF_PY analysis/casf_mutagenesis/scripts/run/08_run_surfdock_variants.py
echo "=== finished $(date -Is) exit=$? ==="
