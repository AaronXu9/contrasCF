#!/bin/bash
# fix_and_run_surfdock.sh
# Apply SurfDock pocket_center shape fix and run SurfDock analysis
# Usage: bash analysis/scripts/fix_and_run_surfdock.sh [--scope all|casf|ligand]

set -e

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

echo "[fix_and_run_surfdock] Applying SurfDock pocket_center shape fix..."
conda run -n rdkit_env python3 analysis/scripts/apply_surfdock_fix.py

if [ $? -eq 0 ]; then
    echo "[fix_and_run_surfdock] Fix applied successfully. Running SurfDock..."
    SCOPE="${1:---scope=all}"
    if [[ "$SCOPE" == --* ]]; then
        conda run -n rdkit_env python3 analysis/scripts/13_run_surfdock.py "$SCOPE"
    else
        conda run -n rdkit_env python3 analysis/scripts/13_run_surfdock.py
    fi
else
    echo "[fix_and_run_surfdock] Failed to apply SurfDock fix." >&2
    exit 1
fi
