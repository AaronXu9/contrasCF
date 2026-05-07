#!/usr/bin/env bash
# Recreate raw/ and crystal_ligands/ symlinks on CARC by re-targeting the
# lab pdbbind_cleansplit symlinks at the CARC HiQBind source.
#
# Usage (from the lab workstation):
#     bash env/setup_casf_data_carc.sh
#
# What it does:
#   1. Enumerates lab symlinks under
#      /home/aoxu/projects/VLS-Benchmark-Dataset/data/pdbbind_cleansplit/
#      {raw,crystal_ligands}/ (for the 20-PDB subset by default; pass
#      `full` as $1 to do all 285 CASF-2016 entries).
#   2. Reads each symlink's target, substitutes the lab prefix
#      `/mnt/katritch_lab2/aoxu/data/hiqbind/` → CARC's
#      `/project2/katritch_223/aoxu/data/hiqbind/`.
#   3. Streams `ln -sfn <new_target> <new_link>` commands over SSH into a
#      bash on discovery.usc.edu, populating
#      /project2/katritch_223/aoxu/projects/VLS-Benchmark-Dataset/data/
#      pdbbind_cleansplit/{raw,crystal_ligands}/.

set -euo pipefail

LAB_PDBB=/home/aoxu/projects/VLS-Benchmark-Dataset/data/pdbbind_cleansplit
LAB_HIQB_PREFIX=/mnt/katritch_lab2/aoxu/data/hiqbind
CARC_HIQB_PREFIX=/project2/katritch_223/aoxu/data/hiqbind
CARC_PDBB=/project2/katritch_223/aoxu/projects/VLS-Benchmark-Dataset/data/pdbbind_cleansplit

# subset selection
SCOPE=${1:-subset20}
case "$SCOPE" in
    subset20)
        IDS=$(python3 -c "import json; print('\n'.join(json.load(open('${LAB_PDBB}/labels/PDBbind_casf2016_subset20.json'))['casf2016']))")
        ;;
    full)
        IDS=$(python3 -c "import json; print('\n'.join(json.load(open('${LAB_PDBB}/labels/PDBbind_data_split_cleansplit.json'))['casf2016']))")
        ;;
    *) echo "Usage: $0 [subset20|full]"; exit 1 ;;
esac
N_IDS=$(echo "$IDS" | wc -l)
echo "Setting up $N_IDS PDBs on CARC ($SCOPE)"

# Build the remote shell commands
{
    echo "set -e"
    echo "mkdir -p ${CARC_PDBB}/raw ${CARC_PDBB}/crystal_ligands"
    while IFS= read -r pdbid; do
        # raw/<pdbid>/{*_protein.pdb, *_ligand.sdf}
        for f in "${LAB_PDBB}/raw/${pdbid}/${pdbid}_protein.pdb" \
                 "${LAB_PDBB}/raw/${pdbid}/${pdbid}_ligand.sdf"; do
            [ -L "$f" ] || { echo "echo 'WARN: not a symlink: $f'" >&2; continue; }
            tgt=$(readlink "$f")
            new_tgt=${tgt//$LAB_HIQB_PREFIX/$CARC_HIQB_PREFIX}
            base=$(basename "$f")
            echo "mkdir -p ${CARC_PDBB}/raw/${pdbid}"
            echo "ln -sfn '${new_tgt}' '${CARC_PDBB}/raw/${pdbid}/${base}'"
        done
        # crystal_ligands/<pdbid>_ligand.sdf
        f="${LAB_PDBB}/crystal_ligands/${pdbid}_ligand.sdf"
        if [ -L "$f" ]; then
            tgt=$(readlink "$f")
            new_tgt=${tgt//$LAB_HIQB_PREFIX/$CARC_HIQB_PREFIX}
            echo "ln -sfn '${new_tgt}' '${CARC_PDBB}/crystal_ligands/${pdbid}_ligand.sdf'"
        fi
    done <<< "$IDS"
    echo "echo done"
} | ssh aoxu@discovery.usc.edu bash
