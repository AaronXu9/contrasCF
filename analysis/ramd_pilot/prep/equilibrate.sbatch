#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --account=katritch_223
#SBATCH --job-name=ramd_equil
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#
# 5-stage equilibration ladder for one pilot system.
# ~3 hr on V100, ~2 hr on A100; 6 hr walltime gives buffer for slow startup.
#
# Submit:
#     sbatch analysis/ramd_pilot/prep/equilibrate.sbatch <SYSTEM_DIR>
#     sbatch analysis/ramd_pilot/prep/equilibrate.sbatch \
#         /scratch1/aoxu/contrasCF/ramd_pilot/bindingsite_wt
#
# To request a specific GPU type, override at submit time:
#     sbatch --gres=gpu:v100:1 .../equilibrate.sbatch <SYSTEM_DIR>
#     sbatch --gres=gpu:a100:1 .../equilibrate.sbatch <SYSTEM_DIR>
#
# NOTE: cannot use `set -u` here — GROMACS's GMXRC references $GMXLDLIB
# and $shell before defining them, which trips the -u check. Use only
# -e and -o pipefail.
set -eo pipefail

SYSTEM_DIR="${1:-${SYSTEM_DIR:-}}"
if [[ -z "$SYSTEM_DIR" ]]; then
    echo "[equilibrate.sbatch] usage: sbatch equilibrate.sbatch <SYSTEM_DIR>"
    echo "                     or:    sbatch --export=SYSTEM_DIR=... equilibrate.sbatch"
    exit 1
fi
SYSTEM_DIR="$(readlink -f "$SYSTEM_DIR")"
[[ -d "$SYSTEM_DIR" ]] || { echo "[equilibrate] not a dir: $SYSTEM_DIR"; exit 1; }

# Need solv_ions.gro + topol.top + index.ndx; bail if missing.
for f in solv_ions.gro topol.top index.ndx; do
    [[ -f "$SYSTEM_DIR/$f" ]] || { echo "[equilibrate] missing $SYSTEM_DIR/$f"; exit 1; }
done

CONTRASCF_ROOT=${CONTRASCF_ROOT:-/project2/katritch_223/aoxu/contrasCF}
source "$CONTRASCF_ROOT/env/carc.sh"

module purge
module load gcc/13.3.0 cuda/12.6.3
source "$CONTRASCF_GMX_RAMD/bin/GMXRC"

# AmberTools env hosts the analysis Python (numpy, scipy, etc).
source /spack/conda/miniforge3/24.3.0/etc/profile.d/conda.sh
conda activate "$CONTRASCF_AMBERTOOLS_ENV"

# OMP env (mirrors slurm_template.sbatch for consistency)
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OMP_DYNAMIC=false
export OMP_NESTED=false
export OMP_SCHEDULE=STATIC
export OMP_PROC_BIND=true
export OMP_PLACES=sockets
export OMP_WAIT_POLICY=ACTIVE
export GMX_MAXBACKUP=1000

date
echo "[equilibrate] job $SLURM_JOBID  on $SLURMD_NODENAME  GPU=$CUDA_VISIBLE_DEVICES"
echo "[equilibrate] system: $SYSTEM_DIR"

MDP_DIR="$CONTRASCF_ROOT/analysis/ramd_pilot/prep/mdp"

# Run the 5-stage ladder. Equilibrate.equilibrate() handles grompp + mdrun
# for each stage and asserts the final .gro lands.
"$CONTRASCF_AMBERTOOLS_ENV/bin/python" -c "
import sys
sys.path.insert(0, '$CONTRASCF_ROOT')
from pathlib import Path
from analysis.ramd_pilot.prep.equilibrate import equilibrate
equilibrate(Path('$SYSTEM_DIR'), Path('$MDP_DIR'))
"

date
echo "[equilibrate] DONE  $SYSTEM_DIR/equil/step6.5_npt_free.gro"
