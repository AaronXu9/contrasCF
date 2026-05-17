#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --account=katritch_223
#SBATCH --job-name="{JOB_NAME}"
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:{GPU_GRES_SPEC}
#SBATCH --time={WALLTIME}
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user={EMAIL}
# To pin a GPU type, override at submit:
#   sbatch --gres=gpu:v100:1 ...   # plentiful, ~50 ns/day
#   sbatch --gres=gpu:a100:1 ...   # faster, more contested
# For full 100 ns production, override walltime: --time=30:00:00
#
# Per-replica τ-RAMD production. One job runs equilibration→production for
# a single (system, replica). Mirrors the labmate VisitorPackage A100 setup
# but runs the patched gromacs-ramd build and stops on ligand exit.
#
# Substituted at submit time by ramd/replica_runner.py:
#   {JOB_NAME}, {EMAIL}, {SYSTEM_DIR}, {REPLICA_IDX}, {RAMD_FORCE},
#   {RAMD_SEED}, {RAMD_RECEPTOR}, {RAMD_LIGAND}, {GMX_RAMD_ROOT}.
#
# Cannot use `set -u`: GMXRC references $GMXLDLIB and $shell before
# defining them. Use only -e and -o pipefail.
set -eo pipefail

module purge
module load gcc/13.3.0 cuda/12.6.3
# gromacs-ramd lives in user space; it provides its own gmx_mpi.
# CUDA must be loaded BEFORE GMXRC for libcufft.so.11 to resolve.
source {GMX_RAMD_ROOT}/bin/GMXRC

# Reset GPU compute mode to allow shared (multi-process) use, matching the
# labmate VisitorPackage convention. Harmless if NVIDIA persistence mode
# already set this; will WARN if user lacks privilege but does not abort.
nvidia-smi -c DEFAULT 2>/dev/null || true

# OMP env (from labmate VisitorPackage; tuned for A100 sockets).
export OMP_DYNAMIC=false
export OMP_NESTED=false
export OMP_SCHEDULE='STATIC'
export OMP_PROC_BIND=true
export OMP_PLACES=sockets
export OMP_WAIT_POLICY=ACTIVE
export OMP_THREAD_LIMIT=4294967295
export OMP_MAX_ACTIVE_LEVELS=2147483647
export GMX_MAXBACKUP=1000

date
echo "[ramd] job $SLURM_JOBID  system={SYSTEM_DIR}  replica={REPLICA_IDX}  force={RAMD_FORCE} kJ/mol/nm"

cd {SYSTEM_DIR}/replicas/r{REPLICA_IDX}

# prod.mdp was pre-rendered by ramd/replica_runner.py at sbatch-render
# time with all RAMD placeholders already substituted. Just grompp it.
gmx_mpi grompp \
    -f prod.mdp \
    -c {SYSTEM_DIR}/equil/step6.5_npt_free.gro \
    -p {SYSTEM_DIR}/topol.top \
    -n {SYSTEM_DIR}/index.ndx \
    -o prod.tpr -maxwarn 2

# ---- run ----
# gromacs-ramd terminates mdrun via SIGTERM when the ligand crosses
# max-dist. mpi/gmx_mpi then exits non-zero, which `set -e` would treat as
# a script failure. Wrap mdrun with explicit exit-code handling: if
# prod.log shows a clean RAMD exit, treat as success regardless of mdrun rc.
set +e
gmx_mpi mdrun -v -deffnm prod -dlb yes
MDRUN_RC=$?
set -e
if grep -q "RAMD group .* has exited" prod.log 2>/dev/null; then
    echo "[ramd] clean RAMD exit detected (mdrun rc=$MDRUN_RC, treating as success)"
    MDRUN_RC=0
fi

# Per-replica replica_result.json is written post-hoc by
# analysis/ramd_pilot/scripts/03_collect_results.py (which calls
# exit_detection.parse_log on prod.log).

date
echo "[ramd] DONE replica={REPLICA_IDX}  mdrun_rc=$MDRUN_RC"
exit $MDRUN_RC
