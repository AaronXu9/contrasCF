#!/bin/bash
# Build the HITS-MCM gromacs-ramd patched fork on USC CARC discovery.
#
# Why a custom build: CARC's gromacs-gpu/2024.3 module is vanilla GROMACS;
# RAMD requires a patched mdrun that applies the random force on the
# ligand COM. HITS-MCM ships gromacs-ramd as a fork (https://github.com/HITS-MCM/gromacs-ramd).
#
# This script must run from a CARC interactive shell (or as a sbatch
# build job). It will NOT submit anything; it only compiles into
# $CONTRASCF_GMX_RAMD/.
#
# Usage:
#     ssh aoxu@discovery.usc.edu
#     bash $CONTRASCF_ROOT/analysis/ramd_pilot/carc_setup/install_gromacs_ramd.sh
#
# Verify after build:
#     source $CONTRASCF_GMX_RAMD/bin/GMXRC
#     gmx_mpi mdrun -h | grep -i ramd
set -euo pipefail

INSTALL_DIR=${CONTRASCF_GMX_RAMD:-/project2/katritch_223/aoxu/gromacs-ramd-2024.3}
SRC_DIR=${INSTALL_DIR}.src
BUILD_DIR=${INSTALL_DIR}.build
GMX_RAMD_REPO=https://github.com/HITS-MCM/gromacs-ramd.git
# IMPORTANT: use the *tag*, not the release-2024 branch. The branch HEAD
# (as of 2026-05) contains RAMD source files but never wires them into
# the CMake graph (no `add_subdirectory(ramd)`), so a build off the branch
# silently produces vanilla GROMACS without RAMD. The tag does wire it.
GMX_RAMD_TAG=gromacs-2024.1-ramd-2.1

echo "[install_gromacs_ramd] target install dir: $INSTALL_DIR"
echo "[install_gromacs_ramd] CARC modules: gcc/13.3.0 cmake/3.29.4 cuda/12.6.3 openmpi/5.0.5"

module purge
module load gcc/13.3.0
module load cmake/3.29.4
module load cuda/12.6.3
module load openmpi/5.0.5
# Verify the toolchain is on PATH; bail loudly otherwise so we don't
# discover a missing tool 30 minutes into the build.
for tool in cmake mpicxx nvcc g++; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "[install_gromacs_ramd] missing tool: $tool. module load chain failed?"
        exit 1
    fi
done

if [[ -d "$INSTALL_DIR" ]]; then
    echo "[install_gromacs_ramd] $INSTALL_DIR already exists; skipping clone+build."
    echo "                       to rebuild: rm -rf $INSTALL_DIR $SRC_DIR $BUILD_DIR"
    exit 0
fi

mkdir -p "$(dirname "$INSTALL_DIR")"
git clone --depth 1 --branch "$GMX_RAMD_TAG" "$GMX_RAMD_REPO" "$SRC_DIR"

mkdir -p "$BUILD_DIR" && cd "$BUILD_DIR"
cmake "$SRC_DIR" \
    -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
    -DGMX_BUILD_OWN_FFTW=ON \
    -DGMX_GPU=CUDA \
    -DGMX_MPI=ON \
    -DGMX_DOUBLE=OFF \
    -DGMX_SIMD=AVX2_256 \
    -DREGRESSIONTEST_DOWNLOAD=OFF \
    -DGMX_DEFAULT_SUFFIX=ON \
    -DGMX_BINARY_SUFFIX=_mpi

make -j8
make install

echo "[install_gromacs_ramd] DONE. To use:"
echo "    source $INSTALL_DIR/bin/GMXRC"
echo "    gmx_mpi mdrun -h | grep -i ramd  # should show RAMD options"
