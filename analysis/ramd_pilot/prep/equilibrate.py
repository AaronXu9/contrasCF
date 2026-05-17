"""Run the 5-stage equilibration ladder for one pilot system.

Stages (driven by mdp files in `prep/mdp/`):
  step6.0_minimization
  step6.1_nvt_heat        (annealed 0 → 300 K, BB=1000 SC=500)
  step6.2_npt_release1    (BB=500 SC=200)
  step6.3_npt_release2    (BB=100 SC=50)
  step6.4_npt_release3    (BB=10 SC=10)
  step6.5_npt_free        (5 ns unrestrained)

Each stage produces .gro/.cpt that feeds the next. Final acceptance gate
is left for a follow-up RMSD analyzer (gmx rms over last 2 ns of step6.5).
Run on the same node where production will run.
"""
from __future__ import annotations
import argparse
import shutil
import subprocess
import sys
from pathlib import Path


GMX = "gmx_mpi"
STAGES = (
    "step6.0_minimization",
    "step6.1_nvt_heat",
    "step6.2_npt_release1",
    "step6.3_npt_release2",
    "step6.4_npt_release3",
    "step6.5_npt_free",
)


def _run(cmd: list[str], cwd: Path) -> None:
    print(f"  $ {' '.join(cmd)}")
    subprocess.run(cmd, cwd=cwd, check=True)


def equilibrate(
    system_dir: Path,
    mdp_dir: Path,
    gmx: str = GMX,
) -> Path:
    """Run all 6 stages in `system_dir/equil/`. Returns final .gro path."""
    system_dir = Path(system_dir).resolve()
    mdp_dir = Path(mdp_dir).resolve()
    equil = system_dir / "equil"
    equil.mkdir(parents=True, exist_ok=True)

    # Stage 0 (minimization) takes solv_ions.gro from system_dir; later
    # stages take the previous stage's .gro.
    prev_gro = system_dir / "solv_ions.gro"
    prev_cpt = None
    for i, stage in enumerate(STAGES):
        mdp = mdp_dir / f"{stage}.mdp"
        if not mdp.exists():
            raise FileNotFoundError(mdp)
        shutil.copy(mdp, equil / f"{stage}.mdp")

        ref_gro = system_dir / "solv_ions.gro"  # always the original for -r
        grompp_cmd = [
            gmx, "grompp", "-f", f"{stage}.mdp", "-o", f"{stage}.tpr",
            "-c", str(prev_gro), "-r", str(ref_gro),
            "-p", str(system_dir / "topol.top"),
            "-n", str(system_dir / "index.ndx"),
            "-maxwarn", "3",
        ]
        # Pass -t prev.cpt only if it actually exists. Steepest-descent
        # minimization (step6.0) writes no .cpt; step6.1 (heating) has
        # gen-vel=yes / continuation=no so it doesn't need one either.
        # From step6.2 onward the prior NPT stage's .cpt carries thermostat
        # + barostat state.
        if prev_cpt is not None and prev_cpt.exists():
            grompp_cmd.extend(["-t", str(prev_cpt)])

        _run(grompp_cmd, cwd=equil)
        _run([gmx, "mdrun", "-v", "-deffnm", stage, "-dlb", "yes"], cwd=equil)

        prev_gro = equil / f"{stage}.gro"
        prev_cpt = equil / f"{stage}.cpt"
        if not prev_gro.exists():
            raise RuntimeError(f"stage {stage} did not produce {prev_gro}")

    print(f"[equilibrate] DONE -> {prev_gro}")
    return prev_gro


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--system-dir", type=Path, required=True,
                   help="dir containing solv_ions.gro, topol.top, index.ndx")
    p.add_argument("--mdp-dir", type=Path, required=True,
                   help="dir containing step6.0..step6.5 .mdp files")
    p.add_argument("--gmx", default=GMX)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    equilibrate(args.system_dir, args.mdp_dir, gmx=args.gmx)
    return 0


if __name__ == "__main__":
    sys.exit(main())
