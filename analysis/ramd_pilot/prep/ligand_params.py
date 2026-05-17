"""GAFF2 + AM1-BCC parameters for ATP via antechamber + acpype.

Output (per ligand):
  <out>/ATP/
    ATP.mol2          input geometry
    ATP.frcmod        parmchk2 missing-parameter patches
    ATP.lib           tleap-emitted residue library
    ATP_GMX.itp       GROMACS topology (acpype output)
    ATP_GMX.gro       GROMACS coordinates
    ATP_posre.itp     position restraints

Requires `$CONTRASCF_AMBERTOOLS_ENV` activated (antechamber, parmchk2, acpype).

Note: for ATP specifically, AMBER-published parameters (Meagher/Ren/Cheatham)
are better-validated than GAFF2 + AM1-BCC. This module is the GAFF2 fallback
and the canonical path for the perturbed variants (methylated glucoses,
charge-modified ATP) where no published parameters exist.
"""
from __future__ import annotations
import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def _check_tool(name: str) -> str:
    p = shutil.which(name)
    if p is None:
        raise FileNotFoundError(
            f"{name!r} not on PATH. Activate $CONTRASCF_AMBERTOOLS_ENV first."
        )
    return p


def parameterize(
    mol2_in: Path,
    out_dir: Path,
    resname: str = "LIG",
    net_charge: int = 0,
    multiplicity: int = 1,
    charge_method: str = "bcc",
    atom_type: str = "gaff2",
) -> dict[str, Path]:
    """Run antechamber → parmchk2 → acpype. Returns dict of output file paths."""
    for t in ("antechamber", "parmchk2", "acpype"):
        _check_tool(t)
    mol2_in = Path(mol2_in).resolve()
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    work = out_dir / "_acpype_work"
    work.mkdir(exist_ok=True)
    shutil.copy(mol2_in, work / f"{resname}.mol2")

    # antechamber: assign GAFF2 atom types, derive AM1-BCC charges.
    subprocess.run([
        "antechamber",
        "-i", f"{resname}.mol2", "-fi", "mol2",
        "-o", f"{resname}_bcc.mol2", "-fo", "mol2",
        "-c", charge_method, "-nc", str(net_charge),
        "-m", str(multiplicity), "-at", atom_type,
        "-rn", resname, "-pf", "yes",
    ], cwd=work, check=True)

    # parmchk2: identify any missing GAFF2 parameters.
    subprocess.run([
        "parmchk2",
        "-i", f"{resname}_bcc.mol2", "-f", "mol2",
        "-o", f"{resname}.frcmod",
    ], cwd=work, check=True)

    # acpype: convert AMBER topology to GROMACS .itp/.gro.
    subprocess.run([
        "acpype", "-i", f"{resname}_bcc.mol2", "-c", "user",
        "-n", str(net_charge), "-a", atom_type, "-o", "gmx",
        "-b", resname,
    ], cwd=work, check=True)

    out = {}
    for src, dst_name in (
        (f"{resname}_bcc.mol2", f"{resname}.mol2"),
        (f"{resname}.frcmod", f"{resname}.frcmod"),
        (f"{resname}.acpype/{resname}_GMX.itp", f"{resname}_GMX.itp"),
        (f"{resname}.acpype/{resname}_GMX.gro", f"{resname}_GMX.gro"),
        (f"{resname}.acpype/posre_{resname}.itp", f"{resname}_posre.itp"),
    ):
        src_p = work / src
        if src_p.exists():
            dst_p = out_dir / dst_name
            shutil.copy(src_p, dst_p)
            out[dst_name] = dst_p
    return out


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--mol2", type=Path, required=True, help="input ligand mol2")
    p.add_argument("--out", type=Path, required=True, help="output dir")
    p.add_argument("--resname", default="LIG")
    p.add_argument("--net-charge", type=int, default=0)
    p.add_argument("--charge-method", default="bcc")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    out = parameterize(
        mol2_in=args.mol2, out_dir=args.out, resname=args.resname,
        net_charge=args.net_charge, charge_method=args.charge_method,
    )
    for k, v in out.items():
        print(f"  {k:24s}  {v}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
