"""End-to-end system builder for one pilot complex.

Steps:
  1. Stage protein PDB + ligand SDF/mol2 → system_dir/inputs/
  2. Parameterize ligand via prep.ligand_params.parameterize
  3. Build solvated, ionized GROMACS system via prep.system_build.build
  4. Run 5-stage equilibration via prep.equilibrate.equilibrate
  5. Acceptance gate (RMSD plateau check) — TODO when on CARC.

Usage on CARC (after `source env/carc.sh` and after both
`carc_setup/install_gromacs_ramd.sh` and
`carc_setup/setup_ambertools_env.sh` have completed):

    conda activate $CONTRASCF_AMBERTOOLS_ENV
    source $CONTRASCF_GMX_RAMD/bin/GMXRC
    $CONTRASCF_PY analysis/ramd_pilot/scripts/00_build_systems.py \\
        --system bindingsite_wt \\
        --protein-pdb <path/to/1B38_chainA.pdb> \\
        --ligand-mol2 <path/to/ATP.mol2> \\
        --ligand-resname ATP --net-charge -4 \\
        --out $CONTRASCF_RAMD_OUT/bindingsite_wt
"""
from __future__ import annotations
import argparse
import shutil
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
_REPO_ROOT = _HERE.parents[3]
sys.path.insert(0, str(_REPO_ROOT))

from analysis.ramd_pilot.prep import equilibrate, ligand_params, system_build  # noqa: E402


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--system", required=True,
                   help="system_id, e.g. bindingsite_wt")
    p.add_argument("--protein-pdb", type=Path, required=True)
    p.add_argument("--ligand-mol2", type=Path, required=True)
    p.add_argument("--ligand-resname", default="ATP")
    p.add_argument("--net-charge", type=int, default=-4,
                   help="ATP at physiological pH = -4")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--gmx", default="gmx_mpi")
    p.add_argument("--ff", default=None,
                   help="protein force field (default: amber99sb-ildn). "
                        "Set to e.g. amber14sb_OL15 if a ported .ff dir is in $GMXLIB.")
    p.add_argument("--water", default=None, help="water model (default: tip3p)")
    p.add_argument("--skip-equil", action="store_true",
                   help="Stop after system build (do not run equilibration).")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)

    # --- 1. Stage inputs ---
    inp = out / "inputs"
    inp.mkdir(exist_ok=True)
    shutil.copy(args.protein_pdb, inp / "protein.pdb")
    shutil.copy(args.ligand_mol2, inp / "ligand.mol2")

    # --- 2. Parameterize ligand ---
    print(f"[00] parameterizing {args.ligand_resname} (net charge {args.net_charge})")
    lig_out = ligand_params.parameterize(
        mol2_in=inp / "ligand.mol2", out_dir=out / "ligand_params",
        resname=args.ligand_resname, net_charge=args.net_charge,
    )
    lig_gro = lig_out[f"{args.ligand_resname}_GMX.gro"]
    lig_itp = lig_out[f"{args.ligand_resname}_GMX.itp"]

    # --- 3. Build solvated, ionized system ---
    print(f"[00] building solvated system in {out}")
    build_kw = dict(
        protein_pdb=inp / "protein.pdb", ligand_gro=lig_gro,
        ligand_itp=lig_itp, out_dir=out,
        ligand_resname=args.ligand_resname, gmx=args.gmx,
    )
    if args.ff:
        build_kw["ff"] = args.ff
    if args.water:
        build_kw["water"] = args.water
    built = system_build.build(**build_kw)
    print(f"[00] system: {built}")

    if args.skip_equil:
        print("[00] --skip-equil set; stopping before equilibration")
        return 0

    # --- 4. Equilibrate ---
    print("[00] running 5-stage equilibration")
    mdp_dir = _REPO_ROOT / "analysis" / "ramd_pilot" / "prep" / "mdp"
    final_gro = equilibrate.equilibrate(out, mdp_dir, gmx=args.gmx)
    print(f"[00] DONE  final equilibrated structure: {final_gro}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
