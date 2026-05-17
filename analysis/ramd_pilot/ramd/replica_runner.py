"""Submit τ-RAMD replicas as independent SLURM jobs on CARC.

Per-system structure produced:

    <system_dir>/
      replicas/
        r00/  prod.tpr prod.gro prod.log replica_result.json
        r01/  ...
        ...
        r14/

Each replica gets a fresh velocity seed (gen_seed in the production mdp).
The slurm template is rendered with system-specific paths and submitted via
`sbatch`. No SLURM dependency chains — each replica fits in 24 hr.
"""
from __future__ import annotations
import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

from .pbcatom import compute_ramd_pbcatoms

# kJ/mol/nm per kcal/mol/Å (conversion factor: 4.184 * 10 = 41.84)
KJ_MOL_NM_PER_KCAL_MOL_A = 41.84


_SLURM_KEYS = (
    "{JOB_NAME}", "{EMAIL}", "{SYSTEM_DIR}", "{REPLICA_IDX}", "{WALLTIME}",
    "{GMX_RAMD_ROOT}",
)
_MDP_KEYS = (
    "{RAMD_FORCE}", "{RAMD_SEED}", "{RAMD_RECEPTOR}", "{RAMD_LIGAND}",
    "{RAMD_RECEPTOR_PBCATOM}", "{RAMD_LIGAND_PBCATOM}", "{NSTEPS}",
)


def _render_mdp(template: Path, out_path: Path, repl: dict[str, str]) -> Path:
    text = Path(template).read_text()
    for k, v in repl.items():
        text = text.replace(k, v)
    out_path.write_text(text)
    return out_path


def render_slurm(
    template: Path,
    out_path: Path,
    *,
    job_name: str,
    email: str,
    system_dir: Path,
    replica_idx: int,
    ramd_force_kcalmol_A: float,
    ramd_seed: int,
    ramd_receptor: str,
    ramd_ligand: str,
    gmx_ramd_root: Path,
    receptor_pbcatom: int,
    ligand_pbcatom: int,
    nsteps: int,
    walltime: str,
    gpu_gres: str,
    mdp_template: Path,
    mdp_out: Path,
) -> Path:
    """Render the slurm script AND write the per-replica prod.mdp with values inlined.

    The slurm script and mdp share `{X}` placeholder names — substituting the
    slurm script first would clobber the mdp-side `sed s/{X}/{X}/g` patterns.
    Splitting the substitution across two files keeps it idempotent: slurm
    keys go into the slurm script, mdp keys go into prod.mdp.
    """
    force_kjmol_nm = ramd_force_kcalmol_A * KJ_MOL_NM_PER_KCAL_MOL_A
    mdp_repl = {
        "{RAMD_FORCE}": f"{force_kjmol_nm:.4f}",
        "{RAMD_SEED}": str(ramd_seed),
        "{RAMD_RECEPTOR}": ramd_receptor,
        "{RAMD_LIGAND}": ramd_ligand,
        "{RAMD_RECEPTOR_PBCATOM}": str(receptor_pbcatom),
        "{RAMD_LIGAND_PBCATOM}": str(ligand_pbcatom),
        "{NSTEPS}": str(nsteps),
    }
    _render_mdp(mdp_template, mdp_out, mdp_repl)

    text = Path(template).read_text()
    slurm_repl = {
        "{JOB_NAME}": job_name,
        "{EMAIL}": email,
        "{SYSTEM_DIR}": str(system_dir),
        "{REPLICA_IDX}": f"{replica_idx:02d}",
        "{GMX_RAMD_ROOT}": str(gmx_ramd_root),
        "{WALLTIME}": walltime,
        "{GPU_GRES_SPEC}": gpu_gres,
    }
    for k, v in slurm_repl.items():
        text = text.replace(k, v)
    out_path.write_text(text)
    return out_path


def submit_replicas(
    system_dir: Path,
    *,
    n_replicas: int,
    ramd_force_kcalmol_A: float,
    email: str,
    max_ns: float = 100.0,
    ramd_receptor: str = "Protein",
    ramd_ligand: str = "ATP",
    gpu_type: str | None = "v100",
    gmx_ramd_root: Path | None = None,
    template: Path | None = None,
    seed_base: int = 12345,
    dry_run: bool = False,
) -> list[int]:
    system_dir = Path(system_dir).resolve()
    repl_root = system_dir / "replicas"
    repl_root.mkdir(parents=True, exist_ok=True)

    if template is None:
        template = Path(__file__).resolve().parent / "slurm_template.sbatch"
    if gmx_ramd_root is None:
        gmx_ramd_root = Path(os.environ.get("CONTRASCF_GMX_RAMD", "/project2/katritch_223/aoxu/gromacs-ramd-2024.3"))

    # Mirror the per-system mdp dir under <system>/mdp so the template's
    # sed substitution finds it.
    sys_mdp = system_dir / "mdp"
    sys_mdp.mkdir(exist_ok=True)
    src_mdp_dir = Path(__file__).resolve().parents[1] / "prep" / "mdp"
    for mdp in src_mdp_dir.glob("*.mdp"):
        shutil.copy(mdp, sys_mdp / mdp.name)

    # Compute pbcatoms once per system (same for every replica).
    pbc = compute_ramd_pbcatoms(system_dir, receptor=ramd_receptor, ligand=ramd_ligand)
    receptor_pbcatom = pbc[ramd_receptor]
    ligand_pbcatom = pbc[ramd_ligand]
    nsteps = int(max_ns * 1000.0 / 0.002)  # dt = 2 fs
    # Pick walltime: ~80 ns/day on V100 at our system size; 2× safety buffer.
    # Cap at 24 hr to fit CARC partition limits; for full 100 ns production
    # the RAMD plugin's max-dist exit usually terminates earlier than the cap.
    walltime_hr = min(24, max(3, int(round(max_ns / 80.0 * 24.0 * 2))))
    walltime = f"{walltime_hr:02d}:00:00"
    # Build the --gres=gpu:... value. v100 is plentiful and ~50 ns/day on
    # this system; a100 is faster but more contested. None = any GPU
    # (risks P100 landing, ~3× slower than V100).
    gpu_gres = f"{gpu_type}:1" if gpu_type else "1"
    print(f"[replica_runner] {system_dir.name}: pbcatom {ramd_receptor}={receptor_pbcatom} "
          f"{ramd_ligand}={ligand_pbcatom}  max_ns={max_ns} (nsteps={nsteps}, walltime={walltime})")

    job_ids: list[int] = []
    for i in range(n_replicas):
        rdir = repl_root / f"r{i:02d}"
        rdir.mkdir(exist_ok=True)
        slurm_path = rdir / "submit.sbatch"
        mdp_template_path = src_mdp_dir / "step7_production_ramd.mdp"
        render_slurm(
            template=template, out_path=slurm_path,
            job_name=f"ramd_{system_dir.name}_r{i:02d}",
            email=email, system_dir=system_dir, replica_idx=i,
            ramd_force_kcalmol_A=ramd_force_kcalmol_A,
            ramd_seed=seed_base + i,
            ramd_receptor=ramd_receptor, ramd_ligand=ramd_ligand,
            receptor_pbcatom=receptor_pbcatom,
            ligand_pbcatom=ligand_pbcatom,
            gmx_ramd_root=gmx_ramd_root,
            nsteps=nsteps,
            walltime=walltime,
            gpu_gres=gpu_gres,
            mdp_template=mdp_template_path,
            mdp_out=rdir / "prod.mdp",
        )
        if dry_run:
            print(f"[dry-run] would submit {slurm_path}")
            job_ids.append(-1)
            continue
        out = subprocess.check_output(["sbatch", str(slurm_path)], text=True).strip()
        # "Submitted batch job NNNN"
        jid = int(out.split()[-1])
        print(f"  r{i:02d}  job {jid}")
        job_ids.append(jid)
    return job_ids


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--system-dir", type=Path, required=True)
    p.add_argument("--n-replicas", type=int, default=15)
    p.add_argument("--ramd-force-kcalmol-A", type=float, default=14.0)
    p.add_argument("--max-ns", type=float, default=100.0,
                   help="trajectory length cap per replica (5 = calibration, 100 = production)")
    p.add_argument("--gpu-type", default="v100",
                   help="GPU type for --gres (v100, a100, l40s, a40); empty string for any")
    p.add_argument("--email", required=True)
    p.add_argument("--ligand-resname", default="ATP")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    submit_replicas(
        system_dir=args.system_dir, n_replicas=args.n_replicas,
        ramd_force_kcalmol_A=args.ramd_force_kcalmol_A,
        max_ns=args.max_ns,
        email=args.email, ramd_ligand=args.ligand_resname,
        gpu_type=(args.gpu_type or None),
        dry_run=args.dry_run,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
