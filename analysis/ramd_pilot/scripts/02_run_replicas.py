"""Submit 15 RAMD replicas per pilot system on CARC.

Reads the calibrated force from `pilot_force_calibration.json` (written
by `01_calibrate_force.py`), then launches 15 independent slurm jobs per
system. Each job runs grompp + 100 ns RAMD on a fresh velocity seed.

Usage on CARC:
    $CONTRASCF_PY analysis/ramd_pilot/scripts/02_run_replicas.py \\
        --system-dir $CONTRASCF_RAMD_OUT/bindingsite_wt \\
        --system-dir $CONTRASCF_RAMD_OUT/bindingsite_pack \\
        --email aoxu@usc.edu
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
_REPO_ROOT = _HERE.parents[3]
sys.path.insert(0, str(_REPO_ROOT))

from analysis.ramd_pilot.config import N_REPLICAS, OUTPUT_ROOT  # noqa: E402
from analysis.ramd_pilot.ramd import replica_runner  # noqa: E402


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--system-dir", type=Path, action="append", required=True,
                   help="repeat for each pilot system")
    p.add_argument("--email", required=True)
    p.add_argument("--ligand-resname", default="ATP")
    p.add_argument("--n-replicas", type=int, default=N_REPLICAS)
    p.add_argument("--force-json", type=Path,
                   default=OUTPUT_ROOT / "pilot_force_calibration.json")
    p.add_argument("--force-kcalmol-A", type=float, default=None,
                   help="override force-json")
    p.add_argument("--max-ns", type=float, default=100.0,
                   help="trajectory length cap per replica (5 = calibration, 100 = production)")
    p.add_argument("--gpu-type", default="v100",
                   help="GPU type to pin (v100, a100, l40s, a40); empty for any (risks P100)")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if args.force_kcalmol_A is not None:
        force = float(args.force_kcalmol_A)
    else:
        if not args.force_json.exists():
            sys.exit(f"force calibration JSON not found at {args.force_json}; "
                     f"run 01_calibrate_force.py first or pass --force-kcalmol-A")
        force = float(json.loads(args.force_json.read_text())["selected_force_kcalmol_A"])
    print(f"[02] using RAMD force = {force:.2f} kcal/mol/Å")

    all_jobs = {}
    for sd in args.system_dir:
        print(f"[02] submitting {args.n_replicas} replicas for {sd.name}  "
              f"(max_ns={args.max_ns}, gpu={args.gpu_type or 'any'})")
        jobs = replica_runner.submit_replicas(
            system_dir=sd, n_replicas=args.n_replicas,
            ramd_force_kcalmol_A=force, email=args.email,
            max_ns=args.max_ns,
            ramd_ligand=args.ligand_resname,
            gpu_type=(args.gpu_type or None),
            dry_run=args.dry_run,
        )
        all_jobs[sd.name] = jobs

    summary_path = OUTPUT_ROOT / "submitted_jobs.json"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(all_jobs, indent=2))
    print(f"[02] submitted {sum(len(v) for v in all_jobs.values())} jobs total")
    print(f"[02] summary: {summary_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
