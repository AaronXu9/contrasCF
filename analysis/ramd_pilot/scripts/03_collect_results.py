"""Walk per-replica output dirs and write per-system pilot_results.json.

Per system, expects:
    <system_dir>/replicas/r{NN}/replica_result.json   (written by ramd plugin or post-hoc)
    <system_dir>/replicas/r{NN}/prod.log              (used as fallback)
    <system_dir>/equil/.equil_passed                  (optional EE-flag override)

Output (consumed by analysis/scripts/04_analyze.py):
    <out>/<system>.json   matching the schema in analysis.ramd_pilot.analysis.ingest

Usage on CARC (after all replicas finish):
    $CONTRASCF_PY analysis/ramd_pilot/scripts/03_collect_results.py \\
        --system-dir $CONTRASCF_RAMD_OUT/bindingsite_wt \\
        --system-dir $CONTRASCF_RAMD_OUT/bindingsite_pack \\
        --out $CONTRASCF_RAMD_OUT/pilot_results
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
_REPO_ROOT = _HERE.parents[3]
sys.path.insert(0, str(_REPO_ROOT))

from analysis.ramd_pilot.config import N_REPLICAS, OUTPUT_ROOT, T_MAX_NS  # noqa: E402
from analysis.ramd_pilot.ramd import exit_detection  # noqa: E402


def collect(system_dir: Path, n_replicas: int, t_max_ns: float) -> dict:
    rep_root = Path(system_dir) / "replicas"
    exit_times: list[float] = []
    exited:    list[int]   = []
    ee_flags:  list[int]   = []
    force = None
    for i in range(n_replicas):
        rdir = rep_root / f"r{i:02d}"
        rj = rdir / "replica_result.json"
        if not rj.exists():
            # Fallback: parse prod.log to (re)build replica_result.json.
            if (rdir / "prod.log").exists():
                exit_detection.write_replica_result(
                    rdir, replica_idx=i, t_max_ns=t_max_ns,
                )
            else:
                # Replica missing entirely; treat as failed (censored at 0)
                exit_times.append(0.0); exited.append(0); ee_flags.append(1)
                continue
        blob = json.loads(rj.read_text())
        exit_times.append(float(blob["exit_time_ns"]))
        exited.append(int(blob["exited"]))
        ee_flags.append(int(blob.get("ee_flag", 0)))
        if force is None and blob.get("force_kcalmolA") is not None:
            force = float(blob["force_kcalmolA"])
    return {
        "system_id": system_dir.name,
        "n_replicas": n_replicas,
        "exit_times_ns": exit_times,
        "exited": exited,
        "ee_flags": ee_flags,
        "force_kcalmolA": force,
    }


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--system-dir", type=Path, action="append", required=True)
    p.add_argument("--n-replicas", type=int, default=N_REPLICAS)
    p.add_argument("--t-max-ns", type=float, default=T_MAX_NS)
    p.add_argument("--out", type=Path,
                   default=OUTPUT_ROOT / "pilot_results")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for sd in args.system_dir:
        blob = collect(sd, args.n_replicas, args.t_max_ns)
        out_path = args.out / f"{sd.name}.json"
        out_path.write_text(json.dumps(blob, indent=2))
        n_done = sum(blob["exited"])
        print(f"[03] {sd.name:24s}  {n_done}/{args.n_replicas} exited  -> {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
