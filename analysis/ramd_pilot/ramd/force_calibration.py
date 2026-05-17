"""Calibrate the RAMD random-force magnitude on the WT system.

Procedure (plan §3.6):
  - Submit 3 short (5 ns) RAMD replicas at force F = 14 kcal/mol/Å.
  - Parse exit times from prod.log.
  - If geom-mean exit time falls in [5, 50] ns → lock F.
  - If < 1 ns → step F down by 2 kcal/mol/Å, retry.
  - If > 50 ns → step F up by 2 kcal/mol/Å, retry.

This module orchestrates the calibration loop. It does NOT run the sims
directly — it submits via `replica_runner.submit_replicas` with a short
override (5 ns instead of 100 ns) and waits/polls for results.

For the first pass on CARC, run the calibration manually:
  scripts/01_calibrate_force.py --system-dir <wt-dir> --start-force 14
"""
from __future__ import annotations
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np


def parse_exit_time_ns(replica_dir: Path, t_max_ns: float) -> float | None:
    """Read replica_result.json (written by slurm template). None if not done."""
    f = replica_dir / "replica_result.json"
    if not f.exists():
        return None
    blob = json.loads(f.read_text())
    return float(blob["exit_time_ns"])


def calibrate(
    system_dir: Path,
    *,
    start_force_kcalmol_A: float,
    target_lo_ns: float = 5.0,
    target_hi_ns: float = 50.0,
    step_force: float = 2.0,
    n_replicas: int = 3,
    t_max_calib_ns: float = 5.0,
    poll_every_s: int = 300,
    timeout_s: int = 6 * 3600,
    submit_fn=None,           # injectable for tests
) -> dict:
    """Run the calibration loop. Returns the chosen force + per-iter log."""
    if submit_fn is None:
        from .replica_runner import submit_replicas as submit_fn
    log = []
    force = float(start_force_kcalmol_A)
    while True:
        rep_dir = Path(system_dir) / f"_calib_F{force:.1f}"
        rep_dir.mkdir(parents=True, exist_ok=True)
        # Submit 3 replicas with t_max overridden to t_max_calib_ns (caller
        # must use a custom mdp; the runner here only stages slurm files).
        # The caller is expected to substitute the production mdp's nsteps
        # if t_max_calib_ns differs from the default.
        job_ids = submit_fn(
            system_dir=rep_dir, n_replicas=n_replicas,
            ramd_force_kcalmol_A=force, email="aoxu@usc.edu",
        )
        # Poll for completion.
        deadline = time.time() + timeout_s
        exit_times: list[float] = []
        while time.time() < deadline:
            rep_root = rep_dir / "replicas"
            results = []
            for i in range(n_replicas):
                t = parse_exit_time_ns(rep_root / f"r{i:02d}", t_max_calib_ns)
                if t is not None:
                    results.append(t)
            if len(results) == n_replicas:
                exit_times = results
                break
            time.sleep(poll_every_s)
        if not exit_times:
            raise TimeoutError(f"calibration replicas didn't finish for F={force}")

        gmean = float(np.exp(np.mean(np.log(np.clip(exit_times, 1e-3, None)))))
        log.append({"force": force, "exit_times_ns": exit_times, "geom_mean_ns": gmean})
        print(f"[calibrate] F={force:.1f}  exits={exit_times}  geom_mean={gmean:.2f} ns")
        if target_lo_ns <= gmean <= target_hi_ns:
            return {"selected_force_kcalmol_A": force, "log": log}
        if gmean < target_lo_ns:
            force -= step_force                 # too fast → drop
        else:
            force += step_force                 # too slow → raise
        if force <= 4 or force >= 30:
            raise RuntimeError(f"force calibration diverged at F={force}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--system-dir", type=Path, required=True)
    p.add_argument("--start-force", type=float, default=14.0)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    out = calibrate(args.system_dir, start_force_kcalmol_A=args.start_force)
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
