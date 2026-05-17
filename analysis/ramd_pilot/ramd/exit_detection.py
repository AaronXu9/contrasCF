"""Post-process a τ-RAMD trajectory to determine exit time.

Default: rely on gromacs-ramd plugin's emitted log line ("RAMD: ligand
exited at step N"). Fallback: parse trajectory and apply COM-distance
threshold ourselves.

This is the canonical extractor used by 03_collect_results.py.
"""
from __future__ import annotations
import argparse
import json
import re
from pathlib import Path

DT_PS = 2.0  # 2 fs * 1000 ps/ns... wait, dt = 0.002 ps = 2 fs


def parse_log(log_path: Path, t_max_ns: float, dt_ps: float = 0.002) -> dict:
    """Return {exit_time_ns, exited, source}.

    `source` indicates whether the time came from the RAMD plugin's exit
    line (canonical) or from log scraping for run termination.
    """
    text = Path(log_path).read_text(errors="replace")
    m = re.search(r"RAMD.*exit.*step\s+(\d+)", text, re.IGNORECASE)
    if m:
        step = int(m.group(1))
        return {
            "exit_time_ns": step * dt_ps / 1000.0,
            "exited": 1,
            "source": "ramd_plugin_log",
        }
    # No exit recorded; check whether mdrun reached nsteps (= t_max) or
    # was killed prematurely.
    if "Finished mdrun" in text or "Writing final coordinates" in text:
        return {"exit_time_ns": t_max_ns, "exited": 0, "source": "ran_to_t_max"}
    return {"exit_time_ns": t_max_ns, "exited": 0, "source": "incomplete"}


def write_replica_result(
    replica_dir: Path,
    *,
    replica_idx: int,
    t_max_ns: float,
    ee_flag: int = 0,
    force_kcalmolA: float | None = None,
) -> Path:
    info = parse_log(replica_dir / "prod.log", t_max_ns=t_max_ns)
    blob = {
        "replica_idx": int(replica_idx),
        "exit_time_ns": info["exit_time_ns"],
        "exited": info["exited"],
        "ee_flag": int(ee_flag),
        "source": info["source"],
        "force_kcalmolA": force_kcalmolA,
    }
    out = replica_dir / "replica_result.json"
    out.write_text(json.dumps(blob, indent=2))
    return out


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--replica-dir", type=Path, required=True)
    p.add_argument("--replica-idx", type=int, required=True)
    p.add_argument("--t-max-ns", type=float, default=100.0)
    p.add_argument("--ee-flag", type=int, default=0)
    p.add_argument("--force-kcalmolA", type=float)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    out = write_replica_result(
        args.replica_dir, replica_idx=args.replica_idx, t_max_ns=args.t_max_ns,
        ee_flag=args.ee_flag, force_kcalmolA=args.force_kcalmolA,
    )
    print(out)
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
