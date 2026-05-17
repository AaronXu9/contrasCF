"""Generate synthetic τ-RAMD pilot data for the analysis-pipeline smoke test.

Models:
  * bindingsite_wt   ~ Exponential(rate=1/30 ns); censor at t_max=100 ns.
  * bindingsite_pack ~ Exponential(rate=1/2  ns).

Optional EE injection: 1 random replica marked f_ee=1 in the WT system to
exercise the EE filter without breaking the PASS criterion.
"""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np


def synth_system(
    system_id: str,
    mean_exit_ns: float,
    n_replicas: int = 15,
    t_max: float = 100.0,
    n_ee: int = 0,
    force_kcalmolA: float = 14.0,
    seed: int = 0,
) -> dict:
    rng = np.random.default_rng(seed)
    raw_times = rng.exponential(scale=mean_exit_ns, size=n_replicas)
    exited = (raw_times <= t_max).astype(int)
    times = np.where(exited == 1, raw_times, t_max)
    ee_flags = np.zeros(n_replicas, dtype=int)
    if n_ee > 0:
        ee_idx = rng.choice(n_replicas, size=min(n_ee, n_replicas), replace=False)
        ee_flags[ee_idx] = 1
        # EE replicas get times set to a small value (departed during equil)
        times[ee_idx] = 0.0
        exited[ee_idx] = 1
    return {
        "system_id": system_id,
        "n_replicas": int(n_replicas),
        "exit_times_ns": [float(x) for x in times],
        "exited":       [int(x) for x in exited],
        "ee_flags":     [int(x) for x in ee_flags],
        "force_kcalmolA": float(force_kcalmolA),
    }


def write_synthetic_fixture(out_dir: Path, seed_wt: int = 1, seed_pack: int = 2) -> Path:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    wt = synth_system("bindingsite_wt",   mean_exit_ns=30.0, n_ee=1, seed=seed_wt)
    pk = synth_system("bindingsite_pack", mean_exit_ns=2.0,  n_ee=0, seed=seed_pack)
    (out_dir / "bindingsite_wt.json").write_text(json.dumps(wt, indent=2))
    (out_dir / "bindingsite_pack.json").write_text(json.dumps(pk, indent=2))
    return out_dir


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    out = write_synthetic_fixture(args.out)
    print(f"wrote synthetic pilot fixture to {out}")
