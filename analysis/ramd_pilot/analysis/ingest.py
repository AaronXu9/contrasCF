"""Ingest per-system pilot_results.json files into a per-replica DataFrame.

Schema (one file per system):
    {
      "system_id": str,
      "n_replicas": int,
      "exit_times_ns": [float, ...],   # length n_replicas
      "exited":        [int, ...],     # 1 = event observed, 0 = right-censored at t_max
      "ee_flags":      [int, ...],     # 1 = early-equilibration; excluded from τ
      "force_kcalmolA": float,
    }
"""
from __future__ import annotations
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


REQUIRED_KEYS = ("system_id", "n_replicas", "exit_times_ns", "exited", "ee_flags")


@dataclass
class PilotResult:
    system_id: str
    exit_times_ns: np.ndarray
    exited: np.ndarray
    ee_flags: np.ndarray
    force_kcalmolA: float | None

    @property
    def n_replicas(self) -> int:
        return len(self.exit_times_ns)


def load_pilot_result(path: Path) -> PilotResult:
    blob = json.loads(Path(path).read_text())
    missing = [k for k in REQUIRED_KEYS if k not in blob]
    if missing:
        raise ValueError(f"{path}: missing keys {missing}")
    n = int(blob["n_replicas"])
    arrs = {k: np.asarray(blob[k]) for k in ("exit_times_ns", "exited", "ee_flags")}
    bad = {k: a.shape for k, a in arrs.items() if a.shape != (n,)}
    if bad:
        raise ValueError(f"{path}: array shapes != ({n},): {bad}")
    fc = blob.get("force_kcalmolA")
    return PilotResult(
        system_id=str(blob["system_id"]),
        exit_times_ns=arrs["exit_times_ns"].astype(float),
        exited=arrs["exited"].astype(int),
        ee_flags=arrs["ee_flags"].astype(int),
        force_kcalmolA=float(fc) if fc is not None else None,
    )


def to_dataframe(results: Iterable[PilotResult]) -> pd.DataFrame:
    rows = []
    for r in results:
        for i in range(r.n_replicas):
            rows.append({
                "system_id": r.system_id,
                "replica_idx": i,
                "exit_time_ns": float(r.exit_times_ns[i]),
                "exited": int(r.exited[i]),
                "ee_flag": int(r.ee_flags[i]),
                "force_kcalmolA": r.force_kcalmolA,
            })
    return pd.DataFrame(rows)


def load_directory(results_dir: Path) -> list[PilotResult]:
    results_dir = Path(results_dir)
    files = sorted(results_dir.glob("*.json"))
    if not files:
        raise FileNotFoundError(f"no *.json under {results_dir}")
    return [load_pilot_result(f) for f in files]
