"""Per-system summary statistics for τ-RAMD pilot output.

Implements:
  * Kaplan-Meier survival + median exit time (right-censored at t_max).
  * Geometric mean over uncensored, non-EE replicas (censored pinned at t_max).
  * f_stable (fraction of non-EE replicas with exit > t_stable_thresh).
  * f_ee (fraction of replicas flagged early-equilibration).
  * Wilson 95% CIs for binomial proportions.
  * Replica-level bootstrap CI on the KM median.

All functions are dependency-light: numpy + scipy.stats only.
"""
from __future__ import annotations
from dataclasses import asdict, dataclass

import numpy as np
from scipy.stats import norm


# --- KM ---------------------------------------------------------------------

def kaplan_meier(times: np.ndarray, events: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Product-limit estimate. Returns (t_grid, S_grid) starting at (0, 1)."""
    times = np.asarray(times, dtype=float)
    events = np.asarray(events, dtype=int)
    n = len(times)
    if n == 0:
        return np.array([0.0]), np.array([1.0])

    order = np.argsort(times, kind="stable")
    t = times[order]
    e = events[order]

    out_t = [0.0]
    out_S = [1.0]
    surv = 1.0
    n_at_risk = n
    i = 0
    while i < n:
        ti = t[i]
        j = i
        events_here = 0
        while j < n and t[j] == ti:
            if e[j] == 1:
                events_here += 1
            j += 1
        if events_here > 0 and n_at_risk > 0:
            surv *= (n_at_risk - events_here) / n_at_risk
            out_t.append(ti)
            out_S.append(surv)
        n_at_risk -= (j - i)
        i = j
    return np.array(out_t), np.array(out_S)


def km_median(t_grid: np.ndarray, S_grid: np.ndarray, t_max: float) -> float:
    """Earliest time at which S(t) ≤ 0.5. Returns t_max if undefined."""
    below = np.where(S_grid <= 0.5)[0]
    if len(below) == 0:
        return float(t_max)
    return float(t_grid[below[0]])


# --- per-replica filters ---------------------------------------------------

def _filter_non_ee(times: np.ndarray, events: np.ndarray, ee_flags: np.ndarray):
    mask = np.asarray(ee_flags) == 0
    return np.asarray(times)[mask], np.asarray(events)[mask]


def tau_KM(times, events, ee_flags, t_max: float) -> float:
    t, e = _filter_non_ee(times, events, ee_flags)
    if len(t) == 0:
        return 0.0
    t_grid, S_grid = kaplan_meier(t, e)
    return km_median(t_grid, S_grid, t_max)


def tau_geom(times, events, ee_flags, t_max: float) -> float:
    """Geometric mean. Censored replicas pinned at t_max (lower-bound estimate)."""
    t, e = _filter_non_ee(times, events, ee_flags)
    if len(t) == 0:
        return float("nan")
    t_eff = np.where(np.asarray(e) == 1, t, t_max)
    t_eff = np.clip(t_eff, 1e-9, None)  # log safety
    return float(np.exp(np.mean(np.log(t_eff))))


def f_stable(times, events, ee_flags, t_stable_thresh: float) -> float:
    t, _ = _filter_non_ee(times, events, ee_flags)
    if len(t) == 0:
        return 0.0
    return float(np.mean(t > t_stable_thresh))


def f_ee(ee_flags) -> float:
    ee_flags = np.asarray(ee_flags)
    if len(ee_flags) == 0:
        return 0.0
    return float(np.mean(ee_flags == 1))


# --- uncertainty ------------------------------------------------------------

def wilson_ci(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    """Wilson 95% CI for a binomial proportion."""
    if n == 0:
        return 0.0, 1.0
    z = norm.ppf(1 - alpha / 2)
    phat = k / n
    denom = 1 + z * z / n
    centre = (phat + z * z / (2 * n)) / denom
    half = z * np.sqrt(phat * (1 - phat) / n + z * z / (4 * n * n)) / denom
    return float(max(0.0, centre - half)), float(min(1.0, centre + half))


def bootstrap_tau_KM(times, events, ee_flags, t_max: float, n_boot: int = 10_000,
                     seed: int | None = 0) -> tuple[float, float]:
    """Replica-level bootstrap 95% CI for tau_KM."""
    rng = np.random.default_rng(seed)
    times = np.asarray(times); events = np.asarray(events); ee_flags = np.asarray(ee_flags)
    n = len(times)
    if n == 0:
        return 0.0, 0.0
    medians = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        medians[b] = tau_KM(times[idx], events[idx], ee_flags[idx], t_max)
    lo, hi = np.percentile(medians, [2.5, 97.5])
    return float(lo), float(hi)


# --- aggregate per-system stats --------------------------------------------

@dataclass
class SystemStats:
    system_id: str
    n_replicas: int
    n_eff: int                # non-EE replicas
    tau_KM: float             # ns
    tau_KM_ci: tuple[float, float]   # bootstrap 95% CI
    tau_geom: float           # ns; nan if all EE
    f_stable: float
    f_stable_ci: tuple[float, float]
    f_ee: float
    f_ee_ci: tuple[float, float]
    n_exited: int
    force_kcalmolA: float | None

    def to_dict(self) -> dict:
        return asdict(self)


def compute_system_stats(
    times, events, ee_flags,
    t_max: float, t_stable_thresh: float,
    system_id: str,
    force_kcalmolA: float | None = None,
    n_boot: int = 10_000,
    seed: int | None = 0,
) -> SystemStats:
    times = np.asarray(times); events = np.asarray(events); ee_flags = np.asarray(ee_flags)
    n = len(times)
    n_ee = int((ee_flags == 1).sum())
    n_eff = n - n_ee
    n_exited = int(((events == 1) & (ee_flags == 0)).sum())

    tau_km = tau_KM(times, events, ee_flags, t_max)
    tau_km_ci = bootstrap_tau_KM(times, events, ee_flags, t_max, n_boot=n_boot, seed=seed)
    tau_g = tau_geom(times, events, ee_flags, t_max)
    fs = f_stable(times, events, ee_flags, t_stable_thresh)
    fs_ci = wilson_ci(int(round(fs * n_eff)), n_eff)
    fe = f_ee(ee_flags)
    fe_ci = wilson_ci(n_ee, n)
    return SystemStats(
        system_id=system_id, n_replicas=n, n_eff=n_eff,
        tau_KM=tau_km, tau_KM_ci=tau_km_ci, tau_geom=tau_g,
        f_stable=fs, f_stable_ci=fs_ci,
        f_ee=fe, f_ee_ci=fe_ci,
        n_exited=n_exited, force_kcalmolA=force_kcalmolA,
    )
