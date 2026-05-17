"""Two-panel pilot figure: Kaplan-Meier curves + per-replica strip plot.

Panel A — KM survival of ligand-in-pocket vs time, one curve per system.
Panel B — Strip plot of per-replica exit times, censored replicas marked.
"""
from __future__ import annotations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # noqa: E402
import matplotlib.pyplot as plt
import numpy as np

from .ingest import PilotResult
from .stats import kaplan_meier


_COLORS = {"binder": "#2c7fb8", "non-binder": "#d7301f"}


def _plot_km(ax, results: list[PilotResult], labels: dict[str, str], t_max: float):
    for r in results:
        # Exclude EE replicas
        mask = r.ee_flags == 0
        t = r.exit_times_ns[mask]
        e = r.exited[mask]
        if len(t) == 0:
            continue
        t_grid, S_grid = kaplan_meier(t, e)
        # Step plot up to t_max
        t_plot = np.append(t_grid, t_max)
        S_plot = np.append(S_grid, S_grid[-1])
        label = labels.get(r.system_id, r.system_id)
        ax.step(t_plot, S_plot, where="post",
                label=f"{r.system_id}  ({label})",
                color=_COLORS.get(label, None), linewidth=2)
    ax.axhline(0.5, linestyle="--", color="grey", linewidth=0.8)
    ax.set_xlabel("time (ns)")
    ax.set_ylabel("P(ligand still in pocket)")
    ax.set_xlim(0, t_max)
    ax.set_ylim(-0.02, 1.02)
    ax.set_title("A. Kaplan-Meier ligand-in-pocket survival")
    ax.legend(loc="lower left", fontsize=9)


def _plot_strip(ax, results: list[PilotResult], labels: dict[str, str], t_max: float):
    rng = np.random.default_rng(0)
    for i, r in enumerate(results):
        x = i + rng.uniform(-0.15, 0.15, size=r.n_replicas)
        # Censored = exited == 0; mark differently
        exited = r.exited == 1
        ee = r.ee_flags == 1
        non_ee_exited = exited & ~ee
        non_ee_cens = ~exited & ~ee
        label = labels.get(r.system_id, r.system_id)
        c = _COLORS.get(label, "grey")
        ax.scatter(x[non_ee_exited], r.exit_times_ns[non_ee_exited],
                   marker="o", s=50, color=c, edgecolor="black", linewidth=0.5,
                   label=f"{label} (exited)" if i == 0 else None)
        ax.scatter(x[non_ee_cens], r.exit_times_ns[non_ee_cens],
                   marker="^", s=50, color=c, edgecolor="black", linewidth=0.5,
                   label="censored at t_max" if i == 0 else None)
        ax.scatter(x[ee], r.exit_times_ns[ee],
                   marker="x", s=50, color="black",
                   label="EE flag" if i == 0 and ee.any() else None)
    ax.set_xticks(range(len(results)))
    ax.set_xticklabels([r.system_id for r in results], rotation=15, ha="right")
    ax.set_yscale("log")
    ax.set_ylabel("exit time (ns, log)")
    ax.set_ylim(0.05, t_max * 1.5)
    ax.axhline(t_max, color="grey", linestyle=":", linewidth=0.8)
    ax.set_title("B. Per-replica exit times")
    ax.legend(loc="lower right", fontsize=8)


def make_pilot_figure(
    results: list[PilotResult],
    fm_labels: dict[str, str],
    out_path: Path,
    t_max: float,
) -> Path:
    """Render the 2-panel pilot figure to {out_path}.pdf and .png."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    _plot_km(axes[0], results, fm_labels, t_max)
    _plot_strip(axes[1], results, fm_labels, t_max)
    fig.tight_layout()
    pdf = out_path.with_suffix(".pdf")
    png = out_path.with_suffix(".png")
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, bbox_inches="tight", dpi=150)
    plt.close(fig)
    return out_path
