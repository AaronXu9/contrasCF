"""Constants for the τ-RAMD pilot module.

Pilot scope: 2 systems (one binder, one non-binder) to de-risk the τ-RAMD
protocol before scaling to the full 13-system FM calibration cohort. See
[docs/ramd_pilot.md] (TBD) and the plan at
`/home/aoxu/.claude/plans/from-the-paper-we-harmonic-snowflake.md`.

Paths are env-aware via the same `_path_env` pattern used in
`analysis/casf_mutagenesis/config.py`.
"""
from __future__ import annotations
import os
from pathlib import Path


def _path_env(name: str, default: str | Path) -> Path:
    val = os.environ.get(name)
    return Path(val) if val else Path(default)


REPO_ROOT = _path_env("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF")
MODULE_ROOT = REPO_ROOT / "analysis" / "ramd_pilot"
OUTPUT_ROOT = _path_env("CONTRASCF_RAMD_OUT", MODULE_ROOT / "outputs")

PILOT_SYSTEMS = ("bindingsite_wt", "bindingsite_pack")
WT_SYSTEM = "bindingsite_wt"
PACK_SYSTEM = "bindingsite_pack"

T_MAX_NS = 100.0
T_STABLE_THRESH_NS = 20.0
N_REPLICAS = 15

R_PASS = 5.0
R_MARGINAL = 2.0
TAU_WT_PASS_NS = 10.0
TAU_WT_MIN_NS = 5.0
F_EE_PASS_MAX = 0.2
F_EE_MARGINAL_MAX = 0.4

B_BOOT_REPL = 10_000

FM_GROUND_TRUTH = {
    "bindingsite_wt":   {"P_bound": 0.61, "SD": 0.05, "label": "binder"},
    "bindingsite_pack": {"P_bound": 0.00, "SD": 0.00, "label": "non-binder"},
}
