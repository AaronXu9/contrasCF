"""Calibrate the RAMD random-force magnitude on the WT pilot system.

Submits 3 short (5 ns) RAMD replicas, polls for completion, picks the
force whose geometric-mean exit time falls in [5, 50] ns. Writes the
selected force to `pilot_force_calibration.json` so downstream scripts
can read it.

Usage on CARC:
    $CONTRASCF_PY analysis/ramd_pilot/scripts/01_calibrate_force.py \\
        --wt-dir $CONTRASCF_RAMD_OUT/bindingsite_wt \\
        --start-force 14.0
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
_REPO_ROOT = _HERE.parents[3]
sys.path.insert(0, str(_REPO_ROOT))

from analysis.ramd_pilot.config import OUTPUT_ROOT  # noqa: E402
from analysis.ramd_pilot.ramd import force_calibration  # noqa: E402


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--wt-dir", type=Path, required=True,
                   help="equilibrated WT system_dir")
    p.add_argument("--start-force", type=float, default=14.0)
    p.add_argument("--out", type=Path, default=None,
                   help="default $CONTRASCF_RAMD_OUT/pilot_force_calibration.json")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    out_json = args.out or (OUTPUT_ROOT / "pilot_force_calibration.json")
    out_json.parent.mkdir(parents=True, exist_ok=True)

    result = force_calibration.calibrate(
        system_dir=args.wt_dir, start_force_kcalmol_A=args.start_force,
    )
    out_json.write_text(json.dumps(result, indent=2))
    print(f"[01] selected force = {result['selected_force_kcalmol_A']:.1f} kcal/mol/Å")
    print(f"[01] log: {out_json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
