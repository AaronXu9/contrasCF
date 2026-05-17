"""Run the τ-RAMD pilot analysis on a directory of pilot_results.json files.

Inputs: a directory containing one .json per system (schema in
`analysis/ramd_pilot/analysis/ingest.py`).

Outputs (written under <out>/):
  - master_table.csv         per-replica records
  - system_stats.json        per-system summary
  - pilot_decision.json      PASS / MARGINAL / FAIL
  - figures/pilot.{pdf,png}  2-panel calibration figure

Usage:
    $CONTRASCF_PY analysis/ramd_pilot/scripts/04_analyze.py \\
        --results-dir <dir-of-pilot_results> --out <out-dir>
"""
from __future__ import annotations
import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

_HERE = Path(__file__).resolve()
_REPO_ROOT = _HERE.parents[3]
sys.path.insert(0, str(_REPO_ROOT))

from analysis.ramd_pilot.config import (  # noqa: E402
    B_BOOT_REPL, FM_GROUND_TRUTH, OUTPUT_ROOT, PILOT_SYSTEMS,
    T_MAX_NS, T_STABLE_THRESH_NS,
)
from analysis.ramd_pilot.analysis import decision as dec  # noqa: E402
from analysis.ramd_pilot.analysis.ingest import load_directory, to_dataframe  # noqa: E402
from analysis.ramd_pilot.analysis.plots import make_pilot_figure  # noqa: E402
from analysis.ramd_pilot.analysis.stats import compute_system_stats  # noqa: E402


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--results-dir", type=Path, required=True,
                   help="directory containing one *.json per system")
    p.add_argument("--out", type=Path, default=None,
                   help="output dir (default: $CONTRASCF_RAMD_OUT/<timestamp>)")
    p.add_argument("--n-boot", type=int, default=B_BOOT_REPL)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    out_dir = args.out or OUTPUT_ROOT / datetime.now().strftime("%Y%m%d-%H%M")
    out_dir.mkdir(parents=True, exist_ok=True)

    results = load_directory(args.results_dir)
    df = to_dataframe(results)
    df.to_csv(out_dir / "master_table.csv", index=False)

    stats_by_system = {}
    for r in results:
        s = compute_system_stats(
            times=r.exit_times_ns, events=r.exited, ee_flags=r.ee_flags,
            t_max=T_MAX_NS, t_stable_thresh=T_STABLE_THRESH_NS,
            system_id=r.system_id, force_kcalmolA=r.force_kcalmolA,
            n_boot=args.n_boot, seed=0,
        )
        stats_by_system[r.system_id] = s
    (out_dir / "system_stats.json").write_text(
        json.dumps({k: v.to_dict() for k, v in stats_by_system.items()}, indent=2)
    )

    decision = dec.evaluate(stats_by_system)
    (out_dir / "pilot_decision.json").write_text(json.dumps(decision.to_dict(), indent=2))

    fm_labels = {sid: FM_GROUND_TRUTH[sid]["label"] for sid in PILOT_SYSTEMS
                 if sid in FM_GROUND_TRUTH}
    make_pilot_figure(
        results=results, fm_labels=fm_labels,
        out_path=out_dir / "figures" / "pilot",
        t_max=T_MAX_NS,
    )

    print(f"[ramd_pilot] decision: {decision.call}  R={decision.R:.2f}  "
          f"τ_KM(WT)={decision.tau_KM_wt:.2f} ns  τ_KM(pack)={decision.tau_KM_pack:.2f} ns")
    print(f"[ramd_pilot] outputs at {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
