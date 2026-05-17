"""End-to-end smoke test: synthetic fixture → ingest → stats → decision → plots.

Runs without any real τ-RAMD data. Exercises the full analysis stack.

Usage (as a script):
    $CONTRASCF_PY -m analysis.ramd_pilot.tests.test_smoke

Asserts:
  * 30 per-replica rows ingested.
  * tau_KM(WT) ≥ 10 ns; tau_KM(pack) ≤ 5 ns; R ≥ 5.
  * f_ee(WT) ≤ 0.2 (we inject 1/15 EE flags = 0.067).
  * decision.call == "PASS".
  * pilot.{pdf,png} created.
"""
from __future__ import annotations
import json
import sys
import tempfile
from pathlib import Path

# Allow running as script from repo root: ensure `analysis` package is importable.
_HERE = Path(__file__).resolve()
_REPO_ROOT = _HERE.parents[3]  # repo/analysis/ramd_pilot/tests/test_smoke.py
sys.path.insert(0, str(_REPO_ROOT))

from analysis.ramd_pilot.config import (
    FM_GROUND_TRUTH, OUTPUT_ROOT, PILOT_SYSTEMS, R_PASS,
    T_MAX_NS, T_STABLE_THRESH_NS, TAU_WT_PASS_NS,
)
from analysis.ramd_pilot.analysis import decision as dec
from analysis.ramd_pilot.analysis.ingest import load_directory, to_dataframe
from analysis.ramd_pilot.analysis.plots import make_pilot_figure
from analysis.ramd_pilot.analysis.stats import compute_system_stats
from analysis.ramd_pilot.tests.synthetic_fixture import write_synthetic_fixture


def run_smoke(out_dir: Path) -> dict:
    out_dir = Path(out_dir)
    fixture_dir = out_dir / "synthetic_fixture"
    write_synthetic_fixture(fixture_dir)

    results = load_directory(fixture_dir)
    df = to_dataframe(results)
    assert len(df) == 30, f"expected 30 rows, got {len(df)}"

    stats_by_system = {}
    for r in results:
        s = compute_system_stats(
            times=r.exit_times_ns, events=r.exited, ee_flags=r.ee_flags,
            t_max=T_MAX_NS, t_stable_thresh=T_STABLE_THRESH_NS,
            system_id=r.system_id, force_kcalmolA=r.force_kcalmolA,
            n_boot=2000, seed=0,  # smaller for speed in smoke test
        )
        stats_by_system[r.system_id] = s

    decision = dec.evaluate(stats_by_system)

    fm_labels = {sid: FM_GROUND_TRUTH[sid]["label"] for sid in PILOT_SYSTEMS
                 if sid in FM_GROUND_TRUTH}
    fig_path = make_pilot_figure(
        results=results, fm_labels=fm_labels,
        out_path=out_dir / "figures" / "pilot",
        t_max=T_MAX_NS,
    )

    df.to_csv(out_dir / "master_table.csv", index=False)
    (out_dir / "system_stats.json").write_text(
        json.dumps({k: v.to_dict() for k, v in stats_by_system.items()}, indent=2)
    )
    (out_dir / "pilot_decision.json").write_text(json.dumps(decision.to_dict(), indent=2))

    # ---- assertions ----
    wt = stats_by_system["bindingsite_wt"]
    pk = stats_by_system["bindingsite_pack"]
    assert wt.tau_KM >= TAU_WT_PASS_NS, f"WT tau_KM={wt.tau_KM} < {TAU_WT_PASS_NS}"
    assert pk.tau_KM <= 5.0, f"pack tau_KM={pk.tau_KM} unexpectedly large"
    assert decision.R >= R_PASS, f"R={decision.R} < {R_PASS}"
    assert wt.f_ee <= 0.2, f"f_ee(WT)={wt.f_ee} > 0.2"
    assert decision.call == "PASS", f"decision={decision.call} reasons={decision.reasons}"
    assert (fig_path.with_suffix(".pdf")).exists()
    assert (fig_path.with_suffix(".png")).exists()

    return {
        "tau_KM_wt": wt.tau_KM, "tau_KM_pack": pk.tau_KM,
        "R": decision.R, "f_ee_wt": wt.f_ee, "call": decision.call,
        "out_dir": str(out_dir),
    }


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="ramd_pilot_smoke_") as td:
        # Persist outputs to module outputs/ so the artifacts survive the test
        persistent = OUTPUT_ROOT / "smoke_test_latest"
        persistent.mkdir(parents=True, exist_ok=True)
        summary = run_smoke(persistent)
    print("[smoke OK]", json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
