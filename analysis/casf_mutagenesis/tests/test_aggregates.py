"""Tests for paired_stats and bootstrap_memorization_ci.

Run from repo root:
    LD_LIBRARY_PATH=/home/aoxu/miniconda3/envs/rdkit_env/lib:$LD_LIBRARY_PATH \
        /home/aoxu/miniconda3/envs/rdkit_env/bin/python \
        -m pytest analysis/casf_mutagenesis/tests/test_aggregates.py -v
"""
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "analysis"))

from casf_mutagenesis.analysis import (  # noqa: E402
    PredictionRecord, paired_stats, bootstrap_memorization_ci,
)


def _rec(pdbid, variant, model, rmsd, pose=0):
    return PredictionRecord(
        pdbid=pdbid, variant=variant, model=model, pose_idx=pose,
        status="ok", ligand_rmsd_a=rmsd,
    )


def test_paired_stats_basic():
    recs = [
        _rec("1abc", "wt",   "Boltz2", 1.0),
        _rec("1abc", "rem",  "Boltz2", 1.5),
        _rec("1abc", "inv",  "Boltz2", 5.0),
        _rec("1abc", "pack", "Boltz2", 3.0),
    ]
    paired = paired_stats(recs)
    by_variant = {p.variant: p for p in paired}
    assert by_variant["rem"].delta_rmsd_a == 0.5
    assert by_variant["rem"].wt_correct_2A is True
    assert by_variant["rem"].memorized_given_wt is True
    assert by_variant["inv"].memorized_given_wt is False
    assert "wt" not in by_variant


def test_paired_stats_no_wt_record():
    """If WT is missing, paired record reports None deltas, not a crash."""
    recs = [_rec("1abc", "rem", "Boltz2", 1.5)]
    paired = paired_stats(recs)
    assert len(paired) == 1
    assert paired[0].wt_rmsd_a is None
    assert paired[0].delta_rmsd_a is None
    assert paired[0].memorized_given_wt is None


def test_paired_stats_oracle_picks_min_rmsd():
    recs = [
        _rec("1abc", "wt",  "Boltz2", 1.0, pose=0),
        _rec("1abc", "rem", "Boltz2", 5.0, pose=0),
        _rec("1abc", "rem", "Boltz2", 1.5, pose=1),  # oracle should pick this
    ]
    p_top1 = paired_stats(recs, pose_selector="top1")
    p_oracle = paired_stats(recs, pose_selector="oracle")
    assert p_top1[0].adv_rmsd_a == 5.0
    assert p_oracle[0].adv_rmsd_a == 1.5


def test_bootstrap_ci_brackets_point_estimate():
    """Trivial property: lo <= point <= hi, and CI bounded by [0, 1]."""
    recs = [
        _rec(f"sys{i:02d}", "rem", "Boltz2", 1.5 if i < 7 else 4.0)
        for i in range(20)
    ]
    ci = bootstrap_memorization_ci(recs, threshold_a=2.0, n_boot=200, seed=0)
    point, lo, hi = ci[("Boltz2", "rem")]
    assert abs(point - 7/20) < 1e-9
    assert lo <= point <= hi
    assert 0.0 <= lo and hi <= 1.0


def test_bootstrap_resamples_systems_not_cells():
    """Three variants per system collapse to one resampling unit per pdbid."""
    # Two systems; 1abc memorises all variants, 1xyz none.
    recs = []
    for v in ("rem", "pack", "inv"):
        recs.append(_rec("1abc", v, "Boltz2", 1.0))
        recs.append(_rec("1xyz", v, "Boltz2", 5.0))
    ci = bootstrap_memorization_ci(recs, threshold_a=2.0, n_boot=500, seed=0)
    for v in ("rem", "pack", "inv"):
        point, lo, hi = ci[("Boltz2", v)]
        assert point == 0.5  # 1 of 2 pdbids


def test_bootstrap_ignores_nonzero_poses():
    """Even if pose 1 has very different RMSD, only pose 0 enters the CI."""
    recs = [
        _rec("1abc", "rem", "Boltz2", 1.0, pose=0),
        _rec("1abc", "rem", "Boltz2", 50.0, pose=1),  # ignored
        _rec("1xyz", "rem", "Boltz2", 5.0, pose=0),
    ]
    ci = bootstrap_memorization_ci(recs, threshold_a=2.0, n_boot=200, seed=0)
    point, _, _ = ci[("Boltz2", "rem")]
    assert point == 0.5  # only 1abc-pose0 (1.0 Å) counts as memorized
