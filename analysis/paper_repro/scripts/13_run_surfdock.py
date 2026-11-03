"""Run SurfDock on the 16 prepared cases.

For each case:
  1. Invoke SurfDock's 4-step pipeline (surface → CSV → ESM → diffusion) on
     docking/inputs/<case>/{receptor.pdb, ligand.sdf} — the same inputs
     UniDock2 and GNINA use (AF3-predicted receptor).
  2. Concatenate the ranked SDFs into contrasCF/data/SurfDock/<case>/poses.sdf,
     embedding each record with SDF tags (rank, confidence, rmsd, sample idx)
     parsed from SurfDock's original filenames.
  3. Build combined_top1.pdb = receptor.pdb + top-1 pose HETATM block via
     docking_io.write_combined_top1 (shared with 11_run_docking.py).

Unlike dockstrat's run_single (which uses a tempdir and loses original
filenames), we drive the 4 steps explicitly with a persistent working
directory so confidence scores are preserved in the SDF tags.

Run:
    /home/aoxu/miniconda3/envs/rdkit_env/bin/python \
        analysis/paper_repro/scripts/13_run_surfdock.py
"""
from __future__ import annotations
import glob
import json
import os
import re
import shutil
import sys
import tempfile
import time
from pathlib import Path


sys.path.insert(0, str(Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF")) / "analysis" / "paper_repro" / "lib"))
sys.path.insert(0, str(Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF")) / "analysis"))
from core.surfdock_engine import (  # noqa: E402,F401
    BATCH_SIZE, DOCKSTRAT_ROOT, INPUTS_ROOT, NUM_POSES, REPO_ROOT, SAMPLES_PER_COMPLEX, SURFDOCK_DIR,
    SURFDOCK_ENV_PREFIX, SURFDOCK_PRECOMPUTED_ARRAYS, SURFDOCK_WEIGHTS, _collect_poses_sdf,
    _load_surfdock_module, _parse_surfdock_score, _run_inference_with_pocket_center,
    _run_surfdock_pipeline, _sdf_record_with_tags, _surfdock_config, _translate_ligand_to_pocket,
)


from config import CASES, DATA_ROOT, SCOPES, cases_in_scope  # noqa: E402
from docking_io import write_combined_top1  # noqa: E402



SURFDOCK_OUT_ROOT = DATA_ROOT / "SurfDock"



















def run_case(case: str) -> dict:
    out_dir = SURFDOCK_OUT_ROOT / case
    out_dir.mkdir(parents=True, exist_ok=True)
    receptor = INPUTS_ROOT / case / "receptor.pdb"
    ligand = INPUTS_ROOT / case / "ligand.sdf"
    combined = out_dir / "combined_top1.pdb"
    poses_sdf = out_dir / "poses.sdf"

    if combined.exists() and poses_sdf.exists():
        return {"case": case, "ok": True, "sec": 0.0, "skipped": True,
                "combined": str(combined)}

    t0 = time.time()
    try:
        work_dir = out_dir / "surfdock_work"
        source_sdfs = _run_surfdock_pipeline(case, receptor, ligand, work_dir)
        if not source_sdfs:
            raise RuntimeError(f"SurfDock produced no poses for {case}")

        poses_sdf = _collect_poses_sdf(source_sdfs, out_dir)
        combined = write_combined_top1(receptor, out_dir / "rank1.sdf", out_dir)
        return {
            "case": case,
            "ok": True,
            "sec": time.time() - t0,
            "poses_sdf": str(poses_sdf),
            "combined": str(combined),
            "n_poses": len(source_sdfs),
        }
    except Exception as e:
        return {
            "case": case,
            "ok": False,
            "sec": time.time() - t0,
            "error": f"{type(e).__name__}: {e}",
        }


def main() -> None:
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--scope", choices=sorted(SCOPES), default="all")
    args = p.parse_args()
    cases = cases_in_scope(args.scope)
    print(f"[surfdock] scope={args.scope}  n_cases={len(cases)}", flush=True)

    os.environ["PROJECT_ROOT"] = str(DOCKSTRAT_ROOT)
    os.environ["SURFDOCK_DIR"] = SURFDOCK_DIR
    os.environ["SURFDOCK_PRECOMPUTED_ARRAYS"] = SURFDOCK_PRECOMPUTED_ARRAYS
    os.environ["precomputed_arrays"] = SURFDOCK_PRECOMPUTED_ARRAYS
    os.environ["SURFDOCK_ENV_PREFIX"] = SURFDOCK_ENV_PREFIX

    SURFDOCK_OUT_ROOT.mkdir(parents=True, exist_ok=True)
    results = []
    for case in cases:
        print(f"[surfdock] {case} ...", flush=True)
        info = run_case(case)
        tag = "ok" if info.get("ok") else "FAIL"
        print(
            f"[surfdock]   {tag} ({info.get('sec', 0):.1f}s)"
            + (f"  err={info['error']}" if not info.get("ok") else "")
            + ("  [skipped]" if info.get("skipped") else "")
        )
        results.append(info)

    (DATA_ROOT / "surfdock_runs.json").write_text(json.dumps(results, indent=2))
    ok = sum(1 for r in results if r.get("ok"))
    print(f"[surfdock] done. ok={ok}/{len(results)}")


if __name__ == "__main__":
    main()
