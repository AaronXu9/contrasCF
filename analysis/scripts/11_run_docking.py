"""Run UniDock2 and GNINA on the 16 prepared cases.

For each case, invoke both tools against the prepared inputs and save:
  contrasCF/data/UniDock2/<case>/{poses.sdf, run.log, combined_top1.pdb}
  contrasCF/data/GNINA/<case>/{poses.sdf,   run.log, combined_top1.pdb}

`combined_top1.pdb` = receptor (protein) + top-1 docked pose embedded as a
HETATM block (resname LIG, chain Z). This is what the existing analysis
pipeline (`select_target_ligand`) consumes unchanged.

Run:
    LD_LIBRARY_PATH=/home/aoxu/miniconda3/envs/rdkit_env/lib:$LD_LIBRARY_PATH \
        /home/aoxu/miniconda3/envs/rdkit_env/bin/python \
        analysis/scripts/11_run_docking.py
"""
from __future__ import annotations
import json
import os
import subprocess
import sys
import time
from pathlib import Path

REPO_ROOT = Path("/mnt/katritch_lab2/aoxu/contrasCF")
sys.path.insert(0, str(REPO_ROOT / "analysis" / "src"))

from config import CASES, DATA_ROOT  # noqa: E402
from docking_io import write_combined_top1  # noqa: E402

INPUTS_ROOT = REPO_ROOT / "docking" / "inputs"
UNIDOCK2_OUT_ROOT = DATA_ROOT / "UniDock2"
GNINA_OUT_ROOT = DATA_ROOT / "GNINA"

GNINA_BIN = "/home/aoxu/projects/PoseBench/forks/GNINA/gnina"
UNIDOCK2_CONDA_ENV = "unidock2"

# GNINA defaults.
GNINA_NUM_MODES = 9
GNINA_EXHAUSTIVENESS = 8
GNINA_SEED = 42
GNINA_CNN = "rescore"


def _load_box(case: str) -> dict:
    return json.loads((INPUTS_ROOT / case / "box.json").read_text())


def run_unidock2(case: str) -> tuple[Path, str]:
    box = _load_box(case)
    cx, cy, cz = box["center"]
    sx, sy, sz = box["size"]
    receptor = INPUTS_ROOT / case / "receptor.pdb"
    ligand = INPUTS_ROOT / case / "ligand.sdf"
    out_dir = UNIDOCK2_OUT_ROOT / case
    out_dir.mkdir(parents=True, exist_ok=True)
    out_sdf = out_dir / "poses.sdf"
    log_path = out_dir / "run.log"

    # UniDock2 size is set via YAML config; pass via -cf if non-default.
    # The CLI only supports --center; size is inherited from YAML defaults of
    # [30,30,30]. We override here by writing a tiny YAML config.
    cfg_path = out_dir / "unidock2_config.yaml"
    cfg_path.write_text(
        "Settings:\n"
        f"    size: [{sx}, {sy}, {sz}]\n"
        "Advanced:\n"
        f"    seed: {GNINA_SEED}\n"
        f"    num_pose: 9\n"
    )
    cmd = [
        "conda", "run", "-n", UNIDOCK2_CONDA_ENV,
        "unidock2", "docking",
        "-r", str(receptor),
        "-l", str(ligand),
        "-c", f"{cx}", f"{cy}", f"{cz}",
        "-cf", str(cfg_path),
        "-o", str(out_sdf),
    ]
    log = subprocess.run(cmd, capture_output=True, text=True)
    log_path.write_text("STDOUT:\n" + log.stdout + "\nSTDERR:\n" + log.stderr)
    if log.returncode != 0 or not out_sdf.exists():
        raise RuntimeError(f"UniDock2 failed for {case}: rc={log.returncode}; see {log_path}")
    return out_sdf, log.stdout


def run_gnina(case: str) -> tuple[Path, str]:
    box = _load_box(case)
    cx, cy, cz = box["center"]
    sx, sy, sz = box["size"]
    receptor = INPUTS_ROOT / case / "receptor.pdb"
    ligand = INPUTS_ROOT / case / "ligand.sdf"
    out_dir = GNINA_OUT_ROOT / case
    out_dir.mkdir(parents=True, exist_ok=True)
    out_sdf = out_dir / "poses.sdf"
    log_path = out_dir / "run.log"

    cmd = [
        GNINA_BIN,
        "-r", str(receptor),
        "-l", str(ligand),
        "--center_x", f"{cx}", "--center_y", f"{cy}", "--center_z", f"{cz}",
        "--size_x", f"{sx}", "--size_y", f"{sy}", "--size_z", f"{sz}",
        "--cnn_scoring", GNINA_CNN,
        "--num_modes", str(GNINA_NUM_MODES),
        "--exhaustiveness", str(GNINA_EXHAUSTIVENESS),
        "--seed", str(GNINA_SEED),
        "-o", str(out_sdf),
    ]
    log = subprocess.run(cmd, capture_output=True, text=True)
    log_path.write_text("STDOUT:\n" + log.stdout + "\nSTDERR:\n" + log.stderr)
    if log.returncode != 0 or not out_sdf.exists():
        raise RuntimeError(f"GNINA failed for {case}: rc={log.returncode}; see {log_path}")
    return out_sdf, log.stdout


def run_case(case: str) -> dict:
    info = {"case": case, "unidock2": {}, "gnina": {}}
    receptor = INPUTS_ROOT / case / "receptor.pdb"

    t0 = time.time()
    try:
        sdf, _ = run_unidock2(case)
        combined = write_combined_top1(receptor, sdf, UNIDOCK2_OUT_ROOT / case)
        info["unidock2"] = {"ok": True, "sdf": str(sdf), "combined": str(combined), "sec": time.time() - t0}
    except Exception as e:
        info["unidock2"] = {"ok": False, "error": f"{type(e).__name__}: {e}", "sec": time.time() - t0}

    t0 = time.time()
    try:
        sdf, _ = run_gnina(case)
        combined = write_combined_top1(receptor, sdf, GNINA_OUT_ROOT / case)
        info["gnina"] = {"ok": True, "sdf": str(sdf), "combined": str(combined), "sec": time.time() - t0}
    except Exception as e:
        info["gnina"] = {"ok": False, "error": f"{type(e).__name__}: {e}", "sec": time.time() - t0}

    return info


def main():
    results = []
    for case in CASES:
        print(f"[dock] {case} ...", flush=True)
        info = run_case(case)
        u = info["unidock2"]; g = info["gnina"]
        print(f"[dock]   UniDock2: {'ok' if u.get('ok') else 'FAIL'} ({u.get('sec', 0):.1f}s)"
              + (f"  err={u.get('error')}" if not u.get('ok') else ""))
        print(f"[dock]   GNINA   : {'ok' if g.get('ok') else 'FAIL'} ({g.get('sec', 0):.1f}s)"
              + (f"  err={g.get('error')}" if not g.get('ok') else ""))
        results.append(info)

    (DATA_ROOT / "docking_runs.json").write_text(json.dumps(results, indent=2))
    ok_u = sum(1 for r in results if r["unidock2"].get("ok"))
    ok_g = sum(1 for r in results if r["gnina"].get("ok"))
    print(f"[dock] done. UniDock2 ok={ok_u}/{len(results)}; GNINA ok={ok_g}/{len(results)}")


if __name__ == "__main__":
    main()
