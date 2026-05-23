"""Run UniDock2 on every (system, variant) cell with docking/ inputs.

Sibling of 08_run_gnina_variants.py — same `discover_cells` logic, same
docking/ contract (receptor.pdb + ligand.sdf + box.json), but the
engine is UniDock2 (Vina-style scoring, GPU-accelerated by AutoDock-GPU).

For each cell:
  - Write a tiny per-cell UniDock2 YAML (Settings.size + Advanced.seed +
    num_pose). UniDock2 only takes `--center` on the CLI; the box size
    must come from the YAML config.
  - Invoke `conda run -n unidock2 unidock2 docking -r ... -l ... -c x y z -cf cfg.yaml -o poses.sdf`.
  - Save log + poses to `<v_dir>/unidock2/{poses.sdf, run.log}`.

Idempotent: skips cells whose `unidock2/poses.sdf` already exists.

Env vars (same family as 08_):
  CONTRASCF_OUTPUTS_ROOT — outputs/ root to walk
  CONTRASCF_VARIANT_FILTER — comma-separated variants (default: all)

Run:
    source env/lab.sh
    CONTRASCF_OUTPUTS_ROOT=analysis/casf_mutagenesis/outputs \\
        $CONTRASCF_PY analysis/casf_mutagenesis/scripts/11_run_unidock2_variants.py
"""
from __future__ import annotations
import json
import os
import subprocess
import sys
import time
from pathlib import Path

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
sys.path.insert(0, str(REPO_ROOT / "analysis"))

# UniDock2 lives in a dedicated conda env on lab.
UNIDOCK2_CONDA_ENV = os.environ.get("CONTRASCF_UNIDOCK2_ENV", "unidock2")
UNIDOCK2_BIN = os.environ.get(
    "CONTRASCF_UNIDOCK2_BIN",
    f"/home/aoxu/miniconda3/envs/{UNIDOCK2_CONDA_ENV}/bin/unidock2",
)
SEED = 42
NUM_POSE = 9
# Per-cell hard timeout. Long UniDock2 runs (>5 min) on outlier systems
# (e.g. 3d4z/meth_3, 3n86/chrg_neu_methyl) appear to trip the NVIDIA
# kernel module into a NVML mismatch state — every subsequent cell then
# errors at ~1 s for the rest of the run. Cap at 120 s to skip those
# rather than crashing the driver. Cells slower than this are too rare /
# too risky to be worth chasing here.
TIMEOUT_S = int(os.environ.get("CONTRASCF_UNIDOCK2_TIMEOUT_S", "120"))

# Abort the script entirely after this many consecutive fast (<5 s)
# errors — strong signal that CUDA is broken and we'd just be wasting
# cycles running 500+ dud cells. Restart after fixing the driver.
FAST_ERROR_ABORT = int(os.environ.get("CONTRASCF_FAST_ERROR_ABORT", "5"))

OUTPUTS_ROOT = Path(
    os.environ.get(
        "CONTRASCF_OUTPUTS_ROOT",
        str(REPO_ROOT / "analysis" / "casf_mutagenesis" / "outputs"),
    )
)
_filter_env = os.environ.get("CONTRASCF_VARIANT_FILTER")
VARIANT_FILTER = (
    set(v.strip() for v in _filter_env.split(",")) if _filter_env else None
)


def discover_cells(root: Path) -> list[tuple[str, str, Path]]:
    """Same as 08_run_gnina_variants — find all <root>/<sys>/<var>/docking/ cells."""
    out: list[tuple[str, str, Path]] = []
    for sys_dir in sorted(root.iterdir()):
        if not sys_dir.is_dir() or len(sys_dir.name) != 4:
            continue
        for var_dir in sorted(sys_dir.iterdir()):
            if not var_dir.is_dir():
                continue
            docking = var_dir / "docking"
            if not docking.is_dir():
                continue
            if not all((docking / f).exists()
                       for f in ("receptor.pdb", "ligand.sdf", "box.json")):
                continue
            variant = var_dir.name
            if VARIANT_FILTER is not None and variant not in VARIANT_FILTER:
                continue
            out.append((sys_dir.name, variant, docking))
    return out


def run_unidock2_one(system: str, variant: str, docking: Path) -> dict:
    entry = {"system": system, "variant": variant}
    v_dir = docking.parent
    out_dir = v_dir / "unidock2"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_sdf = out_dir / "poses.sdf"
    if out_sdf.exists():
        entry["status"] = "skip_existing"
        return entry

    box = json.loads((docking / "box.json").read_text())
    cx, cy, cz = box["center"]
    sx, sy, sz = box["size"]

    cfg_path = out_dir / "unidock2_config.yaml"
    cfg_path.write_text(
        "Settings:\n"
        f"    size: [{sx}, {sy}, {sz}]\n"
        "Advanced:\n"
        f"    seed: {SEED}\n"
        f"    num_pose: {NUM_POSE}\n"
    )
    cmd = [
        "conda", "run", "-n", UNIDOCK2_CONDA_ENV,
        "unidock2", "docking",
        "-r", str(docking / "receptor.pdb"),
        "-l", str(docking / "ligand.sdf"),
        "-c", f"{cx}", f"{cy}", f"{cz}",
        "-cf", str(cfg_path),
        "-o", str(out_sdf),
    ]
    log_path = out_dir / "run.log"
    t0 = time.time()
    try:
        log = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT_S)
        log_path.write_text(
            f"CMD: {' '.join(cmd)}\n\nSTDOUT:\n{log.stdout}\n\nSTDERR:\n{log.stderr}"
        )
        entry["wallclock_s"] = round(time.time() - t0, 1)
        if log.returncode != 0 or not out_sdf.exists():
            entry["status"] = "error"
            entry["error"] = f"unidock2 exit {log.returncode}; see {log_path}"
        else:
            entry["status"] = "ok"
            entry["sdf"] = str(out_sdf)
    except subprocess.TimeoutExpired:
        entry["status"] = "timeout"
        entry["wallclock_s"] = TIMEOUT_S
    return entry


def main() -> int:
    print(f"UniDock2 variants — root: {OUTPUTS_ROOT}")
    if VARIANT_FILTER:
        print(f"variant filter: {sorted(VARIANT_FILTER)}")
    cells = discover_cells(OUTPUTS_ROOT)
    print(f"Discovered {len(cells)} cells with docking/ inputs")
    print(f"UniDock2 env: {UNIDOCK2_CONDA_ENV}")

    log_name = f"unidock2_variants_{OUTPUTS_ROOT.parent.name}_run_log.json"
    log_path = OUTPUTS_ROOT.parent / log_name
    runs: list[dict] = []
    if log_path.exists():
        try:
            runs = json.loads(log_path.read_text()).get("runs", [])
        except Exception:
            runs = []
    seen = {(r["system"], r["variant"]) for r in runs if r.get("status") == "ok"}

    n_total = len(cells)
    n_done_ok = 0
    n_skip = 0
    n_fail = 0
    n_consec_fast_errors = 0   # detect GPU-died pattern (~1 s per error)
    for i, (system, variant, docking) in enumerate(cells, 1):
        if (system, variant) in seen:
            continue
        print(f"  [{i:4d}/{n_total}] {system}/{variant} ...", flush=True)
        entry = run_unidock2_one(system, variant, docking)
        runs.append(entry)
        status = entry.get("status", "?")
        wc = entry.get("wallclock_s", "?")
        print(f"      {status} (wallclock={wc}s)")
        if status == "ok":
            n_done_ok += 1
            n_consec_fast_errors = 0
        elif status == "skip_existing":
            n_skip += 1
        else:
            n_fail += 1
            # Heuristic: an error that returns in <5 s is almost certainly
            # a CUDA-init failure (GPU is dead). Multiple in a row means
            # the GPU went down — abort rather than burn through hundreds
            # of cells producing 1-s fake errors.
            try:
                if isinstance(wc, (int, float)) and wc < 5.0:
                    n_consec_fast_errors += 1
                else:
                    n_consec_fast_errors = 0
            except Exception:
                n_consec_fast_errors = 0
            if n_consec_fast_errors >= FAST_ERROR_ABORT:
                log_path.write_text(json.dumps({"runs": runs}, indent=2))
                print(
                    f"\nABORT: {FAST_ERROR_ABORT} consecutive <5 s errors — "
                    f"likely GPU driver desync. Run `nvidia-smi` and recover "
                    f"the driver before retrying.",
                    flush=True,
                )
                return 2
        if (n_done_ok + n_fail) % 25 == 0:
            log_path.write_text(json.dumps({"runs": runs}, indent=2))

    log_path.write_text(json.dumps({"runs": runs}, indent=2))
    n_ok_total = sum(1 for r in runs if r.get("status") == "ok")
    print(f"\nTotal cells {n_total} | new ok={n_done_ok} skip={n_skip} fail={n_fail}")
    print(f"Cumulative ok={n_ok_total} in log {log_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
