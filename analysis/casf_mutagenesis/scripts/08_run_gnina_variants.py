"""Run GNINA on every (system, variant) cell with docking/ inputs.

Generalized version of 07_run_gnina_wt.py. Discovers all dirs of shape
    <outputs_root>/<system>/<variant>/docking/{receptor.pdb, ligand.sdf, box.json}
and runs GNINA on each, saving poses to
    <outputs_root>/<system>/<variant>/gnina/{poses.sdf, run.log}

Idempotent: re-running skips cells whose poses.sdf already exists.

Works for both `casf_mutagenesis/outputs/` (binding-site rem/pack/inv,
once mutant docking inputs land) AND `ligand_mutagenesis/outputs/`
(halogenation, methylation, charge-swap).

Env vars:
  CONTRASCF_OUTPUTS_ROOT — where to look. Default
    $CONTRASCF_ROOT/analysis/ligand_mutagenesis/outputs (since
    casf_mutagenesis rem/pack/inv docking inputs don't exist yet).
  CONTRASCF_VARIANT_FILTER — comma-separated list of variants to keep
    (e.g. `halo_Br_1,halo_Cl_1,halo_F_1`). Default: all non-wt variants
    (wt already covered by 07_run_gnina_wt.py).
  CONTRASCF_GNINA_BIN, CONTRASCF_CUDA_DEVICE — same as 07_.

Run:
    source env/lab.sh
    CONTRASCF_OUTPUTS_ROOT=analysis/ligand_mutagenesis/outputs \\
        $CONTRASCF_PY analysis/casf_mutagenesis/scripts/08_run_gnina_variants.py
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

GNINA_BIN = os.environ.get(
    "CONTRASCF_GNINA_BIN",
    "/home/aoxu/projects/PoseBench/forks/GNINA/gnina",
)
GNINA_NUM_MODES = 9
GNINA_EXHAUSTIVENESS = 8
GNINA_SEED = 42
GNINA_CNN = "rescore"
TIMEOUT_S = 1800   # 30 min per cell

OUTPUTS_ROOT = Path(
    os.environ.get(
        "CONTRASCF_OUTPUTS_ROOT",
        str(REPO_ROOT / "analysis" / "ligand_mutagenesis" / "outputs"),
    )
)
# Default variant filter: all non-wt (wt is covered by 07_run_gnina_wt).
DEFAULT_FILTER = None  # None = no filter (all variants including wt)
_filter_env = os.environ.get("CONTRASCF_VARIANT_FILTER")
VARIANT_FILTER = (
    set(v.strip() for v in _filter_env.split(",")) if _filter_env else None
)


def discover_cells(root: Path) -> list[tuple[str, str, Path]]:
    """Return [(system, variant, docking_dir), ...] for every directory of
    shape <root>/<system>/<variant>/docking/ that has the three needed files.
    """
    out: list[tuple[str, str, Path]] = []
    for sys_dir in sorted(root.iterdir()):
        if not sys_dir.is_dir():
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


def run_gnina_one(system: str, variant: str, docking: Path) -> dict:
    """Run GNINA on one (system, variant) cell."""
    entry = {"system": system, "variant": variant}
    v_dir = docking.parent
    out_dir = v_dir / "gnina"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_sdf = out_dir / "poses.sdf"
    if out_sdf.exists():
        entry["status"] = "skip_existing"
        return entry

    box = json.loads((docking / "box.json").read_text())
    cx, cy, cz = box["center"]
    sx, sy, sz = box["size"]
    cmd = [
        GNINA_BIN,
        "-r", str(docking / "receptor.pdb"),
        "-l", str(docking / "ligand.sdf"),
        "--center_x", f"{cx}", "--center_y", f"{cy}", "--center_z", f"{cz}",
        "--size_x", f"{sx}", "--size_y", f"{sy}", "--size_z", f"{sz}",
        "--cnn_scoring", GNINA_CNN,
        "--num_modes", str(GNINA_NUM_MODES),
        "--exhaustiveness", str(GNINA_EXHAUSTIVENESS),
        "--seed", str(GNINA_SEED),
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
            entry["error"] = f"gnina exit {log.returncode}; see {log_path}"
        else:
            entry["status"] = "ok"
            entry["sdf"] = str(out_sdf)
    except subprocess.TimeoutExpired:
        entry["status"] = "timeout"
        entry["wallclock_s"] = TIMEOUT_S
    return entry


def main() -> int:
    print(f"GNINA variants — root: {OUTPUTS_ROOT}")
    if VARIANT_FILTER:
        print(f"variant filter: {sorted(VARIANT_FILTER)}")
    cells = discover_cells(OUTPUTS_ROOT)
    print(f"Discovered {len(cells)} cells with docking/ inputs")
    print(f"GNINA bin: {GNINA_BIN}")

    # Per-module run log so different roots get separate logs
    log_name = f"gnina_variants_{OUTPUTS_ROOT.parent.name}_run_log.json"
    log_path = OUTPUTS_ROOT.parent / "outputs" / log_name if False else OUTPUTS_ROOT / f"../{log_name}"
    # Simpler: put it alongside the outputs root
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
    for i, (system, variant, docking) in enumerate(cells, 1):
        if (system, variant) in seen:
            continue
        print(f"  [{i:4d}/{n_total}] {system}/{variant} ...", flush=True)
        entry = run_gnina_one(system, variant, docking)
        runs.append(entry)
        status = entry.get("status", "?")
        wc = entry.get("wallclock_s", "?")
        print(f"      {status} (wallclock={wc}s)")
        if status == "ok":
            n_done_ok += 1
        elif status in ("skip_existing",):
            n_skip += 1
        else:
            n_fail += 1
        if (n_done_ok + n_fail) % 25 == 0:
            log_path.write_text(json.dumps({"runs": runs}, indent=2))

    log_path.write_text(json.dumps({"runs": runs}, indent=2))
    n_ok_total = sum(1 for r in runs if r.get("status") == "ok")
    print(f"\nTotal cells {n_total} | new ok={n_done_ok} skip={n_skip} fail={n_fail}")
    print(f"Cumulative ok={n_ok_total} in log {log_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
