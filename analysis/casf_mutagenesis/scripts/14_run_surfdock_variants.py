"""Run SurfDock on every (system, variant) cell with `docking/` inputs.

Sibling of 08_run_gnina_variants.py and 11_run_unidock2_variants.py for the
SurfDock engine. Walks `<outputs_root>/<system>/<variant>/docking/` cells
and writes poses to `<outputs_root>/<system>/<variant>/surfdock/poses.sdf`.

Heavy lifting (4-step SurfDock pipeline: surface → CSV → ESM → diffusion)
is delegated to `_run_surfdock_pipeline` from the existing 16-case runner
`analysis/scripts/13_run_surfdock.py`, loaded via importlib because
`scripts/` isn't a Python package.

SurfDock is LAB-BOX-ONLY today — CARC doesn't have the SurfDock conda env,
the model weights, or the precomputed arrays. Don't try to submit this
via SLURM until those are provisioned.

Env vars:
  CONTRASCF_OUTPUTS_ROOT — which outputs dir to walk. Default:
    $CONTRASCF_ROOT/analysis/casf_mutagenesis/outputs
  CONTRASCF_VARIANT_FILTER — comma-separated variants to keep
    (e.g. "wt,rem,pack,inv"). Default: all variants found.
  CONTRASCF_SYSTEM_LIMIT — int, run only this many systems
    (alphabetical order). Default: all. Use for smoke tests.

Run:
    source env/lab.sh
    # subset20-flavoured smoke test (just 1 system, all variants):
    CONTRASCF_OUTPUTS_ROOT=analysis/casf_mutagenesis/outputs \\
        CONTRASCF_SYSTEM_LIMIT=1 \\
        $CONTRASCF_PY analysis/casf_mutagenesis/scripts/14_run_surfdock_variants.py
"""
from __future__ import annotations
import importlib.util
import json
import os
import shutil
import sys
import time
import types
from pathlib import Path

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
DOCKSTRAT_ROOT = Path("/mnt/katritch_lab2/aoxu/CogLigandBench")
SURFDOCK_DIR = "/home/aoxu/projects/SurfDock"
SURFDOCK_PRECOMPUTED_ARRAYS = "/home/aoxu/projects/precomputed/precomputed_arrays"

sys.path.insert(0, str(REPO_ROOT / "analysis" / "src"))


def _load_module(path: Path, mod_name: str):
    """Load a Python file as a module, registering stub parents if needed."""
    parts = mod_name.split(".")
    for i in range(1, len(parts)):
        name = ".".join(parts[:i])
        if name not in sys.modules:
            stub = types.ModuleType(name)
            sys.modules[name] = stub
    spec = importlib.util.spec_from_file_location(mod_name, str(path))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[mod_name] = mod
    spec.loader.exec_module(mod)
    return mod


def _load_runner():
    """Load 13_run_surfdock.py as a module so we can reuse its helpers."""
    script_path = REPO_ROOT / "analysis" / "scripts" / "13_run_surfdock.py"
    return _load_module(script_path, "contrasCF_surfdock_runner_16case")


# ---------------------------------------------------------------------------
# Cell discovery (same contract as 08_/11_)
# ---------------------------------------------------------------------------

def discover_cells(outputs_root: Path,
                   variant_filter: set[str] | None,
                   system_limit: int | None,
                   pdbid_filter: set[str] | None = None) -> list[tuple[str, str, Path]]:
    """Return sorted (pdbid, variant, docking_dir) for every cell that has
    receptor.pdb + ligand.sdf + box.json under
    <outputs_root>/<pdbid>/<variant>/docking/."""
    out: list[tuple[str, str, Path]] = []
    systems = sorted(d for d in outputs_root.iterdir()
                     if d.is_dir() and not d.name.startswith("_"))
    if pdbid_filter is not None:
        systems = [d for d in systems if d.name in pdbid_filter]
    if system_limit is not None:
        systems = systems[:system_limit]
    for sys_dir in systems:
        for v_dir in sorted(sys_dir.iterdir()):
            if not v_dir.is_dir():
                continue
            if variant_filter is not None and v_dir.name not in variant_filter:
                continue
            dock_dir = v_dir / "docking"
            if not all((dock_dir / f).exists()
                       for f in ("receptor.pdb", "ligand.sdf", "box.json")):
                continue
            out.append((sys_dir.name, v_dir.name, dock_dir))
    return out


# ---------------------------------------------------------------------------
# Per-cell runner
# ---------------------------------------------------------------------------

def run_cell(pdbid: str, variant: str, dock_dir: Path, runner_mod) -> dict:
    """Run SurfDock on one cell. Returns a small status dict.

    Bridges the 16-case runner's hardcoded `INPUTS_ROOT/<case>/box.json` path
    to our per-cell `dock_dir/box.json` via a symlink in `docking/inputs/`.
    The receptor + ligand are passed in directly, so they don't need symlinks.
    """
    case = f"{pdbid}_{variant}"
    out_dir = dock_dir.parent / "surfdock"
    out_dir.mkdir(parents=True, exist_ok=True)
    poses_sdf = out_dir / "poses.sdf"
    if poses_sdf.exists():
        return {"pdbid": pdbid, "variant": variant, "status": "skip_existing",
                "poses_sdf": str(poses_sdf)}

    receptor = dock_dir / "receptor.pdb"
    ligand = dock_dir / "ligand.sdf"
    work_dir = out_dir / "surfdock_work"

    # Bridge: ensure docking/inputs/<case>/box.json points at our box.json
    # (16-case runner's _run_surfdock_pipeline reads box.json from that hardcoded path).
    inputs_root_case = REPO_ROOT / "docking" / "inputs" / case
    inputs_root_case.mkdir(parents=True, exist_ok=True)
    bridge_box = inputs_root_case / "box.json"
    if bridge_box.is_symlink() or bridge_box.exists():
        bridge_box.unlink()
    bridge_box.symlink_to((dock_dir / "box.json").resolve())

    t0 = time.time()
    try:
        source_sdfs = runner_mod._run_surfdock_pipeline(
            case, receptor, ligand, work_dir,
        )
        if not source_sdfs:
            raise RuntimeError("SurfDock pipeline produced no poses")
        # Collect ranked SDFs into poses.sdf with score tags.
        runner_mod._collect_poses_sdf(source_sdfs, out_dir)
        shutil.rmtree(work_dir, ignore_errors=True)
        return {"pdbid": pdbid, "variant": variant, "status": "ok",
                "poses_sdf": str(poses_sdf),
                "wallclock_s": round(time.time() - t0, 1)}
    except Exception as exc:
        return {"pdbid": pdbid, "variant": variant, "status": "error",
                "error": f"{type(exc).__name__}: {exc}",
                "wallclock_s": round(time.time() - t0, 1)}


# ---------------------------------------------------------------------------

def main() -> int:
    # 13_run_surfdock.py sets these in its own main(), which we bypass by
    # calling _run_surfdock_pipeline directly -- so set them here too or the
    # SurfDock subprocesses inherit an environment without them.
    os.environ["PROJECT_ROOT"] = str(DOCKSTRAT_ROOT)
    os.environ["SURFDOCK_DIR"] = SURFDOCK_DIR
    os.environ["SURFDOCK_PRECOMPUTED_ARRAYS"] = SURFDOCK_PRECOMPUTED_ARRAYS
    os.environ["precomputed_arrays"] = SURFDOCK_PRECOMPUTED_ARRAYS

    outputs_root = Path(os.environ.get(
        "CONTRASCF_OUTPUTS_ROOT",
        REPO_ROOT / "analysis" / "casf_mutagenesis" / "outputs",
    ))
    outputs_root = outputs_root.resolve() if not outputs_root.is_absolute() else outputs_root
    if not outputs_root.exists():
        print(f"outputs root not found: {outputs_root}", file=sys.stderr)
        return 1

    variant_filter_env = os.environ.get("CONTRASCF_VARIANT_FILTER", "").strip()
    variant_filter: set[str] | None = (
        set(s.strip() for s in variant_filter_env.split(",") if s.strip())
        if variant_filter_env else None
    )
    pdbid_filter_env = os.environ.get("CONTRASCF_PDBID_FILTER", "").strip()
    pdbid_filter: set[str] | None = (
        set(s.strip() for s in pdbid_filter_env.split(",") if s.strip())
        if pdbid_filter_env else None
    )
    system_limit = os.environ.get("CONTRASCF_SYSTEM_LIMIT")
    system_limit = int(system_limit) if system_limit else None

    print(f"SurfDock variant runner")
    print(f"  outputs_root:    {outputs_root}")
    print(f"  variant_filter:  {variant_filter or '(all)'}")
    print(f"  pdbid_filter:    {pdbid_filter or '(all)'}")
    print(f"  system_limit:    {system_limit or '(unlimited)'}")

    cells = discover_cells(outputs_root, variant_filter, system_limit, pdbid_filter)
    print(f"  discovered cells: {len(cells)}")
    if not cells:
        print("nothing to do.")
        return 0

    runner = _load_runner()
    log_path = outputs_root / "surfdock_run_log.json"
    runs: list[dict] = []
    n_ok = n_skip = n_fail = 0

    for i, (pdbid, variant, dock_dir) in enumerate(cells, 1):
        print(f"\n[{i}/{len(cells)}] {pdbid}/{variant} ...", flush=True)
        entry = run_cell(pdbid, variant, dock_dir, runner)
        if entry["status"] == "skip_existing":
            n_skip += 1
        elif entry["status"] == "ok":
            n_ok += 1
            print(f"   done in {entry.get('wallclock_s', 0)}s → {entry['poses_sdf']}")
        else:
            n_fail += 1
            print(f"   FAILED: {entry.get('error')}")
        runs.append(entry)
        log_path.write_text(json.dumps({"runs": runs}, indent=2))

    print(f"\nTotal {len(cells)} | ok={n_ok} skip={n_skip} fail={n_fail}")
    print(f"Run log: {log_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
