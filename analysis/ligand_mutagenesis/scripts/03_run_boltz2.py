"""Run Boltz-2 (with binding-affinity head) on every ligand_mutagenesis cell.

Sibling of `analysis/casf_mutagenesis/scripts/03_run_boltz2_subset20.py` but
specialized for the ligand-mutation module: variant names are not the fixed
{wt, rem, pack, inv} set — they're per-system (halo_F_1, chrg_pos_2, meth_3, …).
So we DISCOVER cells by walking `analysis/ligand_mutagenesis/outputs/<pdbid>/<variant>/boltz.yaml`
instead of iterating a hard-coded (pdbid, variant) list.

For chunked CARC array execution:
  CONTRASCF_START / CONTRASCF_END — slice the sorted pdbid list (not the
  full cell list). Variants for a given pdbid are always processed together.

Boltz-2 parameters and helpers are imported from the casf_mutagenesis module
to keep them in one place.

Run:
    LD_LIBRARY_PATH=/home/aoxu/miniconda3/envs/rdkit_env/lib:$LD_LIBRARY_PATH \\
        /home/aoxu/miniconda3/envs/rdkit_env/bin/python \\
        analysis/ligand_mutagenesis/scripts/03_run_boltz2.py
"""
from __future__ import annotations
import json
import os
import shutil
import sys
import time
from pathlib import Path

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
sys.path.insert(0, str(REPO_ROOT / "analysis"))

# Boltz invocation helpers are vendored below from
# `casf_mutagenesis/scripts/03_run_boltz2_subset20.py` — scripts/ isn't a
# package so direct cross-script import isn't ergonomic, and the helpers
# are small enough to duplicate.
import subprocess  # noqa: E402

from casf_mutagenesis.config import (  # noqa: E402
    BOLTZ_BIN as _BOLTZ_BIN_PATH, CUDA_DEVICE,
    assert_boltz2_binary,
)

BOLTZ_BIN = str(_BOLTZ_BIN_PATH)
_v = assert_boltz2_binary(BOLTZ_BIN)
print(f"Boltz binary: {BOLTZ_BIN} (v{_v})")

DIFFUSION_SAMPLES = 5
RECYCLING_STEPS = 3
SAMPLING_STEPS = 200
SEED = 42
MAX_TOTAL_LENGTH = 800

OUTPUTS_ROOT = REPO_ROOT / "analysis" / "ligand_mutagenesis" / "outputs"


# --------------------------- vendored from casf -----------------------------

def _seq_total_len(yaml_path: Path) -> int:
    total = 0
    for line in yaml_path.read_text().splitlines():
        line = line.strip()
        if line.startswith("sequence:"):
            total += len(line.split(":", 1)[1].strip())
    return total


def _run_boltz_one(yaml_path: Path, prefix: str, work_dir: Path) -> Path:
    work_dir.mkdir(parents=True, exist_ok=True)
    yaml_renamed = work_dir / f"{prefix}.yaml"
    shutil.copy(yaml_path, yaml_renamed)
    cmd = [
        BOLTZ_BIN, "predict", str(yaml_renamed),
        "--out_dir", str(work_dir),
        "--model", "boltz2",
        "--output_format", "mmcif",
        "--diffusion_samples", str(DIFFUSION_SAMPLES),
        "--recycling_steps", str(RECYCLING_STEPS),
        "--sampling_steps", str(SAMPLING_STEPS),
        "--seed", str(SEED),
    ]
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = CUDA_DEVICE
    result = subprocess.run(cmd, capture_output=True, text=True, env=env)
    (work_dir / "boltz_stdout.log").write_text(result.stdout)
    (work_dir / "boltz_stderr.log").write_text(result.stderr)
    if result.returncode != 0:
        raise RuntimeError(
            f"boltz predict failed (exit {result.returncode}); "
            f"see {work_dir}/boltz_stderr.log"
        )
    pred_dir = work_dir / f"boltz_results_{prefix}" / "predictions" / prefix
    if not pred_dir.exists():
        raise FileNotFoundError(f"missing predictions dir: {pred_dir}")
    return pred_dir


def _copy_all_samples(pred_dir: Path, prefix: str, dst_dir: Path) -> dict:
    cifs = sorted(pred_dir.glob(f"{prefix}_model_*.cif"))
    if not cifs:
        raise FileNotFoundError(f"no model cifs in {pred_dir}")
    out: dict = {"cifs": [], "confs": [], "affinity": None}
    for src_cif in cifs:
        dst_cif = dst_dir / src_cif.name
        shutil.copy(src_cif, dst_cif)
        out["cifs"].append(str(dst_cif))
        src_conf = pred_dir / f"confidence_{src_cif.stem}.json"
        if src_conf.exists():
            dst_conf = dst_dir / src_conf.name
            shutil.copy(src_conf, dst_conf)
            out["confs"].append(str(dst_conf))
    aff = pred_dir / f"affinity_{prefix}.json"
    if aff.exists():
        dst_aff = dst_dir / f"affinity_{prefix}.json"
        shutil.copy(aff, dst_aff)
        out["affinity"] = str(dst_aff)
    return out


# --------------------------- ligand-side discovery --------------------------

def discover_cells(outputs_root: Path) -> list[tuple[str, str, Path]]:
    """Return sorted (pdbid, variant, yaml_path) for every boltz.yaml under
    `outputs_root/<pdbid>/<variant>/`. Skips `_backup*` subtrees."""
    out: list[tuple[str, str, Path]] = []
    for sys_dir in sorted(outputs_root.iterdir()):
        if not sys_dir.is_dir() or sys_dir.name.startswith("_"):
            continue
        for v_dir in sorted(sys_dir.iterdir()):
            if not v_dir.is_dir():
                continue
            yaml = v_dir / "boltz.yaml"
            if yaml.exists():
                out.append((sys_dir.name, v_dir.name, yaml))
    return out


def main() -> int:
    cells = discover_cells(OUTPUTS_ROOT)
    # Optional pdbid-level slicing for CARC array jobs
    all_pdbids = sorted({pid for pid, _, _ in cells})
    start = int(os.environ.get("CONTRASCF_START", "0"))
    end = int(os.environ.get("CONTRASCF_END", str(len(all_pdbids))))
    slice_pdbids = set(all_pdbids[start:end])
    cells = [c for c in cells if c[0] in slice_pdbids]
    print(f"ligand_mutagenesis Boltz-2: pdbid slice [{start}:{end}] "
          f"= {len(slice_pdbids)} systems, {len(cells)} cells")
    n_total = len(cells)
    n_done = n_skip = n_fail = 0
    log_path = OUTPUTS_ROOT / "boltz2_run_log.json"
    runs: list[dict] = []
    for pdbid, variant, yaml_path in cells:
        n_done += 1
        v_dir = yaml_path.parent
        prefix = f"{pdbid}_{variant}"
        existing_cif = v_dir / f"{prefix}_model_0.cif"
        entry: dict = {"pdbid": pdbid, "variant": variant, "prefix": prefix}
        if existing_cif.exists():
            entry["status"] = "skip_existing"
            n_skip += 1
            runs.append(entry); continue
        seq_len = _seq_total_len(yaml_path)
        entry["seq_len"] = seq_len
        if seq_len > MAX_TOTAL_LENGTH:
            entry["status"] = "skip_too_long"
            n_skip += 1
            print(f"  [{n_done}/{n_total}] {prefix} SKIP (len={seq_len} > {MAX_TOTAL_LENGTH})")
            runs.append(entry); continue
        print(f"  [{n_done}/{n_total}] {prefix} (len={seq_len}) ...", flush=True)
        t0 = time.time()
        try:
            work = v_dir / "_boltz_work"
            if work.exists():
                shutil.rmtree(work)
            pred_dir = _run_boltz_one(yaml_path, prefix, work)
            outs = _copy_all_samples(pred_dir, prefix, v_dir)
            entry.update({"status": "ok", "outputs": outs,
                          "wallclock_s": round(time.time() - t0, 1)})
            shutil.rmtree(work, ignore_errors=True)
            print(f"      done in {entry['wallclock_s']}s")
        except Exception as exc:
            entry.update({"status": "error", "error": str(exc),
                          "wallclock_s": round(time.time() - t0, 1)})
            n_fail += 1
            print(f"      FAILED: {exc}")
        runs.append(entry)
        log_path.write_text(json.dumps({"runs": runs}, indent=2))
    print(f"\nTotal {n_total} | done={n_done - n_skip - n_fail} skip={n_skip} fail={n_fail}")
    print(f"Run log: {log_path}")
    # Always return 0 — per-cell failures shouldn't kill the SLURM script.
    return 0


if __name__ == "__main__":
    sys.exit(main())
