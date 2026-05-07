"""Run Boltz-2 on every (system, variant) under outputs/.

Iterates over outputs/<pdbid>/<variant>/boltz.yaml, runs `boltz predict`, and
copies the top-ranked mmCIF + confidence JSON into the same directory using
Boltz-1 naming (`<pdbid>_<variant>_model_0.cif`,
`confidence_<pdbid>_<variant>_model_0.json`).

Run:
    LD_LIBRARY_PATH=/home/aoxu/miniconda3/envs/rdkit_env/lib:$LD_LIBRARY_PATH \\
        /home/aoxu/miniconda3/envs/rdkit_env/bin/python \\
        analysis/casf_mutagenesis/scripts/03_run_boltz2_subset20.py
"""
from __future__ import annotations
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO_ROOT = Path("/mnt/katritch_lab2/aoxu/contrasCF")
sys.path.insert(0, str(REPO_ROOT / "analysis"))

from casf_mutagenesis.config import OUTPUT_ROOT, SUBSET20_JSON, VARIANTS

BOLTZ_BIN = "/home/aoxu/miniconda3/envs/boltzina_env/bin/boltz"
DIFFUSION_SAMPLES = 5
RECYCLING_STEPS = 3
SAMPLING_STEPS = 200
SEED = 42
CUDA_DEVICE = "0"

# Skip outliers that may OOM on a 24 GB RTX 4090. 3dx2 is 1016 residues
# (homo-tetramer), 4ih5/4ih7 are 562 residues — keep those, only skip the
# truly huge case.
MAX_TOTAL_LENGTH = 800


def _seq_total_len(yaml_path: Path) -> int:
    """Sum of all `sequence:` lengths in the YAML (cheap line scan)."""
    total = 0
    for line in yaml_path.read_text().splitlines():
        line = line.strip()
        if line.startswith("sequence:"):
            total += len(line.split(":", 1)[1].strip())
    return total


def _run_boltz_one(yaml_path: Path, prefix: str, work_dir: Path) -> Path:
    """Boltz uses the input YAML's filename stem as the output prefix, so we
    copy `boltz.yaml` → `<prefix>.yaml` inside the work dir before invoking."""
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


def _copy_top(pred_dir: Path, prefix: str, dst_dir: Path) -> dict:
    src_cif = pred_dir / f"{prefix}_model_0.cif"
    src_conf = pred_dir / f"confidence_{prefix}_model_0.json"
    if not src_cif.exists():
        cifs = sorted(pred_dir.glob(f"{prefix}_model_*.cif"))
        if not cifs:
            raise FileNotFoundError(f"no model cifs in {pred_dir}")
        src_cif = cifs[0]
        src_conf = pred_dir / f"confidence_{src_cif.stem}.json"
    dst_cif = dst_dir / f"{prefix}_model_0.cif"
    shutil.copy(src_cif, dst_cif)
    out = {"cif": str(dst_cif), "conf": None, "affinity": None}
    if src_conf.exists():
        dst_conf = dst_dir / f"confidence_{prefix}_model_0.json"
        shutil.copy(src_conf, dst_conf)
        out["conf"] = str(dst_conf)
    aff = pred_dir / f"affinity_{prefix}.json"
    if aff.exists():
        dst_aff = dst_dir / f"affinity_{prefix}.json"
        shutil.copy(aff, dst_aff)
        out["affinity"] = str(dst_aff)
    return out


def main() -> int:
    ids = json.loads(SUBSET20_JSON.read_text())["casf2016"]
    n_total = len(ids) * len(VARIANTS)
    n_done = 0
    n_skip = 0
    n_fail = 0
    log_path = OUTPUT_ROOT / "boltz2_run_log.json"
    runs: list[dict] = []

    for pdbid in ids:
        for variant in VARIANTS:
            n_done += 1
            v_dir = OUTPUT_ROOT / pdbid / variant
            yaml_path = v_dir / "boltz.yaml"
            prefix = f"{pdbid}_{variant}"
            existing_cif = v_dir / f"{prefix}_model_0.cif"
            entry = {"pdbid": pdbid, "variant": variant, "prefix": prefix}
            if not yaml_path.exists():
                entry["status"] = "missing_input"
                n_skip += 1
                runs.append(entry); continue
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
                outs = _copy_top(pred_dir, prefix, v_dir)
                entry.update({"status": "ok", "outputs": outs,
                              "wallclock_s": round(time.time() - t0, 1)})
                # tidy up: remove the work dir to save disk
                shutil.rmtree(work, ignore_errors=True)
                print(f"      done in {entry['wallclock_s']}s")
            except Exception as exc:
                entry.update({"status": "error", "error": str(exc),
                              "wallclock_s": round(time.time() - t0, 1)})
                n_fail += 1
                print(f"      FAILED: {exc}")
            runs.append(entry)
            # checkpoint after each run so a long batch is recoverable
            log_path.write_text(json.dumps({"runs": runs}, indent=2))

    print(f"\nTotal {n_total} | done={n_done - n_skip - n_fail} "
          f"skip={n_skip} fail={n_fail}")
    print(f"Run log: {log_path}")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
