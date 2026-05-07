"""Run AlphaFold3 (no-MSA / single-sequence mode) on every (system, variant)
under outputs/.

Mirrors `dockstrat.models.alphafold3_inference._run_af3_subprocess`:
  - Uses the conda env at /mnt/katritch_lab2/aoxu/CogLigandBench/envs/alphafold3
  - Calls run_alphafold.py with --norun_data_pipeline --run_inference=true
  - Sets XLA_PYTHON_CLIENT_PREALLOCATE=false and prepends the env's nvidia
    wheel `lib/` directories to LD_LIBRARY_PATH (JAX needs them to load CUDA).

For each (pdbid, variant): reads outputs/<pdbid>/<variant>/af3.json, runs AF3,
copies the rank-0 mmCIF and ranking_scores.csv into the variant dir.

Run:
    LD_LIBRARY_PATH=/home/aoxu/miniconda3/envs/rdkit_env/lib:$LD_LIBRARY_PATH \\
        /home/aoxu/miniconda3/envs/rdkit_env/bin/python \\
        analysis/casf_mutagenesis/scripts/04_run_af3_subset20.py
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

AF3_ENV = Path("/mnt/katritch_lab2/aoxu/CogLigandBench/envs/alphafold3")
AF3_DIR = Path("/mnt/katritch_lab2/aoxu/CogLigandBench/forks/alphafold3/alphafold3")
AF3_MODEL_DIR = Path("/mnt/katritch_lab2/aoxu/CogLigandBench/forks/alphafold3/models")
NUM_DIFFUSION_SAMPLES = 5
NUM_RECYCLES = 10
TIMEOUT_S = 3600
CUDA_DEVICE = "0"
MAX_TOTAL_LENGTH = 800   # skip OOM-risk inputs on a 24 GB GPU


def _build_env() -> dict:
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = CUDA_DEVICE
    env["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
    nvidia_base = AF3_ENV / "lib" / "python3.12" / "site-packages" / "nvidia"
    if nvidia_base.is_dir():
        nvidia_lib_dirs = [
            str(p / "lib") for p in nvidia_base.iterdir()
            if (p / "lib").is_dir()
        ]
        existing = env.get("LD_LIBRARY_PATH", "")
        env["LD_LIBRARY_PATH"] = ":".join(
            nvidia_lib_dirs + ([existing] if existing else [])
        )
    return env


def _seq_total_len(json_path: Path) -> int:
    d = json.loads(json_path.read_text())
    return sum(
        len(s["protein"]["sequence"])
        for s in d.get("sequences", [])
        if "protein" in s
    )


def _run_af3_one(json_path: Path, work_dir: Path) -> Path:
    """Run AF3 once. AF3 writes outputs to work_dir/<name>/ where <name> is
    the `name` field from the JSON. Returns that per-system output dir."""
    work_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(AF3_ENV / "bin" / "python"),
        str(AF3_DIR / "run_alphafold.py"),
        f"--json_path={json_path}",
        f"--output_dir={work_dir}",
        f"--model_dir={AF3_MODEL_DIR}",
        "--norun_data_pipeline",
        "--run_inference=true",
        f"--num_diffusion_samples={NUM_DIFFUSION_SAMPLES}",
        f"--num_recycles={NUM_RECYCLES}",
    ]
    result = subprocess.run(
        cmd, capture_output=True, text=True,
        env=_build_env(), timeout=TIMEOUT_S,
    )
    (work_dir / "af3_stdout.log").write_text(result.stdout)
    (work_dir / "af3_stderr.log").write_text(result.stderr)
    if result.returncode != 0:
        raise RuntimeError(
            f"AF3 failed (exit {result.returncode}); see "
            f"{work_dir}/af3_stderr.log"
        )
    name = json.loads(json_path.read_text())["name"]
    sys_dir = work_dir / name
    if not sys_dir.exists():
        raise FileNotFoundError(f"AF3 output dir missing: {sys_dir}")
    return sys_dir


def _copy_top(af3_sys_dir: Path, name: str, dst_dir: Path) -> dict:
    """Copy the top-ranked AF3 model into dst_dir with an `af3_` prefix so it
    doesn't collide with Boltz outputs in the same dir.

    AF3 v3.0.1 writes:
      <af3_sys_dir>/<name>_ranking_scores.csv
      <af3_sys_dir>/seed-<s>_sample-<n>/<name>_seed-<s>_sample-<n>_model.cif
      <af3_sys_dir>/seed-<s>_sample-<n>/<name>_seed-<s>_sample-<n>_summary_confidences.json
    """
    import csv
    out: dict = {"cif": None, "ranking": None, "conf": None}
    ranking = af3_sys_dir / f"{name}_ranking_scores.csv"
    if not ranking.exists():
        raise FileNotFoundError(f"missing ranking CSV: {ranking}")
    with ranking.open() as f:
        rows = list(csv.DictReader(f))
    rows.sort(key=lambda r: float(r["ranking_score"]), reverse=True)
    top = rows[0]
    seed, sample = int(top["seed"]), int(top["sample"])
    sub = af3_sys_dir / f"seed-{seed}_sample-{sample}"
    src_cif = sub / f"{name}_seed-{seed}_sample-{sample}_model.cif"
    src_conf = sub / f"{name}_seed-{seed}_sample-{sample}_summary_confidences.json"
    if not src_cif.exists():
        raise FileNotFoundError(f"missing model cif: {src_cif}")
    dst_cif = dst_dir / f"af3_{name}_model_0.cif"
    shutil.copy(src_cif, dst_cif)
    out["cif"] = str(dst_cif)
    if src_conf.exists():
        dst_conf = dst_dir / f"af3_summary_confidences_{name}.json"
        shutil.copy(src_conf, dst_conf)
        out["conf"] = str(dst_conf)
    dst_rank = dst_dir / f"af3_ranking_scores_{name}.csv"
    shutil.copy(ranking, dst_rank)
    out["ranking"] = str(dst_rank)
    return out


def main() -> int:
    ids = json.loads(SUBSET20_JSON.read_text())["casf2016"]
    n_total = len(ids) * len(VARIANTS)
    n_done = 0
    n_skip = 0
    n_fail = 0
    log_path = OUTPUT_ROOT / "af3_run_log.json"
    runs: list[dict] = []

    for pdbid in ids:
        for variant in VARIANTS:
            n_done += 1
            v_dir = OUTPUT_ROOT / pdbid / variant
            json_path = v_dir / "af3.json"
            prefix = f"{pdbid}_{variant}"
            af3_cif = v_dir / f"af3_{prefix}_model_0.cif"
            entry = {"pdbid": pdbid, "variant": variant, "prefix": prefix}
            if not json_path.exists():
                entry["status"] = "missing_input"
                n_skip += 1
                runs.append(entry); continue
            if af3_cif.exists():
                entry["status"] = "skip_existing"
                n_skip += 1
                runs.append(entry); continue
            seq_len = _seq_total_len(json_path)
            entry["seq_len"] = seq_len
            if seq_len > MAX_TOTAL_LENGTH:
                entry["status"] = "skip_too_long"
                n_skip += 1
                print(f"  [{n_done}/{n_total}] {prefix} SKIP (len={seq_len} > {MAX_TOTAL_LENGTH})")
                runs.append(entry); continue
            print(f"  [{n_done}/{n_total}] {prefix} (len={seq_len}) ...", flush=True)
            t0 = time.time()
            try:
                work = v_dir / "_af3_work"
                if work.exists():
                    shutil.rmtree(work)
                af3_sys = _run_af3_one(json_path, work)
                outs = _copy_top(af3_sys, prefix, v_dir)
                entry.update({
                    "status": "ok", "outputs": outs,
                    "wallclock_s": round(time.time() - t0, 1),
                })
                shutil.rmtree(work, ignore_errors=True)
                print(f"      done in {entry['wallclock_s']}s")
            except Exception as exc:
                entry.update({
                    "status": "error", "error": str(exc),
                    "wallclock_s": round(time.time() - t0, 1),
                })
                n_fail += 1
                print(f"      FAILED: {exc}")
            runs.append(entry)
            log_path.write_text(json.dumps({"runs": runs}, indent=2))

    print(f"\nTotal {n_total} | done={n_done - n_skip - n_fail} "
          f"skip={n_skip} fail={n_fail}")
    print(f"Run log: {log_path}")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
