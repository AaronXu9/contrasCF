"""Run AF3 on subset20 with MSAs (replaces no-MSA results).

Two-phase pipeline:
  1. For each unique WT sequence: fetch a uniref+bfd merged A3M via Boltz's
     `--use_msa_server` machinery (cached by sha1 in outputs/_msa_cache/).
     ~2 s + ~6 s Boltz minimum-inference per system; ~150 s total for 19.
  2. For each (pdbid, variant): take the WT A3M, rewrite the query row at
     mutation positions, render af3.json with the MSA injected, run AF3,
     copy outputs as `af3msa_<prefix>_*` (the no-MSA `af3_*` outputs from
     `04_run_af3_subset20.py` are preserved as the baseline).

Run:
    LD_LIBRARY_PATH=/home/aoxu/miniconda3/envs/rdkit_env/lib:$LD_LIBRARY_PATH \\
        /home/aoxu/miniconda3/envs/rdkit_env/bin/python \\
        analysis/casf_mutagenesis/scripts/06_run_af3_msa_subset20.py
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
from casf_mutagenesis.inputs_af3 import (
    clean_a3m_for_af3, render_af3, rewrite_a3m_query,
)
from casf_mutagenesis.msa_via_boltz import fetch_msa_via_boltz

AF3_ENV = Path("/mnt/katritch_lab2/aoxu/CogLigandBench/envs/alphafold3")
AF3_DIR = Path("/mnt/katritch_lab2/aoxu/CogLigandBench/forks/alphafold3/alphafold3")
AF3_MODEL_DIR = Path("/mnt/katritch_lab2/aoxu/CogLigandBench/forks/alphafold3/models")
NUM_DIFFUSION_SAMPLES = 5
NUM_RECYCLES = 10
TIMEOUT_S = 3600
CUDA_DEVICE = "0"
MAX_TOTAL_LENGTH = 800

MSA_CACHE = OUTPUT_ROOT / "_msa_cache"


def _build_af3_env() -> dict:
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


def _read_af3_protein_sequence(af3_json_path: Path) -> str:
    """Extract the (single) protein sequence from a no-MSA af3.json."""
    d = json.loads(af3_json_path.read_text())
    for s in d["sequences"]:
        if "protein" in s:
            return s["protein"]["sequence"]
    raise RuntimeError(f"no protein in {af3_json_path}")


def _read_ligand_smiles(af3_json_path: Path) -> str:
    d = json.loads(af3_json_path.read_text())
    for s in d["sequences"]:
        if "ligand" in s:
            return s["ligand"]["smiles"]
    raise RuntimeError(f"no ligand in {af3_json_path}")


def _run_af3(json_path: Path, work_dir: Path) -> Path:
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
        env=_build_af3_env(), timeout=TIMEOUT_S,
    )
    (work_dir / "af3_stdout.log").write_text(result.stdout)
    (work_dir / "af3_stderr.log").write_text(result.stderr)
    if result.returncode != 0:
        raise RuntimeError(
            f"AF3 failed (exit {result.returncode}); "
            f"see {work_dir}/af3_stderr.log"
        )
    name = json.loads(json_path.read_text())["name"]
    sys_dir = work_dir / name
    if not sys_dir.exists():
        raise FileNotFoundError(f"AF3 output dir missing: {sys_dir}")
    return sys_dir


def _copy_top(af3_sys_dir: Path, name: str, dst_dir: Path) -> dict:
    """Copy the top-ranked AF3 model with `af3msa_` prefix."""
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
    dst_cif = dst_dir / f"af3msa_{name}_model_0.cif"
    shutil.copy(src_cif, dst_cif)
    out["cif"] = str(dst_cif)
    if src_conf.exists():
        dst_conf = dst_dir / f"af3msa_summary_confidences_{name}.json"
        shutil.copy(src_conf, dst_conf)
        out["conf"] = str(dst_conf)
    dst_rank = dst_dir / f"af3msa_ranking_scores_{name}.csv"
    shutil.copy(ranking, dst_rank)
    out["ranking"] = str(dst_rank)
    return out


def main() -> int:
    ids = json.loads(SUBSET20_JSON.read_text())["casf2016"]

    # Phase 1: ensure WT MSA is cached for each system
    print("Phase 1: fetch MSAs for WT sequences\n")
    wt_msas: dict[str, str] = {}
    for i, pdbid in enumerate(ids, 1):
        wt_json = OUTPUT_ROOT / pdbid / "wt" / "af3.json"
        if not wt_json.exists():
            print(f"  [{i}/{len(ids)}] {pdbid}: no af3.json — skip"); continue
        wt_seq = _read_af3_protein_sequence(wt_json)
        if len(wt_seq) > MAX_TOTAL_LENGTH:
            print(f"  [{i}/{len(ids)}] {pdbid}: len={len(wt_seq)} > {MAX_TOTAL_LENGTH} — skip")
            continue
        t0 = time.time()
        try:
            a3m = fetch_msa_via_boltz(wt_seq, cache_dir=MSA_CACHE)
            n_seqs = sum(1 for ln in a3m.splitlines() if ln.startswith(">"))
            wt_msas[pdbid] = a3m
            print(f"  [{i}/{len(ids)}] {pdbid}: MSA ok "
                  f"(n_seqs={n_seqs}, len={len(wt_seq)}, {time.time()-t0:.1f}s)")
        except Exception as exc:
            print(f"  [{i}/{len(ids)}] {pdbid}: MSA fetch failed: {exc}")

    # Phase 2: AF3 with MSA for each (pdbid, variant)
    print("\nPhase 2: run AF3 with MSA for each (pdbid, variant)\n")
    n_total = len(ids) * len(VARIANTS)
    n_done = 0
    n_skip = 0
    n_fail = 0
    log_path = OUTPUT_ROOT / "af3_msa_run_log.json"
    runs: list[dict] = []

    for pdbid in ids:
        wt_a3m = wt_msas.get(pdbid)
        if wt_a3m is None:
            for v in VARIANTS:
                runs.append({"pdbid": pdbid, "variant": v, "status": "skip_no_msa"})
                n_skip += 1
                n_done += 1
            continue
        for variant in VARIANTS:
            n_done += 1
            v_dir = OUTPUT_ROOT / pdbid / variant
            af3_json = v_dir / "af3.json"
            prefix = f"{pdbid}_{variant}"
            af3msa_cif = v_dir / f"af3msa_{prefix}_model_0.cif"
            entry = {"pdbid": pdbid, "variant": variant, "prefix": prefix}
            if not af3_json.exists():
                entry["status"] = "missing_input"
                n_skip += 1
                runs.append(entry); continue
            if af3msa_cif.exists():
                entry["status"] = "skip_existing"
                n_skip += 1
                runs.append(entry); continue
            variant_seq = _read_af3_protein_sequence(af3_json)
            smiles = _read_ligand_smiles(af3_json)
            try:
                variant_a3m = rewrite_a3m_query(wt_a3m, variant_seq)
                variant_a3m = clean_a3m_for_af3(variant_a3m, len(variant_seq))
            except Exception as exc:
                entry.update({"status": "error", "error": f"a3m rewrite: {exc}"})
                n_fail += 1
                runs.append(entry); continue

            # Render fresh JSON with MSA
            json_msa = v_dir / "af3_msa.json"
            render_af3(
                name=prefix, chain_seqs=[("A", variant_seq)],
                ligand_smiles=smiles, out_path=json_msa,
                chain_msas={"A": variant_a3m},
            )
            print(f"  [{n_done}/{n_total}] {prefix} ...", flush=True)
            t0 = time.time()
            try:
                work = v_dir / "_af3_msa_work"
                if work.exists():
                    shutil.rmtree(work)
                af3_sys = _run_af3(json_msa, work)
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
