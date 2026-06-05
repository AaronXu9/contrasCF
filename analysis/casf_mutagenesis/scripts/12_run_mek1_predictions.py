"""Run AF3+MSA, Boltz-1, and Boltz-2 predictions for MEK1 wt/rem/pack/inv.

Three phases (one MEK1 system, 4 variants each):
  1. Fetch one WT MSA via Boltz piggyback. Rewrite query row per variant.
  2. AF3+MSA: 4 jobs, ~95 s/job on RTX 4090 → ~6 min.
  3. Boltz-1 + Boltz-2: 4+4 = 8 jobs, ~30-40 s each → ~5 min.

Outputs land under contrasCF/data/{AF3,Boltz,Boltz2}/mek1_<variant>/ with the
canonical Boltz/AF3 file names expected by analysis/src/config.py::MODELS.

Pre-requisite: run 11_build_mek1.py first (writes the input JSON/YAML files).

Run:
    LD_LIBRARY_PATH=/home/aoxu/miniconda3/envs/rdkit_env/lib:$LD_LIBRARY_PATH \\
        /home/aoxu/miniconda3/envs/rdkit_env/bin/python \\
        analysis/casf_mutagenesis/scripts/12_run_mek1_predictions.py [phase]
where phase ∈ {all, msa, af3, boltz1, boltz2} (default: all).
"""
from __future__ import annotations
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
sys.path.insert(0, str(REPO_ROOT / "analysis"))

from casf_mutagenesis.config import (  # noqa: E402
    AF3_DIR, AF3_ENV, AF3_MODEL_DIR, BOLTZ_BIN as BOLTZ2_BIN_PATH, CUDA_DEVICE,
    OUTPUT_ROOT, assert_boltz2_binary,
)
from casf_mutagenesis.inputs_af3 import (  # noqa: E402
    clean_a3m_for_af3, render_af3, rewrite_a3m_query,
)
from casf_mutagenesis.msa_via_boltz import fetch_msa_via_boltz  # noqa: E402

DATA_ROOT = REPO_ROOT / "contrasCF" / "data"
VARIANTS = ("wt", "rem", "pack", "inv")

BOLTZ2_BIN = str(BOLTZ2_BIN_PATH)
BOLTZ1_BIN = "/home/aoxu/miniconda3/envs/rdkit_env/bin/boltz"  # v0.4.1

# AF3
AF3_NUM_DIFFUSION_SAMPLES = 5
AF3_NUM_RECYCLES = 10
AF3_TIMEOUT_S = 3600

# Boltz (both versions)
BOLTZ_DIFFUSION_SAMPLES = 5
BOLTZ_RECYCLING_STEPS = 3
BOLTZ_SAMPLING_STEPS = 200
BOLTZ_SEED = 42

MEK1_OUT = OUTPUT_ROOT / "7xlp_mek1"
MSA_CACHE = OUTPUT_ROOT / "_msa_cache"


# ============================================================================
# MSA helpers
# ============================================================================

def _read_af3_protein_sequence(af3_json_path: Path) -> str:
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


def phase_msa() -> str:
    """Fetch WT MEK1 MSA via Boltz piggyback. Return A3M text."""
    wt_json = MEK1_OUT / "wt" / "af3.json"
    wt_seq = _read_af3_protein_sequence(wt_json)
    print(f"[msa] WT sequence length: {len(wt_seq)}")
    t0 = time.time()
    a3m = fetch_msa_via_boltz(wt_seq, cache_dir=MSA_CACHE)
    n_seqs = sum(1 for ln in a3m.splitlines() if ln.startswith(">"))
    print(f"[msa] fetched WT MSA: n_seqs={n_seqs} in {time.time() - t0:.1f}s")
    return a3m


# ============================================================================
# AF3
# ============================================================================

def _build_af3_env() -> dict:
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = CUDA_DEVICE
    env["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
    xla_flag = "--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"
    existing_xla = env.get("XLA_FLAGS", "")
    env["XLA_FLAGS"] = (
        f"{existing_xla} {xla_flag}".strip() if xla_flag not in existing_xla
        else existing_xla
    )
    nvidia_lib_dirs: list[str] = []
    for sp in AF3_ENV.glob("lib/python*/site-packages/nvidia"):
        if sp.is_dir():
            for child in sp.iterdir():
                if (child / "lib").is_dir():
                    nvidia_lib_dirs.append(str(child / "lib"))
    if nvidia_lib_dirs:
        existing = env.get("LD_LIBRARY_PATH", "")
        env["LD_LIBRARY_PATH"] = ":".join(
            nvidia_lib_dirs + ([existing] if existing else [])
        )
    return env


def _run_af3_one(json_path: Path, work_dir: Path) -> Path:
    work_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(AF3_ENV / "bin" / "python"),
        str(AF3_DIR / "run_alphafold.py"),
        f"--json_path={json_path}",
        f"--output_dir={work_dir}",
        f"--model_dir={AF3_MODEL_DIR}",
        "--norun_data_pipeline",
        "--run_inference=true",
        f"--num_diffusion_samples={AF3_NUM_DIFFUSION_SAMPLES}",
        f"--num_recycles={AF3_NUM_RECYCLES}",
        "--flash_attention_implementation=xla",
    ]
    result = subprocess.run(
        cmd, capture_output=True, text=True,
        env=_build_af3_env(), timeout=AF3_TIMEOUT_S,
    )
    (work_dir / "af3_stdout.log").write_text(result.stdout)
    (work_dir / "af3_stderr.log").write_text(result.stderr)
    if result.returncode != 0:
        raise RuntimeError(
            f"AF3 failed (exit {result.returncode}); see {work_dir}/af3_stderr.log"
        )
    name = json.loads(json_path.read_text())["name"]
    sys_dir = work_dir / name
    if not sys_dir.exists():
        raise FileNotFoundError(f"AF3 output dir missing: {sys_dir}")
    return sys_dir


def _copy_af3_to_data(sys_dir: Path, name: str, dst_dir: Path) -> dict:
    """Copy AF3 outputs to contrasCF/data/AF3/<case>/ with canonical naming.

    The analysis pipeline reads `*_model_0.cif`; we copy ranked samples as
    `<name>_model_<rank>.cif` and `<name>_summary_confidences_<rank>.json`.
    Also writes a `_ranking_scores.csv` for traceability.
    """
    import csv
    ranking = sys_dir / f"{name}_ranking_scores.csv"
    if not ranking.exists():
        raise FileNotFoundError(f"missing ranking CSV: {ranking}")
    with ranking.open() as f:
        rows = list(csv.DictReader(f))
    rows.sort(key=lambda r: float(r["ranking_score"]), reverse=True)
    dst_dir.mkdir(parents=True, exist_ok=True)
    out = {"cifs": [], "confs": []}
    for rank, row in enumerate(rows):
        seed, sample = int(row["seed"]), int(row["sample"])
        sub = sys_dir / f"seed-{seed}_sample-{sample}"
        src_cif = sub / f"{name}_seed-{seed}_sample-{sample}_model.cif"
        src_conf = sub / f"{name}_seed-{seed}_sample-{sample}_summary_confidences.json"
        if not src_cif.exists():
            raise FileNotFoundError(f"missing model cif: {src_cif}")
        dst_cif = dst_dir / f"{name}_model_{rank}.cif"
        shutil.copy(src_cif, dst_cif)
        out["cifs"].append(str(dst_cif))
        if src_conf.exists():
            dst_conf = dst_dir / f"{name}_summary_confidences_{rank}.json"
            shutil.copy(src_conf, dst_conf)
            out["confs"].append(str(dst_conf))
    shutil.copy(ranking, dst_dir / f"{name}_ranking_scores.csv")
    return out


def phase_af3(wt_a3m: str) -> dict:
    print("\n=== Phase AF3+MSA ===")
    log: dict = {"runs": []}
    for variant in VARIANTS:
        v_dir = MEK1_OUT / variant
        af3_json_path = v_dir / "af3.json"
        prefix = f"mek1_{variant}"
        target_dst = DATA_ROOT / "AF3" / f"mek1_{variant}"
        # skip if already done
        if (target_dst / f"{prefix}_model_0.cif").exists():
            print(f"[af3] {variant}: skip (existing model_0.cif)")
            log["runs"].append({"variant": variant, "status": "skip_existing"})
            continue
        variant_seq = _read_af3_protein_sequence(af3_json_path)
        smiles = _read_ligand_smiles(af3_json_path)
        # MSA rewrite + clean
        variant_a3m = rewrite_a3m_query(wt_a3m, variant_seq)
        variant_a3m = clean_a3m_for_af3(variant_a3m, len(variant_seq))
        # Render fresh JSON with MSA
        json_msa = v_dir / "af3_msa.json"
        render_af3(
            name=prefix,
            chain_seqs=[("A", variant_seq)],
            ligand_smiles=smiles,
            out_path=json_msa,
            chain_msas={"A": variant_a3m},
        )
        # Run AF3
        t0 = time.time()
        work = v_dir / "_af3_work"
        if work.exists():
            shutil.rmtree(work)
        print(f"[af3] {variant}: running ...")
        try:
            sys_dir = _run_af3_one(json_msa, work)
            res = _copy_af3_to_data(sys_dir, prefix, target_dst)
            dt = time.time() - t0
            print(f"[af3]   ok ({dt:.0f}s) — {len(res['cifs'])} cifs copied")
            log["runs"].append({"variant": variant, "status": "ok",
                                "elapsed_s": dt, "n_samples": len(res["cifs"])})
        except Exception as exc:
            print(f"[af3]   FAIL: {type(exc).__name__}: {exc}")
            log["runs"].append({"variant": variant, "status": "error",
                                "error": f"{type(exc).__name__}: {exc}"})
    return log


# ============================================================================
# Boltz (both versions)
# ============================================================================

def _run_boltz_one(yaml_path: Path, prefix: str, work_dir: Path,
                   boltz_bin: str, model_kind: str) -> Path:
    """Run a single Boltz prediction.

    `model_kind` is "boltz1" or "boltz2"; only Boltz-2 takes `--model`
    (Boltz-1 doesn't have that flag).
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    yaml_renamed = work_dir / f"{prefix}.yaml"
    shutil.copy(yaml_path, yaml_renamed)
    cmd = [
        boltz_bin, "predict", str(yaml_renamed),
        "--out_dir", str(work_dir),
        "--output_format", "mmcif",
        "--diffusion_samples", str(BOLTZ_DIFFUSION_SAMPLES),
        "--recycling_steps", str(BOLTZ_RECYCLING_STEPS),
        "--sampling_steps", str(BOLTZ_SAMPLING_STEPS),
        "--seed", str(BOLTZ_SEED),
    ]
    if model_kind == "boltz2":
        cmd += ["--model", "boltz2"]
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = CUDA_DEVICE
    result = subprocess.run(cmd, capture_output=True, text=True, env=env,
                            timeout=AF3_TIMEOUT_S)
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


def _make_boltz_yaml_msa_empty(src_yaml: Path, dst_yaml: Path) -> None:
    """Boltz-1 doesn't ship a 2-arg `msa: empty` shorthand the same way; it
    accepts the same YAML so just copy."""
    shutil.copy(src_yaml, dst_yaml)


def _copy_boltz_to_data(pred_dir: Path, prefix: str, dst_dir: Path,
                        kind: str) -> dict:
    """Copy Boltz outputs to contrasCF/data/{Boltz,Boltz2}/<case>/.

    Canonical names expected by analysis/src/config.py::MODELS:
      structure: `*_model_0.cif`
      confidence: `confidence_*_model_0.json`
    """
    cifs = sorted(pred_dir.glob(f"{prefix}_model_*.cif"))
    if not cifs:
        raise FileNotFoundError(f"no model cifs in {pred_dir}")
    dst_dir.mkdir(parents=True, exist_ok=True)
    out = {"cifs": [], "confs": [], "affinity": None}
    for src_cif in cifs:
        dst_cif = dst_dir / src_cif.name
        shutil.copy(src_cif, dst_cif)
        out["cifs"].append(str(dst_cif))
        rank = src_cif.stem.rsplit("_", 1)[-1]
        src_conf = pred_dir / f"confidence_{prefix}_model_{rank}.json"
        if src_conf.exists():
            dst_conf = dst_dir / src_conf.name
            shutil.copy(src_conf, dst_conf)
            out["confs"].append(str(dst_conf))
    # Boltz-2 affinity sidecar
    aff = pred_dir / f"affinity_{prefix}.json"
    if aff.exists():
        dst_aff = dst_dir / aff.name
        shutil.copy(aff, dst_aff)
        out["affinity"] = str(dst_aff)
    return out


def phase_boltz(kind: str) -> dict:
    assert kind in ("boltz1", "boltz2")
    boltz_bin = BOLTZ1_BIN if kind == "boltz1" else BOLTZ2_BIN
    model_dir = "Boltz" if kind == "boltz1" else "Boltz2"
    print(f"\n=== Phase {kind.upper()} ({boltz_bin}) ===")
    log: dict = {"runs": []}
    for variant in VARIANTS:
        v_dir = MEK1_OUT / variant
        yaml_src = v_dir / "boltz.yaml"
        prefix = f"mek1_{variant}"
        target_dst = DATA_ROOT / model_dir / f"mek1_{variant}"
        if (target_dst / f"{prefix}_model_0.cif").exists():
            print(f"[{kind}] {variant}: skip (existing model_0.cif)")
            log["runs"].append({"variant": variant, "status": "skip_existing"})
            continue
        work = v_dir / f"_{kind}_work"
        if work.exists():
            shutil.rmtree(work)
        t0 = time.time()
        print(f"[{kind}] {variant}: running ...")
        try:
            pred_dir = _run_boltz_one(yaml_src, prefix, work, boltz_bin, kind)
            res = _copy_boltz_to_data(pred_dir, prefix, target_dst, kind)
            dt = time.time() - t0
            print(f"[{kind}]   ok ({dt:.0f}s) — {len(res['cifs'])} cifs copied")
            log["runs"].append({"variant": variant, "status": "ok",
                                "elapsed_s": dt, "n_samples": len(res["cifs"])})
        except Exception as exc:
            print(f"[{kind}]   FAIL: {type(exc).__name__}: {exc}")
            log["runs"].append({"variant": variant, "status": "error",
                                "error": f"{type(exc).__name__}: {exc}"})
    return log


# ============================================================================
# Main
# ============================================================================

def main() -> int:
    phase = sys.argv[1] if len(sys.argv) > 1 else "all"
    assert phase in ("all", "msa", "af3", "boltz1", "boltz2"), phase

    overall_log: dict = {"phases": {}}

    if phase in ("all", "msa", "af3"):
        wt_a3m = phase_msa()
    if phase in ("all", "af3"):
        # Boltz version assertion only matters for the MSA phase; do it here
        # to fail early before AF3.
        assert_boltz2_binary(BOLTZ2_BIN)
        overall_log["phases"]["af3"] = phase_af3(wt_a3m)
    if phase in ("all", "boltz1"):
        overall_log["phases"]["boltz1"] = phase_boltz("boltz1")
    if phase in ("all", "boltz2"):
        overall_log["phases"]["boltz2"] = phase_boltz("boltz2")

    log_path = MEK1_OUT / "predictions_run_log.json"
    log_path.write_text(json.dumps(overall_log, indent=2))
    print(f"\nlog: {log_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
