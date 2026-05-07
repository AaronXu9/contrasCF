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
        analysis/scripts/13_run_surfdock.py
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

REPO_ROOT = Path("/mnt/katritch_lab2/aoxu/contrasCF")
DOCKSTRAT_ROOT = Path("/mnt/katritch_lab2/aoxu/CogLigandBench")
SURFDOCK_DIR = "/home/aoxu/projects/SurfDock"
SURFDOCK_PRECOMPUTED_ARRAYS = "/home/aoxu/projects/precomputed/precomputed_arrays"

sys.path.insert(0, str(REPO_ROOT / "analysis" / "src"))

from config import CASES, DATA_ROOT  # noqa: E402
from docking_io import write_combined_top1  # noqa: E402


def _load_surfdock_module():
    """Import dockstrat.models.surfdock_inference without triggering
    dockstrat.__init__.py (which needs rootutils, not installed in rdkit_env).
    """
    import importlib.util
    import types

    # Register stub parent packages so the loaded module's dotted name resolves.
    for name, path in [
        ("dockstrat", DOCKSTRAT_ROOT / "dockstrat"),
        ("dockstrat.models", DOCKSTRAT_ROOT / "dockstrat" / "models"),
    ]:
        if name not in sys.modules:
            m = types.ModuleType(name)
            m.__path__ = [str(path)]
            sys.modules[name] = m

    sf_path = DOCKSTRAT_ROOT / "dockstrat" / "models" / "surfdock_inference.py"
    spec = importlib.util.spec_from_file_location(
        "dockstrat.models.surfdock_inference", str(sf_path)
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules["dockstrat.models.surfdock_inference"] = mod
    spec.loader.exec_module(mod)
    return mod

INPUTS_ROOT = REPO_ROOT / "docking" / "inputs"
SURFDOCK_OUT_ROOT = DATA_ROOT / "SurfDock"

SAMPLES_PER_COMPLEX = 40
NUM_POSES = 10
BATCH_SIZE = 40


def _surfdock_config() -> dict:
    return {
        "surfdock_env": "SurfDock",
        "surfdock_dir": SURFDOCK_DIR,
        "precomputed_arrays": SURFDOCK_PRECOMPUTED_ARRAYS,
        "diffusion_model_dir": str(DOCKSTRAT_ROOT / "forks" / "SurfDock" / "model_weights" / "docking"),
        "confidence_model_dir": str(DOCKSTRAT_ROOT / "forks" / "SurfDock" / "model_weights" / "posepredict"),
        "num_gpus": 1,
        "main_process_port": 29510,
        "batch_size": BATCH_SIZE,
        "samples_per_complex": SAMPLES_PER_COMPLEX,
        "num_poses": NUM_POSES,
        "mdn_dist_threshold": 3,
        "skip_existing": True,
    }


_SDF_FNAME_RE = re.compile(
    r"_sample_idx_(?P<sample>\d+)_rank_(?P<rank>\d+)"
    r"_rmsd_(?P<rmsd>[\d.-]+)_confidence_(?P<conf>-?[\d.]+)\.sdf$"
)


def _parse_surfdock_score(name: str) -> dict:
    m = _SDF_FNAME_RE.search(name)
    if not m:
        return {}
    return {
        "SurfDock_sample_idx": m.group("sample"),
        "SurfDock_rank": m.group("rank"),
        "SurfDock_rmsd": m.group("rmsd"),
        "SurfDock_confidence": m.group("conf"),
    }


def _sdf_record_with_tags(sdf_text: str, tags: dict) -> str:
    """Insert SDF tag block before $$$$ (which terminates a record)."""
    if "$$$$" not in sdf_text:
        parts = [sdf_text.rstrip()]
        for k, v in tags.items():
            parts.append(f">  <{k}>")
            parts.append(str(v))
            parts.append("")
        parts.append("$$$$")
        return "\n".join(parts) + "\n"
    head, _, rest = sdf_text.partition("$$$$")
    head = head.rstrip()
    tag_lines: list[str] = []
    for k, v in tags.items():
        tag_lines.append(f">  <{k}>")
        tag_lines.append(str(v))
        tag_lines.append("")
    return head + "\n" + "\n".join(tag_lines) + "$$$$" + rest


def _translate_ligand_to_pocket(ligand_sdf: Path, pocket_center: tuple[float, float, float],
                                out_sdf: Path) -> None:
    """Translate the ligand's 3D conformer so its centroid sits at pocket_center.

    SurfDock uses the ligand.sdf coordinates to crop an 8 Å pocket around the
    ligand. docking/inputs/<case>/ligand.sdf has RDKit-generated coordinates
    centered near the origin, which puts the SurfDock pocket outside the
    receptor's binding site. This helper shifts the ligand into the binding
    pocket (determined from the docking box.json, which in turn comes from
    the AF3-predicted ligand centroid).
    """
    from rdkit import Chem
    import numpy as np
    suppl = Chem.SDMolSupplier(str(ligand_sdf), removeHs=False)
    mol = next(iter(suppl), None)
    if mol is None:
        raise ValueError(f"cannot read {ligand_sdf}")
    conf = mol.GetConformer()
    xyz = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z]
                    for i in range(mol.GetNumAtoms())])
    centroid = xyz.mean(axis=0)
    shift = np.array(pocket_center) - centroid
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (p.x + shift[0], p.y + shift[1], p.z + shift[2]))
    writer = Chem.SDWriter(str(out_sdf))
    writer.write(mol)
    writer.close()


def _run_inference_with_pocket_center(csv_path: str, esm_pt: str, out_dir: str, config: dict) -> str:
    """Local SurfDock inference launcher, mirroring dockstrat's _run_inference but
    with --ligand_to_pocket_center=True appended.

    SurfDock's diffusion sampler uses randomize_position(); without
    ligand_to_pocket_center=True it adds a Normal(0, tr_sigma_max) translation
    and lets the model learn translations from there. For OOD adversarial
    ligands (e.g. propyl-ATP, penta-methyl-glucose) the diffusion can produce
    pose centroids thousands of Angstroms from the pocket. Anchoring the
    initial position to the predicted pocket center prevents that runaway.
    """
    import subprocess

    surfdock_dir = config["surfdock_dir"]
    os.makedirs(out_dir, exist_ok=True)

    env = os.environ.copy()
    env["precomputed_arrays"] = config.get("precomputed_arrays",
        os.path.join(os.path.dirname(surfdock_dir), "precomputed", "precomputed_arrays"))

    accelerate = f"/home/aoxu/miniconda3/envs/{config.get('surfdock_env', 'SurfDock')}/bin/accelerate"
    cmd = [
        accelerate, "launch",
        "--num_processes", str(config.get("num_gpus", 1)),
        "--main_process_port", str(config.get("main_process_port", 29510)),
        os.path.join(surfdock_dir, "inference_accelerate.py"),
        "--data_csv", csv_path,
        "--model_dir", config["diffusion_model_dir"],
        "--ckpt", "best_ema_inference_epoch_model.pt",
        "--confidence_model_dir", config["confidence_model_dir"],
        "--confidence_ckpt", "best_model.pt",
        "--save_docking_result",
        "--mdn_dist_threshold_test", str(config.get("mdn_dist_threshold", 3)),
        "--esm_embeddings_path", esm_pt,
        "--run_name", "dockstrat_run",
        "--project", "surfdock_inference",
        "--out_dir", out_dir,
        "--batch_size", str(config.get("batch_size", 40)),
        "--batch_size_molecule", "1",
        "--samples_per_complex", str(config.get("samples_per_complex", 40)),
        "--save_docking_result_number", str(config.get("num_poses", 40)),
        "--head_index", "0",
        "--tail_index", "10000",
        "--inference_mode", "evaluate",
        "--wandb_dir", os.path.join(out_dir, "wandb"),
        "--ligand_to_pocket_center",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=surfdock_dir, env=env)
    if result.returncode != 0:
        raise RuntimeError(
            f"SurfDock inference failed (exit {result.returncode}):\n"
            f"  stderr tail: {result.stderr[-2000:]}"
        )
    return os.path.join(out_dir, "SurfDock_docking_result")


def _run_surfdock_pipeline(case: str, receptor: Path, ligand: Path, work_dir: Path) -> list[Path]:
    """Run SurfDock's 4 steps with a persistent work_dir.

    Returns the list of ranked source SDFs (with SurfDock's original filename
    encoding score info). Highest confidence (rank_1) is first.
    """
    sf = _load_surfdock_module()
    _compute_surface = sf._compute_surface
    _build_input_csv = sf._build_input_csv
    _compute_esm_embeddings = sf._compute_esm_embeddings
    cfg = _surfdock_config()
    work_dir.mkdir(parents=True, exist_ok=True)

    data_dir = work_dir / "data"
    sysdir = data_dir / case
    sysdir.mkdir(parents=True, exist_ok=True)
    shutil.copy(receptor, sysdir / f"{case}_protein_processed.pdb")

    # Position ligand at the docking-box center so SurfDock's 8Å pocket crop
    # captures the actual binding site.
    box_json = INPUTS_ROOT / case / "box.json"
    box = json.loads(box_json.read_text())
    pocket_center = tuple(box["center"])
    _translate_ligand_to_pocket(ligand, pocket_center, sysdir / f"{case}_ligand.sdf")

    surface_dir = work_dir / "surface"
    csv_dir = work_dir / "csv"; csv_dir.mkdir(parents=True, exist_ok=True)
    csv_path = csv_dir / "input.csv"
    esm_dir = work_dir / "esm"
    infer_out = work_dir / "inference_out"

    _compute_surface(str(data_dir), str(surface_dir), cfg)
    _build_input_csv(str(data_dir), str(surface_dir), str(csv_path), cfg)

    # Inject pocket_center column (read from box.json) so inference_accelerate.py
    # populates receptor['pocket_center'] for --ligand_to_pocket_center.
    import csv as csv_mod
    cx, cy, cz = pocket_center
    with open(csv_path) as fh:
        rows = list(csv_mod.DictReader(fh))
    if rows and "pocket_center" not in rows[0]:
        for r in rows:
            r["pocket_center"] = f"{cx},{cy},{cz}"
        with open(csv_path, "w", newline="") as fh:
            writer = csv_mod.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    esm_pt = _compute_esm_embeddings(str(csv_path), str(esm_dir), cfg)
    docking_result_dir = _run_inference_with_pocket_center(str(csv_path), esm_pt, str(infer_out), cfg)

    pocket_stem = f"{case}_protein_processed_8A"
    candidates = glob.glob(os.path.join(docking_result_dir, f"{pocket_stem}_*"))
    if not candidates:
        candidates = [d for d in glob.glob(os.path.join(docking_result_dir, "*"))
                      if pocket_stem in os.path.basename(d)]
    if not candidates:
        raise RuntimeError(f"SurfDock produced no output dir under {docking_result_dir}")
    subdir = candidates[0]
    sdf_files = glob.glob(os.path.join(subdir, "*.sdf"))

    def _rank(f: str) -> int:
        m = re.search(r"_rank_(\d+)_", os.path.basename(f))
        return int(m.group(1)) if m else 9999

    sdf_files.sort(key=_rank)
    return [Path(p) for p in sdf_files[:NUM_POSES]]


def _collect_poses_sdf(source_sdfs: list[Path], out_dir: Path) -> Path:
    """Copy rank SDFs and build a single poses.sdf with score tags embedded."""
    chunks: list[str] = []
    for i, src in enumerate(source_sdfs, start=1):
        text = src.read_text()
        tags = _parse_surfdock_score(src.name)
        tags["SurfDock_rank_out"] = str(i)
        chunks.append(_sdf_record_with_tags(text, tags))
        shutil.copy(src, out_dir / f"rank{i}.sdf")
    poses_sdf = out_dir / "poses.sdf"
    poses_sdf.write_text("".join(chunks))
    return poses_sdf


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
    os.environ["PROJECT_ROOT"] = str(DOCKSTRAT_ROOT)
    os.environ["SURFDOCK_DIR"] = SURFDOCK_DIR
    os.environ["SURFDOCK_PRECOMPUTED_ARRAYS"] = SURFDOCK_PRECOMPUTED_ARRAYS
    os.environ["precomputed_arrays"] = SURFDOCK_PRECOMPUTED_ARRAYS

    SURFDOCK_OUT_ROOT.mkdir(parents=True, exist_ok=True)
    results = []
    for case in CASES:
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
