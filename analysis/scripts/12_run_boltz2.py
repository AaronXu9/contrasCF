"""Run Boltz-2 on the 16 prepared cases.

For each case:
  1. Build a Boltz-2 input YAML from
     ``docking/inputs/<case>/{receptor.pdb, ligand.sdf}`` plus the per-case
     cofactor / metal flags in ``analysis/src/config.py``:
         protein:  chains from receptor.pdb CA atoms, msa=empty (no-MSA mode)
         ligand:   target SMILES from ligand.sdf
         ligand:   ccd: NAP   (if has_nadp_input)
         ligand:   ccd: MG    (if has_mg_input)
         ligand:   ccd: ZN    (if has_zn_input)
         properties:
           - affinity: {binder: <target ligand chain>}
  2. Invoke ``boltz predict`` from the boltzina_env (Boltz 2.2.1).
  3. Copy the top-ranked mmCIF + confidence JSON + affinity JSON to
     ``contrasCF/data/Boltz2/<case>/`` with filenames matching the Boltz-1
     convention (``<prefix>_model_0.cif``, ``confidence_<prefix>_model_0.json``).

Run:
    /home/aoxu/miniconda3/envs/rdkit_env/bin/python \
        analysis/scripts/12_run_boltz2.py
"""
from __future__ import annotations
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import yaml

REPO_ROOT = Path("/mnt/katritch_lab2/aoxu/contrasCF")
sys.path.insert(0, str(REPO_ROOT / "analysis" / "src"))

from config import CASES, DATA_ROOT  # noqa: E402

# Inlined from dockstrat.utils.sequence to avoid pulling in dockstrat's
# engine.py import chain (needs rootutils, not installed in rdkit_env).
AA_3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "HSD": "H", "HSE": "H", "HSP": "H", "HIE": "H", "HID": "H",
    "MSE": "M", "SEC": "U",
}


def extract_protein_sequence(pdb_path: str) -> list[str]:
    """Return one amino-acid sequence string per protein chain in the PDB."""
    import gemmi
    st = gemmi.read_structure(pdb_path)
    model = st[0]
    sequences: list[str] = []
    for chain in model:
        seen_resid: set[tuple[int, str]] = set()
        letters: list[str] = []
        for r in chain:
            if r.name not in AA_3TO1:
                continue
            key = (r.seqid.num, r.seqid.icode or "")
            if key in seen_resid:
                continue
            if not any(a.name == "CA" for a in r):
                continue
            seen_resid.add(key)
            letters.append(AA_3TO1[r.name])
        if letters:
            sequences.append("".join(letters))
    return sequences

INPUTS_ROOT = REPO_ROOT / "docking" / "inputs"
BOLTZ2_OUT_ROOT = DATA_ROOT / "Boltz2"

BOLTZ_BIN = "/home/aoxu/miniconda3/envs/boltzina_env/bin/boltz"

# Boltz-2 inference knobs (match dockstrat_config/model/boltz2_inference.yaml).
DIFFUSION_SAMPLES = 5
RECYCLING_STEPS = 3
SAMPLING_STEPS = 200
SEED = 42
CUDA_DEVICE = "0"


def _smiles_from_sdf(sdf_path: Path) -> str:
    from rdkit import Chem
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=True)
    mol = next(iter(suppl), None)
    if mol is None:
        raise ValueError(f"cannot parse molecule from {sdf_path}")
    return Chem.MolToSmiles(mol)


def _build_yaml(case: str) -> tuple[dict, str]:
    """Return (yaml_dict, ligand_chain_id) for the case."""
    spec = CASES[case]
    receptor = INPUTS_ROOT / case / "receptor.pdb"
    ligand = INPUTS_ROOT / case / "ligand.sdf"
    sequences = extract_protein_sequence(str(receptor))
    if not sequences:
        raise RuntimeError(f"no protein chains found in {receptor}")
    smiles = _smiles_from_sdf(ligand)

    entries: list[dict] = []
    # Protein chain(s): A, B, ...
    for i, seq in enumerate(sequences):
        entries.append({
            "protein": {
                "id": chr(ord("A") + i),
                "sequence": seq,
                "msa": "empty",
            }
        })

    # Target ligand (the binder).
    next_id = chr(ord("A") + len(sequences))
    binder_chain = next_id
    entries.append({"ligand": {"id": binder_chain, "smiles": smiles}})

    # Cofactors / metals — appended as CCD ligands, one chain id each.
    next_ord = ord(next_id) + 1
    if spec.has_nadp_input:
        entries.append({"ligand": {"id": chr(next_ord), "ccd": "NAP"}}); next_ord += 1
    if spec.has_mg_input:
        entries.append({"ligand": {"id": chr(next_ord), "ccd": "MG"}}); next_ord += 1
    if spec.has_zn_input:
        entries.append({"ligand": {"id": chr(next_ord), "ccd": "ZN"}}); next_ord += 1

    return (
        {
            "version": 1,
            "sequences": entries,
            "properties": [{"affinity": {"binder": binder_chain}}],
        },
        binder_chain,
    )


def _run_boltz(yaml_path: Path, out_dir: Path) -> None:
    cmd = [
        BOLTZ_BIN, "predict", str(yaml_path),
        "--out_dir", str(out_dir),
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
    # Always write log regardless of success.
    (out_dir / "boltz_stdout.log").write_text(result.stdout)
    (out_dir / "boltz_stderr.log").write_text(result.stderr)
    if result.returncode != 0:
        raise RuntimeError(
            f"boltz predict failed (exit {result.returncode}); see {out_dir}/boltz_stderr.log"
        )


def _collect_outputs(prefix: str, run_dir: Path, final_dir: Path) -> dict:
    """Copy the top-ranked Boltz-2 outputs into ``final_dir`` using Boltz-1 naming."""
    # Boltz writes: <run_dir>/boltz_results_<prefix>/predictions/<prefix>/{*_model_*.cif, confidence_*_model_*.json, affinity_*.json}
    pred_dir = run_dir / f"boltz_results_{prefix}" / "predictions" / prefix
    if not pred_dir.exists():
        raise FileNotFoundError(f"Boltz-2 predictions directory missing: {pred_dir}")

    # Top model = model_0 (Boltz writes all samples; model_0 is highest-ranked).
    src_cif = pred_dir / f"{prefix}_model_0.cif"
    src_conf = pred_dir / f"confidence_{prefix}_model_0.json"
    if not src_cif.exists():
        # Fallback: pick any model file.
        cifs = sorted(pred_dir.glob(f"{prefix}_model_*.cif"))
        if not cifs:
            raise FileNotFoundError(f"no predicted mmCIF in {pred_dir}")
        src_cif = cifs[0]
        src_conf = pred_dir / f"confidence_{src_cif.stem}.json"

    final_dir.mkdir(parents=True, exist_ok=True)
    dst_cif = final_dir / f"{prefix}_model_0.cif"
    shutil.copy(src_cif, dst_cif)

    dst_conf = None
    if src_conf.exists():
        dst_conf = final_dir / f"confidence_{prefix}_model_0.json"
        shutil.copy(src_conf, dst_conf)

    dst_aff = None
    aff_src = pred_dir / f"affinity_{prefix}.json"
    if aff_src.exists():
        dst_aff = final_dir / f"affinity_{prefix}.json"
        shutil.copy(aff_src, dst_aff)

    return {"cif": str(dst_cif), "conf": str(dst_conf), "affinity": str(dst_aff)}


def run_case(case: str) -> dict:
    spec = CASES[case]
    out_dir = BOLTZ2_OUT_ROOT / case
    out_dir.mkdir(parents=True, exist_ok=True)

    # Use the case name as the prefix — Boltz derives input stem from the YAML filename.
    prefix = case
    yaml_path = out_dir / f"{prefix}.yaml"
    run_dir = out_dir / "run"

    # Idempotency: skip if top-ranked CIF already exists.
    if (out_dir / f"{prefix}_model_0.cif").exists():
        return {"case": case, "ok": True, "sec": 0.0, "skipped": True}

    t0 = time.time()
    try:
        yaml_dict, binder_chain = _build_yaml(case)
        with open(yaml_path, "w") as fh:
            yaml.safe_dump(yaml_dict, fh, default_flow_style=False, sort_keys=False)

        # Boltz writes to <run_dir>/boltz_results_<prefix>/...
        run_dir.mkdir(parents=True, exist_ok=True)
        _run_boltz(yaml_path, run_dir)

        copied = _collect_outputs(prefix, run_dir, out_dir)
        return {
            "case": case,
            "ok": True,
            "binder_chain": binder_chain,
            "sec": time.time() - t0,
            **copied,
        }
    except Exception as e:
        return {
            "case": case,
            "ok": False,
            "error": f"{type(e).__name__}: {e}",
            "sec": time.time() - t0,
        }


def main() -> None:
    BOLTZ2_OUT_ROOT.mkdir(parents=True, exist_ok=True)
    results = []
    for case in CASES:
        print(f"[boltz2] {case} ...", flush=True)
        info = run_case(case)
        tag = "ok" if info.get("ok") else "FAIL"
        print(
            f"[boltz2]   {tag} ({info.get('sec', 0):.1f}s)"
            + (f"  err={info['error']}" if not info.get("ok") else "")
            + ("  [skipped — already done]" if info.get("skipped") else "")
        )
        results.append(info)

    (DATA_ROOT / "boltz2_runs.json").write_text(json.dumps(results, indent=2))
    ok = sum(1 for r in results if r.get("ok"))
    print(f"[boltz2] done. ok={ok}/{len(results)}")


if __name__ == "__main__":
    main()
