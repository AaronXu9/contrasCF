"""Run GNINA on the WT variant of every system with docking/ inputs.

For each <pdbid>/wt/docking/{receptor.pdb, ligand.sdf, box.json}, invoke
GNINA and save:
    <pdbid>/wt/gnina/{poses.sdf, run.log}
Top-1 score (CNN affinity, Vina-like score) is parsed from the SDF tags.

Mutant variants (rem/pack/inv) are NOT handled here — they need a
mutant protein structure (from AF3+MSA mutant CIFs, currently being
generated on CARC). Once those land on lab, a follow-on script will
strip them to receptor.pdb and run GNINA against them.

GNINA defaults match analysis/scripts/11_run_docking.py:
  --cnn_scoring rescore --num_modes 9 --exhaustiveness 8 --seed 42

Scope (env-var controlled, defaults to subset20):
  CONTRASCF_SCOPE=subset20 | full

Run:
    source env/lab.sh
    $CONTRASCF_PY analysis/casf_mutagenesis/scripts/07_run_gnina_wt.py
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

from casf_mutagenesis.config import (  # noqa: E402
    OUTPUT_ROOT, SPLIT_JSON, SUBSET20_JSON,
)

GNINA_BIN = os.environ.get(
    "CONTRASCF_GNINA_BIN",
    "/home/aoxu/projects/PoseBench/forks/GNINA/gnina",
)
GNINA_NUM_MODES = 9
GNINA_EXHAUSTIVENESS = 8
GNINA_SEED = 42
GNINA_CNN = "rescore"
TIMEOUT_S = 1800   # 30 min per system; bigger systems can time out


def _resolve_ids() -> tuple[list[str], str]:
    scope = os.environ.get("CONTRASCF_SCOPE", "subset20")
    if scope == "full":
        ids = json.loads(SPLIT_JSON.read_text())["casf2016"]
    else:
        ids = json.loads(SUBSET20_JSON.read_text())["casf2016"]
    start = int(os.environ.get("CONTRASCF_START", "0"))
    end = int(os.environ.get("CONTRASCF_END", str(len(ids))))
    return ids[start:end], f"{scope}[{start}:{end}]"


def run_gnina_one(pdbid: str) -> dict:
    """Run GNINA on <pdbid>/wt/. Returns a dict for the run log."""
    entry = {"pdbid": pdbid, "variant": "wt"}
    v_dir = OUTPUT_ROOT / pdbid / "wt"
    inputs = v_dir / "docking"
    receptor = inputs / "receptor.pdb"
    ligand = inputs / "ligand.sdf"
    box_json = inputs / "box.json"
    if not all(p.exists() for p in (receptor, ligand, box_json)):
        entry["status"] = "missing_input"
        return entry

    out_dir = v_dir / "gnina"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_sdf = out_dir / "poses.sdf"
    if out_sdf.exists():
        entry["status"] = "skip_existing"
        return entry

    box = json.loads(box_json.read_text())
    cx, cy, cz = box["center"]
    sx, sy, sz = box["size"]
    cmd = [
        GNINA_BIN,
        "-r", str(receptor),
        "-l", str(ligand),
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
    ids, scope_label = _resolve_ids()
    print(f"GNINA WT scope: {scope_label}, n_pdb={len(ids)}")
    print(f"GNINA bin: {GNINA_BIN}")
    log_path = OUTPUT_ROOT / "gnina_wt_run_log.json"
    # Resume from existing log if any
    runs: list[dict] = []
    if log_path.exists():
        try:
            runs = json.loads(log_path.read_text()).get("runs", [])
        except Exception:
            runs = []
    seen = {(r["pdbid"], r.get("variant", "wt")) for r in runs}

    n_total = len(ids)
    n_done = 0
    n_skip = 0
    n_fail = 0
    for pdbid in ids:
        n_done += 1
        if (pdbid, "wt") in seen:
            continue
        print(f"  [{n_done:3d}/{n_total}] {pdbid} ...", flush=True)
        entry = run_gnina_one(pdbid)
        runs.append(entry)
        status = entry.get("status", "?")
        wc = entry.get("wallclock_s", "?")
        print(f"      {status} (wallclock={wc}s)")
        if status == "ok":
            pass
        elif status in ("skip_existing", "missing_input"):
            n_skip += 1
        else:
            n_fail += 1
        log_path.write_text(json.dumps({"runs": runs}, indent=2))

    n_ok = sum(1 for r in runs if r.get("status") == "ok")
    print(f"\nTotal {n_total} systems | new this run: ok={n_ok-len(seen)+n_skip} "
          f"skip/missing={n_skip} fail={n_fail}")
    print(f"Cumulative ok={n_ok} in log {log_path}")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
