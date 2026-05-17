"""Build wt/rem/pack/inv inputs for the full CASF-2016 set (n=285).

Reads the 285 PDB ids from
  data/casf2016/labels/PDBbind_data_split_cleansplit.json["casf2016"]
and runs the same per-system pipeline as 01_build_subset20.py. Writes
manifest_full_casf.json.

Optional slicing via env vars (for SLURM job arrays):
  CONTRASCF_START — 0-indexed first PDB id (inclusive)
  CONTRASCF_END   — 0-indexed last PDB id (exclusive)
Slicing only affects which IDs are processed; outputs land in the same
shared `outputs/<pdbid>/<variant>/` tree, so concurrent array tasks
write disjoint subtrees with no contention.

Run:
    source env/lab.sh         # or env/carc.sh
    $CONTRASCF_PY analysis/casf_mutagenesis/scripts/02_build_full_casf.py
"""
from __future__ import annotations
import json
import os
import sys
from pathlib import Path

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
sys.path.insert(0, str(REPO_ROOT / "analysis"))

from casf_mutagenesis.build import build_system, write_manifest
from casf_mutagenesis.config import OUTPUT_ROOT, SPLIT_JSON


def main() -> int:
    ids = json.loads(SPLIT_JSON.read_text())["casf2016"]
    start = int(os.environ.get("CONTRASCF_START", "0"))
    end = int(os.environ.get("CONTRASCF_END", str(len(ids))))
    ids_slice = ids[start:end]
    print(f"full CASF-2016: total={len(ids)} processing [{start}:{end}] = {len(ids_slice)} PDBs")

    entries = []
    for i, pdbid in enumerate(ids_slice, 1):
        entry = build_system(pdbid)
        n_pocket = len(entry.pocket)
        n_warn = len(entry.warnings)
        tag = "OK" if entry.status == "ok" and n_warn == 0 else (
            entry.status.upper() if entry.status != "ok" else "WARN"
        )
        print(
            f"  [{i:3d}/{len(ids_slice)}] {pdbid:6s} {tag:5s} "
            f"chains={entry.n_chains} pocket={n_pocket} "
            f"het={entry.het_code} warnings={n_warn}",
            flush=True,
        )
        if entry.status != "ok":
            print(f"      error: {entry.error}")
        for w in entry.warnings[:3]:
            print(f"      warn: {w}")
        entries.append(entry)

    # Per-slice manifest so array tasks don't stomp each other; merge later
    # if needed.
    suffix = f"_{start}_{end}" if (start, end) != (0, len(ids)) else ""
    manifest = OUTPUT_ROOT / f"manifest_full_casf{suffix}.json"
    write_manifest(entries, manifest)
    print(f"\nManifest written: {manifest}")

    n_ok = sum(1 for e in entries if e.status == "ok")
    n_err = sum(1 for e in entries if e.status == "error")
    n_warn = sum(1 for e in entries if e.status == "ok" and e.warnings)
    print(f"Summary: ok={n_ok} warn={n_warn} error={n_err}")
    # Per-system errors (missing raw/ entries, malformed crystal SDF, etc.)
    # are data issues, not pipeline failures — downstream Boltz / AF3
    # phases skip those cells via missing_input. Always exit 0 unless ALL
    # systems failed (which means something larger is wrong: bad paths,
    # missing dependencies, etc.).
    if entries and n_err == len(entries):
        print("ERROR: all systems failed — likely a configuration problem")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
