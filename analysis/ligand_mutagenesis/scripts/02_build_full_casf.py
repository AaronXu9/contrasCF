"""Build the full CASF-2016 set (~285 systems) of ligand-side variants.

Iterates the PDB-id list from `casf2016` in
`PDBbind_data_split_cleansplit.json`. For chunked re-runs (the run can take
~30 min and you may hit a transient SDF / sanitize failure on 1-2 systems),
use `--start` and `--limit` to restart from a given index.

Run:
    LD_LIBRARY_PATH=/home/aoxu/miniconda3/envs/rdkit_env/lib:$LD_LIBRARY_PATH \\
        /home/aoxu/miniconda3/envs/rdkit_env/bin/python \\
        analysis/ligand_mutagenesis/scripts/02_build_full_casf.py [--limit N] [--start N]
"""
from __future__ import annotations
import argparse
import json
import os
import sys
from pathlib import Path

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
sys.path.insert(0, str(REPO_ROOT / "analysis"))

from ligand_mutagenesis.build import build_system, write_manifest
from ligand_mutagenesis.config import OUTPUT_ROOT, SPLIT_JSON


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--start", type=int, default=0, help="0-based index to start from")
    p.add_argument("--limit", type=int, default=None, help="max systems to build")
    p.add_argument("--manifest", type=str, default=None,
                   help="output manifest path (default: outputs/manifest_full.json)")
    args = p.parse_args()

    ids: list[str] = json.loads(SPLIT_JSON.read_text())["casf2016"]
    total = len(ids)
    end = total if args.limit is None else min(total, args.start + args.limit)
    chunk = ids[args.start:end]
    print(f"casf2016 size: {total} — building [{args.start}:{end}] = {len(chunk)} systems")

    entries = []
    for i, pdbid in enumerate(chunk, args.start + 1):
        entry = build_system(pdbid)
        n_var = len(entry.variants)
        n_warn = len(entry.warnings)
        elig = ",".join(
            r for r, v in entry.rules_eligible.items() if v["eligible"]
        ) or "none"
        tag = "OK" if entry.status == "ok" and n_warn == 0 else (
            entry.status.upper() if entry.status != "ok" else "WARN"
        )
        print(
            f"  [{i:3d}/{total}] {pdbid:6s} {tag:5s} "
            f"chains={entry.n_chains} variants={n_var} "
            f"eligible={elig} warnings={n_warn}"
        )
        if entry.status != "ok":
            print(f"      error: {entry.error}")
        entries.append(entry)

    manifest_path = Path(args.manifest) if args.manifest else (OUTPUT_ROOT / "manifest_full.json")
    write_manifest(entries, manifest_path)
    print(f"\nManifest written: {manifest_path}")

    n_ok = sum(1 for e in entries if e.status == "ok")
    n_err = sum(1 for e in entries if e.status == "error")
    n_warn = sum(1 for e in entries if e.status == "ok" and e.warnings)
    total_variants = sum(len(e.variants) for e in entries)
    print(f"Summary: ok={n_ok} warn={n_warn} error={n_err} "
          f"total_variants={total_variants}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
