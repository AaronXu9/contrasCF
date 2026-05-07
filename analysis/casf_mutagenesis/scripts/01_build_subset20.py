"""Smoke test on the 20-complex CASF-2016 subset.

For each of the 20 PDB ids in
  data/casf2016/labels/PDBbind_casf2016_subset20.json
we build wt/rem/pack/inv variants × {AF3 JSON, Boltz YAML}, plus WT docking
inputs. Manifest written to outputs/manifest_subset20.json.

Run:
    LD_LIBRARY_PATH=/home/aoxu/miniconda3/envs/rdkit_env/lib:$LD_LIBRARY_PATH \\
        /home/aoxu/miniconda3/envs/rdkit_env/bin/python \\
        analysis/casf_mutagenesis/scripts/01_build_subset20.py
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

REPO_ROOT = Path("/mnt/katritch_lab2/aoxu/contrasCF")
sys.path.insert(0, str(REPO_ROOT / "analysis"))

from casf_mutagenesis.build import build_system, write_manifest
from casf_mutagenesis.config import OUTPUT_ROOT, SUBSET20_JSON


def main() -> int:
    ids = json.loads(SUBSET20_JSON.read_text())["casf2016"]
    print(f"subset20 size: {len(ids)}")
    entries = []
    for i, pdbid in enumerate(ids, 1):
        entry = build_system(pdbid)
        n_pocket = len(entry.pocket)
        n_warn = len(entry.warnings)
        tag = "OK" if entry.status == "ok" and n_warn == 0 else (
            entry.status.upper() if entry.status != "ok" else "WARN"
        )
        print(
            f"  [{i:2d}/{len(ids)}] {pdbid:6s} {tag:5s} "
            f"chains={entry.n_chains} pocket={n_pocket} "
            f"het={entry.het_code} warnings={n_warn}"
        )
        if entry.status != "ok":
            print(f"      error: {entry.error}")
        for w in entry.warnings[:3]:
            print(f"      warn: {w}")
        entries.append(entry)

    manifest = OUTPUT_ROOT / "manifest_subset20.json"
    write_manifest(entries, manifest)
    print(f"\nManifest written: {manifest}")

    # summary
    n_ok = sum(1 for e in entries if e.status == "ok")
    n_err = sum(1 for e in entries if e.status == "error")
    n_warn = sum(1 for e in entries if e.status == "ok" and e.warnings)
    print(f"Summary: ok={n_ok} warn={n_warn} error={n_err}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
