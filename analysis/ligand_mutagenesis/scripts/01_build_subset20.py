"""Smoke test on the 20-complex CASF-2016 subset.

For each PDB id in `data/casf2016/labels/PDBbind_casf2016_subset20.json`
we build wt + every chemistry-rule variant × {AF3 JSON, Boltz YAML, docking
inputs}. Manifest written to outputs/manifest_subset20.json.

Run:
    LD_LIBRARY_PATH=/home/aoxu/miniconda3/envs/rdkit_env/lib:$LD_LIBRARY_PATH \\
        /home/aoxu/miniconda3/envs/rdkit_env/bin/python \\
        analysis/ligand_mutagenesis/scripts/01_build_subset20.py
"""
from __future__ import annotations
import json
import os
import sys
from pathlib import Path

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
sys.path.insert(0, str(REPO_ROOT / "analysis"))

from ligand_mutagenesis.build import build_system, write_manifest
from ligand_mutagenesis.config import OUTPUT_ROOT, SUBSET20_JSON


def main() -> int:
    ids = json.loads(SUBSET20_JSON.read_text())["casf2016"]
    print(f"subset20 size: {len(ids)}")
    entries = []
    for i, pdbid in enumerate(ids, 1):
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
            f"  [{i:2d}/{len(ids)}] {pdbid:6s} {tag:5s} "
            f"chains={entry.n_chains} variants={n_var} "
            f"eligible={elig} warnings={n_warn}"
        )
        if entry.status != "ok":
            print(f"      error: {entry.error}")
        for w in entry.warnings[:3]:
            print(f"      warn: {w}")
        entries.append(entry)

    manifest = OUTPUT_ROOT / "manifest_subset20.json"
    write_manifest(entries, manifest)
    print(f"\nManifest written: {manifest}")

    n_ok = sum(1 for e in entries if e.status == "ok")
    n_err = sum(1 for e in entries if e.status == "error")
    n_warn = sum(1 for e in entries if e.status == "ok" and e.warnings)
    total_variants = sum(len(e.variants) for e in entries)
    print(f"Summary: ok={n_ok} warn={n_warn} error={n_err} "
          f"total_variants={total_variants}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
