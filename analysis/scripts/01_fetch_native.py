#!/usr/bin/env python
"""Download native crystal references (1B38 CDK2-ATP, 2VWH GDH-glucose) into analysis/native/."""
from __future__ import annotations
import sys, urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from config import NATIVES, NATIVE_DIR


def fetch(pdb_id: str, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{pdb_id}.cif"
    if out.exists() and out.stat().st_size > 0:
        print(f"  [cached] {out}")
        return out
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    print(f"  [fetch ] {url}")
    urllib.request.urlretrieve(url, out)
    print(f"  [saved ] {out} ({out.stat().st_size:,} bytes)")
    return out


def main() -> None:
    print(f"Native dir: {NATIVE_DIR}")
    for key, nref in NATIVES.items():
        print(f"* {key}: {nref.pdb_id}")
        fetch(nref.pdb_id, NATIVE_DIR)


if __name__ == "__main__":
    main()
