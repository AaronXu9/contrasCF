#!/usr/bin/env python3
"""Protonate AF3-predicted mutant docking receptors to match the WT arm.

WHY. The wt cell docks into the HiQBind-curated crystal receptor, which is
protonated -- `REMARK 1 CREATED WITH OPENMM 8.1.1` and PDB v3 hydrogen naming
(HB2/HB3, HG2/HG3) identify PDBFixer/OpenMM as the tool. Every rem/pack/inv cell
docks into an AF3 prediction with ZERO hydrogens. Vina-family engines are
insensitive to this (affinity is bit-identical, verified), but SurfDock is not:
its MSMS surface takes per-atom radii including H, and its H-bond channel is
literally `if atom_name in polarHydrogens[resname]`
(SurfDock/comp_surface/prepare_target/computeCharges.py:127). Stripping H from a
crystal receptor costs up to 36% of surface vertices and moves the rest 0.57 A.

So the WT and mutant arms hand SurfDock differently-prepared receptors. This
script removes that asymmetry by reproducing the WT arm's own protonation --
same tool, same version (openmm 8.1.1 in the SurfDock env) -- rather than
picking a nominally better one. Consistency between arms is the goal.

Writes `<cell>/docking/receptor_protonated.pdb` next to the original; nothing is
overwritten, so the unprotonated cells stay reproducible.

Run (SurfDock env has pdbfixer + openmm 8.1.1):
  /home/aoxu/miniconda3/envs/SurfDock/bin/python build/06_protonate_receptors.py \
      --systems 3mss,4eo8 --variants rem,pack,inv
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
sys.path.insert(0, str(REPO_ROOT / "analysis"))

from casf_mutagenesis.config import OUTPUT_ROOT, VARIANTS  # noqa: E402

PH = float(os.environ.get("CONTRASCF_PROTONATION_PH", "7.0"))


def protonate(src: Path, dst: Path, ph: float = PH) -> dict:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    fixer = PDBFixer(filename=str(src))
    # Hydrogens ONLY. Do not let PDBFixer build missing loops: an AF3 prediction
    # is already complete, and inventing residues would change the pocket we are
    # trying to hold fixed.
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    n_missing = sum(len(v) for v in fixer.missingAtoms.values())
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    with dst.open("w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)

    def count_h(p: Path) -> tuple[int, int]:
        h = tot = 0
        for line in p.read_text().splitlines():
            if line.startswith(("ATOM", "HETATM")):
                tot += 1
                if line[76:78].strip() == "H" or line[12:16].strip().startswith("H"):
                    h += 1
        return h, tot

    h0, t0 = count_h(src)
    h1, t1 = count_h(dst)
    return {"h_before": h0, "atoms_before": t0, "h_after": h1, "atoms_after": t1,
            "heavy_atoms_added": n_missing}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--systems", default=None, help="comma-separated ids (default: all)")
    ap.add_argument("--variants", default="rem,pack,inv",
                    help="wt is already protonated by HiQBind; default skips it")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    systems = ([s.strip().lower() for s in args.systems.split(",") if s.strip()]
               if args.systems else
               sorted(d.name for d in OUTPUT_ROOT.iterdir()
                      if d.is_dir() and len(d.name) == 4))
    variants = [v.strip() for v in args.variants.split(",") if v.strip() in VARIANTS]

    n_ok = n_skip = n_fail = 0
    for system in systems:
        for variant in variants:
            src = OUTPUT_ROOT / system / variant / "docking" / "receptor.pdb"
            if not src.exists():
                continue
            dst = src.with_name("receptor_protonated.pdb")
            if dst.exists() and not args.force:
                n_skip += 1
                continue
            try:
                st = protonate(src, dst)
                n_ok += 1
                print(f"  {system}/{variant}: {st['atoms_before']} -> "
                      f"{st['atoms_after']} atoms  (+{st['h_after'] - st['h_before']} H, "
                      f"+{st['heavy_atoms_added']} heavy)", flush=True)
            except Exception as exc:
                n_fail += 1
                print(f"  {system}/{variant}: FAILED {type(exc).__name__}: {exc}",
                      flush=True)
    print(f"\nprotonated={n_ok} skipped={n_skip} failed={n_fail}  (pH {PH})")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
