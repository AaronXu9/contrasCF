"""Build docking inputs for the rem/pack/inv mutant variants.

For each (pdbid, variant ∈ {rem, pack, inv}) cell that has an AF3+MSA
mutant CIF at `outputs/<pdbid>/<variant>/af3msa_<prefix>_model_0.cif`:

  1. Strip the predicted complex to protein only (all standard-AA chains
     preserved) → `outputs/<pdbid>/<variant>/docking/receptor.pdb`.
  2. Copy the crystal ligand SDF → `docking/ligand.sdf` (binding-site
     mutations modify the PROTEIN, not the ligand — same crystal ligand
     chemistry for all four variants).
  3. Extract the AF3-predicted ligand heavy-atom centroid → `box.json`
     (25 × 25 × 25 Å). The box at the predicted ligand position gives
     GNINA the model's best guess at the pocket; if the model placed the
     ligand at an alternative site, the box covers that.

Idempotent: skips any cell whose `docking/box.json` already exists.
"""
from __future__ import annotations
import json
import os
import shutil
import sys
from pathlib import Path

import gemmi
import numpy as np

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
sys.path.insert(0, str(REPO_ROOT / "analysis"))

from casf_mutagenesis.config import (  # noqa: E402
    CASF_LIGANDS, DOCKING_BOX_A, OUTPUT_ROOT, AA3_STANDARD,
)
from casf_mutagenesis.inputs_docking import _strip_to_protein_pdb  # noqa: E402

MUTANT_VARIANTS = ("rem", "pack", "inv")


def _predicted_ligand_centroid(st: gemmi.Structure) -> tuple[float, float, float]:
    """Centroid of the AF3+MSA predicted ligand (largest non-protein, non-water
    residue), in the predicted-protein frame. Used to position the GNINA box.
    """
    best_xyz: list[tuple[float, float, float]] | None = None
    best_n = 0
    for chain in st[0]:
        for r in chain:
            if r.name in AA3_STANDARD or r.name in ("HOH", "WAT"):
                continue
            xyz = [(a.pos.x, a.pos.y, a.pos.z) for a in r if a.element.name != "H"]
            if len(xyz) > best_n:
                best_n = len(xyz)
                best_xyz = xyz
    if best_xyz is None:
        raise RuntimeError("no candidate ligand in predicted CIF")
    arr = np.asarray(best_xyz)
    return tuple(arr.mean(axis=0).tolist())


def build_one(pdbid: str, variant: str) -> str:
    """Return a status string."""
    v_dir = OUTPUT_ROOT / pdbid / variant
    cif = v_dir / f"af3msa_{pdbid}_{variant}_model_0.cif"
    if not cif.exists():
        return "missing_cif"
    docking = v_dir / "docking"
    box_path = docking / "box.json"
    if box_path.exists():
        return "skip_existing"
    docking.mkdir(parents=True, exist_ok=True)

    st = gemmi.read_structure(str(cif))

    # Receptor: protein chains stripped from predicted complex.
    (docking / "receptor.pdb").write_text(_strip_to_protein_pdb(st))

    # Ligand: copy the crystal SDF (binding-site mutations don't change
    # the ligand chemistry).
    crystal_lig = CASF_LIGANDS / f"{pdbid}_ligand.sdf"
    if not crystal_lig.exists():
        return "missing_crystal_lig"
    shutil.copy(crystal_lig, docking / "ligand.sdf")

    # Box centered on AF3-predicted ligand position. If the model placed
    # the ligand at the original pocket, GNINA docks there; if it ejected
    # to a different site, GNINA explores that.
    cx, cy, cz = _predicted_ligand_centroid(st)
    box_path.write_text(json.dumps(
        {"center": [cx, cy, cz], "size": list(DOCKING_BOX_A)}, indent=2
    ))
    return "ok"


def main() -> int:
    pdb_dirs = sorted(p for p in OUTPUT_ROOT.iterdir() if p.is_dir() and len(p.name) == 4)
    print(f"Found {len(pdb_dirs)} PDB output directories")
    counts: dict[str, int] = {}
    for pdb_dir in pdb_dirs:
        for variant in MUTANT_VARIANTS:
            status = build_one(pdb_dir.name, variant)
            counts[status] = counts.get(status, 0) + 1
    print("Summary:", counts)
    return 0


if __name__ == "__main__":
    sys.exit(main())
