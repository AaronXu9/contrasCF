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
import argparse
import json
import os
import shutil
import sys
from difflib import SequenceMatcher
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


_AA3to1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M",
}


def _af3_chains(path: Path) -> dict[str, str]:
    """{chain_id: sequence} from an af3.json (scalar or list `id`)."""
    d = json.loads(path.read_text())
    out: dict[str, str] = {}
    for s in d.get("sequences", []):
        pr = s.get("protein")
        if not pr:
            continue
        ids = pr["id"]
        for cid in (ids if isinstance(ids, list) else [ids]):
            out[str(cid)] = pr["sequence"]
    return out


def check_mutations_present(pdbid: str, variant: str, receptor_pdb: Path) -> tuple[int, int]:
    """(n_specified, n_present) — do the spec's mutations survive into the receptor?

    Guard for the defect where a receptor was written from a prediction that
    silently omitted the mutated chain, so docking ran against an effectively
    wild-type pocket (worst case 1bcu: a 26-residue stub with no mutation).
    A residue-count tolerance would NOT catch that — only checking the
    mutated positions themselves does.
    """
    wt_json = OUTPUT_ROOT / pdbid / "wt" / "af3.json"
    var_json = OUTPUT_ROOT / pdbid / variant / "af3.json"
    if not wt_json.exists() or not var_json.exists():
        return (0, 0)
    wt, mut = _af3_chains(wt_json), _af3_chains(var_json)
    muts = {cid: [(i, a, b) for i, (a, b) in enumerate(zip(wt[cid], seq)) if a != b]
            for cid, seq in mut.items()
            if cid in wt and len(wt[cid]) == len(seq)}
    n_spec = sum(len(v) for v in muts.values())
    if n_spec == 0:
        return (0, 0)

    st = gemmi.read_structure(str(receptor_pdb))
    st.remove_ligands_and_waters()
    rec = {ch.name: "".join(_AA3to1[r.name] for r in ch if r.name in _AA3to1)
           for ch in st[0]}
    rec = {k: v for k, v in rec.items() if v}

    n_present = 0
    for cid, mlist in muts.items():
        spec = mut[cid]
        best, best_r = None, 0.0
        for rname, rseq in rec.items():
            r = SequenceMatcher(None, spec, rseq, autojunk=False).ratio()
            if r > best_r:
                best, best_r = rname, r
        if best is None or best_r < 0.50:
            continue
        m = {}
        for a, b, size in SequenceMatcher(
                None, spec, rec[best], autojunk=False).get_matching_blocks():
            for k in range(size):
                m[a + k] = b + k
        for (i, _a, b) in mlist:
            j = m.get(i)
            if j is not None and rec[best][j] == b:
                n_present += 1
    return (n_spec, n_present)


def build_one(pdbid: str, variant: str, *, force: bool = False) -> str:
    """Return a status string."""
    v_dir = OUTPUT_ROOT / pdbid / variant
    cif = v_dir / f"af3msa_{pdbid}_{variant}_model_0.cif"
    if not cif.exists():
        return "missing_cif"
    docking = v_dir / "docking"
    box_path = docking / "box.json"
    if box_path.exists() and not force:
        return "skip_existing"
    docking.mkdir(parents=True, exist_ok=True)

    st = gemmi.read_structure(str(cif))

    # Receptor: protein chains stripped from predicted complex.
    receptor = docking / "receptor.pdb"
    receptor.write_text(_strip_to_protein_pdb(st))

    # Guard: refuse to ship a receptor that lost the mutations it is supposed
    # to carry. Written first so the file can be inspected, then flagged.
    n_spec, n_present = check_mutations_present(pdbid, variant, receptor)
    if n_spec and n_present == 0:
        (docking / "MUTATIONS_ABSENT").write_text(
            f"{n_spec} mutation(s) specified, 0 present in receptor.pdb\n")
        return "no_mutation_in_receptor"
    (docking / "MUTATIONS_ABSENT").unlink(missing_ok=True)
    if n_spec and n_present < n_spec:
        (docking / "MUTATIONS_PARTIAL").write_text(
            f"{n_present}/{n_spec} mutation(s) present in receptor.pdb\n")
    else:
        (docking / "MUTATIONS_PARTIAL").unlink(missing_ok=True)

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


def _backup_docking(pdbid: str, variant: str, backup_root: Path) -> bool:
    """Copy an existing docking/ cell aside once, before it is rebuilt."""
    src = OUTPUT_ROOT / pdbid / variant / "docking"
    dst = backup_root / pdbid / variant / "docking"
    if not src.is_dir() or dst.exists():
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(src, dst)
    return True


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true",
                    help="rebuild cells that already have box.json")
    ap.add_argument("--ids", default="",
                    help="comma-separated pdbids (default: all)")
    ap.add_argument("--backup-dir", default="",
                    help="copy existing docking/ cells here before rebuilding")
    args = ap.parse_args()

    if args.ids:
        names = [s.strip() for s in args.ids.split(",") if s.strip()]
    else:
        names = sorted(p.name for p in OUTPUT_ROOT.iterdir()
                       if p.is_dir() and len(p.name) == 4)
    print(f"{len(names)} PDB output directories; force={args.force}")

    backup_root = Path(args.backup_dir) if args.backup_dir else None
    counts: dict[str, int] = {}
    n_backed = 0
    flagged: list[str] = []
    for name in names:
        for variant in MUTANT_VARIANTS:
            if backup_root is not None and args.force:
                n_backed += int(_backup_docking(name, variant, backup_root))
            status = build_one(name, variant, force=args.force)
            counts[status] = counts.get(status, 0) + 1
            if status == "no_mutation_in_receptor":
                flagged.append(f"{name}/{variant}")
    if backup_root is not None:
        print(f"backed up {n_backed} docking cells -> {backup_root}")
    print("Summary:", counts)
    if flagged:
        print(f"\nMUTATIONS ABSENT in {len(flagged)} cells (marker file written):")
        print("  " + " ".join(flagged))
    return 0


if __name__ == "__main__":
    sys.exit(main())
