"""GNINA pose RMSD analysis.

GNINA dockings are against the WT crystal protein (receptor.pdb) with the
box centered on the crystal ligand centroid. So the predicted pose is
already in the crystal frame — no Cα superposition needed.

For each (system, variant) cell with `<v_dir>/gnina/poses.sdf`:
  - Load top-1 GNINA pose (first conformer in poses.sdf).
  - Load crystal ligand from `data/casf2016/crystal_ligands/<pdbid>_ligand.sdf`.
  - Compute heavy-atom RMSD via maximum-common-substructure (MCS) so that
    modified ligands (halogenation / methylation / charge swap) can still
    be compared to crystal on their shared scaffold. For WT GNINA where
    the ligand is unchanged, this reduces to a standard heavy-atom RMSD.

Returns RMSD plus a `n_matched_heavy` count so a tiny-MCS match (e.g.,
charge_pos_3 vs ATP) can be filtered out from aggregates.
"""
from __future__ import annotations
import os
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFMCS

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
sys.path.insert(0, str(REPO_ROOT / "analysis"))

from casf_mutagenesis.config import CASF_LIGANDS, CASF_RAW  # noqa: E402
from casf_mutagenesis.analysis import (  # noqa: E402
    extract_protein_ca_near, superpose_by_index,
)
import gemmi  # noqa: E402


@dataclass
class GninaRecord:
    system: str
    variant: str
    module: str                  # "casf" | "ligand"
    engine: str = "gnina"        # "gnina" | "unidock2"
    status: str = "ok"
    error: str | None = None
    rmsd_a: float | None = None
    n_matched_heavy: int | None = None
    n_heavy_pose: int | None = None
    n_heavy_crystal: int | None = None


def _read_top_pose(sdf_path: Path) -> Chem.Mol | None:
    """Top-1 pose = first molecule record in an SDF that parses cleanly."""
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=True, sanitize=False)
    for m in suppl:
        if m is None:
            continue
        try:
            Chem.SanitizeMol(m)
        except Exception:
            pass
        return m
    return None


def _heavy_coords(mol: Chem.Mol, indices: list[int]) -> np.ndarray:
    conf = mol.GetConformer()
    return np.array(
        [(conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z)
         for i in indices],
        dtype=float,
    )


def _mcs_match_indices(
    crystal: Chem.Mol, pose: Chem.Mol
) -> tuple[list[int], list[int], int]:
    """Find the maximum common substructure and return (crystal_idx, pose_idx, n)."""
    # Quick path: if pose is a superstructure of crystal (halogenation,
    # methylation add atoms), GetSubstructMatch is fast.
    direct = pose.GetSubstructMatch(crystal)
    if direct:
        return list(range(crystal.GetNumAtoms())), list(direct), crystal.GetNumAtoms()
    # Reverse direction (crystal contains pose — unusual but try).
    direct = crystal.GetSubstructMatch(pose)
    if direct:
        return list(direct), list(range(pose.GetNumAtoms())), pose.GetNumAtoms()
    # Fall back to MCS (handles charge_swap etc. where neither is a
    # substructure of the other).
    res = rdFMCS.FindMCS(
        [crystal, pose],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        timeout=10,
        completeRingsOnly=False,
    )
    if res.numAtoms == 0 or res.smartsString == "":
        return [], [], 0
    pattern = Chem.MolFromSmarts(res.smartsString)
    if pattern is None:
        return [], [], 0
    crystal_match = crystal.GetSubstructMatch(pattern)
    pose_match = pose.GetSubstructMatch(pattern)
    if not crystal_match or not pose_match:
        return [], [], 0
    return list(crystal_match), list(pose_match), len(crystal_match)


def _aligned_rmsd(crystal_pts: np.ndarray, pose_pts: np.ndarray) -> float:
    """Heavy-atom RMSD assuming matched atom-by-atom (no superposition).

    GNINA poses are already in the crystal frame so no superposition is
    needed. RMSD reflects the distance between the docked pose and the
    crystal ligand position.
    """
    if crystal_pts.shape != pose_pts.shape:
        return float("nan")
    diff = crystal_pts - pose_pts
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def _ca_transform(receptor_pdb: Path, crystal_pdb: Path, lig_xyz: np.ndarray):
    """Return a Superposition that maps receptor → crystal frame.

    Used for casf_mutagenesis rem/pack/inv where GNINA docked against an
    AF3-predicted mutant receptor — pose coords are in the AF3 frame and
    must be transformed before comparing to the crystal ligand.

    For WT and ligand_mutagenesis, the docking receptor IS the crystal
    protein → identity transform; we skip this call entirely there.
    """
    st_rec = gemmi.read_structure(str(receptor_pdb))
    st_crys = gemmi.read_structure(str(crystal_pdb))
    rec_ca = extract_protein_ca_near(st_rec, lig_xyz)
    crys_ca = extract_protein_ca_near(st_crys, lig_xyz)
    return superpose_by_index(rec_ca, crys_ca)


def analyze_gnina(
    system: str, variant: str, module: str,
    gnina_sdf: Path, crystal_sdf: Path | None = None,
    *, receptor_pdb: Path | None = None,
    engine: str = "gnina",
) -> GninaRecord:
    rec = GninaRecord(system=system, variant=variant, module=module, engine=engine)
    if crystal_sdf is None:
        crystal_sdf = CASF_LIGANDS / f"{system}_ligand.sdf"
    try:
        if not gnina_sdf.exists():
            rec.status = "missing_pose"; return rec
        if not crystal_sdf.exists():
            rec.status = "missing_crystal"; return rec
        pose = _read_top_pose(gnina_sdf)
        crystal = _read_top_pose(crystal_sdf)
        if pose is None or crystal is None:
            rec.status = "parse_error"; return rec
        rec.n_heavy_pose = pose.GetNumHeavyAtoms()
        rec.n_heavy_crystal = crystal.GetNumHeavyAtoms()
        c_idx, p_idx, n = _mcs_match_indices(crystal, pose)
        if n < 3:
            rec.status = "no_match"; return rec
        rec.n_matched_heavy = n
        c_pts = _heavy_coords(crystal, c_idx)
        p_pts = _heavy_coords(pose, p_idx)
        # If a receptor frame differs from the crystal frame (only happens
        # for casf rem/pack/inv where GNINA docked into AF3-predicted
        # mutants), align receptor → crystal and apply the transform to
        # the pose before computing RMSD.
        needs_align = (module == "casf" and variant != "wt") and receptor_pdb is not None
        if needs_align:
            crystal_pdb = CASF_RAW / system / f"{system}_protein.pdb"
            if crystal_pdb.exists():
                # Use crystal ligand centroid for pocket-aligned chain pick
                # (we know where the ligand sits in the crystal frame).
                crys_centroid = np.asarray(_heavy_coords(crystal, list(range(crystal.GetNumHeavyAtoms()))))
                try:
                    sup = _ca_transform(receptor_pdb, crystal_pdb, crys_centroid)
                    p_pts = sup.apply(p_pts)
                except Exception as exc:
                    rec.status = "align_failed"
                    rec.error = f"{type(exc).__name__}: {exc}"
                    return rec
        rec.rmsd_a = round(_aligned_rmsd(c_pts, p_pts), 3)
    except Exception as exc:
        rec.status = "error"
        rec.error = f"{type(exc).__name__}: {exc}"
    return rec
