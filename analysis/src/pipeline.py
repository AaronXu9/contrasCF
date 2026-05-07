"""Orchestrator: `run_one(model, case) -> dict` producing a results row."""
from __future__ import annotations
from pathlib import Path

import numpy as np

from config import CASES, MODELS, COMMON_SUBSETS, case_dir, find_first
from loaders import read_structure, extract_protein_ca, select_target_ligand, extract_all_protein_heavy
from native import load_native
from align import superpose_ca
from ligand_match import match_all, match_common
from rmsd import paired_rmsd, bestfit_rmsd
from confidence import extract_confidence
from clashes import count_clashes


NAN = float("nan")


def run_one(model: str, case: str) -> dict:
    spec = CASES[case]
    d = case_dir(model, case)
    structure_path = find_first(d, MODELS[model].structure_globs)
    row = dict(
        model=model, case=case, family=spec.family, native_pdb=None,
        structure_path=str(structure_path) if structure_path else None,
        n_ca_aligned=np.nan, ca_rmsd=np.nan,
        ligand_rmsd_all=np.nan, ligand_rmsd_common=np.nan,
        ligand_rmsd_bestfit=np.nan,
        n_atoms_all=np.nan, n_atoms_common=np.nan,
        ptm=np.nan, iptm=np.nan, ligand_iptm=np.nan,
        ligand_plddt_mean=np.nan, ranking_score=np.nan,
        n_clashes=np.nan, min_prot_lig_dist=np.nan,
        has_clash_geom=None, has_clash_reported=None,
        lig_chain=None, lig_resname=None, n_lig_heavy=np.nan,
        # Physics-based docking scores (NaN for co-folding models).
        dock_vina_score=np.nan, dock_cnn_score=np.nan, dock_cnn_affinity=np.nan,
        error=None,
    )
    if structure_path is None:
        row["error"] = "no structure file"
        return row
    try:
        nat = load_native(spec.native_key)
        row["native_pdb"] = nat.pdb_id

        st = read_structure(structure_path)
        pred_ca = extract_protein_ca(st)
        pred_lig, dbg = select_target_ligand(st, spec.target_smiles)
        if pred_lig is None:
            row["error"] = f"target ligand not found: {dbg}"
            return row
        row["lig_chain"] = pred_lig.chain
        row["lig_resname"] = pred_lig.resname
        row["n_lig_heavy"] = len(pred_lig.heavy_coords)

        # Superpose predicted Cα onto native Cα.
        sup = superpose_ca(pred_ca, nat.ca)
        row["n_ca_aligned"] = sup.n_paired
        row["ca_rmsd"] = sup.rmsd

        # Apply transform to predicted ligand heavy atoms and predicted protein heavy atoms.
        pred_lig_xyz = sup.apply(pred_lig.heavy_coords)
        pred_protein_heavy_xyz = sup.apply(
            extract_all_protein_heavy(st, pred_ca.residues[0][0])
        )

        nat_lig_xyz = nat.ligand.heavy_coords
        pred_mol = pred_lig.mol
        nat_mol = nat.ligand.mol

        # RMSD over full-molecule MCS (symmetry-minimizing).
        if pred_mol is not None and nat_mol is not None:
            p_idx, n_idx, rmsd_all = match_all(pred_mol, nat_mol, pred_lig_xyz, nat_lig_xyz)
            if p_idx:
                row["ligand_rmsd_all"] = rmsd_all
                row["n_atoms_all"] = len(p_idx)

        # RMSD over common core SMARTS.
        smarts = COMMON_SUBSETS.get(spec.common_smarts_key)
        if pred_mol is not None and nat_mol is not None and smarts is not None:
            p_idx, n_idx, rmsd_c = match_common(
                pred_mol, nat_mol, pred_lig_xyz, nat_lig_xyz, smarts,
            )
            if p_idx:
                row["ligand_rmsd_common"] = rmsd_c
                row["n_atoms_common"] = len(p_idx)

        # Paper-style RMSD: symmetric, optimal-rigid-body-fit, conformation-only.
        # Input coords can be in any frame (GetBestRMS does its own Kabsch fit),
        # but we pass aligned coords anyway for consistency with other metrics.
        row["ligand_rmsd_bestfit"] = bestfit_rmsd(
            pred_mol, nat_mol, pred_lig_xyz, nat_lig_xyz,
        )

        # Clashes (post-alignment, on full heavy atoms of predicted protein).
        n_cl, min_d = count_clashes(pred_protein_heavy_xyz, pred_lig_xyz, cutoff=2.0)
        row["n_clashes"] = n_cl
        row["min_prot_lig_dist"] = min_d
        row["has_clash_geom"] = bool(n_cl > 0)

        # Confidence.
        conf = extract_confidence(model, case, pred_lig)
        row.update(conf)

    except Exception as e:
        row["error"] = f"{type(e).__name__}: {e}"
    return row
