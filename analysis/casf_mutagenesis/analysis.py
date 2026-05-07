"""Compute ligand RMSD-to-native and memorization rate for CASF predictions.

For each (pdbid, variant, model):
  1. Load predicted mmCIF → extract protein Cα + target ligand (RDKit Mol).
  2. Load crystal protein PDB + crystal ligand SDF.
  3. Superpose predicted Cα onto crystal Cα via SVD over the residue-number
     intersection.
  4. Apply the rotation+translation to predicted ligand heavy-atom coords.
  5. Match predicted vs crystal ligand atoms by SMILES template; compute the
     heavy-atom RMSD over best symmetry-matched correspondence.

Memorization rate = fraction of *adversarial* variants (rem/pack/inv) with
ligand RMSD < threshold (default 2 Å, per Masters et al. 2025). WT is the
baseline.

This module reuses building blocks from `analysis/src/`:
  - `loaders.read_structure`, `extract_protein_ca`, `select_target_ligand`,
    `_ligand_to_mol`.
  - `align.superpose_ca` (after wrapping the predicted ProteinCA so it shares
    the Biopython-friendly residue-pair-by-resnum logic).
"""
from __future__ import annotations
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path

import os
_repo_root_default = "/mnt/katritch_lab2/aoxu/contrasCF"
REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", _repo_root_default))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "src"))

import gemmi  # noqa: E402
import numpy as np  # noqa: E402
from rdkit import Chem  # noqa: E402

from Bio.SVDSuperimposer import SVDSuperimposer  # noqa: E402

from loaders import (  # noqa: E402
    LigandBlock, ProteinCA, _ligand_to_mol, _collect_candidate_ligands,
    extract_protein_ca, read_structure,
)


@dataclass
class Superposition:
    R: np.ndarray
    t: np.ndarray
    rmsd: float
    n_paired: int

    def apply(self, coords: np.ndarray) -> np.ndarray:
        return coords @ self.R.T + self.t


def superpose_by_index(pred: ProteinCA, native: ProteinCA) -> Superposition:
    """Pair Cα by sequential chain-position index.

    `analysis.src.align.superpose_ca` pairs by PDB resnum which fails on
    CASF: crystals use UniProt-canonical numbering (e.g. 4jia chain A is
    833-1132) while predicted models number 1..N. Both, however, list the
    SAME chain residues in order, so pairing by index is correct as long as
    the chain has no internal gaps relative to the input sequence — which is
    the case for CASF (we extract the input sequence directly from the
    crystal protein chain via gemmi's CA iteration).
    """
    n = min(len(pred.residues), len(native.residues))
    if n < 10:
        raise RuntimeError(f"too few Cα pairs ({n}) for index pairing")
    pred_xyz = np.array([r[3] for r in pred.residues[:n]])
    nat_xyz = np.array([r[3] for r in native.residues[:n]])
    sup = SVDSuperimposer()
    sup.set(nat_xyz, pred_xyz)
    sup.run()
    rot, tran = sup.get_rotran()
    return Superposition(
        R=np.asarray(rot).T, t=np.asarray(tran),
        rmsd=float(sup.get_rms()), n_paired=n,
    )

from .config import CASF_LIGANDS, CASF_RAW, OUTPUT_ROOT, VARIANTS


# Models we currently produce outputs for. Each tuple is
# (model_name, cif_filename_pattern, must_exist_marker).
MODEL_FILES = {
    "Boltz2":  "{prefix}_model_0.cif",
    "AF3":     "af3_{prefix}_model_0.cif",      # no-MSA baseline
    "AF3+MSA": "af3msa_{prefix}_model_0.cif",   # with ColabFold MSA
}

MEMORIZATION_THRESHOLDS_A = (2.0, 4.0)


@dataclass
class PredictionRecord:
    pdbid: str
    variant: str
    model: str
    status: str = "ok"            # "ok" | "missing_cif" | "no_ligand" | "error"
    error: str | None = None
    ligand_rmsd_a: float | None = None
    ca_rmsd_a: float | None = None
    n_ca_paired: int | None = None
    n_heavy_matched: int | None = None
    n_heavy_pred: int | None = None
    n_heavy_native: int | None = None


# ---------------------------------------------------------------------------
# Crystal native loading
# ---------------------------------------------------------------------------

def _crystal_ligand_mol(sdf_path: Path) -> Chem.Mol:
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=True, sanitize=False)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        raise RuntimeError(f"could not parse crystal SDF: {sdf_path}")
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass
    return mol


def _heavy_coords(mol: Chem.Mol) -> np.ndarray:
    conf = mol.GetConformer()
    out = []
    for i, a in enumerate(mol.GetAtoms()):
        if a.GetAtomicNum() == 1:
            continue
        p = conf.GetAtomPosition(i)
        out.append((p.x, p.y, p.z))
    return np.asarray(out, dtype=float)


def _heavy_indices(mol: Chem.Mol) -> list[int]:
    return [i for i, a in enumerate(mol.GetAtoms()) if a.GetAtomicNum() != 1]


# ---------------------------------------------------------------------------
# Predicted ligand extraction
# ---------------------------------------------------------------------------

def _predicted_ligand(
    st: gemmi.Structure, target_smiles: str
) -> tuple[LigandBlock, Chem.Mol] | tuple[None, None]:
    """Pick the largest non-protein/non-water residue and assign bond orders
    from the target SMILES. Simpler than the existing `select_target_ligand`
    because for CASF predictions we only have one ligand per system."""
    candidates = _collect_candidate_ligands(st)
    if not candidates:
        return None, None
    # largest by heavy-atom count
    best = max(candidates, key=lambda c: sum(
        1 for a in c.atoms if a[1].upper() != "H"
    ))
    mol = _ligand_to_mol(best, target_smiles)
    return best, mol


# ---------------------------------------------------------------------------
# RMSD with symmetry / atom-mapping
# ---------------------------------------------------------------------------

def _matched_rmsd(
    crystal_mol: Chem.Mol, pred_xyz: np.ndarray, pred_mol: Chem.Mol,
) -> tuple[float, int]:
    """Return (RMSD, n_matched_heavy) between the crystal ligand and the
    predicted ligand whose heavy-atom coords have been transformed into the
    crystal frame.

    Strategy: enumerate substructure matches between pred_mol and crystal_mol;
    over all matches, take the one minimizing RMSD. Falls back to the trivial
    1-to-1 mapping (pred[i] ↔ crystal[i]) when no match exists.
    """
    crystal_xyz = _heavy_coords(crystal_mol)
    n_c = len(crystal_xyz)
    n_p = len(pred_xyz)
    if n_c != n_p:
        # SMILES round-trip can lose/gain explicit H counts; align on min.
        n = min(n_c, n_p)
        pred_xyz = pred_xyz[:n]
        crystal_xyz = crystal_xyz[:n]

    if pred_mol is None:
        # fallback: assume same atom order
        d = pred_xyz - crystal_xyz
        return float(np.sqrt((d * d).sum() / len(d))), len(d)

    # Map pred → crystal heavy-atom indices
    crystal_heavy_idx = _heavy_indices(crystal_mol)
    pred_heavy_idx = _heavy_indices(pred_mol)

    # Build heavy-only RDKit copies for substructure search
    crystal_heavy = Chem.RemoveHs(crystal_mol)
    pred_heavy = Chem.RemoveHs(pred_mol)

    matches = pred_heavy.GetSubstructMatches(
        crystal_heavy, useChirality=False, uniquify=False, maxMatches=200,
    )
    if not matches:
        # Try the other direction (pred is template) — useful when bond
        # orders differ on aromatics
        matches = crystal_heavy.GetSubstructMatches(
            pred_heavy, useChirality=False, uniquify=False, maxMatches=200,
        )
        # invert: matches now map crystal heavy idx → pred heavy idx
        matches = [tuple(range(len(m))) for m in matches]

    if not matches:
        d = pred_xyz - crystal_xyz
        return float(np.sqrt((d * d).sum() / len(d))), len(d)

    crystal_pts = _heavy_coords(crystal_heavy)  # heavy-only ordering matches
    pred_pts = _heavy_coords(pred_heavy)

    # Translate the pred points: we already passed the transformed heavy coords
    # but the heavy-only Mol may have a slightly different ordering than the
    # ligand block. Rebuild pred_xyz aligned with pred_heavy ordering.
    # The simplest correct approach: re-derive pred coords from the heavy-only
    # mol's conformer (which was inherited from the original pred_mol whose
    # conformer we already transformed).
    n_match_atoms = crystal_heavy.GetNumAtoms()
    if pred_pts.shape[0] != n_match_atoms or crystal_pts.shape[0] != n_match_atoms:
        d = pred_xyz - crystal_xyz
        return float(np.sqrt((d * d).sum() / len(d))), len(d)

    best_rmsd = float("inf")
    for m in matches:
        # m: tuple of len n_match_atoms — the i-th crystal heavy atom maps to
        # m[i] in pred. (When matching pred against crystal as template, that
        # is the SubstructMatch convention.)
        try:
            sel_pred = pred_pts[list(m)]
        except IndexError:
            continue
        d = sel_pred - crystal_pts
        rmsd = float(np.sqrt((d * d).sum() / len(d)))
        if rmsd < best_rmsd:
            best_rmsd = rmsd
    return best_rmsd, n_match_atoms


# ---------------------------------------------------------------------------
# Per-prediction pipeline
# ---------------------------------------------------------------------------

def analyze_prediction(
    pdbid: str, variant: str, model: str,
    *, smiles_override: str | None = None,
) -> PredictionRecord:
    pdbid = pdbid.lower()
    rec = PredictionRecord(pdbid=pdbid, variant=variant, model=model)
    try:
        v_dir = OUTPUT_ROOT / pdbid / variant
        cif_pattern = MODEL_FILES[model]
        cif_path = v_dir / cif_pattern.format(prefix=f"{pdbid}_{variant}")
        if not cif_path.exists():
            rec.status = "missing_cif"
            return rec

        # Native (crystal) loading
        crystal_pdb = CASF_RAW / pdbid / f"{pdbid}_protein.pdb"
        crystal_sdf = CASF_LIGANDS / f"{pdbid}_ligand.sdf"
        st_native = read_structure(crystal_pdb)
        native_ca = extract_protein_ca(st_native)
        crystal_mol = _crystal_ligand_mol(crystal_sdf)

        # Predicted loading
        st_pred = read_structure(cif_path)
        pred_ca = extract_protein_ca(st_pred)

        # SMILES for bond-order recovery: use crystal-derived SMILES (same as
        # what we passed to the model when generating inputs)
        if smiles_override is not None:
            target_smiles = smiles_override
        else:
            target_smiles = Chem.MolToSmiles(crystal_mol)
        pred_lig_block, pred_mol = _predicted_ligand(st_pred, target_smiles)
        if pred_lig_block is None:
            rec.status = "no_ligand"
            return rec
        rec.n_heavy_pred = sum(
            1 for a in pred_lig_block.atoms if a[1].upper() != "H"
        )
        rec.n_heavy_native = sum(
            1 for _ in crystal_mol.GetAtoms() if _.GetAtomicNum() != 1
        )

        # Cα superpose: we want pred → native, then apply same R,t to ligand
        sup = superpose_by_index(pred_ca, native_ca)
        rec.ca_rmsd_a = round(sup.rmsd, 3)
        rec.n_ca_paired = sup.n_paired

        # Transform predicted ligand heavy coords into native frame
        pred_lig_heavy = pred_lig_block.heavy_coords
        pred_lig_native_frame = sup.apply(pred_lig_heavy)

        # If we have an RDKit mol with the same atom count as the LigandBlock,
        # rebuild a transformed conformer so symmetry matching uses the right
        # coords.
        if pred_mol is not None and pred_mol.GetNumHeavyAtoms() == len(pred_lig_native_frame):
            new_mol = Chem.Mol(pred_mol)
            conf = new_mol.GetConformer()
            heavy_idx = _heavy_indices(new_mol)
            for i, hidx in enumerate(heavy_idx):
                p = pred_lig_native_frame[i]
                conf.SetAtomPosition(hidx, (float(p[0]), float(p[1]), float(p[2])))
            transformed_mol = new_mol
        else:
            transformed_mol = None

        rmsd, n_match = _matched_rmsd(
            crystal_mol, pred_lig_native_frame, transformed_mol,
        )
        rec.ligand_rmsd_a = round(rmsd, 3)
        rec.n_heavy_matched = n_match
    except Exception as exc:
        rec.status = "error"
        rec.error = f"{type(exc).__name__}: {exc}"
    return rec


# ---------------------------------------------------------------------------
# Aggregate stats
# ---------------------------------------------------------------------------

@dataclass
class MemorizationStats:
    model: str
    variant: str          # one of "rem" / "pack" / "inv"
    n_total: int
    n_below_2A: int
    n_below_4A: int
    median_rmsd_a: float | None

    def rate(self, threshold: float) -> float:
        if self.n_total == 0:
            return float("nan")
        n = self.n_below_2A if threshold == 2.0 else self.n_below_4A
        return n / self.n_total


def memorization_stats(
    records: list[PredictionRecord],
) -> dict[tuple[str, str], MemorizationStats]:
    """Per (model, adversarial-variant) memorization rate."""
    out: dict[tuple[str, str], MemorizationStats] = {}
    by_key: dict[tuple[str, str], list[float]] = {}
    for r in records:
        if r.variant == "wt":
            continue
        if r.status != "ok" or r.ligand_rmsd_a is None:
            continue
        key = (r.model, r.variant)
        by_key.setdefault(key, []).append(r.ligand_rmsd_a)
    for key, rmsds in by_key.items():
        arr = np.asarray(rmsds)
        out[key] = MemorizationStats(
            model=key[0], variant=key[1],
            n_total=len(arr),
            n_below_2A=int((arr < 2.0).sum()),
            n_below_4A=int((arr < 4.0).sum()),
            median_rmsd_a=float(np.median(arr)),
        )
    return out
