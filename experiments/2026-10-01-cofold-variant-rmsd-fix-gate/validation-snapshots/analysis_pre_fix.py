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
import json
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
# Each entry is a per-model file spec. `{prefix}` is substituted with
# f"{pdbid}_{variant}" and `{rank}` enumerates 0..N-1 of diffusion samples.
# If no `_model_<rank>` files exist for rank > 0, we fall back to the legacy
# single-pose layout where only `_model_0` exists on disk.
MODEL_FILES = {
    "Boltz2": {
        "cif_glob":    "{prefix}_model_*.cif",
        "conf":        "confidence_{prefix}_model_{rank}.json",
        # Boltz-2 also emits binding-affinity predictions (one per system,
        # not per pose) when the YAML requests `properties: - affinity:`.
        "affinity":    "affinity_{prefix}.json",
    },
    "AF3": {  # no-MSA baseline
        "cif_glob":    "af3_{prefix}_model_*.cif",
        "conf":        "af3_summary_confidences_{prefix}_{rank}.json",
        "conf_legacy": "af3_summary_confidences_{prefix}.json",  # rank-0 fallback
    },
    "AF3+MSA": {  # with ColabFold MSA
        "cif_glob":    "af3msa_{prefix}_model_*.cif",
        "conf":        "af3msa_summary_confidences_{prefix}_{rank}.json",
        "conf_legacy": "af3msa_summary_confidences_{prefix}.json",
    },
}

MEMORIZATION_THRESHOLDS_A = (2.0, 4.0)


@dataclass
class PredictionRecord:
    pdbid: str
    variant: str
    model: str
    pose_idx: int = 0             # 0 = top-ranked by model confidence
    status: str = "ok"            # "ok" | "missing_cif" | "no_ligand" | "error"
    error: str | None = None
    # Three RMSD variants, all in Å, all symmetry-corrected on heavy atoms:
    #   ligand_rmsd_a         : Cα-superpose on the *ligand-near chain only*,
    #                           apply R/t to the ligand, then heavy-atom RMSD
    #                           (canonical "did the model place it in the right
    #                           pocket?" metric).
    #   ligand_rmsd_fullca_a  : Cα-superpose on ALL chains' Cα (sequential pair
    #                           per chain), apply R/t to the ligand, then heavy
    #                           RMSD. Larger than ligand_rmsd_a on multi-chain
    #                           systems when AF3 misorients inter-chain.
    #   bestfit_rmsd_a        : symmetry-corrected Kabsch rigid fit of the
    #                           ligand onto the crystal ligand — pocket-blind.
    #                           Answers "is the ligand's internal geometry
    #                           correct?", which is the paper's quoted metric.
    ligand_rmsd_a: float | None = None
    ligand_rmsd_fullca_a: float | None = None
    bestfit_rmsd_a: float | None = None
    ca_rmsd_a: float | None = None
    ca_rmsd_fullca_a: float | None = None
    n_ca_paired: int | None = None
    n_ca_paired_fullca: int | None = None
    n_heavy_matched: int | None = None
    n_heavy_pred: int | None = None
    n_heavy_native: int | None = None
    # Confidence sidecar fields — None if JSON missing / key absent.
    # Boltz-2: confidence_score, iptm, ptm, ligand_iptm, complex_plddt
    # AF3 / AF3+MSA: ranking_score, iptm, ptm (no ligand_iptm / plddt)
    confidence_score: float | None = None
    iptm: float | None = None
    ptm: float | None = None
    ligand_iptm: float | None = None
    complex_plddt: float | None = None
    ranking_score: float | None = None
    # Boltz-2 binding-affinity sidecar (one value per system, broadcast to
    # every pose of that cell). None for AF3 / AF3+MSA which don't predict
    # affinity, and None if the affinity JSON is absent.
    #   affinity_pred_value:        log[IC50] in µM (lower = tighter)
    #   affinity_probability_binary: P(binder) ∈ [0, 1]
    affinity_pred_value: float | None = None
    affinity_probability_binary: float | None = None


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
# Pocket-aligned Cα extraction (chain closest to ligand)
# ---------------------------------------------------------------------------

def extract_protein_ca_near(
    st: gemmi.Structure, lig_xyz: np.ndarray,
) -> ProteinCA:
    """Like `loaders.extract_protein_ca` but pick the protein chain whose
    nearest Cα is closest to the ligand centroid — not the largest chain.

    Rationale: CASF systems include homo-multimers where the largest chain
    is not necessarily the chain binding the ligand. 4w9l is a homo-trimer
    (chains C/K/L); crystal ligand binds chain C, predicted-Boltz ligand
    binds chain L. Picking "largest" → both pick chain L → superpose L→L →
    ligand RMSD reflects inter-chain offset rather than placement accuracy.
    Picking "closest to ligand" on each side gives the biophysically
    meaningful pocket-aligned superposition.

    Falls back to largest-chain when no chain has > 20 polymer residues
    (single-chain systems).
    """
    from loaders import STANDARD_AA, ProteinCA
    centroid = lig_xyz.mean(axis=0)
    model = st[0]
    candidates: list[tuple[float, "gemmi.Chain"]] = []
    for ch in model:
        # Need both polymer residues AND a Cα to measure proximity.
        ca_coords = []
        for r in ch:
            if r.name not in STANDARD_AA:
                continue
            for a in r:
                if a.name == "CA":
                    ca_coords.append((a.pos.x, a.pos.y, a.pos.z))
                    break
        if len(ca_coords) < 20:
            continue
        d_min = float(
            np.min(np.linalg.norm(np.asarray(ca_coords) - centroid, axis=1))
        )
        candidates.append((d_min, ch))
    if not candidates:
        # No protein chain with >20 residues — fall back to library helper.
        from loaders import extract_protein_ca
        return extract_protein_ca(st)
    candidates.sort(key=lambda x: x[0])
    chain = candidates[0][1]

    residues: list[tuple[str, int, str, np.ndarray]] = []
    for r in chain:
        if r.name not in STANDARD_AA:
            continue
        ca = next((a for a in r if a.name == "CA"), None)
        if ca is None:
            continue
        pos = np.array([ca.pos.x, ca.pos.y, ca.pos.z], dtype=float)
        residues.append((chain.name, r.seqid.num, r.name, pos))
    return ProteinCA(residues=residues)


def extract_protein_ca_all(st: gemmi.Structure) -> "ProteinCAByChain":
    """Return Cα coordinates of every polymer chain (>= 20 residues), grouped
    by chain name and preserving residue order within each chain.

    Used by `superpose_all_chains` to build a multi-chain alignment that
    answers "did the model place both proteins correctly?", in contrast to
    `extract_protein_ca_near` which picks just the ligand-binding chain.
    """
    from loaders import STANDARD_AA
    out: dict[str, list[tuple[str, int, str, np.ndarray]]] = {}
    for ch in st[0]:
        residues: list[tuple[str, int, str, np.ndarray]] = []
        for r in ch:
            if r.name not in STANDARD_AA:
                continue
            ca = next((a for a in r if a.name == "CA"), None)
            if ca is None:
                continue
            residues.append((
                ch.name, r.seqid.num, r.name,
                np.array([ca.pos.x, ca.pos.y, ca.pos.z], dtype=float),
            ))
        if len(residues) >= 20:
            out[ch.name] = residues
    return ProteinCAByChain(by_chain=out)


@dataclass
class ProteinCAByChain:
    by_chain: dict[str, list[tuple[str, int, str, np.ndarray]]]

    @property
    def chain_names(self) -> list[str]:
        return sorted(self.by_chain.keys())


def superpose_all_chains(
    pred: "ProteinCAByChain", native: "ProteinCAByChain",
) -> Superposition:
    """Multi-chain Cα superposition.

    Pair chains by sorted name (predicted-`A` ↔ native-`A`, etc.). Within each
    paired chain, pair Cα by sequential index. Concatenate all pairs and run
    one global SVD fit. This is the "align two proteins" comparison metric:
    homo-multimer ligand placements that look fine under single-chain pocket
    alignment may diverge here when inter-chain orientation is wrong.

    If `pred` and `native` share no chain names (chain IDs renumbered), fall
    back to pairing by sorted index — pred chain[0] ↔ nat chain[0], etc.
    """
    pred_names = pred.chain_names
    nat_names = native.chain_names
    shared = [c for c in pred_names if c in nat_names]
    if shared:
        pairs = [(c, c) for c in shared]
    else:
        # Different naming: pair by sorted index, limited to min length.
        k = min(len(pred_names), len(nat_names))
        pairs = list(zip(pred_names[:k], nat_names[:k]))

    pred_xyz_all: list[np.ndarray] = []
    nat_xyz_all: list[np.ndarray] = []
    for pc, nc in pairs:
        pr = pred.by_chain[pc]
        nr = native.by_chain[nc]
        n = min(len(pr), len(nr))
        pred_xyz_all.extend(r[3] for r in pr[:n])
        nat_xyz_all.extend(r[3] for r in nr[:n])

    if len(pred_xyz_all) < 10:
        raise RuntimeError(f"too few Cα pairs ({len(pred_xyz_all)}) across all chains")
    pred_xyz = np.asarray(pred_xyz_all)
    nat_xyz = np.asarray(nat_xyz_all)
    sup = SVDSuperimposer()
    sup.set(nat_xyz, pred_xyz)
    sup.run()
    rot, tran = sup.get_rotran()
    return Superposition(
        R=np.asarray(rot).T, t=np.asarray(tran),
        rmsd=float(sup.get_rms()), n_paired=len(pred_xyz_all),
    )


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


def _bestfit_rmsd(
    crystal_mol: Chem.Mol, pred_heavy_xyz: np.ndarray, pred_mol: Chem.Mol,
) -> float:
    """Symmetry-corrected Kabsch RMSD between predicted and crystal ligand.

    Uses the SAME symmetry-aware matcher as `_matched_rmsd` (rdkit substructure
    enumeration), but for each candidate atom mapping performs an INDEPENDENT
    rigid-body Kabsch superposition of the predicted heavy coords onto the
    crystal heavy coords. Returns the minimum RMSD over correspondences.

    This is the paper-style "conformation-only" RMSD: ignores where the ligand
    sits in the protein, asks only "is the internal heavy-atom geometry right?"
    Returns NaN if matching fails on both directions.
    """
    crystal_heavy = Chem.RemoveHs(crystal_mol)
    crystal_xyz = _heavy_coords(crystal_heavy)
    if pred_mol is None:
        n = min(len(crystal_xyz), len(pred_heavy_xyz))
        return _kabsch_rmsd(pred_heavy_xyz[:n], crystal_xyz[:n])

    pred_heavy = Chem.RemoveHs(pred_mol)
    pred_xyz = _heavy_coords(pred_heavy)
    n_match_atoms = crystal_heavy.GetNumAtoms()
    matches = pred_heavy.GetSubstructMatches(
        crystal_heavy, useChirality=False, uniquify=False, maxMatches=200,
    )
    if not matches:
        # Bond-order asymmetry: try the other direction; this returns matches
        # in crystal indexing, so the trivial mapping is the identity.
        if crystal_heavy.GetSubstructMatches(
            pred_heavy, useChirality=False, uniquify=False, maxMatches=1
        ):
            matches = [tuple(range(n_match_atoms))]
    if not matches or pred_xyz.shape[0] != n_match_atoms or crystal_xyz.shape[0] != n_match_atoms:
        n = min(len(crystal_xyz), len(pred_heavy_xyz))
        return _kabsch_rmsd(pred_heavy_xyz[:n], crystal_xyz[:n])

    best = float("inf")
    for m in matches:
        try:
            sel_pred = pred_xyz[list(m)]
        except IndexError:
            continue
        r = _kabsch_rmsd(sel_pred, crystal_xyz)
        if r < best:
            best = r
    return best if best != float("inf") else float("nan")


def _kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """Minimum RMSD of P onto Q under a rigid-body (rotation + translation)."""
    if P.shape != Q.shape or len(P) < 3:
        return float("nan")
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    H = Pc.T @ Qc
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    fitted = Pc @ R.T
    diff = fitted - Qc
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


# ---------------------------------------------------------------------------
# Per-prediction pipeline
# ---------------------------------------------------------------------------

def _read_affinity(spec: dict, prefix: str, sys_dir: Path) -> dict:
    """Return Boltz-2 binding-affinity fields for this cell. Empty/None if the
    model doesn't predict affinity (AF3) or the sidecar is absent.

    Affinity is a per-system quantity in Boltz-2 — it does NOT vary across
    diffusion samples, so the same values are broadcast to every pose record.
    """
    out: dict = {"affinity_pred_value": None,
                 "affinity_probability_binary": None}
    if "affinity" not in spec:
        return out
    aff_path = sys_dir / spec["affinity"].format(prefix=prefix)
    if not aff_path.exists():
        return out
    try:
        d = json.loads(aff_path.read_text())
    except Exception:
        return out
    for k in ("affinity_pred_value", "affinity_probability_binary"):
        if k in d:
            out[k] = float(d[k])
    return out


def _read_confidence(spec: dict, prefix: str, rank: int, sys_dir: Path) -> dict:
    """Return a dict of confidence fields for one pose. All keys may be None.

    Looks for spec["conf"] formatted with prefix+rank first, then
    spec.get("conf_legacy") formatted with prefix only (rank-0 only).
    """
    out: dict = {"confidence_score": None, "iptm": None, "ptm": None,
                 "ligand_iptm": None, "complex_plddt": None,
                 "ranking_score": None}
    conf_path = sys_dir / spec["conf"].format(prefix=prefix, rank=rank)
    if not conf_path.exists() and rank == 0 and "conf_legacy" in spec:
        conf_path = sys_dir / spec["conf_legacy"].format(prefix=prefix)
    if not conf_path.exists():
        return out
    try:
        d = json.loads(conf_path.read_text())
    except Exception:
        return out
    for k in ("confidence_score", "iptm", "ptm", "ligand_iptm", "complex_plddt"):
        if k in d:
            out[k] = float(d[k])
    if "ranking_score" in d:
        out["ranking_score"] = float(d["ranking_score"])
    return out


def _analyze_single_pose(
    pdbid: str, variant: str, model: str, pose_idx: int,
    cif_path: Path, spec: dict, sys_dir: Path,
    *, smiles_override: str | None = None,
) -> PredictionRecord:
    pdbid = pdbid.lower()
    prefix = f"{pdbid}_{variant}"
    rec = PredictionRecord(pdbid=pdbid, variant=variant, model=model,
                           pose_idx=pose_idx)
    try:
        if not cif_path.exists():
            rec.status = "missing_cif"
            return rec

        # Native (crystal) loading
        crystal_pdb = CASF_RAW / pdbid / f"{pdbid}_protein.pdb"
        crystal_sdf = CASF_LIGANDS / f"{pdbid}_ligand.sdf"
        st_native = read_structure(crystal_pdb)
        crystal_mol = _crystal_ligand_mol(crystal_sdf)
        crystal_lig_xyz = _heavy_coords(crystal_mol)

        # Predicted loading
        st_pred = read_structure(cif_path)

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

        # Pick the chain CLOSEST to the ligand on each side, not the largest
        # chain. CASF includes homo-multimers (e.g. 4w9l C/K/L trimer) where
        # crystal vs predicted may put the ligand on a different homologous
        # chain. Without this, 4w9l-WT shows 47 Å ligand RMSD on a fold that
        # was actually predicted correctly (Cα 0.66 Å).
        pred_lig_heavy = pred_lig_block.heavy_coords
        native_ca = extract_protein_ca_near(st_native, crystal_lig_xyz)
        pred_ca = extract_protein_ca_near(st_pred, pred_lig_heavy)

        # Cα superpose: we want pred → native, then apply same R,t to ligand
        sup = superpose_by_index(pred_ca, native_ca)
        rec.ca_rmsd_a = round(sup.rmsd, 3)
        rec.n_ca_paired = sup.n_paired

        # Transform predicted ligand heavy coords into native frame
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

        # Variant 2: full-protein (all-chains) Cα alignment then ligand RMSD.
        try:
            pred_all = extract_protein_ca_all(st_pred)
            nat_all = extract_protein_ca_all(st_native)
            sup_all = superpose_all_chains(pred_all, nat_all)
            rec.ca_rmsd_fullca_a = round(sup_all.rmsd, 3)
            rec.n_ca_paired_fullca = sup_all.n_paired
            pred_lig_all_frame = sup_all.apply(pred_lig_heavy)
            if (pred_mol is not None
                    and pred_mol.GetNumHeavyAtoms() == len(pred_lig_all_frame)):
                m2 = Chem.Mol(pred_mol)
                cf2 = m2.GetConformer()
                for i, hidx in enumerate(_heavy_indices(m2)):
                    p = pred_lig_all_frame[i]
                    cf2.SetAtomPosition(hidx, (float(p[0]), float(p[1]), float(p[2])))
                transformed_all = m2
            else:
                transformed_all = None
            rmsd_all, _ = _matched_rmsd(
                crystal_mol, pred_lig_all_frame, transformed_all,
            )
            rec.ligand_rmsd_fullca_a = round(rmsd_all, 3)
        except Exception:
            pass  # leave the columns at None

        # Variant 3: paper-style symmetry-Kabsch RMSD on the ligand alone.
        # Uses untransformed pred_lig_heavy — bestfit is pocket-blind.
        try:
            rec.bestfit_rmsd_a = round(
                _bestfit_rmsd(crystal_mol, pred_lig_heavy, pred_mol), 3
            )
        except Exception:
            pass
    except Exception as exc:
        rec.status = "error"
        rec.error = f"{type(exc).__name__}: {exc}"

    # Merge confidence sidecar values (independent of RMSD success — even a
    # failed-RMSD cell can report its model's reported confidence).
    conf = _read_confidence(spec, prefix, pose_idx, sys_dir)
    for k, v in conf.items():
        setattr(rec, k, v)
    # Merge Boltz-2 affinity sidecar (None for AF3 — _read_affinity short-
    # circuits when spec has no "affinity" key).
    aff = _read_affinity(spec, prefix, sys_dir)
    for k, v in aff.items():
        setattr(rec, k, v)
    return rec


def analyze_predictions(
    pdbid: str, variant: str, model: str,
    *, smiles_override: str | None = None,
) -> list[PredictionRecord]:
    """Enumerate every available diffusion sample for this cell. Returns a
    list with at least one record (missing_cif if no files at all)."""
    pdbid = pdbid.lower()
    spec = MODEL_FILES[model]
    sys_dir = OUTPUT_ROOT / pdbid / variant
    prefix = f"{pdbid}_{variant}"
    cifs = sorted(sys_dir.glob(spec["cif_glob"].format(prefix=prefix)))
    if not cifs:
        return [PredictionRecord(pdbid=pdbid, variant=variant, model=model,
                                 pose_idx=0, status="missing_cif")]
    records: list[PredictionRecord] = []
    for rank, cif in enumerate(cifs):
        rec = _analyze_single_pose(
            pdbid, variant, model, rank, cif, spec, sys_dir,
            smiles_override=smiles_override,
        )
        records.append(rec)
    return records


def analyze_prediction(
    pdbid: str, variant: str, model: str,
    *, smiles_override: str | None = None,
) -> PredictionRecord:
    """Back-compat wrapper: return only the top-ranked (rank-0) pose."""
    return analyze_predictions(
        pdbid, variant, model, smiles_override=smiles_override,
    )[0]


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


# ---------------------------------------------------------------------------
# Paired WT-vs-adversarial framing
# ---------------------------------------------------------------------------

@dataclass
class PairedRecord:
    pdbid: str
    model: str
    variant: str          # rem | pack | inv
    pose_idx: int         # rank within the adversarial-variant cell
    wt_rmsd_a: float | None
    adv_rmsd_a: float | None
    delta_rmsd_a: float | None
    wt_correct_2A: bool | None
    memorized_given_wt: bool | None  # True iff wt_correct_2A and adv_rmsd < 2 Å


def paired_stats(
    records: list[PredictionRecord],
    pose_selector: str = "top1",      # "top1" = rank-0; "oracle" = min-rmsd
) -> list[PairedRecord]:
    """Per (pdbid, model) join WT ↔ each adversarial variant; report Δ RMSD
    and conditional memorization. Only `ok`-status records contribute.
    `pose_selector` picks which pose per cell to use as the comparison
    point — "top1" is the model's own rank-0, "oracle" picks the
    minimum-RMSD pose (an upper bound on what best-of-N could achieve).
    """
    def _pick(cell: list[PredictionRecord]) -> PredictionRecord | None:
        ok = [r for r in cell if r.status == "ok" and r.ligand_rmsd_a is not None]
        if not ok:
            return None
        if pose_selector == "oracle":
            return min(ok, key=lambda r: r.ligand_rmsd_a)
        return min(ok, key=lambda r: r.pose_idx)  # top1 = rank 0

    cells: dict[tuple[str, str, str], list[PredictionRecord]] = {}
    for r in records:
        cells.setdefault((r.pdbid, r.model, r.variant), []).append(r)
    selected: dict[tuple[str, str, str], PredictionRecord | None] = {
        k: _pick(v) for k, v in cells.items()
    }

    out: list[PairedRecord] = []
    for (pdbid, model, variant), adv in selected.items():
        if variant == "wt":
            continue
        wt = selected.get((pdbid, model, "wt"))
        wt_rmsd = wt.ligand_rmsd_a if wt is not None else None
        adv_rmsd = adv.ligand_rmsd_a if adv is not None else None
        delta = (adv_rmsd - wt_rmsd) if (wt_rmsd is not None and adv_rmsd is not None) else None
        wt_ok = (wt_rmsd < 2.0) if wt_rmsd is not None else None
        mem = (wt_ok and adv_rmsd < 2.0) if (wt_ok is not None and adv_rmsd is not None) else None
        out.append(PairedRecord(
            pdbid=pdbid, model=model, variant=variant,
            pose_idx=(adv.pose_idx if adv is not None else 0),
            wt_rmsd_a=wt_rmsd, adv_rmsd_a=adv_rmsd, delta_rmsd_a=delta,
            wt_correct_2A=wt_ok, memorized_given_wt=mem,
        ))
    return out


# ---------------------------------------------------------------------------
# Bootstrap CIs on memorization rate
# ---------------------------------------------------------------------------

@dataclass
class AffinityPairedRecord:
    pdbid: str
    model: str               # "Boltz2" (the only model predicting affinity)
    variant: str             # rem | pack | inv
    wt_affinity: float | None        # log[IC50] µM
    adv_affinity: float | None
    delta_affinity: float | None     # adv − wt: positive = recognized
    wt_probability: float | None     # P(binder)
    adv_probability: float | None
    delta_probability: float | None  # adv − wt: negative = recognized


def affinity_paired_stats(
    records: list[PredictionRecord],
) -> list[AffinityPairedRecord]:
    """Per (pdbid, model) join WT ↔ each adversarial variant on affinity
    fields. Only records where `affinity_pred_value is not None` contribute
    (i.e. Boltz-2 cells with the affinity sidecar present).

    Affinity is a per-system quantity — we use pose_idx==0 to dedupe
    (every pose record carries the same affinity values for that cell).

    Interpretation:
      - WT (true binder): low affinity_pred_value, high probability.
      - Adversarial (broken pocket): biophysically should have HIGHER
        affinity_pred_value and LOWER probability.
      - delta_affinity > 0 ⇒ model recognized the perturbation.
      - delta_affinity ≈ 0 ⇒ model memorized the affinity.
    """
    # Filter to one record per (pdbid, model, variant) on pose 0 with an
    # affinity value.
    by_key: dict[tuple[str, str, str], PredictionRecord] = {}
    for r in records:
        if r.pose_idx != 0 or r.affinity_pred_value is None:
            continue
        by_key[(r.pdbid, r.model, r.variant)] = r

    out: list[AffinityPairedRecord] = []
    seen_pairs: set[tuple[str, str, str]] = set()
    for (pdbid, model, variant), adv in by_key.items():
        if variant == "wt":
            continue
        wt = by_key.get((pdbid, model, "wt"))
        adv_aff = adv.affinity_pred_value
        adv_prob = adv.affinity_probability_binary
        wt_aff = wt.affinity_pred_value if wt is not None else None
        wt_prob = wt.affinity_probability_binary if wt is not None else None
        d_aff = (adv_aff - wt_aff) if (wt_aff is not None and adv_aff is not None) else None
        d_prob = (adv_prob - wt_prob) if (wt_prob is not None and adv_prob is not None) else None
        out.append(AffinityPairedRecord(
            pdbid=pdbid, model=model, variant=variant,
            wt_affinity=wt_aff, adv_affinity=adv_aff, delta_affinity=d_aff,
            wt_probability=wt_prob, adv_probability=adv_prob,
            delta_probability=d_prob,
        ))
        seen_pairs.add((pdbid, model, variant))
    return out


def bootstrap_memorization_ci(
    records: list[PredictionRecord],
    threshold_a: float = 2.0,
    n_boot: int = 1000,
    seed: int = 42,
) -> dict[tuple[str, str], tuple[float, float, float]]:
    """For each (model, adversarial-variant), return (point_rate, lo95, hi95).

    Bootstrap resamples *pdbids* (not individual cells) so the unit of
    independence matches the experimental design — three variants per system
    are not independent. Only rank-0 (top-1-by-confidence) poses are used.
    """
    keyed: dict[tuple[str, str], dict[str, float]] = {}
    for r in records:
        if r.pose_idx != 0 or r.status != "ok" or r.ligand_rmsd_a is None:
            continue
        if r.variant == "wt":
            continue
        keyed.setdefault((r.model, r.variant), {})[r.pdbid] = r.ligand_rmsd_a

    rng = np.random.default_rng(seed)
    out: dict[tuple[str, str], tuple[float, float, float]] = {}
    for key, by_pdb in keyed.items():
        pdbids = list(by_pdb.keys())
        rmsds = np.asarray([by_pdb[p] for p in pdbids])
        if len(rmsds) == 0:
            continue
        point = float((rmsds < threshold_a).mean())
        n = len(pdbids)
        if n < 2:
            out[key] = (point, point, point)
            continue
        idx = rng.integers(0, n, size=(n_boot, n))
        boot_rates = (rmsds[idx] < threshold_a).mean(axis=1)
        lo, hi = np.percentile(boot_rates, [2.5, 97.5])
        out[key] = (point, float(lo), float(hi))
    return out
