"""Per-model confidence extraction.

Fields returned (use NaN when not applicable):
  ptm, iptm, ligand_iptm, ligand_plddt_mean, ranking_score, has_clash_reported

AF3 summary JSON fields observed:
  ptm, iptm, ranking_score, chain_pair_iptm (2D), chain_pair_pae_min, has_clash,
  fraction_disordered
Boltz confidence JSON fields:
  confidence_score, ptm, iptm, ligand_iptm, protein_iptm, complex_plddt, ...
Chai scores.rank_0.json fields:
  aggregate_score, ptm, iptm, per_chain_ptm, per_chain_pair_iptm,
  has_inter_chain_clashes, chain_chain_clashes
RFAA: no JSON -> NaN everywhere except ligand_plddt_mean (from ligand B-factors in PDB).
"""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np

from config import MODELS, case_dir, find_first
from loaders import LigandBlock


NAN = float("nan")


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text())


def _conf_paths(model: str, case: str) -> tuple[Path | None, Path | None, Path | None]:
    """Return (summary_json, data_json, plddt_npz) where applicable."""
    d = case_dir(model, case)
    globs = MODELS[model].confidence_globs
    summary = None
    for g in globs:
        hits = sorted(d.glob(g))
        if hits:
            summary = hits[0]; break
    # AF3 also has per-atom confidences JSON
    data_json = None
    if model == "AF3":
        # Try: <prefix>_confidences.json  OR  <prefix>_confidences_0.json
        for g in ("*_confidences_0.json", "*_confidences.json"):
            hits = [p for p in d.glob(g) if "summary" not in p.name]
            if hits:
                data_json = hits[0]; break
    plddt_npz = None
    if model in ("Boltz", "Boltz2"):
        hits = sorted(d.glob("plddt_*_model_0.npz"))
        if hits:
            plddt_npz = hits[0]
    return summary, data_json, plddt_npz


def _af3_ligand_plddt(af3_conf_json: Path, lig_chain: str) -> float:
    j = _load_json(af3_conf_json)
    atom_chain = j.get("atom_chain_ids") or []
    atom_plddt = j.get("atom_plddts") or []
    if not atom_chain or not atom_plddt or len(atom_chain) != len(atom_plddt):
        return NAN
    vals = [p for ch, p in zip(atom_chain, atom_plddt) if ch == lig_chain]
    return float(np.mean(vals)) if vals else NAN


def _af3_ligand_iptm(summary: dict, lig_chain: str) -> float:
    cp = summary.get("chain_pair_iptm")
    if not cp:
        return NAN
    # chain_pair_iptm is a 2D array; index -> chain letter mapping assumed A,B,C,...
    idx = {chr(ord("A") + i): i for i in range(len(cp))}
    i_lig = idx.get(lig_chain)
    i_prot = idx.get("A")
    if i_lig is None or i_prot is None:
        return NAN
    try:
        return float((cp[i_prot][i_lig] + cp[i_lig][i_prot]) / 2)
    except Exception:
        return NAN


def extract_af3(case: str, lig: LigandBlock) -> dict:
    summary_p, conf_p, _ = _conf_paths("AF3", case)
    if summary_p is None:
        return {k: NAN for k in ("ptm","iptm","ligand_iptm","ligand_plddt_mean","ranking_score","has_clash_reported")}
    s = _load_json(summary_p)
    ligand_plddt = _af3_ligand_plddt(conf_p, lig.chain) if conf_p else NAN
    # fall back to ligand B-factor mean if no per-atom plddt JSON
    if np.isnan(ligand_plddt):
        ligand_plddt = lig.mean_bfactor
    return dict(
        ptm=float(s.get("ptm", NAN)) if s.get("ptm") is not None else NAN,
        iptm=float(s.get("iptm", NAN)) if s.get("iptm") is not None else NAN,
        ligand_iptm=_af3_ligand_iptm(s, lig.chain),
        ligand_plddt_mean=ligand_plddt,
        ranking_score=float(s.get("ranking_score", NAN)) if s.get("ranking_score") is not None else NAN,
        has_clash_reported=bool(s.get("has_clash")) if s.get("has_clash") is not None else None,
    )


def extract_boltz(case: str, lig: LigandBlock) -> dict:
    summary_p, _, plddt_p = _conf_paths("Boltz", case)
    if summary_p is None:
        return {k: NAN for k in ("ptm","iptm","ligand_iptm","ligand_plddt_mean","ranking_score","has_clash_reported")}
    s = _load_json(summary_p)
    # ligand pLDDT: prefer ligand B-factor from CIF (Boltz writes pLDDT there).
    ligand_plddt = lig.mean_bfactor
    return dict(
        ptm=float(s.get("ptm", NAN)),
        iptm=float(s.get("iptm", NAN)),
        ligand_iptm=float(s.get("ligand_iptm", NAN)) if s.get("ligand_iptm") is not None else NAN,
        ligand_plddt_mean=ligand_plddt,
        ranking_score=float(s.get("confidence_score", NAN)),
        has_clash_reported=None,
    )


def extract_boltz2(case: str, lig: LigandBlock) -> dict:
    """Boltz-2 confidence — identical JSON shape to Boltz-1, plus an adjacent
    ``affinity_*.json`` with ``affinity_pred_value`` and
    ``affinity_probability_binary``. The latter are surfaced as extra columns
    so downstream plots can show them alongside ranking_score.
    """
    summary_p, _, _ = _conf_paths("Boltz2", case)
    out = {
        "ptm": NAN, "iptm": NAN, "ligand_iptm": NAN,
        "ligand_plddt_mean": lig.mean_bfactor, "ranking_score": NAN,
        "has_clash_reported": None,
        "boltz2_affinity_pred_value": NAN,
        "boltz2_affinity_probability_binary": NAN,
    }
    if summary_p is None:
        return out
    s = _load_json(summary_p)
    out["ptm"] = float(s.get("ptm", NAN))
    out["iptm"] = float(s.get("iptm", NAN))
    out["ligand_iptm"] = (
        float(s.get("ligand_iptm", NAN)) if s.get("ligand_iptm") is not None else NAN
    )
    out["ranking_score"] = float(s.get("confidence_score", NAN))

    # Affinity JSON lives next to the confidence JSON: affinity_<prefix>.json
    d = case_dir("Boltz2", case)
    for aff_p in sorted(d.glob("affinity_*.json")):
        a = _load_json(aff_p)
        try:
            out["boltz2_affinity_pred_value"] = float(a.get("affinity_pred_value", NAN))
        except (TypeError, ValueError):
            pass
        try:
            out["boltz2_affinity_probability_binary"] = float(
                a.get("affinity_probability_binary", NAN)
            )
        except (TypeError, ValueError):
            pass
        break
    return out


def extract_chai(case: str, lig: LigandBlock) -> dict:
    summary_p, _, _ = _conf_paths("Chai", case)
    if summary_p is None:
        return {k: NAN for k in ("ptm","iptm","ligand_iptm","ligand_plddt_mean","ranking_score","has_clash_reported")}
    s = _load_json(summary_p)
    # per_chain_pair_iptm is nested 4D; Chai writes it as [[[AA,AB],[BA,BB]]].
    ligand_iptm = NAN
    pcp = s.get("per_chain_pair_iptm")
    try:
        # Try assuming [0][0][1] is A->B
        a_b = pcp[0][0][1]
        b_a = pcp[0][1][0]
        ligand_iptm = float((a_b + b_a) / 2)
    except Exception:
        pass
    return dict(
        ptm=float(s.get("ptm", NAN)),
        iptm=float(s.get("iptm", NAN)),
        ligand_iptm=ligand_iptm,
        ligand_plddt_mean=lig.mean_bfactor,  # Chai writes pLDDT into CIF B-factor
        ranking_score=float(s.get("aggregate_score", NAN)),
        has_clash_reported=bool(s.get("has_inter_chain_clashes")) if s.get("has_inter_chain_clashes") is not None else None,
    )


def extract_rfaa(case: str, lig: LigandBlock) -> dict:
    return dict(
        ptm=NAN, iptm=NAN, ligand_iptm=NAN,
        ligand_plddt_mean=lig.mean_bfactor,     # RFAA writes pLDDT in B-factor column
        ranking_score=NAN,
        has_clash_reported=None,
    )


def _read_sdf_top_tags(sdf_path: Path) -> dict:
    """Parse the tag block of the FIRST record in an SDF without rdkit."""
    text = sdf_path.read_text()
    # First record = everything up to first $$$$
    first = text.split("\n$$$$", 1)[0]
    tags: dict[str, str] = {}
    lines = first.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("> <") and ">" in line[3:]:
            key = line[3:line.index(">", 3)]
            if i + 1 < len(lines):
                tags[key] = lines[i + 1].strip()
            i += 2
        elif line.startswith(">  <") and ">" in line[4:]:
            key = line[4:line.index(">", 4)]
            if i + 1 < len(lines):
                tags[key] = lines[i + 1].strip()
            i += 2
        else:
            i += 1
    return tags


def extract_unidock2(case: str, lig: LigandBlock) -> dict:
    d = case_dir("UniDock2", case)
    sdf = d / "poses.sdf"
    out = dict(
        ptm=NAN, iptm=NAN, ligand_iptm=NAN,
        ligand_plddt_mean=NAN, ranking_score=NAN,
        has_clash_reported=None,
        dock_vina_score=NAN, dock_cnn_score=NAN, dock_cnn_affinity=NAN,
    )
    if not sdf.exists():
        return out
    tags = _read_sdf_top_tags(sdf)
    try:
        vina = float(tags.get("vina_binding_free_energy", "nan"))
    except ValueError:
        vina = NAN
    out["dock_vina_score"] = vina
    # UniDock2 has no CNN scoring. Mirror vina score into ranking_score so plots
    # that consume ranking_score see something for this model.
    out["ranking_score"] = vina
    return out


def extract_gnina(case: str, lig: LigandBlock) -> dict:
    d = case_dir("GNINA", case)
    sdf = d / "poses.sdf"
    out = dict(
        ptm=NAN, iptm=NAN, ligand_iptm=NAN,
        ligand_plddt_mean=NAN, ranking_score=NAN,
        has_clash_reported=None,
        dock_vina_score=NAN, dock_cnn_score=NAN, dock_cnn_affinity=NAN,
    )
    if not sdf.exists():
        return out
    tags = _read_sdf_top_tags(sdf)
    def _f(k: str) -> float:
        try:
            return float(tags.get(k, "nan"))
        except ValueError:
            return NAN
    out["dock_vina_score"] = _f("minimizedAffinity")
    out["dock_cnn_score"] = _f("CNNscore")
    out["dock_cnn_affinity"] = _f("CNNaffinity")
    # Use CNNaffinity as the "ranking score" for GNINA (pK units, higher = better).
    out["ranking_score"] = out["dock_cnn_affinity"]
    return out


def extract_surfdock(case: str, lig: LigandBlock) -> dict:
    """SurfDock confidence — read from SDF tags of poses.sdf top record.

    13_run_surfdock.py writes each pose record with tags:
      SurfDock_rank       (int, 1-based)
      SurfDock_confidence (float; higher = more confident per SurfScore)
      SurfDock_rmsd       (float; SurfDock's self-reported RMSD vs co-crystal
                            if one was provided — not meaningful here)
    """
    d = case_dir("SurfDock", case)
    sdf = d / "poses.sdf"
    out = dict(
        ptm=NAN, iptm=NAN, ligand_iptm=NAN,
        ligand_plddt_mean=NAN, ranking_score=NAN,
        has_clash_reported=None,
        dock_vina_score=NAN, dock_cnn_score=NAN, dock_cnn_affinity=NAN,
        surfdock_confidence=NAN,
    )
    if not sdf.exists():
        return out
    tags = _read_sdf_top_tags(sdf)
    try:
        conf = float(tags.get("SurfDock_confidence", "nan"))
    except ValueError:
        conf = NAN
    out["surfdock_confidence"] = conf
    # Mirror SurfDock confidence into ranking_score so cross-model plots that
    # consume ranking_score see something (higher = better, same convention
    # as AF3 ranking_score / Boltz confidence_score / GNINA CNNaffinity).
    out["ranking_score"] = conf
    return out


EXTRACTORS = {
    "AF3": extract_af3,
    "Boltz": extract_boltz,
    "Boltz2": extract_boltz2,
    "Chai": extract_chai,
    "RFAA": extract_rfaa,
    "UniDock2": extract_unidock2,
    "GNINA": extract_gnina,
    "SurfDock": extract_surfdock,
}


def extract_confidence(model: str, case: str, lig: LigandBlock) -> dict:
    return EXTRACTORS[model](case, lig)
