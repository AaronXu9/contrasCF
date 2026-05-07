"""Per-system build pipeline for ligand-side perturbations.

For one CASF PDB id:
  1. Load the crystal protein + ligand.
  2. Resolve the WT ligand SMILES (canonical, RDKit-clean).
  3. Run each chemistry rule (methylation / charge_swap / halogenation),
     recording per-rule eligibility and generated variants.
  4. Emit AF3 JSON + Boltz-2 YAML + docking inputs for the WT and every
     generated variant. The protein is constant across variants in this
     module, so receptor.pdb / box.json are written once for `wt` and
     copied (hardlinked) for the rest.
"""
from __future__ import annotations
import json
import shutil
import traceback
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import gemmi

from casf_mutagenesis.inputs_af3 import render_af3
from casf_mutagenesis.inputs_boltz import render_boltz
from casf_mutagenesis.inputs_docking import _embed_ligand, _strip_to_protein_pdb
from casf_mutagenesis.pocket import load_ligand_heavy_coords
from casf_mutagenesis.sequence import extract_chain_sequences

from .config import CASF_LIGANDS, CASF_RAW, DOCKING_BOX_A, OUTPUT_ROOT
from .ligand_io import data_dict, resolve_wt_smiles
from .rules.base import LigandRule, LigandVariant
from .rules.charge_swap import ChargeSwapRule
from .rules.halogenation import HalogenationRule
from .rules.methylation import MethylationRule


_RULES: list[LigandRule] = [
    MethylationRule(),
    ChargeSwapRule(),
    HalogenationRule(),
]


@dataclass
class SystemEntry:
    pdbid: str
    status: str = "ok"               # "ok" | "error"
    error: str | None = None
    het_code: str | None = None
    wt_smiles: str | None = None
    wt_smiles_source: str | None = None
    n_chains: int = 0
    seq_lengths: list[int] = field(default_factory=list)
    rules_eligible: dict[str, dict] = field(default_factory=dict)
    variants: dict[str, dict] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)


def build_system(pdbid: str, output_root: Path = OUTPUT_ROOT) -> SystemEntry:
    pdbid = pdbid.lower()
    entry = SystemEntry(pdbid=pdbid)
    out_base = output_root / pdbid

    try:
        protein_pdb = CASF_RAW / pdbid / f"{pdbid}_protein.pdb"
        ligand_sdf = CASF_LIGANDS / f"{pdbid}_ligand.sdf"
        if not protein_pdb.exists():
            raise FileNotFoundError(f"missing {protein_pdb}")
        if not ligand_sdf.exists():
            raise FileNotFoundError(f"missing {ligand_sdf}")

        meta = data_dict().get(pdbid, {})
        het = meta.get("ligand_name")
        entry.het_code = het

        wt_smiles, src = resolve_wt_smiles(ligand_sdf, het)
        if wt_smiles is None:
            entry.warnings.append(
                f"no SMILES from crystal SDF or het cache (HET '{het}')"
            )
            entry.status = "error"
            entry.error = "missing WT SMILES"
            return entry
        entry.wt_smiles = wt_smiles
        entry.wt_smiles_source = src

        chains = extract_chain_sequences(protein_pdb)
        entry.n_chains = len(chains)
        entry.seq_lengths = [len(recs) for _, recs in chains]
        wt_seqs = [(cid, "".join(r.aa1 for r in recs)) for cid, recs in chains]

        # WT variant is always present.
        wt_variant = LigandVariant(
            name="wt",
            rule="",                # empty — wt is not the product of any rule
            smiles=wt_smiles,
            applied=[],
            paper_faithful=True,
        )
        all_variants: list[LigandVariant] = [wt_variant]

        # Run rules.
        from rdkit import Chem
        wt_mol = Chem.MolFromSmiles(wt_smiles)
        if wt_mol is None:
            entry.warnings.append(f"RDKit could not parse WT SMILES: {wt_smiles}")
            entry.status = "error"
            entry.error = "RDKit parse failure on WT SMILES"
            return entry
        for rule in _RULES:
            eligible, reason, n_matches = rule.is_eligible(wt_mol)
            entry.rules_eligible[rule.name] = {
                "eligible": eligible,
                "reason": reason,
                "n_matches": n_matches,
            }
            if not eligible:
                continue
            try:
                rule_variants = rule.generate(wt_mol)
            except Exception as exc:
                entry.warnings.append(f"{rule.name} generate failed: {exc}")
                continue
            all_variants.extend(rule_variants)

        # Render WT docking inputs once (receptor + box are constant per system).
        wt_dock_dir = out_base / "wt" / "docking"
        wt_dock_dir.mkdir(parents=True, exist_ok=True)
        try:
            st = gemmi.read_structure(str(protein_pdb))
            (wt_dock_dir / "receptor.pdb").write_text(_strip_to_protein_pdb(st))
            lig_xyz = load_ligand_heavy_coords(ligand_sdf)
            centre = lig_xyz.mean(axis=0)
            (wt_dock_dir / "box.json").write_text(json.dumps({
                "center": [float(centre[0]), float(centre[1]), float(centre[2])],
                "size": list(DOCKING_BOX_A),
            }, indent=2))
        except Exception as exc:
            entry.warnings.append(f"WT docking shared inputs failed: {exc}")
            wt_dock_dir = None

        # Per-variant render.
        for v in all_variants:
            v_dir = out_base / v.name
            v_dir.mkdir(parents=True, exist_ok=True)
            try:
                render_af3(
                    name=f"{pdbid}_{v.name}",
                    chain_seqs=wt_seqs,
                    ligand_smiles=v.smiles,
                    out_path=v_dir / "af3.json",
                )
            except Exception as exc:
                entry.warnings.append(f"{v.name} af3 render failed: {exc}")
            try:
                render_boltz(
                    chain_seqs=wt_seqs,
                    ligand_smiles=v.smiles,
                    out_path=v_dir / "boltz.yaml",
                )
            except Exception as exc:
                entry.warnings.append(f"{v.name} boltz render failed: {exc}")

            v_dock_dir = v_dir / "docking"
            v_dock_dir.mkdir(parents=True, exist_ok=True)
            if wt_dock_dir is not None:
                # Reuse receptor + box from WT.
                _link_or_copy(wt_dock_dir / "receptor.pdb", v_dock_dir / "receptor.pdb")
                _link_or_copy(wt_dock_dir / "box.json", v_dock_dir / "box.json")
            try:
                _embed_ligand(v.smiles, v_dock_dir / "ligand.sdf")
            except Exception as exc:
                entry.warnings.append(f"{v.name} ligand embed failed: {exc}")

            entry.variants[v.name] = {
                "rule": v.rule,
                "smiles": v.smiles,
                "applied": list(v.applied),
                "paper_faithful": v.paper_faithful,
                "out_dir": str(v_dir),
            }
    except Exception as exc:
        entry.status = "error"
        entry.error = f"{type(exc).__name__}: {exc}"
        entry.warnings.append(traceback.format_exc().splitlines()[-1])
    return entry


def _link_or_copy(src: Path, dst: Path) -> None:
    if dst.exists():
        return
    try:
        dst.hardlink_to(src)
    except (OSError, AttributeError):
        shutil.copyfile(src, dst)


def write_manifest(entries: list[SystemEntry], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {"systems": [asdict(e) for e in entries]}
    path.write_text(json.dumps(payload, indent=2))
