"""Per-system build pipeline: take a CASF PDB id, produce four variants
(wt/rem/pack/inv) of AF3 + Boltz inputs, plus WT-only docking inputs.

Returns a manifest entry dict suitable for serialisation.
"""
from __future__ import annotations
import json
import traceback
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from rdkit import Chem

from .config import (
    CASF_LIGANDS, CASF_RAW, DATA_DICT_JSON, OUTPUT_ROOT, SMILES_CACHE,
    VARIANTS,
)
from .inputs_af3 import render_af3
from .inputs_boltz import render_boltz
from .inputs_docking import render_docking_wt
from .mutate import apply_mutations
from .pocket import detect_pocket
from .sequence import extract_chain_sequences


def _smiles_from_crystal_sdf(sdf_path: Path) -> str | None:
    """Canonical SMILES from the crystal ligand SDF — preferred over the
    cached HET SMILES because it's guaranteed RDKit-parseable.

    The cached SMILES from `het_smiles_cache.json` uses lowercase aromatic
    notation that fails RDKit kekulization for some fused N-heterocycles
    (e.g. 4ih5's pyrazolopyrimidine, 4de1's tetrazole-phthalazine, 4ivc's
    cyano-fused-heterocycle). Round-tripping through the SDF (which has
    explicit bonds) → RDKit `MolToSmiles` produces a canonical, parseable
    string both Boltz-2 and AF3 accept.
    """
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=True, sanitize=False)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        # SanitizeMol may fail on exotic ligands; canonical SMILES often
        # still works without full sanitisation.
        pass
    try:
        return Chem.MolToSmiles(mol)
    except Exception:
        return None

_DATA_DICT_CACHE: dict[str, Any] | None = None
_SMILES_CACHE: dict[str, str] | None = None


def _data_dict() -> dict[str, Any]:
    global _DATA_DICT_CACHE
    if _DATA_DICT_CACHE is None:
        _DATA_DICT_CACHE = json.loads(DATA_DICT_JSON.read_text())
    return _DATA_DICT_CACHE


def _smiles_cache() -> dict[str, str]:
    global _SMILES_CACHE
    if _SMILES_CACHE is None:
        _SMILES_CACHE = json.loads(SMILES_CACHE.read_text())
    return _SMILES_CACHE


@dataclass
class SystemEntry:
    pdbid: str
    status: str = "ok"           # "ok" or "error"
    error: str | None = None
    het_code: str | None = None
    smiles: str | None = None
    smiles_source: str | None = None  # "crystal_sdf" | "het_cache" | None
    n_chains: int = 0
    seq_lengths: list[int] = field(default_factory=list)
    pocket: list[dict] = field(default_factory=list)
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

        meta = _data_dict().get(pdbid, {})
        het = meta.get("ligand_name")
        entry.het_code = het

        # Prefer SMILES derived from the crystal SDF — it's RDKit-canonical
        # and guaranteed to round-trip through RDKit (which both Boltz-2 and
        # AF3 use internally to add hydrogens). Fall back to the cached HET
        # SMILES only if the SDF can't be parsed.
        smiles = _smiles_from_crystal_sdf(ligand_sdf)
        smiles_source = "crystal_sdf"
        if smiles is None:
            smiles = _smiles_cache().get(het) if het else None
            smiles_source = "het_cache"
        if smiles is None:
            entry.warnings.append(
                f"no SMILES from crystal SDF or het cache (HET '{het}')"
            )
        entry.smiles = smiles
        entry.smiles_source = smiles_source if smiles else None

        chains = extract_chain_sequences(protein_pdb)
        entry.n_chains = len(chains)
        entry.seq_lengths = [len(recs) for _, recs in chains]

        pocket = detect_pocket(protein_pdb, ligand_sdf)
        entry.pocket = [
            {"chain": p.chain, "resnum": p.resnum, "ins": p.ins_code,
             "wt_aa3": p.aa3}
            for p in pocket
        ]
        if not pocket:
            entry.warnings.append("empty pocket — no residues within 3.5 Å sc")

        wt_seqs = [(cid, "".join(r.aa1 for r in recs)) for cid, recs in chains]
        wt_len = sum(len(s) for _, s in wt_seqs)

        for variant in VARIANTS:
            seqs, applied = apply_mutations(chains, pocket, variant)
            # invariant: all sequences same length
            assert sum(len(s) for _, s in seqs) == wt_len, "length changed"
            v_dir = out_base / variant
            v_dir.mkdir(parents=True, exist_ok=True)

            # AF3 + Boltz need a SMILES; if missing, skip those but still
            # produce the protein-only manifest entry
            if smiles is not None:
                render_af3(
                    name=f"{pdbid}_{variant}",
                    chain_seqs=seqs,
                    ligand_smiles=smiles,
                    out_path=v_dir / "af3.json",
                )
                render_boltz(
                    chain_seqs=seqs,
                    ligand_smiles=smiles,
                    out_path=v_dir / "boltz.yaml",
                )

            entry.variants[variant] = {
                "applied": [m.code() for m in applied],
                "n_mutations": len(applied),
                "out_dir": str(v_dir),
            }

        # WT-only docking inputs (mutant-protein docking is downstream)
        if smiles is not None:
            try:
                docking_dir = out_base / "wt" / "docking"
                render_docking_wt(
                    pdb_path=protein_pdb,
                    crystal_ligand_sdf=ligand_sdf,
                    target_smiles=smiles,
                    out_dir=docking_dir,
                )
                entry.variants["wt"]["docking_dir"] = str(docking_dir)
            except Exception as exc:
                entry.warnings.append(f"docking inputs failed: {exc}")
    except Exception as exc:
        entry.status = "error"
        entry.error = f"{type(exc).__name__}: {exc}"
        entry.warnings.append(traceback.format_exc().splitlines()[-1])
    return entry


def write_manifest(entries: list[SystemEntry], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {"systems": [asdict(e) for e in entries]}
    path.write_text(json.dumps(payload, indent=2))
