"""Extract per-chain protein sequence + PDB residue-number mapping.

Returns a list of (chain_id, [ResidueRecord, ...]) preserving PDB numbering and
insertion codes. Pattern follows analysis/scripts/12_run_boltz2.py:51-72 but
keeps PDB resnums alongside the 1-letter sequence (which that script discards).
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path

import gemmi

from .config import AA3_TO_AA1


@dataclass(frozen=True)
class ResidueRecord:
    resnum: int        # PDB seqid number
    ins_code: str      # PDB insertion code, '' if none
    aa3: str           # 3-letter (may be non-standard, e.g. 'MSE')
    aa1: str           # 1-letter; 'X' if unknown


def extract_chain_sequences(
    pdb_path: Path | str,
) -> list[tuple[str, list[ResidueRecord]]]:
    """Return [(chain_id, [ResidueRecord, ...])] for every protein chain.

    A residue is recorded once per (resnum, ins_code) — duplicate alt-locs are
    deduplicated by first occurrence (matches gemmi's CA-iteration behaviour).
    Non-protein chains (no Cα atoms) are skipped.
    """
    st = gemmi.read_structure(str(pdb_path))
    out: list[tuple[str, list[ResidueRecord]]] = []
    model = st[0]
    for chain in model:
        seen = set()
        recs: list[ResidueRecord] = []
        for res in chain:
            # protein-only filter: must have a CA atom
            if not any(a.name == "CA" for a in res):
                continue
            key = (res.seqid.num, res.seqid.icode.strip())
            if key in seen:
                continue
            seen.add(key)
            aa3 = res.name.strip().upper()
            aa1 = AA3_TO_AA1.get(aa3, "X")
            recs.append(ResidueRecord(
                resnum=res.seqid.num,
                ins_code=res.seqid.icode.strip(),
                aa3=aa3,
                aa1=aa1,
            ))
        if recs:
            out.append((chain.name, recs))
    return out


def sequence_string(records: list[ResidueRecord]) -> str:
    """Join 1-letter codes into a contiguous sequence string."""
    return "".join(r.aa1 for r in records)
