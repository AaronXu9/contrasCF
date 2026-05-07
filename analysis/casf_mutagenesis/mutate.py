"""Apply pocket residue mutations to per-chain sequences."""
from __future__ import annotations
from dataclasses import dataclass

from .config import AA3_TO_AA1, INV_TABLE, PACK_AA1, REM_AA1
from .pocket import PocketResidue
from .sequence import ResidueRecord


@dataclass(frozen=True)
class MutationApplied:
    chain: str
    resnum: int
    ins_code: str
    wt_aa1: str
    mut_aa1: str

    def code(self) -> str:
        ins = self.ins_code or ""
        return f"{self.wt_aa1}{self.resnum}{ins}{self.mut_aa1}"


def _new_aa(variant: str, wt_aa3: str) -> str:
    if variant == "rem":
        return REM_AA1
    if variant == "pack":
        return PACK_AA1
    if variant == "inv":
        return INV_TABLE[wt_aa3]
    raise ValueError(f"unknown variant: {variant}")


def apply_mutations(
    chains: list[tuple[str, list[ResidueRecord]]],
    pocket: list[PocketResidue],
    variant: str,
) -> tuple[list[tuple[str, str]], list[MutationApplied]]:
    """Return ([(chain_id, mutated_sequence_string)], [MutationApplied, ...]).

    `variant` is one of {"wt", "rem", "pack", "inv"}.
    For "wt" the chain sequences are returned unchanged and `applied` is empty.

    Pocket entries that point at a residue absent from the chain records are
    silently skipped — this should not happen if `pocket` came from the same
    structure, but is defensive against external inputs.
    """
    chains_dict = {cid: list(recs) for cid, recs in chains}
    applied: list[MutationApplied] = []
    if variant != "wt":
        # index recs per chain by (resnum, ins_code) for O(1) lookup
        for p in pocket:
            recs = chains_dict.get(p.chain)
            if recs is None:
                continue
            for i, r in enumerate(recs):
                if r.resnum == p.resnum and r.ins_code == p.ins_code:
                    wt_aa1 = AA3_TO_AA1.get(p.aa3, "X")
                    new_aa1 = _new_aa(variant, p.aa3)
                    if new_aa1 == wt_aa1:
                        # identity substitution (shouldn't happen with current
                        # tables since INV_TABLE has no fixed points)
                        break
                    recs[i] = ResidueRecord(
                        resnum=r.resnum, ins_code=r.ins_code,
                        aa3=p.aa3, aa1=new_aa1,
                    )
                    applied.append(MutationApplied(
                        chain=p.chain, resnum=p.resnum, ins_code=p.ins_code,
                        wt_aa1=wt_aa1, mut_aa1=new_aa1,
                    ))
                    break
    seqs = [(cid, "".join(r.aa1 for r in chains_dict[cid]))
            for cid, _ in chains]
    return seqs, applied
