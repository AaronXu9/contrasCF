"""Common types for ligand-mutation rules."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Protocol, runtime_checkable

from rdkit.Chem import Mol


@dataclass(frozen=True)
class LigandVariant:
    """One generated ligand variant ready to render into AF3/Boltz/docking inputs."""
    name: str                   # e.g. "meth_3", "chrg_pos_2", "halo_F_1"
    rule: str                   # "methylation" | "charge_swap" | "halogenation"
    smiles: str                 # canonical isomeric SMILES of the variant
    applied: list[str] = field(default_factory=list)   # human-readable change log
    paper_faithful: bool = True


@runtime_checkable
class LigandRule(Protocol):
    """Protocol every chemistry rule implements."""
    name: str

    def is_eligible(self, mol: Mol) -> tuple[bool, str, int]:
        """Return (eligible, reason_if_not, n_matches).

        n_matches is the number of SMARTS hits (0 if not eligible) — used by
        the orchestrator to record per-rule eligibility statistics in the
        manifest without re-doing the SMARTS match.
        """
        ...

    def generate(self, mol: Mol) -> list[LigandVariant]:
        """Produce all ladder variants this rule yields for this molecule."""
        ...
