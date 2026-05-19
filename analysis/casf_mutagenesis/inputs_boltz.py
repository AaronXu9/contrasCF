"""Render Boltz YAML input."""
from __future__ import annotations
from pathlib import Path

_LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def _next_id(used: set[str]) -> str:
    for c in _LETTERS:
        if c not in used:
            return c
    raise RuntimeError("ran out of single-letter chain IDs")


def render_boltz(
    chain_seqs: list[tuple[str, str]],   # [(chain_id, sequence), ...]
    ligand_smiles: str,
    out_path: Path | str,
) -> None:
    """Write a Boltz YAML matching the Cofolding-Tools template format.

    YAML written by hand (no PyYAML dependency); all sequence/SMILES strings
    are emitted unquoted-safe (no special chars in 1-letter AA codes; SMILES
    is wrapped on a single line and contains no YAML-significant chars).
    """
    used = {cid for cid, _ in chain_seqs}
    lig_id = _next_id(used)
    lines = ["version: 1", "sequences:"]
    for cid, seq in chain_seqs:
        lines += [
            "  - protein:",
            f"      id: [{cid}]",
            f"      sequence: {seq}",
            "      msa: empty",
        ]
    lines += [
        "  - ligand:",
        f"      id: [{lig_id}]",
        f"      smiles: {ligand_smiles}",
    ]
    # Request Boltz-2 binding-affinity prediction for the ligand. This produces
    # `affinity_<prefix>.json` alongside the structure outputs, containing
    # `affinity_pred_value` (log[IC50] µM; lower = tighter) and
    # `affinity_probability_binary` (P(binder) ∈ [0, 1]). The runner already
    # copies the affinity sidecar into the variant directory.
    lines += [
        "properties:",
        "  - affinity:",
        f"      binder: {lig_id}",
    ]
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n")
