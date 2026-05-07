"""Render AlphaFold3 fold_input.json (local 'alphafold3' dialect)."""
from __future__ import annotations
import json
from pathlib import Path

# Chain ID alphabet for ligand IDs that don't collide with protein chains.
_LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def _next_id(used: set[str]) -> str:
    for c in _LETTERS:
        if c not in used:
            return c
    raise RuntimeError("ran out of single-letter chain IDs")


def render_af3(
    name: str,
    chain_seqs: list[tuple[str, str]],   # [(chain_id, sequence), ...]
    ligand_smiles: str,
    out_path: Path | str,
    seed: int = 1,
    *,
    chain_msas: dict[str, str] | None = None,
) -> None:
    """Write an AF3 local-dialect job input.

    Format matches contrasCF/Cofolding-Tools-main/examples/af_input/fold_input.json
    (top-level dict with `dialect: alphafold3`, `version: 2`).

    Args:
        chain_msas: optional {chain_id → A3M string} mapping. The query row of
            each A3M MUST match the corresponding protein sequence; the caller
            is responsible for rewriting the query row when the sequence has
            been mutated relative to the MSA's original target. When `None`
            (or a chain is missing from the dict), single-sequence mode is
            used for that chain — the query becomes the only row.
    """
    used = {cid for cid, _ in chain_seqs}
    lig_id = _next_id(used)
    proteins = []
    for cid, seq in chain_seqs:
        a3m = (chain_msas or {}).get(cid)
        if a3m is None:
            a3m = f">query\n{seq}\n"
        proteins.append({"protein": {
            "id": cid,
            "sequence": seq,
            "unpairedMsa": a3m,
            "pairedMsa": "",
            "templates": [],
        }})
    payload = {
        "name": name,
        "sequences": [
            *proteins,
            {"ligand": {"id": lig_id, "smiles": ligand_smiles}},
        ],
        "modelSeeds": [seed],
        "dialect": "alphafold3",
        "version": 2,
    }
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2))


def clean_a3m_for_af3(a3m: str, expected_query_len: int) -> str:
    """Drop MSA records whose concatenated sequence has the wrong AF3 shape.

    AF3's `MSA.featurize` requires every record's sequence to have exactly
    `query_length` non-lowercase characters (uppercase letters + '-' + 'X';
    lowercase letters are A3M insertions and don't count). Boltz/MMseqs2
    sometimes returns hits whose non-lowercase count exceeds the query
    length (e.g. a hit against a tandem-duplicated protein where two regions
    of the hit align to the query, encoded as one record of length ~2*Q).
    Such rows trigger:
        Error: Invalid shape — all strings must have the same number of
        non-lowercase characters

    Implementation: parse A3M as records (header + multi-line sequence;
    sequence terminates at the next `>` or EOF), strip null bytes from each
    sequence, drop records whose non-lowercase count != Q. Always keep the
    first record (the query). The output emits each record on two lines:
    header, then sequence — which is also a valid A3M.
    """
    a3m = a3m.replace("\x00", "")
    records: list[tuple[str, str]] = []
    cur_hdr: str | None = None
    cur_seq_parts: list[str] = []
    for ln in a3m.splitlines():
        if ln.startswith(">"):
            if cur_hdr is not None:
                records.append((cur_hdr, "".join(cur_seq_parts)))
            cur_hdr = ln
            cur_seq_parts = []
        elif ln:
            cur_seq_parts.append(ln)
    if cur_hdr is not None:
        records.append((cur_hdr, "".join(cur_seq_parts)))

    out_lines: list[str] = []
    for i, (hdr, seq) in enumerate(records):
        if i == 0:
            # Always keep the query, regardless of shape (it should match Q).
            out_lines.append(hdr)
            out_lines.append(seq)
            continue
        non_lower = sum(1 for c in seq if not c.islower())
        if non_lower == expected_query_len:
            out_lines.append(hdr)
            out_lines.append(seq)
    return "\n".join(out_lines) + ("\n" if a3m.endswith("\n") else "")


def rewrite_a3m_query(a3m: str, new_query_seq: str) -> str:
    """Replace the FIRST sequence row of an A3M with `new_query_seq`.

    The MSA columns are preserved (we do NOT re-align the query); we only
    swap the query letters. This is appropriate for variants of the WT
    sequence that differ at a few positions — the MSA's column structure
    still reflects the conserved core, and AF3 sees the variant residues as
    non-conserved positions in the alignment.

    The query row's length must equal `new_query_seq` length; if not, raises.
    """
    lines = a3m.splitlines()
    out_lines = list(lines)
    # find first header row (must exist by A3M spec)
    for i, ln in enumerate(lines):
        if ln.startswith(">"):
            # the next non-empty, non-header line is the query
            for j in range(i + 1, len(lines)):
                if lines[j] and not lines[j].startswith(">"):
                    if len(lines[j]) != len(new_query_seq):
                        raise ValueError(
                            f"a3m query length {len(lines[j])} != "
                            f"new sequence length {len(new_query_seq)}"
                        )
                    out_lines[j] = new_query_seq
                    return "\n".join(out_lines) + ("\n" if a3m.endswith("\n") else "")
            break
    raise ValueError("a3m has no records")
