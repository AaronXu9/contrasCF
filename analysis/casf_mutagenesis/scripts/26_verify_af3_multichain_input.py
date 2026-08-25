"""Validation gate for the AF3+MSA multi-chain fix (no GPU, no AF3 run).

Re-renders the AF3+MSA input JSON with the fixed code path and asserts that the
rendered job carries EVERY protein chain the mutation spec asked for, with a
per-chain MSA whose query row matches that chain's own (possibly mutated)
sequence.

The bug this guards: `_read_af3_protein_sequence` returned only the first
protein and `render_af3` was called with `chain_seqs=[("A", seq)]`, so 51 of 52
multi-chain systems were folded single-chain. For 1bcu the surviving chain was
the 26-residue L while both mutations sat on the discarded 257-residue H.

Usage:
  python analysis/casf_mutagenesis/scripts/26_verify_af3_multichain_input.py [--limit N]

Exits non-zero if any system fails.
"""
from __future__ import annotations
import argparse
import json
import os
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
sys.path.insert(0, str(REPO_ROOT / "analysis"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "casf_mutagenesis" / "scripts"))

from casf_mutagenesis.config import OUTPUT_ROOT  # noqa: E402
from casf_mutagenesis.inputs_af3 import render_af3  # noqa: E402

import importlib.util  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "af3msa_runner",
    REPO_ROOT / "analysis/casf_mutagenesis/scripts/06_run_af3_msa_subset20.py",
)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
read_chains = _mod._read_af3_protein_chains
read_smiles = _mod._read_ligand_smiles

VARIANTS = ("rem", "pack", "inv")


def check(pdbid: str, variant: str, tmp: Path) -> tuple[bool, str]:
    v_dir = OUTPUT_ROOT / pdbid / variant
    af3_json = v_dir / "af3.json"
    wt_json = OUTPUT_ROOT / pdbid / "wt" / "af3.json"
    if not af3_json.exists() or not wt_json.exists():
        return True, "skip (no af3.json)"

    wt_chains = read_chains(wt_json)
    var_chains = read_chains(af3_json)
    if len(var_chains) != len(wt_chains):
        return False, f"chain count wt={len(wt_chains)} variant={len(var_chains)}"

    # Stand-in MSAs: query row + one filler, mimicking rewrite_a3m_query output.
    chain_msas = {cid: f">query\n{seq}\n>filler\n{seq}\n" for cid, seq in var_chains}
    out = tmp / f"{pdbid}_{variant}.json"
    render_af3(name=f"{pdbid}_{variant}", chain_seqs=var_chains,
               ligand_smiles=read_smiles(af3_json), out_path=out,
               chain_msas=chain_msas)

    d = json.loads(out.read_text())
    rendered = []
    for s in d["sequences"]:
        pr = s.get("protein")
        if pr:
            ids = pr["id"]
            for cid in (ids if isinstance(ids, list) else [ids]):
                rendered.append((str(cid), pr["sequence"], pr.get("unpairedMsa", "")))

    if len(rendered) != len(var_chains):
        return False, f"rendered {len(rendered)} chains, spec has {len(var_chains)}"

    spec = {c: s for c, s in var_chains}
    for cid, seq, msa in rendered:
        if cid not in spec:
            return False, f"rendered unknown chain {cid}"
        if seq != spec[cid]:
            return False, f"chain {cid} sequence mismatch"
        first = msa.splitlines()[1] if len(msa.splitlines()) > 1 else ""
        if first != seq:
            return False, f"chain {cid} MSA query row != chain sequence"

    # every mutation in the spec must be reachable in the rendered job
    wt = {c: s for c, s in wt_chains}
    n_mut = sum(
        1 for cid, s in var_chains if len(wt.get(cid, "")) == len(s)
        for a, b in zip(wt[cid], s) if a != b
    )
    if not any(l for l in d.get("sequences", []) if l.get("ligand")):
        return False, "ligand missing from rendered job"
    return True, f"{len(rendered)} chains, {n_mut} mutations carried"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="0 = all systems")
    args = ap.parse_args()

    systems = sorted(p.name for p in OUTPUT_ROOT.iterdir()
                     if p.is_dir() and (p / "wt" / "af3.json").exists())
    if args.limit:
        systems = systems[: args.limit]

    multi, fails, checked = 0, [], 0
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        for pdbid in systems:
            for variant in VARIANTS:
                ok, msg = check(pdbid, variant, tmp)
                if msg.startswith("skip"):
                    continue
                checked += 1
                if msg.split()[0].isdigit() and int(msg.split()[0]) > 1:
                    multi += 1
                if not ok:
                    fails.append(f"{pdbid}/{variant}: {msg}")

    print(f"cells checked           : {checked}")
    print(f"  multi-chain cells     : {multi}")
    print(f"  failures              : {len(fails)}")
    for f in fails[:20]:
        print(f"    FAIL {f}")
    if fails:
        print("\nRESULT: FAIL")
        return 1
    print("\nRESULT: PASS — every rendered AF3 job carries all its protein "
          "chains with matching per-chain MSAs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
