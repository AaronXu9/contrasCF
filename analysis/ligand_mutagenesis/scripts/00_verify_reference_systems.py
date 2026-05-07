"""Gating test: run the methylation and charge-swap generators against the
two ligand systems Masters et al. 2025 hand-built — α-D-glucose (paper Fig.
4) and ATP (paper Fig. 5) — and assert that the auto-generated SMILES match
the paper's literal SMILES already pinned in `analysis/src/config.py`.

Failure here means the rule logic is wrong — do not scale to CASF until both
ladders pass. The comparison is done on RDKit-canonicalized SMILES so that
permutation-equivalent SMILES writings of the same molecule compare equal.

Run:
    LD_LIBRARY_PATH=/home/aoxu/miniconda3/envs/rdkit_env/lib:$LD_LIBRARY_PATH \\
        /home/aoxu/miniconda3/envs/rdkit_env/bin/python \\
        analysis/ligand_mutagenesis/scripts/00_verify_reference_systems.py
"""
from __future__ import annotations
import os
import sys
from pathlib import Path

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
sys.path.insert(0, str(REPO_ROOT / "analysis"))

from rdkit import Chem

# Paper's literal SMILES (already pinned in the analysis pipeline).
from src.config import ATP_SMILES, ATP_CHARGE_SMILES, GLUCOSE_SMILES  # noqa: E402

from ligand_mutagenesis.rules.methylation import MethylationRule  # noqa: E402
from ligand_mutagenesis.rules.charge_swap import ChargeSwapRule  # noqa: E402


def _canon(smi: str) -> str:
    """Canonical isomeric SMILES — same molecule → same string."""
    m = Chem.MolFromSmiles(smi)
    if m is None:
        raise ValueError(f"unparseable SMILES: {smi}")
    try:
        Chem.SanitizeMol(m)
    except Exception:
        pass
    return Chem.MolToSmiles(m, canonical=True, isomericSmiles=True)


def _diff_smiles(label: str, got: str, expected: str) -> bool:
    g = _canon(got)
    e = _canon(expected)
    if g == e:
        print(f"  ✓ {label}: SMILES match ({g})")
        return True
    print(f"  ✗ {label}: SMILES MISMATCH")
    print(f"    got:      {g}")
    print(f"    expected: {e}")
    return False


def verify_glucose_methylation() -> bool:
    """Methylation generator on glucose_0 must reproduce glucose_1..5."""
    print("\n=== glucose hydroxyl methylation (paper Fig. 4) ===")
    parent_smi = GLUCOSE_SMILES["glucose_0"]
    parent = Chem.MolFromSmiles(parent_smi)
    if parent is None:
        print(f"  FAIL: could not parse glucose_0 SMILES: {parent_smi}")
        return False

    rule = MethylationRule()
    eligible, reason, n_matches = rule.is_eligible(parent)
    print(f"  eligibility: {eligible} (n_matches={n_matches}) {reason}")
    if not eligible:
        return False

    variants = rule.generate(parent)
    print(f"  generated {len(variants)} variants: "
          f"{[v.name for v in variants]}")

    by_name = {v.name: v for v in variants}
    all_ok = True
    for k in range(1, 6):
        name = f"meth_{k}"
        if name not in by_name:
            print(f"  ✗ meth_{k}: variant not generated")
            all_ok = False
            continue
        ok = _diff_smiles(
            f"meth_{k} vs glucose_{k}",
            by_name[name].smiles,
            GLUCOSE_SMILES[f"glucose_{k}"],
        )
        all_ok = all_ok and ok
    return all_ok


def verify_atp_charge_swap() -> bool:
    """Charge-swap generator on ATP must reproduce atp_charge_*."""
    print("\n=== ATP triphosphate charge swap (paper Fig. 5) ===")
    parent = Chem.MolFromSmiles(ATP_SMILES)
    if parent is None:
        print(f"  FAIL: could not parse ATP_SMILES")
        return False

    rule = ChargeSwapRule()
    eligible, reason, n_matches = rule.is_eligible(parent)
    print(f"  eligibility: {eligible} (n_matches={n_matches}) {reason}")
    if not eligible:
        return False

    variants = rule.generate(parent)
    print(f"  generated {len(variants)} variants: "
          f"{[v.name for v in variants]}")
    by_name = {v.name: v for v in variants}

    expected_pairs = [
        ("chrg_neu_methyl", "atp_charge_methyl"),
        ("chrg_neu_ethyl", "atp_charge_ethyl"),
        ("chrg_neu_propyl", "atp_charge_propyl"),
        ("chrg_pos_1", "atp_charge_1"),
        ("chrg_pos_2", "atp_charge_2"),
        ("chrg_pos_3", "atp_charge_3"),
    ]
    all_ok = True
    for our_name, paper_name in expected_pairs:
        if our_name not in by_name:
            print(f"  ✗ {our_name}: variant not generated")
            all_ok = False
            continue
        ok = _diff_smiles(
            f"{our_name} vs {paper_name}",
            by_name[our_name].smiles,
            ATP_CHARGE_SMILES[paper_name],
        )
        all_ok = all_ok and ok
    return all_ok


def verify_negative_eligibility() -> bool:
    """Sanity: glucose has no phosphate; charge_swap should refuse it."""
    print("\n=== negative eligibility checks ===")
    glc = Chem.MolFromSmiles(GLUCOSE_SMILES["glucose_0"])
    cs = ChargeSwapRule()
    eligible, reason, _ = cs.is_eligible(glc)
    if eligible:
        print(f"  ✗ charge_swap accepted glucose (should refuse): {reason}")
        return False
    print(f"  ✓ charge_swap refuses glucose: {reason}")
    return True


def main() -> int:
    ok_glucose = verify_glucose_methylation()
    ok_atp = verify_atp_charge_swap()
    ok_neg = verify_negative_eligibility()

    print("\n=== summary ===")
    print(f"  glucose methylation : {'PASS' if ok_glucose else 'FAIL'}")
    print(f"  ATP charge swap     : {'PASS' if ok_atp else 'FAIL'}")
    print(f"  negative eligibility: {'PASS' if ok_neg else 'FAIL'}")
    return 0 if (ok_glucose and ok_atp and ok_neg) else 1


if __name__ == "__main__":
    sys.exit(main())
