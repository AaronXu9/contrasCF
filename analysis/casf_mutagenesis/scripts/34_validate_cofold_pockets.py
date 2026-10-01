#!/usr/bin/env python3
"""Is the co-folded mutant complex PHYSICALLY VALID, or did the model strain the
structure to preserve a memorized ligand pose?

Two independent tests on the Boltz-2 rank-0 predictions, per variant:

A. PoseBusters (`dock` config) on (predicted ligand, predicted protein).
   Covers sanitization, connectivity, bond lengths/angles, internal clash, ring
   flatness, internal energy, minimum protein-ligand distance and volume
   overlap. If mutant complexes pass at the WT rate, the model is producing
   genuinely legal structures around the remembered pose rather than clashing
   ones -- which makes the memorization claim stronger, not weaker.

B. Side-chain rotamer plausibility for the INTRODUCED residues (the ones the
   mutation actually created), plus their orientation relative to the ligand.
   `pack` is the sharp case: filling a pocket with Phe should occlude it, so if
   there are no clashes, the Phe must have gone somewhere. Two possibilities:
     - they adopt normal rotamers and the site was never really occluded, or
     - they adopt strained/outlier rotamers pointing away from the ligand,
       i.e. the model deformed the protein to protect the memorized pose.
   chi1 outlier = >40 deg from the nearest canonical value (-60 / 60 / 180).
   Baselines: the same residue type in the crystal structure, and in the
   model's own WT prediction.

Run (PoseBench env has posebusters + gemmi + rdkit):
  /home/aoxu/miniconda3/envs/PoseBench/bin/python 34_validate_cofold_pockets.py --limit 60
  ... --out outputs/cofold_validity.csv
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import tempfile
import warnings
from collections import defaultdict
from difflib import SequenceMatcher
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
sys.path.insert(0, str(REPO_ROOT / "analysis"))

import gemmi  # noqa: E402
from rdkit import Chem, RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")

from casf_mutagenesis.config import CASF_LIGANDS, CASF_RAW, OUTPUT_ROOT, VARIANTS  # noqa: E402

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M",
}

# chi1 = N-CA-CB-XG ; chi2 = CA-CB-XG-XD  (first branch atom)
CHI1 = {"N", "CA", "CB"}
CHI_G = {"PHE": "CG", "TRP": "CG", "TYR": "CG", "LEU": "CG", "ASP": "CG",
         "GLU": "CG", "ASN": "CG", "GLN": "CG", "HIS": "CG", "ARG": "CG",
         "LYS": "CG", "MET": "CG", "ILE": "CG1", "VAL": "CG1", "THR": "OG1",
         "SER": "OG", "CYS": "SG", "PRO": "CG"}
CHI_D = {"PHE": "CD1", "TRP": "CD1", "TYR": "CD1", "LEU": "CD1", "ASP": "OD1",
         "ASN": "OD1", "GLU": "CD", "GLN": "CD", "HIS": "ND1", "ARG": "CD",
         "LYS": "CD", "MET": "SD", "ILE": "CD1"}
CANONICAL_CHI1 = (-60.0, 60.0, 180.0)

# PoseBusters "dock" checks that are MEANINGFUL here. The cofactor/water checks
# are vacuous: co-folding predictions contain no waters or cofactors, so the
# conditioning PDB has none and those checks fail 100% of the time for every
# variant including wt. Reporting them would make the comparison meaningless.
CORE_CHECKS = (
    "sanitization", "all_atoms_connected", "bond_lengths", "bond_angles",
    "internal_steric_clash", "aromatic_ring_flatness", "double_bond_flatness",
    "internal_energy", "protein-ligand_maximum_distance",
    "minimum_distance_to_protein", "volume_overlap_with_protein",
)
CHI1_TOL = 40.0


# --------------------------------------------------------------------------

def dihedral(p0, p1, p2, p3) -> float:
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return float(np.degrees(np.arctan2(
        np.dot(np.cross(b1, v), w), np.dot(v, w))))


def chi1_offset(chi1: float) -> float:
    """Angular distance to the nearest canonical chi1 rotamer well."""
    return min(abs((chi1 - c + 180) % 360 - 180) for c in CANONICAL_CHI1)


def af3_chains(path: Path) -> dict[str, str]:
    d = json.loads(path.read_text())
    out = {}
    for s in d.get("sequences", []):
        pr = s.get("protein")
        if not pr:
            continue
        ids = pr["id"]
        for cid in (ids if isinstance(ids, list) else [ids]):
            out[str(cid)] = pr["sequence"]
    return out


def mutated_positions(system: str, variant: str) -> dict[str, list[tuple[int, str, str]]]:
    """{chain_id: [(seq_index, wt_aa, mut_aa), ...]} from the af3.json specs."""
    wt_j = OUTPUT_ROOT / system / "wt" / "af3.json"
    var_j = OUTPUT_ROOT / system / variant / "af3.json"
    if not wt_j.exists() or not var_j.exists():
        return {}
    wt, mut = af3_chains(wt_j), af3_chains(var_j)
    return {cid: [(i, a, b) for i, (a, b) in enumerate(zip(wt[cid], seq)) if a != b]
            for cid, seq in mut.items()
            if cid in wt and len(wt[cid]) == len(seq)}


def split_complex(cif: Path):
    """(gemmi Structure of protein only, ligand residue) from a prediction CIF."""
    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    st.remove_hydrogens()
    lig, lig_n = None, 0
    for model in st:
        for ch in model:
            for res in ch:
                if res.name in AA3 or res.name in ("HOH", "WAT"):
                    continue
                n = sum(1 for a in res if a.element.name != "H")
                if n > lig_n:
                    lig, lig_n = res, n
    return st, lig


def protein_pdb(st: gemmi.Structure, path: Path) -> None:
    s2 = st.clone()
    s2.setup_entities()
    s2.remove_ligands_and_waters()   # gemmi 0.7: Chain has no remove_residue()
    s2.remove_empty_chains()
    s2.write_pdb(str(path))


def ligand_mol(lig, smiles: str) -> Chem.Mol | None:
    """RDKit mol for the predicted ligand with bond orders from the crystal SMILES."""
    blk = ["HETATM%5d %-4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s" %
           (i + 1, a.name[:4], "LIG", 1, a.pos.x, a.pos.y, a.pos.z,
            a.element.name.upper())
           for i, a in enumerate(lig) if a.element.name != "H"]
    m = Chem.MolFromPDBBlock("\n".join(blk) + "\nEND\n", removeHs=True, sanitize=False)
    if m is None:
        return None
    tmpl = Chem.MolFromSmiles(smiles)
    if tmpl is not None:
        try:
            from rdkit.Chem import AllChem
            return AllChem.AssignBondOrdersFromTemplate(tmpl, m)
        except Exception:
            pass
    try:
        Chem.SanitizeMol(m)
    except Exception:
        pass
    return m


# --------------------------------------------------------------------------
# B. rotamer + orientation of the introduced side chains
# --------------------------------------------------------------------------

def residue_chis(res) -> tuple[float | None, float | None]:
    pos = {a.name: np.array([a.pos.x, a.pos.y, a.pos.z]) for a in res}
    g, d = CHI_G.get(res.name), CHI_D.get(res.name)
    chi1 = chi2 = None
    if g and CHI1 <= set(pos) and g in pos:
        chi1 = dihedral(pos["N"], pos["CA"], pos["CB"], pos[g])
        if d and d in pos:
            chi2 = dihedral(pos["CA"], pos["CB"], pos[g], pos[d])
    return chi1, chi2


def sidechain_orientation(res, lig_centroid) -> float | None:
    """Angle (deg) between CB->sidechain-centroid and CB->ligand-centroid.

    ~0 deg = side chain points AT the ligand; ~180 deg = points away.
    """
    pos = {a.name: np.array([a.pos.x, a.pos.y, a.pos.z]) for a in res}
    if "CB" not in pos:
        return None
    sc = [v for k, v in pos.items() if k not in ("N", "CA", "C", "O", "OXT", "CB")]
    if not sc:
        return None
    v1 = np.mean(sc, axis=0) - pos["CB"]
    v2 = lig_centroid - pos["CB"]
    n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
    if n1 < 1e-6 or n2 < 1e-6:
        return None
    return float(np.degrees(math.acos(np.clip(np.dot(v1, v2) / (n1 * n2), -1, 1))))


def analyze_introduced(st, lig, system, variant) -> list[dict]:
    """chi1/chi2 + orientation for every residue the mutation actually created."""
    muts = mutated_positions(system, variant)
    if not muts or lig is None:
        return []
    lig_xyz = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in lig
                        if a.element.name != "H"])
    lig_c = lig_xyz.mean(axis=0)

    chains = {}
    for ch in st[0]:
        rs = [r for r in ch if r.name in AA3]
        if rs:
            chains[ch.name] = rs

    out = []
    for cid, mlist in muts.items():
        spec = "".join(a if False else b for (_, a, b) in [])  # placeholder
        # Build the variant sequence string for this chain from the spec file.
        var_seq = af3_chains(OUTPUT_ROOT / system / variant / "af3.json").get(cid)
        if var_seq is None:
            continue
        # Match the spec sequence to whichever predicted chain it corresponds to.
        best, best_r = None, 0.0
        for name, rs in chains.items():
            seq = "".join(AA3[r.name] for r in rs)
            r = SequenceMatcher(None, var_seq, seq, autojunk=False).ratio()
            if r > best_r:
                best, best_r = name, r
        if best is None or best_r < 0.50:
            continue
        rs = chains[best]
        seq = "".join(AA3[r.name] for r in rs)
        idx = {}
        for a, b, size in SequenceMatcher(
                None, var_seq, seq, autojunk=False).get_matching_blocks():
            for k in range(size):
                idx[a + k] = b + k
        for (i, wt_aa, mut_aa) in mlist:
            j = idx.get(i)
            if j is None:
                continue
            res = rs[j]
            if AA3.get(res.name) != mut_aa:
                continue          # mutation did not survive into this structure
            chi1, chi2 = residue_chis(res)
            ca = next((a for a in res if a.name == "CA"), None)
            if ca is None:
                continue
            d_lig = float(np.linalg.norm(
                lig_xyz - np.array([ca.pos.x, ca.pos.y, ca.pos.z]), axis=1).min())
            out.append({
                "system": system, "variant": variant, "chain": best,
                "wt_aa": wt_aa, "mut_aa": mut_aa, "resname": res.name,
                "chi1": chi1, "chi2": chi2,
                "chi1_offset": chi1_offset(chi1) if chi1 is not None else None,
                "orientation_deg": sidechain_orientation(res, lig_c),
                "ca_to_lig_a": round(d_lig, 2),
            })
    return out


def crystal_baseline(system: str, resnames: set[str]) -> list[dict]:
    """chi1 offsets for the same residue types in the CRYSTAL structure."""
    pdb = CASF_RAW / system / f"{system}_protein.pdb"
    if not pdb.exists():
        return []
    st = gemmi.read_structure(str(pdb))
    st.remove_hydrogens()
    out = []
    for ch in st[0]:
        for res in ch:
            if res.name not in resnames:
                continue
            chi1, chi2 = residue_chis(res)
            if chi1 is None:
                continue
            out.append({"system": system, "resname": res.name, "chi1": chi1,
                        "chi1_offset": chi1_offset(chi1)})
    return out


# --------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--systems", default=None)
    ap.add_argument("--out", default=str(OUTPUT_ROOT / "cofold_validity.csv"))
    ap.add_argument("--rot-out", default=str(OUTPUT_ROOT / "cofold_rotamers.csv"))
    ap.add_argument("--skip-pb", action="store_true")
    args = ap.parse_args()

    from posebusters import PoseBusters

    if args.systems:
        systems = [s.strip().lower() for s in args.systems.split(",") if s.strip()]
    else:
        systems = sorted(d.name for d in OUTPUT_ROOT.iterdir()
                         if d.is_dir() and len(d.name) == 4)
    if args.limit:
        systems = systems[:args.limit]
    print(f"{len(systems)} system(s)")

    buster = PoseBusters(config="dock")
    pb_rows: list[dict] = []
    rot_rows: list[dict] = []
    base_rows: list[dict] = []
    tmp = Path(tempfile.mkdtemp(prefix="cofold_val_"))

    for n, system in enumerate(systems, 1):
        crystal_sdf = CASF_LIGANDS / f"{system}_ligand.sdf"
        if not crystal_sdf.exists():
            continue
        cm = next((m for m in Chem.SDMolSupplier(str(crystal_sdf), removeHs=True,
                                                 sanitize=False) if m is not None), None)
        if cm is None:
            continue
        try:
            Chem.SanitizeMol(cm)
            smiles = Chem.MolToSmiles(cm)
        except Exception:
            smiles = ""

        seen_resnames: set[str] = set()
        for variant in VARIANTS:
            cifs = sorted((OUTPUT_ROOT / system / variant).glob(
                f"{system}_{variant}_model_*.cif"))
            if not cifs:
                continue
            st, lig = split_complex(cifs[0])
            if lig is None:
                continue

            if variant != "wt":
                r = analyze_introduced(st, lig, system, variant)
                rot_rows.extend(r)
                seen_resnames.update(x["resname"] for x in r)

            if args.skip_pb:
                continue
            mol = ligand_mol(lig, smiles)
            if mol is None:
                continue
            pp = tmp / f"{system}_{variant}_prot.pdb"
            ls = tmp / f"{system}_{variant}_lig.sdf"
            protein_pdb(st, pp)
            w = Chem.SDWriter(str(ls)); w.write(mol); w.close()
            try:
                df = buster.bust([ls], None, pp)
            except Exception as exc:
                print(f"  {system}/{variant}: PB failed ({type(exc).__name__})")
                continue
            rec = {"system": system, "variant": variant}
            for col in df.columns:
                v = df.iloc[0][col]
                if isinstance(v, (bool, np.bool_)):
                    rec[col] = bool(v)
            rec["pb_valid"] = all(rec.get(c, True) for c in CORE_CHECKS)
            pb_rows.append(rec)

        if seen_resnames:
            base_rows.extend(crystal_baseline(system, seen_resnames))
        if n % 10 == 0:
            print(f"  [{n}/{len(systems)}] pb={len(pb_rows)} rot={len(rot_rows)}")

    if pb_rows:
        keys = list(dict.fromkeys(k for r in pb_rows for k in r))
        with open(args.out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
            for r in pb_rows: w.writerow(r)
        print(f"\nPoseBusters per-cell: {args.out}  ({len(pb_rows)} cells)")

    if rot_rows:
        keys = list(rot_rows[0].keys())
        with open(args.rot_out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
            for r in rot_rows: w.writerow(r)
        print(f"Rotamers per-residue: {args.rot_out}  ({len(rot_rows)} residues)")
    if base_rows:
        p = Path(args.rot_out).with_name("crystal_rotamer_baseline.csv")
        with open(p, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(base_rows[0].keys())); w.writeheader()
            for r in base_rows: w.writerow(r)
        print(f"Crystal baseline:     {p}  ({len(base_rows)} residues)")

    # ---- summaries ----
    if pb_rows:
        print("\nA. PoseBusters validity of the co-folded complex")
        by = defaultdict(list)
        for r in pb_rows: by[r["variant"]].append(r)
        checks = [c for c in CORE_CHECKS if c in pb_rows[0]]
        print(f"   {'variant':<6s} {'n':>4s} {'all-pass':>9s}   worst failing checks")
        for v in VARIANTS:
            rs = by.get(v)
            if not rs: continue
            fails = sorted(((sum(1 for r in rs if not r.get(c, True)) / len(rs), c)
                            for c in checks), reverse=True)[:3]
            ftxt = ", ".join(f"{c} {p:.0%}" for p, c in fails if p > 0) or "none"
            print(f"   {v:<6s} {len(rs):>4d} "
                  f"{sum(1 for r in rs if r['pb_valid']) / len(rs):>8.1%}   {ftxt}")

    if rot_rows:
        print("\nB. Introduced side chains: rotamer strain + orientation")
        by = defaultdict(list)
        for r in rot_rows:
            if r["chi1_offset"] is not None: by[r["variant"]].append(r)
        print(f"   {'variant':<6s} {'n_res':>6s} {'chi1 outlier':>13s} "
              f"{'med offset':>11s} {'pointing away':>14s}")
        for v in ("rem", "pack", "inv"):
            rs = by.get(v)
            if not rs: continue
            off = [r["chi1_offset"] for r in rs]
            ori = [r["orientation_deg"] for r in rs if r["orientation_deg"] is not None]
            print(f"   {v:<6s} {len(rs):>6d} "
                  f"{sum(1 for o in off if o > CHI1_TOL) / len(off):>12.1%} "
                  f"{np.median(off):>10.1f}° "
                  f"{(sum(1 for o in ori if o > 90) / len(ori) if ori else float('nan')):>13.1%}")
        if base_rows:
            off = [r["chi1_offset"] for r in base_rows]
            print(f"   {'CRYSTAL':<6s} {len(off):>6d} "
                  f"{sum(1 for o in off if o > CHI1_TOL) / len(off):>12.1%} "
                  f"{np.median(off):>10.1f}°   (experimental baseline)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
