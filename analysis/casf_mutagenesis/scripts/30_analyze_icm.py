"""Score the ICM docking poses with the SAME metric as the other engines.

ICM was run on CARC by a collaborator (2026-09-01) against
`docking/receptor_aligned.pdb` — the corrected multi-chain receptors
(2026-08-25) pre-transformed into the CRYSTAL frame. Verified: mean Cα
distance to the crystal receptor is 0.80 Å (1bcu) / 1.61 Å (2qnq) for
`receptor_aligned.pdb` versus 65.95 Å / 29.75 Å for the untransformed
`receptor.pdb`. Poses are therefore already in the crystal frame and need
NO superposition — the same situation as the WT GNINA cells.

Deliberately reuses `gnina_analysis._mcs_match_indices` + `_aligned_rmsd` so
ICM is directly comparable to GNINA / UniDock2 / SurfDock. That matcher is
known to take a single substructure match rather than minimising over
symmetry-equivalent ones (docs/rmsd_and_failure_handling.md §4), which
understates rates by ~2.4 points; the bias applies equally to all four
engines, so the comparison stays internally consistent. A symmetry-minimised
value is emitted alongside as `rmsd_symmin_a` so the effect is visible.

Input:  outputs/_icm_poses/<sys>/<variant>/ICM/D_<sys>_<variant>_pose1.sdf
Output: outputs/icm_results.csv  (schema-compatible with docking_results.csv)

Usage:
  python analysis/casf_mutagenesis/scripts/30_analyze_icm.py
"""
from __future__ import annotations
import csv
import os
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
sys.path.insert(0, str(REPO_ROOT / "analysis"))

from casf_mutagenesis.config import CASF_LIGANDS, OUTPUT_ROOT  # noqa: E402
from casf_mutagenesis.gnina_analysis import (  # noqa: E402
    _aligned_rmsd, _mcs_match_indices,
)

ICM_ROOT = OUTPUT_ROOT / "_icm_poses"
OUT_CSV = OUTPUT_ROOT / "icm_results.csv"
VARIANTS = ("wt", "rem", "pack", "inv")


def _heavy(mol: Chem.Mol, sanitize: bool = True) -> tuple[Chem.Mol, np.ndarray]:
    m = Chem.RemoveHs(mol, sanitize=sanitize)
    c = m.GetConformer()
    return m, np.array([list(c.GetAtomPosition(i)) for i in range(m.GetNumAtoms())])


def _read_icm_pose(path: Path) -> tuple[Chem.Mol | None, str]:
    """Read an ICM pose SDF, falling back to a lenient parse.

    ICM writes valences RDKit rejects — e.g. `AtomValenceException: Explicit
    valence for atom # 12 P, 7, is greater than permitted` on phosphates. A
    strict read drops 228 of 948 cells (24 %), and the failure is entirely on
    the pose side (the crystal SDFs parse fine).

    Fallback re-reads unsanitized and then runs every sanitization step EXCEPT
    `SANITIZE_PROPERTIES` (the valence check). Ring perception and aromaticity
    still run, so MCS matching behaves normally, and only heavy-atom
    COORDINATES are used downstream. Verified: 60/60 sampled cells sanitize
    this way and 25/25 keep a heavy-atom count equal to the crystal ligand.
    """
    m = next(iter(Chem.SDMolSupplier(str(path), sanitize=True)), None)
    if m is not None:
        return m, "strict"
    m = next(iter(Chem.SDMolSupplier(str(path), sanitize=False)), None)
    if m is None:
        return None, "failed"
    try:
        m = Chem.Mol(m)
        Chem.SanitizeMol(
            m, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        return m, "lenient"
    except Exception:
        return None, "failed"


def score_cell(sysid: str, variant: str) -> dict:
    rec = {"system": sysid, "variant": variant, "module": "casf", "engine": "icm",
           "status": "ok", "error": "", "rmsd_a": "", "rmsd_symmin_a": "",
           "n_matched_heavy": "", "icm_self_rmsd": "", "finished": "",
           "parse_mode": ""}
    cell = ICM_ROOT / sysid / variant / "ICM"
    pose = cell / f"D_{sysid}_{variant}_pose1.sdf"
    rec["finished"] = "yes" if (cell / "FINISHED").exists() else "no"

    # ICM's own reported rank-1 RMSD, for cross-checking only
    rr = cell / "redock_rmsd_results.csv"
    if rr.exists():
        try:
            first = rr.read_text().strip().splitlines()[0].split(",")
            rec["icm_self_rmsd"] = round(float(first[2]), 3)
        except Exception:
            pass

    if not pose.exists():
        rec["status"] = "missing_pose"
        return rec
    cl = CASF_LIGANDS / f"{sysid}_ligand.sdf"
    if not cl.exists():
        rec["status"] = "missing_crystal"
        return rec
    try:
        cm = Chem.MolFromMolFile(str(cl), sanitize=True)
        pm, mode = _read_icm_pose(pose)
        rec["parse_mode"] = mode
        if cm is None or pm is None:
            rec["status"] = "parse_error"
            return rec
        cm, cx = _heavy(cm)
        pm, px = _heavy(pm, sanitize=(mode == "strict"))

        c_idx, p_idx, n = _mcs_match_indices(cm, pm)
        if n == 0:
            rec["status"] = "no_match"
            return rec
        rec["n_matched_heavy"] = n
        # No superposition: receptor_aligned.pdb put the poses in the crystal frame.
        rec["rmsd_a"] = round(_aligned_rmsd(cx[c_idx], px[p_idx]), 3)

        # Symmetry-minimised companion (see module docstring).
        try:
            best = None
            for m in pm.GetSubstructMatches(cm, uniquify=False, maxMatches=200):
                r = float(np.sqrt(((px[list(m)] - cx) ** 2).sum(1).mean()))
                best = r if best is None else min(best, r)
            if best is not None:
                rec["rmsd_symmin_a"] = round(best, 3)
        except Exception:
            pass
    except Exception as exc:
        rec["status"] = "error"
        rec["error"] = f"{type(exc).__name__}: {exc}"[:200]
    return rec


def main() -> int:
    if not ICM_ROOT.is_dir():
        print(f"missing {ICM_ROOT}")
        return 1
    systems = sorted(p.name for p in ICM_ROOT.iterdir() if p.is_dir() and len(p.name) == 4)
    rows = []
    for s in systems:
        for v in VARIANTS:
            if not (ICM_ROOT / s / v / "ICM").is_dir():
                continue
            rows.append(score_cell(s, v))

    with OUT_CSV.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    from collections import Counter
    print(f"cells scored: {len(rows)} -> {OUT_CSV}")
    print("  status:", dict(Counter(r["status"] for r in rows)))
    print("  FINISHED marker:", dict(Counter(r["finished"] for r in rows)))
    print(f"\n{'variant':<7}{'n':>5}{'<2A':>8}{'median':>9}{'symmin<2A':>11}{'ICM self':>10}")
    for v in VARIANTS:
        ok = [r for r in rows if r["variant"] == v and r["status"] == "ok"]
        if not ok:
            continue
        a = np.array([r["rmsd_a"] for r in ok])
        sm = np.array([r["rmsd_symmin_a"] for r in ok if r["rmsd_symmin_a"] != ""])
        self_ = np.array([r["icm_self_rmsd"] for r in ok if r["icm_self_rmsd"] != ""])
        print(f"{v:<7}{len(a):>5}{(a<2).mean():>8.3f}{np.median(a):>9.2f}"
              f"{(sm<2).mean() if len(sm) else float('nan'):>11.3f}"
              f"{(self_<2).mean() if len(self_) else float('nan'):>10.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
