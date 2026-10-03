"""Five-method retention table for the ligand-mutagenesis arm.

Writes, in analysis/ligand_mutagenesis/outputs/:
  memorization_ligand.csv   one row per (method, variant)
  retention_by_dheavy.csv   one row per (method, Δheavy stratum)

READ THIS BEFORE QUOTING A NUMBER. On the ligand arm the metric is RETENTION:
P(top-1 ligand RMSD < 2 Å | the method solved this system's WT < 2 Å). Unlike
the protein arm, high retention is not evidence of memorization — a single
halogen swap plausibly SHOULD leave the pose where it was. The file keeps the
`memorization_` name only to sit beside `memorization_full.csv`.

Always quote Δheavy with a rate: perturbation size is exact per variant here
(halo_* = 1, meth_k = k, chrg_* = 2..10) and retention is strongly
dose-dependent. Never pool across families.

Sources (top-1 everywhere):
  co-folding : results_ligand.csv, pose_idx == 0           (Boltz2, AF3+MSA)
  docking    : casf_mutagenesis/outputs/docking_results.csv, module == ligand
               (GNINA, UniDock2, SurfDock; rank-1 pose)
Failure cells (missing_cif, docking error) are counted in n_fail, never dropped.
"""
from __future__ import annotations
import json
import math
import os
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")
REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
LIG = REPO_ROOT / "analysis" / "ligand_mutagenesis" / "outputs"
CASF = REPO_ROOT / "analysis" / "casf_mutagenesis" / "outputs"
THRESH = 2.0
ENGINE_NAME = {"gnina": "GNINA", "unidock2": "UniDock2", "surfdock": "SurfDock"}


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (math.nan, math.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (round(c - h, 3), round(c + h, 3))


def delta_heavy() -> pd.DataFrame:
    man = json.loads((LIG / "manifest_full.json").read_text())
    rows = []
    for s in man["systems"]:
        if s.get("status") != "ok":
            continue
        v = s["variants"]
        wt = Chem.MolFromSmiles(v["wt"]["smiles"])
        if wt is None:
            continue
        for name, info in v.items():
            m = Chem.MolFromSmiles(info.get("smiles") or "")
            if m is None:
                continue
            rows.append({"pdbid": s["pdbid"], "variant": name,
                         "delta_heavy": m.GetNumHeavyAtoms() - wt.GetNumHeavyAtoms()})
    return pd.DataFrame(rows)


def load_cells() -> pd.DataFrame:
    """One row per (method, pdbid, variant): top-1 rmsd or NaN + ok flag."""
    co = pd.read_csv(LIG / "results_ligand.csv")
    co = co[co.pose_idx == 0]
    co = pd.DataFrame({"method": co.model, "pdbid": co.pdbid, "variant": co.variant,
                       "ok": co.status == "ok", "rmsd": co.ligand_rmsd_a})
    dk = pd.read_csv(CASF / "docking_results.csv")
    dk = dk[dk.module == "ligand"]
    dk = pd.DataFrame({"method": dk.engine.map(ENGINE_NAME), "pdbid": dk.system,
                       "variant": dk.variant, "ok": dk.status == "ok", "rmsd": dk.rmsd_a})
    cells = pd.concat([co, dk], ignore_index=True)
    cells.loc[~cells.ok, "rmsd"] = math.nan
    return cells


def main() -> int:
    dh = delta_heavy()
    cells = load_cells().merge(dh, on=["pdbid", "variant"], how="left")
    cells["family"] = cells.variant.str.replace(r"_\d+$", "", regex=True)

    wt = cells[cells.variant == "wt"]
    solved = {(m, p) for m, p, r in zip(wt.method, wt.pdbid, wt.rmsd) if r < THRESH}
    cells["wtok"] = [(m, p) in solved for m, p in zip(cells.method, cells.pdbid)]

    out = []
    for (method, variant), g in cells.groupby(["method", "variant"]):
        ok = g[g.ok]
        c = g[g.wtok & g.ok]
        k = int((c.rmsd < THRESH).sum())
        lo, hi = wilson(k, len(c))
        wt_m = wt[wt.method == method]
        out.append({
            "method": method, "family": g.family.iloc[0], "variant": variant,
            "delta_heavy_median": g.delta_heavy.median(),
            "n_total": len(g), "n_fail": int((~g.ok).sum()),
            "n_wtok": len(c), "n_fail_wtok": int((g.wtok & ~g.ok).sum()),
            "retention_2A_wtok": round(k / len(c), 3) if len(c) else math.nan,
            "ci_lo_2A_wtok": lo, "ci_hi_2A_wtok": hi,
            "retention_4A_wtok": round(float((c.rmsd < 4).mean()), 3) if len(c) else math.nan,
            "median_rmsd_A_wtok": round(float(c.rmsd.median()), 2) if len(c) else math.nan,
            "retention_2A": round(float((ok.rmsd < THRESH).mean()), 3) if len(ok) else math.nan,
            "wt_correct_systems": int((wt_m.rmsd < THRESH).sum()),
            "wt_scored_systems": int(wt_m.ok.sum()),
        })
    tab = pd.DataFrame(out).sort_values(["method", "delta_heavy_median", "variant"])
    tab.to_csv(LIG / "memorization_ligand.csv", index=False)

    adv = cells[(cells.variant != "wt") & cells.wtok & cells.ok].copy()
    adv["stratum"] = pd.cut(adv.delta_heavy, [-100, 0, 1, 2, 100],
                            labels=["<=0", "1", "2", ">=3"])
    strat = []
    for (method, st), g in adv.groupby(["method", "stratum"], observed=True):
        k = int((g.rmsd < THRESH).sum())
        lo, hi = wilson(k, len(g))
        strat.append({"method": method, "delta_heavy_stratum": st, "n_wtok": len(g),
                      "retention_2A_wtok": round(k / len(g), 3), "ci_lo": lo, "ci_hi": hi})
    pd.DataFrame(strat).to_csv(LIG / "retention_by_dheavy.csv", index=False)

    print(f"wrote {LIG/'memorization_ligand.csv'} ({len(tab)} rows)")
    print(f"wrote {LIG/'retention_by_dheavy.csv'}")
    print(pd.DataFrame(strat).pivot(index="method", columns="delta_heavy_stratum",
                                    values="retention_2A_wtok").to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
