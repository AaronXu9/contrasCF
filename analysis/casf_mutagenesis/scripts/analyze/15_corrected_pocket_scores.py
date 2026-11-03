"""Corrected pocket-arm scores for all six methods (audit follow-up 1).

Generalizes the 2026-10-21 audit scorer (experiments/2026-10-21-pipeline-audit/
pocket/audit_pockets.py: Set B, pack/inv, three docking engines) to every system,
variant (wt/rem/pack/inv) and method (GNINA, UniDock2, ICM, SurfDock, Boltz-2,
AF3+MSA). Protocol, identical for every method:

1. Chain correspondence by EXACT sequence: each crystal chain must equal its
   WT input sequence; each target chain (docking receptor, ICM-aligned receptor,
   or co-folded structure) must equal the variant input sequence after removing
   non-standard residues. No nearest-chain fallback — the legacy fallback put
   multimer poses onto a differently mutated site (audit §1).
2. Frame: Kabsch fit of the target's Cα onto the crystal's Cα over the
   crystal-defined NATIVE-chain 5 Å residue shell (ligand never fitted).
3. Ligand RMSD in that frame, symmetry-aware (core.ligand_rmsd), top-1 pose.
4. Flags: no-op (zero actual substitutions), number of substitutions, pose
   contacts with changed side chains, native ligand inside the docking box.

GNINA cells with the uppercase CL/BR input defect use the audit's corrected
replays (gnina_replay/replacement_manifest.csv). Output:
outputs/corrected_pocket_scores.csv (one row per system x variant x method) and
outputs/corrected_pocket_table.csv (summary). Read-only on predictions.
"""
from __future__ import annotations
import csv, json, os, sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem, RDLogger

REPO = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
sys.path.insert(0, str(REPO / "analysis"))
from casf_mutagenesis.config import AA3_TO_AA1, AA3_STANDARD, CASF_LIGANDS, CASF_RAW  # noqa: E402
from casf_mutagenesis.gnina_analysis import _read_top_pose  # noqa: E402
from core.ligand_rmsd import _matched_rmsd  # noqa: E402

RDLogger.DisableLog("rdApp.*")
OUT = REPO / "analysis/casf_mutagenesis/outputs"
ICM_ALIGNED = REPO / "experiments/2026-11-03-icm-aligned-receptors"
GNINA_FIX = REPO / "experiments/2026-10-21-pipeline-audit/gnina_replay/replacement_manifest.csv"
VARIANTS = ("wt", "rem", "pack", "inv")
METHODS = ("GNINA", "UniDock2", "ICM", "SurfDock", "Boltz-2", "AF3+MSA")
BACKBONE = {"N", "CA", "C", "O", "OXT"}
STD = set("ARNDCQEGHILKMFPSTWYV")


def residues(st):
    out = {}
    for ch in st[0]:
        seen, rr = set(), []
        for r in ch:
            key = (r.seqid.num, r.seqid.icode.strip())
            if not any(a.name == "CA" for a in r) or key in seen:
                continue
            seen.add(key); rr.append(r)
        if rr:
            out[ch.name] = rr
    return out


def seq(rr): return "".join(AA3_TO_AA1.get(r.name, "X") for r in rr)
def ca(r): a = next(x for x in r if x.name == "CA"); return np.array([a.pos.x, a.pos.y, a.pos.z])
def heavy_atoms(r, sc=False):
    return [np.array([a.pos.x, a.pos.y, a.pos.z]) for a in r
            if a.element.name not in ("H", "D") and not (sc and a.name in BACKBONE)]


def fit(p, q):
    pc, qc = p.mean(0), q.mean(0); u, s, vt = np.linalg.svd((p - pc).T @ (q - qc)); R = u @ vt
    if np.linalg.det(R) < 0: u[:, -1] *= -1; R = u @ vt
    t = qc - pc @ R
    return R, t, float(np.sqrt(np.mean(np.sum((p @ R + t - q) ** 2, axis=1))))


def inputseq(p):
    d = json.loads(Path(p).read_text()); out = {}
    for s in d["sequences"]:
        if "protein" in s:
            v = s["protein"]
            for c in (v["id"] if isinstance(v["id"], list) else [v["id"]]): out[c] = v["sequence"]
    return out


def crystal_mol(pdbid):
    m = Chem.MolFromMolFile(str(CASF_LIGANDS / f"{pdbid}_ligand.sdf"))
    if m is None:
        m = Chem.MolFromMolFile(str(CASF_LIGANDS / f"{pdbid}_ligand.sdf"), sanitize=False)
        m.UpdatePropertyCache(strict=False); Chem.GetSymmSSSR(m)
    return m


def _embed(want, have, rightmost=False):
    """Map every position of `have` onto `want` as an exact-identity subsequence
    (`have` may only DELETE residues relative to `want`). Returns list of want
    indices, or None. Leftmost and rightmost embeddings bracket gap ambiguity."""
    if rightmost:
        r = _embed(want[::-1], have[::-1])
        return None if r is None else [len(want) - 1 - i for i in r[::-1]]
    out, j = [], 0
    for c in have:
        while j < len(want) and want[j] != c: j += 1
        if j == len(want): return None
        out.append(j); j += 1
    return out


def correspondence(crystal_res, target_res, wt, mut):
    """Crystal chain -> target chain map with residue index maps.

    mode "exact": target chain sequence equals the variant input (non-standard
    residues removed) — the audit's rule. mode "gapped_exact": the target only
    LACKS residues (prepared receptors drop unresolvable ones); every mapped
    residue is identical. Ambiguous gap placement is resolved later, against
    the pocket shell, by comparing leftmost and rightmost embeddings.
    Never maps a mismatched residue; never falls back to a nearest chain."""
    cmap, imap, changed, modes, amb = {}, {}, [], {}, {}
    for cid, w in wt.items():
        if cid not in crystal_res or seq(crystal_res[cid]) != w:
            raise ValueError(f"crystal chain {cid} != WT input")
        m = mut[cid]
        if len(m) != len(w): raise ValueError("input length change")
        keep_in = [i for i, a in enumerate(m) if a in STD]
        want = "".join(m[i] for i in keep_in)
        best = None
        for k, rr in target_res.items():
            tseq = seq(rr); tkeep = [j for j, a in enumerate(tseq) if a in STD]
            have = "".join(tseq[j] for j in tkeep)
            if have == want:
                cand = ("exact", k, {keep_in[a]: tkeep[a] for a in range(len(want))}, set())
            else:
                L = _embed(want, have); Rr = _embed(want, have, rightmost=True)
                if L is None or len(have) < 0.8 * len(want): continue
                ambiguous = {keep_in[a] for a, b in zip(L, Rr) if a != b} | {keep_in[b] for a, b in zip(L, Rr) if a != b}
                cand = ("gapped_exact", k, {keep_in[a]: tkeep[n] for n, a in enumerate(L)}, ambiguous)
            if best is None or (cand[0] == "exact" and best[0] != "exact") or (cand[1] == cid and best[1] != cid and cand[0] == best[0]):
                best = cand
        if best is None: raise ValueError(f"no exact or gapped-exact target chain for {cid}")
        modes[cid], cmap[cid], imap[cid], amb[cid] = best[0], best[1], best[2], best[3]
        changed += [(cid, i) for i, (a, b) in enumerate(zip(w, m)) if a != b]
    missing = [(c, i) for c, i in changed if i not in imap[c]]
    if missing: raise ValueError(f"mutated site absent from target: {missing[:3]}")
    return cmap, imap, changed, modes, amb


def frame(crystal_res, target_res, cmap, imap, native, lig, amb=frozenset()):
    shell_all = [i for i, r in enumerate(crystal_res[native])
                 if min(np.linalg.norm(x - lig, axis=1).min() for x in heavy_atoms(r)) <= 5.0]
    if set(shell_all) & set(amb): raise ValueError("gap placement ambiguous inside the pocket shell")
    shell = [i for i in shell_all if i in imap[native]]
    if len(shell) < 6: raise ValueError(f"too few shell CAs ({len(shell)})")
    q = np.array([ca(crystal_res[native][i]) for i in shell])
    p = np.array([ca(target_res[cmap[native]][imap[native][i]]) for i in shell])
    R, t, r = fit(p, q)
    return R, t, r, len(shell)


def pose_source(pdbid, variant, method, gfix):
    v = OUT / pdbid / variant
    if method == "GNINA":
        p = gfix.get((pdbid, variant), v / "gnina/poses.sdf"); return "sdf", Path(p), v / "docking/receptor.pdb"
    if method == "UniDock2": return "sdf", v / "unidock2/poses.sdf", v / "docking/receptor.pdb"
    if method == "SurfDock": return "sdf", v / "surfdock/poses.sdf", v / "docking/receptor.pdb"
    if method == "ICM":
        rec = ICM_ALIGNED / pdbid / variant / "docking/receptor_aligned.pdb"
        if variant == "wt" or not rec.exists(): rec = v / "docking/receptor.pdb"
        return "sdf", OUT / "_icm_poses" / pdbid / variant / "ICM" / f"D_{pdbid}_{variant}_pose1.sdf", rec
    if method == "Boltz-2": c = v / f"{pdbid}_{variant}_model_0.cif"; return "cif", c, c
    if method == "AF3+MSA": c = v / f"af3msa_{pdbid}_{variant}_model_0.cif"; return "cif", c, c


def cofold_ligand(cif, smiles):
    from casf_mutagenesis.analysis import _predicted_ligand
    block, mol = _predicted_ligand(gemmi.read_structure(str(cif)), smiles)
    if block is None or mol is None: raise ValueError("no predicted ligand")
    return mol


def score_system(pdbid):
    rows = []
    gfix = {}
    if GNINA_FIX.exists():
        for r in csv.DictReader(GNINA_FIX.open()):
            if r["chemistry_audit_passed"] == "True": gfix[(r["system"], r["variant"])] = r["corrected_output"]
    try:
        cst = gemmi.read_structure(str(CASF_RAW / pdbid / f"{pdbid}_protein.pdb")); cres = residues(cst)
        xmol = crystal_mol(pdbid); lig = Chem.RemoveHs(xmol, sanitize=False).GetConformer().GetPositions()
        wt = inputseq(OUT / pdbid / "wt/af3.json")
        native = max(cres, key=lambda c: sum(1 for r in cres[c] for x in heavy_atoms(r)
                                              if np.linalg.norm(x - lig, axis=1).min() <= 5.0))
    except Exception as e:
        return [dict(system=pdbid, status="system_error", error=repr(e)[:200])]
    smiles = Chem.MolToSmiles(Chem.RemoveHs(xmol, sanitize=False))
    for variant in VARIANTS:
        af3 = OUT / pdbid / variant / "af3.json"
        if not af3.exists(): continue
        mut = inputseq(af3)
        for method in METHODS:
            row = dict(system=pdbid, variant=variant, method=method, native_chain=native)
            try:
                kind, pose, target = pose_source(pdbid, variant, method, gfix)
                if not pose.exists(): row["status"] = "missing_pose"; rows.append(row); continue
                tst = gemmi.read_structure(str(target)); tres = residues(tst)
                cmap, imap, changed, modes, amb = correspondence(cres, tres, wt, mut)
                R, t, carms, nsh = frame(cres, tres, cmap, imap, native, lig, amb[native])
                row["chain_match"] = "exact" if all(v == "exact" for v in modes.values()) else "gapped_exact"
                m = _read_top_pose(pose) if kind == "sdf" else cofold_ligand(pose, smiles)
                if m is None: raise ValueError("pose parse failed")
                if any(a.GetAtomicNum() == 0 for a in m.GetAtoms()): raise ValueError("dummy atom in pose")
                mh = Chem.RemoveHs(m, sanitize=False)
                xyz = mh.GetConformer().GetPositions() @ R + t
                conf = mh.GetConformer()
                for i, p in enumerate(xyz): conf.SetAtomPosition(i, p.tolist())
                try:
                    rmsd, nmatch = _matched_rmsd(xmol, xyz, mh)
                except Exception:
                    # Some engine SDFs (notably ICM) carry bond orders RDKit rejects;
                    # take bond orders from the crystal ligand, coordinates unchanged.
                    from rdkit.Chem import AllChem
                    mh = AllChem.AssignBondOrdersFromTemplate(Chem.RemoveHs(xmol, sanitize=False), mh)
                    rmsd, nmatch = _matched_rmsd(xmol, xyz, mh)
                    row["bond_orders_from_crystal"] = True
                # changed side chains, in the crystal frame
                ch_atoms = [x @ R + t for (cid, i) in changed if i in imap[cid]
                            for x in heavy_atoms(tres[cmap[cid]][imap[cid][i]], sc=True)]
                dmin = float(np.min(np.linalg.norm(np.array(ch_atoms)[:, None] - xyz[None], axis=2))) if ch_atoms else None
                row.update(status="ok", ligand_rmsd_a=round(rmsd, 3), n_matched=nmatch,
                           n_mutations=len(changed), no_op=len(changed) == 0,
                           shell_ca=nsh, shell_ca_rmsd_a=round(carms, 3),
                           target_chain=cmap[native], chain_remapped=cmap[native] != native,
                           pose_min_dist_changed_sc_a=None if dmin is None else round(dmin, 2),
                           pose_centroid_to_native_a=round(float(np.linalg.norm(xyz.mean(0) - lig.mean(0))), 2),
                           gnina_corrected=method == "GNINA" and (pdbid, variant) in gfix)
                if kind == "sdf" and method != "ICM":
                    box = json.loads((OUT / pdbid / variant / "docking/box.json").read_text())
                    nat = (lig - t) @ R.T   # crystal ligand transplanted into the docking-receptor frame
                    row["native_atoms_in_box_frac"] = round(float((np.abs(nat - box["center"]) <= np.array(box["size"]) / 2).all(1).mean()), 3)
            except Exception as e:
                row.update(status="error", error=repr(e)[:200])
            rows.append(row)
    return rows


def summarize(rows):
    ok = [r for r in rows if r.get("status") == "ok"]
    wt_ok = {(r["method"], r["system"]) for r in ok if r["variant"] == "wt" and r["ligand_rmsd_a"] < 2}
    out = []
    for method in METHODS:
        wt = [r for r in ok if r["method"] == method and r["variant"] == "wt"]
        out.append(dict(method=method, variant="wt", n=len(wt), retained=sum(r["ligand_rmsd_a"] < 2 for r in wt),
                        rate=round(sum(r["ligand_rmsd_a"] < 2 for r in wt) / len(wt), 3) if wt else None))
        for v in VARIANTS[1:]:
            cells = [r for r in ok if r["method"] == method and r["variant"] == v and not r["no_op"]]
            cond = [r for r in cells if (method, r["system"]) in wt_ok]
            out.append(dict(method=method, variant=v, n=len(cells), retained=sum(r["ligand_rmsd_a"] < 2 for r in cells),
                            rate=round(sum(r["ligand_rmsd_a"] < 2 for r in cells) / len(cells), 3) if cells else None,
                            n_wtok=len(cond), retained_wtok=sum(r["ligand_rmsd_a"] < 2 for r in cond),
                            rate_wtok=round(sum(r["ligand_rmsd_a"] < 2 for r in cond) / len(cond), 3) if cond else None,
                            n_noop_excluded=sum(1 for r in ok if r["method"] == method and r["variant"] == v and r["no_op"])))
    # matched set: systems with an ok non-no-op row for every method in this variant and wt
    for v in VARIANTS[1:]:
        by = defaultdict(dict)
        for r in ok:
            if r["variant"] in ("wt", v): by[r["system"]][(r["method"], r["variant"])] = r
        common = [s for s, d in by.items() if all((m, x) in d for m in METHODS for x in ("wt", v)) and not d[("GNINA", v)]["no_op"]]
        for method in METHODS:
            ret = sum(by[s][(method, v)]["ligand_rmsd_a"] < 2 for s in common)
            out.append(dict(method=method, variant=v, view="matched_six_method", n=len(common), retained=ret,
                            rate=round(ret / len(common), 3) if common else None))
    return out


def main(argv):
    systems = argv or sorted(p.parent.parent.name for p in OUT.glob("*/wt/af3.json") if not p.parent.parent.name.startswith("_"))
    rows = []
    with ProcessPoolExecutor(max_workers=int(os.environ.get("CONTRASCF_WORKERS", "16"))) as ex:
        for i, rr in enumerate(ex.map(score_system, systems), 1):
            rows += rr
            if i % 25 == 0: print(f"  {i}/{len(systems)} systems", flush=True)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with (OUT / "corrected_pocket_scores.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
    summ = summarize(rows)
    keys = list(dict.fromkeys(k for r in summ for k in r))
    with (OUT / "corrected_pocket_table.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(summ)
    from collections import Counter
    print("row status:", Counter((r.get("method"), r.get("status")) for r in rows))
    print(f"wrote {OUT/'corrected_pocket_scores.csv'} and corrected_pocket_table.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
