#!/usr/bin/env python
"""Recover the CASF-2016 core receptors that HiQBind does not carry.

Item 1 of `docs/data_prep_todo.md`. 34 of the 285 CASF-2016 core ids have no
entry under `$CASF/raw/`, because our receptor source -- HiQBind -- lacks
them. They split cleanly in two:

  * **31 are ordinary small-molecule complexes** whose GT crystal ligand is
    ALREADY on disk at `$CASF/crystal_ligands/<id>_ligand.sdf` (those were
    fetched from RCSB ModelServer, not from HiQBind). Only the receptor is
    missing. This script rebuilds it directly from the RCSB mmCIF.
  * **3 have PEPTIDE ligands** and are deliberately NOT recovered:
    `1a30` (3-mer), `3bv9` (6-mer), `3uri` (8-mer) -- `3uri` has zero
    non-polymer entities at all. HiQBind filing two of them under
    `raw_data_hiq_poly` was a categorical call, not a quality veto. The
    ligand machinery here (SMILES cache, RDKit MCS matching, the
    ligand-mutagenesis rules) assumes small molecules.

The 31 are absent from HiQBind's 32,275-row metadata entirely, so recovering
them does not override a HiQBind quality judgement -- HiQBind never assessed
them.

**Receptor construction** mirrors HiQBind's own protein definition -- "chains
within 10 Angstrom of the ligand structure" -- so chain selection is identical
to that of the 251. It differs in exactly one documented way: HiQBind ships
`*_protein_refined.pdb` (PDBBind-Opt -- hydrogens added, missing atoms and
residues rebuilt) whereas this writes deposited coordinates, unrefined.
`--controls` measures how much that difference actually costs, by running the
identical recovery on ids HiQBind DOES have and diffing the result against the
refined receptor (residue count and the detected 3.5 A pocket set).

Usage:
    source env/lab.sh
    $CONTRASCF_PY .../25_recover_missing_receptors.py --controls 15   # validate
    $CONTRASCF_PY .../25_recover_missing_receptors.py --recover       # write the 31
    $CONTRASCF_PY .../25_recover_missing_receptors.py                 # dry run
"""
from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import gemmi
import numpy as np

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
sys.path.insert(0, str(REPO_ROOT / "analysis"))

from casf_mutagenesis.config import (  # noqa: E402
    CASF_LIGANDS, CASF_RAW, MODULE_ROOT, POCKET_CUTOFF_A,
)
from casf_mutagenesis.pocket import detect_pocket, load_ligand_heavy_coords  # noqa: E402

# --- the 34, and why each is where it is ----------------------------------

PEPTIDE_LIGAND = {"1a30": "3-mer", "3bv9": "6-mer", "3uri": "8-mer"}

MISSING_34 = [
    "1a30", "1c5z", "2p15", "2vvn", "2wca", "2zcr", "2zda", "2zy1", "3acw",
    "3arp", "3arq", "3aru", "3arv", "3ary", "3b65", "3bv9", "3dxg", "3f3a",
    "3kgp", "3myg", "3o9i", "3prs", "3pww", "3qqs", "3r88", "3twp", "3uex",
    "3uri", "3zsx", "4abg", "4kzq", "4kzu", "4owm", "4pcs",
]
RECOVERABLE = [i for i in MISSING_34 if i not in PEPTIDE_LIGAND]  # 31

# --- construction parameters ----------------------------------------------

CHAIN_CUTOFF_A = 10.0          # HiQBind's protein definition
LIGAND_MATCH_TOL_A = 2.0       # centroid tolerance when locating the ligand copy
RCSB_URL = "https://files.rcsb.org/download/{}.cif"
CACHE_DIR = MODULE_ROOT / "_rcsb_cache"
PROVENANCE_JSON = MODULE_ROOT / "recovered_receptors.json"
STATUS_JSON = CASF_LIGANDS / "download_status.json"
IDEAL_BACKUP = CACHE_DIR / "ideal_ligand_backup"
MODELSERVER_URL = ("https://models.rcsb.org/v1/{pdb}/ligand"
                   "?auth_comp_id={het}&encoding=sdf&copy_all_categories=false")


def ideal_sourced_ids() -> set:
    """Ids whose crystal_ligands SDF is an IDEALISED CCD template, not a pose.

    `download_crystal_ligands.py` queries RCSB ModelServer by PDBbind's
    `ligand_name`; when that disagrees with the deposited het code the query
    misses and it silently falls back to `<het>_ideal.sdf`, whose coordinates
    are origin-centred and carry no crystal frame. 168 entries corpus-wide are
    like this. None are among the 251 currently in `raw/` -- but `1c5z` is one,
    so recovering it naively would put the first garbage GT pose into the
    benchmark.
    """
    if not STATUS_JSON.exists():
        return set()
    d = json.loads(STATUS_JSON.read_text())
    return {k for k, v in d.items()
            if isinstance(v, dict) and v.get("source") == "ideal"}


def fetch_cif(pdbid: str, retries: int = 3) -> Path:
    """Download (and cache) the deposited mmCIF for `pdbid`."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    dest = CACHE_DIR / f"{pdbid}.cif"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    url = RCSB_URL.format(pdbid.upper())
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                data = r.read()
            if not data:
                raise RuntimeError("empty response")
            dest.write_bytes(data)
            return dest
        except (urllib.error.URLError, RuntimeError, TimeoutError) as exc:
            last = exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"RCSB fetch failed for {pdbid}: {last}")


def _heavy_coords(res: gemmi.Residue) -> np.ndarray:
    pts = [(a.pos.x, a.pos.y, a.pos.z) for a in res if a.element.name != "H"]
    return np.asarray(pts, dtype=float) if pts else np.empty((0, 3))


def nonpolymer_residues(model: gemmi.Model):
    """Yield (chain_name, residue) for residues outside any polymer.

    Membership in the polymer -- not residue chemistry -- is the right test.
    `3f3a` binds a FREE TRYPTOPHAN as its ligand; filtering on
    `is_amino_acid()` would discard exactly the molecule we are looking for.
    """
    for chain in model:
        poly_ids = {(r.seqid.num, r.seqid.icode, r.name)
                    for r in chain.get_polymer()}
        for res in chain:
            if res.is_water():
                continue
            if (res.seqid.num, res.seqid.icode, res.name) in poly_ids:
                continue
            yield chain.name, res


def repair_ideal_ligand(pdbid: str, model: gemmi.Model, n_heavy: int,
                        dest_sdf: Path, write: bool) -> dict:
    """Re-fetch a real crystal pose for an idealised-SDF system.

    Identifies the ligand by heavy-atom count (the ideal SDF's coordinates are
    unusable for matching) and re-queries ModelServer with the het code that is
    actually deposited.
    """
    names = {}
    for _, res in nonpolymer_residues(model):
        if _heavy_coords(res).shape[0] == n_heavy:
            names[res.name] = names.get(res.name, 0) + 1
    if len(names) != 1:
        raise RuntimeError(
            f"cannot disambiguate the ligand by heavy-atom count ({n_heavy}): "
            f"candidates {sorted(names) or 'none'}"
        )
    het = next(iter(names))
    url = MODELSERVER_URL.format(pdb=pdbid.upper(), het=het.upper())
    with urllib.request.urlopen(url, timeout=60) as r:
        data = r.read()
    if not data or b"$$$$" not in data:
        raise RuntimeError(f"ModelServer returned no pose for {pdbid}/{het}")
    if write:
        IDEAL_BACKUP.mkdir(parents=True, exist_ok=True)
        if dest_sdf.exists():
            (IDEAL_BACKUP / dest_sdf.name).write_bytes(dest_sdf.read_bytes())
        dest_sdf.write_bytes(data)
    else:
        dest_sdf = CACHE_DIR / "_dryrun" / dest_sdf.name
        dest_sdf.parent.mkdir(parents=True, exist_ok=True)
        dest_sdf.write_bytes(data)
    return {"het_code": het, "repaired_sdf": str(dest_sdf),
            "backup": str(IDEAL_BACKUP / dest_sdf.name) if write else None}


def locate_ligand(model: gemmi.Model, lig_xyz: np.ndarray) -> tuple:
    """Find the residue in `model` that IS the ligand in the reference SDF.

    The crystal_ligands SDF carries deposited coordinates, so the true copy
    sits essentially on top of it. Match on centroid, break ties on
    heavy-atom count -- this pins the receptor to the same ligand copy the
    pipeline will later use for pocket detection.
    """
    target_cen = lig_xyz.mean(axis=0)
    n_target = len(lig_xyz)
    best = None
    for cname, res in nonpolymer_residues(model):
        hv = _heavy_coords(res)
        if hv.shape[0] == 0:
            continue
        dist = float(np.linalg.norm(hv.mean(axis=0) - target_cen))
        key = (dist, abs(hv.shape[0] - n_target))
        if best is None or key < best[0]:
            best = (key, cname, res.name, str(res.seqid.num), hv)
    if best is None:
        raise RuntimeError("no candidate ligand residue found")
    (dist, dn), cname, rname, rnum, hv = best
    if dist > LIGAND_MATCH_TOL_A:
        raise RuntimeError(
            f"no ligand within {LIGAND_MATCH_TOL_A} A of the reference SDF "
            f"centroid (closest {rname} {cname}{rnum} at {dist:.2f} A)"
        )
    return cname, rname, rnum, hv, dist, dn


def build_receptor(pdbid: str, ligand_sdf: Path, out_pdb: Path,
                   write_store: bool = False) -> dict:
    """Write a protein-only receptor for `pdbid` next to its reference ligand.

    Chains are selected by HiQBind's rule: any polymer chain with a heavy atom
    within CHAIN_CUTOFF_A of the ligand. Waters, hydrogens, alternate
    conformations and all non-polymer residues are dropped.
    """
    lig_xyz = load_ligand_heavy_coords(ligand_sdf)
    cif = fetch_cif(pdbid)

    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    model = st[0]

    # An idealised SDF has no crystal frame, so it cannot anchor the receptor
    # (and would be a worthless RMSD reference downstream). Repair it first.
    repair = None
    if pdbid in ideal_sourced_ids():
        repair = repair_ideal_ligand(pdbid, model, len(lig_xyz), ligand_sdf,
                                     write=write_store)
        ligand_sdf = Path(repair["repaired_sdf"])
        lig_xyz = load_ligand_heavy_coords(ligand_sdf)

    cname, rname, rnum, lig_hv, match_dist, count_delta = locate_ligand(model, lig_xyz)

    # Which polymer chains touch that ligand copy?
    keep: list[str] = []
    for chain in model:
        poly = chain.get_polymer()
        if len(poly) == 0:
            continue
        pts = [(a.pos.x, a.pos.y, a.pos.z)
               for res in poly for a in res if a.element.name != "H"]
        if not pts:
            continue
        arr = np.asarray(pts, dtype=float)
        d2 = ((arr[:, None, :] - lig_hv[None, :, :]) ** 2).sum(-1)
        if float(np.sqrt(d2.min())) <= CHAIN_CUTOFF_A:
            keep.append(chain.name)

    if not keep:
        raise RuntimeError(f"no polymer chain within {CHAIN_CUTOFF_A} A of the ligand")

    out = gemmi.Structure()
    out.spacegroup_hm = st.spacegroup_hm
    out.cell = st.cell
    out.name = pdbid.upper()
    m = gemmi.Model("1")
    n_res = 0
    for chain in model:
        if chain.name not in keep:
            continue
        nc = gemmi.Chain(chain.name)
        for res in chain.get_polymer():
            nc.add_residue(res)
            n_res += 1
        if len(nc):
            m.add_chain(nc)
    out.add_model(m)
    out.setup_entities()

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    out.write_pdb(str(out_pdb))

    return {
        "pdbid": pdbid,
        "source": "RCSB",
        "cif_url": RCSB_URL.format(pdbid.upper()),
        "ligand_resname": rname,
        "ligand_chain": cname,
        "ligand_seqid": rnum,
        "ligand_centroid_offset_a": round(match_dist, 3),
        "ligand_heavy_count_delta": count_delta,
        "chains_kept": keep,
        "n_residues": n_res,
        "refined": False,
        "chain_cutoff_a": CHAIN_CUTOFF_A,
        "ligand_sdf_used": str(ligand_sdf),
        "ideal_sdf_repaired": repair,
    }


# --- validation gate -------------------------------------------------------

def _pocket_set(pdb: Path, sdf: Path) -> set:
    return {(r.chain, r.resnum, r.ins_code, r.aa3)
            for r in detect_pocket(pdb, sdf, cutoff=POCKET_CUTOFF_A)}


def _n_residues(pdb: Path) -> int:
    st = gemmi.read_structure(str(pdb))
    st.remove_hydrogens()
    return sum(1 for ch in st[0] for _ in ch)


def run_controls(n: int, seed: int = 42) -> int:
    """Rebuild `n` receptors that HiQBind DOES have, and diff against it.

    This is the gate: if the RCSB-derived receptor reproduces HiQBind's
    detected pocket, the 31 recovered systems are comparable to the 251.
    """
    have = sorted(p.name for p in CASF_RAW.iterdir()
                  if (p / f"{p.name}_protein.pdb").exists()
                  and (CASF_LIGANDS / f"{p.name}_ligand.sdf").exists())
    random.Random(seed).shuffle(have)
    picked = have[:n]
    tmp = CACHE_DIR / "_controls"
    tmp.mkdir(parents=True, exist_ok=True)

    print(f"CONTROLS -- rebuilding {len(picked)} HiQBind-backed receptors from RCSB")
    print(f"{'id':6s} {'ref_res':>7s} {'new_res':>7s} {'ref_pkt':>7s} {'new_pkt':>7s} "
          f"{'shared':>6s} {'jacc':>5s}  verdict")

    exact = n_close = 0
    rows = []
    for pdbid in picked:
        sdf = CASF_LIGANDS / f"{pdbid}_ligand.sdf"
        ref = CASF_RAW / pdbid / f"{pdbid}_protein.pdb"
        new = tmp / f"{pdbid}_protein.pdb"
        try:
            build_receptor(pdbid, sdf, new)
            ref_pkt, new_pkt = _pocket_set(ref, sdf), _pocket_set(new, sdf)
            inter = ref_pkt & new_pkt
            union = ref_pkt | new_pkt
            jac = len(inter) / len(union) if union else 1.0
            ok = ref_pkt == new_pkt
            exact += ok
            n_close += (jac >= 0.8)
            verdict = "EXACT" if ok else ("close" if jac >= 0.8 else "DIVERGES")
            rows.append({"pdbid": pdbid, "jaccard": round(jac, 3),
                         "ref_pocket": len(ref_pkt), "new_pocket": len(new_pkt),
                         "exact": ok})
            print(f"{pdbid:6s} {_n_residues(ref):7d} {_n_residues(new):7d} "
                  f"{len(ref_pkt):7d} {len(new_pkt):7d} {len(inter):6d} "
                  f"{jac:5.2f}  {verdict}")
        except Exception as exc:  # noqa: BLE001 -- report, don't abort the sweep
            rows.append({"pdbid": pdbid, "error": str(exc)})
            print(f"{pdbid:6s} {'-':>7s} {'-':>7s} {'-':>7s} {'-':>7s} {'-':>6s} "
                  f"{'-':>5s}  ERROR: {exc}")

    n_ok = sum(1 for r in rows if "error" not in r)
    print()
    print(f"built {n_ok}/{len(picked)}; pocket EXACT {exact}/{n_ok}, "
          f"Jaccard>=0.8 {n_close}/{n_ok}")
    (CACHE_DIR / "controls_report.json").write_text(json.dumps(rows, indent=2))
    print(f"wrote {CACHE_DIR / 'controls_report.json'}")
    return 0 if n_ok == len(picked) else 1


def run_recover(dry_run: bool) -> int:
    print(f"{'RECOVER' if not dry_run else 'DRY RUN'} -- {len(RECOVERABLE)} systems "
          f"(excluding {len(PEPTIDE_LIGAND)} peptide-ligand: "
          f"{', '.join(f'{k} [{v}]' for k, v in PEPTIDE_LIGAND.items())})")
    records, failures = [], []
    for i, pdbid in enumerate(RECOVERABLE, 1):
        sdf = CASF_LIGANDS / f"{pdbid}_ligand.sdf"
        if not sdf.exists():
            failures.append((pdbid, "no reference ligand SDF"))
            print(f"  [{i:2d}/{len(RECOVERABLE)}] {pdbid} SKIP -- no reference ligand")
            continue
        dest = CASF_RAW / pdbid / f"{pdbid}_protein.pdb"
        target = dest if not dry_run else CACHE_DIR / "_dryrun" / f"{pdbid}_protein.pdb"
        try:
            rec = build_receptor(pdbid, sdf, target, write_store=not dry_run)
            pkt = _pocket_set(target, Path(rec["ligand_sdf_used"]))
            rec["n_pocket_residues"] = len(pkt)
            records.append(rec)
            print(f"  [{i:2d}/{len(RECOVERABLE)}] {pdbid} OK  "
                  f"chains={','.join(rec['chains_kept'])} res={rec['n_residues']:4d} "
                  f"pocket={len(pkt):2d} lig={rec['ligand_resname']} "
                  f"(offset {rec['ligand_centroid_offset_a']} A)"
                  + ("  [ideal SDF REPAIRED -> "
                     f"{rec['ideal_sdf_repaired']['het_code']}]"
                     if rec["ideal_sdf_repaired"] else ""))
        except Exception as exc:  # noqa: BLE001
            failures.append((pdbid, str(exc)))
            print(f"  [{i:2d}/{len(RECOVERABLE)}] {pdbid} FAIL -- {exc}")

    print()
    print(f"recovered {len(records)}/{len(RECOVERABLE)}; failures {len(failures)}")
    for pdbid, err in failures:
        print(f"  {pdbid}: {err}")
    if not dry_run and records:
        PROVENANCE_JSON.write_text(json.dumps(
            {"generated_by": Path(__file__).name,
             "note": "RCSB-derived, UNREFINED receptors. The other 251 are "
                     "HiQBind PDBBind-Opt refined. Stratify on `source` before "
                     "pooling.",
             "excluded_peptide_ligand": PEPTIDE_LIGAND,
             "records": records}, indent=2))
        print(f"wrote {PROVENANCE_JSON}")
    return 0 if not failures else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--controls", type=int, metavar="N",
                    help="validate: rebuild N HiQBind-backed receptors and diff")
    ap.add_argument("--recover", action="store_true",
                    help="write the 31 recovered receptors into $CASF/raw/")
    args = ap.parse_args()
    if args.controls:
        return run_controls(args.controls)
    return run_recover(dry_run=not args.recover)


if __name__ == "__main__":
    raise SystemExit(main())
