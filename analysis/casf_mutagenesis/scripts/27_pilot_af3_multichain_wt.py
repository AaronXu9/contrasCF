"""WT-only AF3+MSA pilot for the multi-chain systems (issue 4e).

Before the 2026-08-23 fix, `06_run_af3_msa_subset20.py` handed AF3 a
single-chain job for every system, so all 55 multi-chain systems were folded
with only their first protein chain. This re-runs **WT only** for those systems
with the fixed multi-chain input, so the size of the correction can be measured
before committing GPU time to all four variants.

Safety: the existing `af3msa_*` artifacts are **moved** (not copied) into
`outputs/_af3_mc_pilot_backup/<pdbid>/wt/` first. That preserves the
single-chain baseline for the before/after comparison AND clears the runner's
skip-existing guard. `--restore` puts them back.

Usage:
  # 1. see what would run
  python .../27_pilot_af3_multichain_wt.py --dry-run
  # 2. smoke-test the smallest system, with timing
  python .../27_pilot_af3_multichain_wt.py --limit 1
  # 3. full pilot
  python .../27_pilot_af3_multichain_wt.py
  # 4. compare single-chain baseline vs multi-chain re-run
  python .../27_pilot_af3_multichain_wt.py --compare
  # rollback
  python .../27_pilot_af3_multichain_wt.py --restore
"""
from __future__ import annotations
import argparse
import importlib.util
import json
import os
import shutil
import sys
import time
from pathlib import Path

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
sys.path.insert(0, str(REPO_ROOT / "analysis"))

from casf_mutagenesis.config import OUTPUT_ROOT  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "af3msa_runner", REPO_ROOT / "analysis/casf_mutagenesis/scripts/06_run_af3_msa_subset20.py")
_R = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_R)

BACKUP = OUTPUT_ROOT / "_af3_mc_pilot_backup"
LOG = OUTPUT_ROOT / "af3_mc_pilot_log.json"
MAX_TOTAL_LENGTH = _R.MAX_TOTAL_LENGTH


def multichain_systems() -> list[tuple[str, int, int]]:
    """[(pdbid, n_chains, total_len)] for systems with >1 protein chain."""
    out = []
    for p in sorted(OUTPUT_ROOT.glob("*/wt/af3.json")):
        pdbid = p.parts[-3]
        try:
            ch = _R._read_af3_protein_chains(p)
        except Exception:
            continue
        if len(ch) > 1:
            out.append((pdbid, len(ch), sum(len(s) for _, s in ch)))
    return sorted(out, key=lambda r: r[2])


def backup_variant(pdbid: str, variant: str) -> int:
    """Move the single-chain artifacts aside. Returns -1 if already backed up.

    A backup that already exists is the ORIGINAL single-chain baseline. Moving
    on top of it would overwrite that baseline with the multi-chain re-run and
    silently destroy the comparison, so refuse instead.
    """
    dst = BACKUP / pdbid / variant
    if dst.is_dir() and any(dst.glob("af3msa_*")):
        return -1
    src = OUTPUT_ROOT / pdbid / variant
    dst.mkdir(parents=True, exist_ok=True)
    n = 0
    for f in list(src.glob("af3msa_*")) + list(src.glob("af3_msa.json")):
        shutil.move(str(f), str(dst / f.name))
        n += 1
    return n


def restore_variant(pdbid: str, variant: str) -> int:
    src = BACKUP / pdbid / variant
    if not src.is_dir():
        return 0
    dst = OUTPUT_ROOT / pdbid / variant
    dst.mkdir(parents=True, exist_ok=True)
    n = 0
    for f in list(src.iterdir()):
        shutil.move(str(f), str(dst / f.name))
        n += 1
    return n


def run_one(pdbid: str, variant: str = "wt") -> dict:
    """Fetch per-chain MSAs, render the multi-chain job, run AF3.

    The MSA is always fetched for the **WT** chain sequence and its query row
    rewritten onto this variant's (possibly mutated) sequence — the same
    contract as the fixed runner. Fetching an MSA for a mutated sequence would
    query the server with a non-natural sequence.
    """
    v_dir = OUTPUT_ROOT / pdbid / variant
    af3_json = v_dir / "af3.json"
    wt_json = OUTPUT_ROOT / pdbid / "wt" / "af3.json"
    prefix = f"{pdbid}_{variant}"
    rec = {"pdbid": pdbid, "variant": variant, "prefix": prefix}

    if not af3_json.exists() or not wt_json.exists():
        rec["status"] = "missing_input"
        return rec
    chains = _R._read_af3_protein_chains(af3_json)
    wt_chains = dict(_R._read_af3_protein_chains(wt_json))
    total_len = sum(len(s) for _, s in chains)
    rec.update(n_chains=len(chains), total_len=total_len)
    if total_len > MAX_TOTAL_LENGTH:
        rec["status"] = "skip_too_long"
        return rec

    t0 = time.time()
    try:
        by_seq: dict[str, str] = {}
        chain_msas: dict[str, str] = {}
        for cid, seq in chains:
            wt_seq = wt_chains.get(cid)
            if wt_seq is None:
                raise RuntimeError(f"chain {cid} not in WT af3.json")
            if wt_seq not in by_seq:
                by_seq[wt_seq] = _R.fetch_msa_via_boltz(wt_seq, cache_dir=_R.MSA_CACHE)
            a3m = _R.rewrite_a3m_query(by_seq[wt_seq], seq)
            chain_msas[cid] = _R.clean_a3m_for_af3(a3m, len(seq))
        rec["msa_s"] = round(time.time() - t0, 1)
        rec["distinct_seqs"] = len(by_seq)
    except Exception as exc:
        rec.update(status="error", error=f"msa: {exc}")
        return rec

    json_msa = v_dir / "af3_msa.json"
    _R.render_af3(name=prefix, chain_seqs=chains,
                  ligand_smiles=_R._read_ligand_smiles(af3_json),
                  out_path=json_msa, chain_msas=chain_msas)

    t1 = time.time()
    work = v_dir / "_af3_msa_work"
    if work.exists():
        shutil.rmtree(work)
    try:
        af3_sys = _R._run_af3(json_msa, work)
        _R._copy_all_samples(af3_sys, prefix, v_dir, "af3msa")
        rec.update(status="ok", af3_s=round(time.time() - t1, 1))
    except Exception as exc:
        rec.update(status="error", error=f"af3: {str(exc)[:300]}",
                   af3_s=round(time.time() - t1, 1))
    return rec


def compare() -> int:
    """Single-chain baseline (backup) vs multi-chain re-run, same RMSD code."""
    from casf_mutagenesis.analysis import analyze_prediction
    import numpy as np

    rows = []
    for pdbid_dir in sorted(BACKUP.iterdir()):
        pdbid = pdbid_dir.name
        new_cif = OUTPUT_ROOT / pdbid / "wt" / f"af3msa_{pdbid}_wt_model_0.cif"
        if not new_cif.exists():
            continue
        try:
            r = analyze_prediction(pdbid, "wt", "AF3+MSA")
            rows.append({"pdbid": pdbid, "new_rmsd": r.ligand_rmsd_a,
                         "status": r.status})
        except Exception as exc:
            rows.append({"pdbid": pdbid, "new_rmsd": None, "status": str(exc)[:80]})

    import csv
    old = {}
    with open(OUTPUT_ROOT / "results_full.csv") as f:
        for r in csv.DictReader(f):
            if r["model"] == "AF3+MSA" and r["variant"] == "wt" and r["pose_idx"] == "0":
                try:
                    old[r["pdbid"]] = float(r["ligand_rmsd_a"])
                except (TypeError, ValueError):
                    pass

    ok = [r for r in rows if r["new_rmsd"] is not None and r["pdbid"] in old]
    if not ok:
        print("no comparable systems yet — run the pilot first")
        return 1
    o = np.array([old[r["pdbid"]] for r in ok])
    n = np.array([r["new_rmsd"] for r in ok])
    print(f"systems compared: {len(ok)}\n")
    print(f"{'':22}{'single-chain':>14}{'multi-chain':>14}")
    print(f"{'rate < 2 A':22}{(o<2).mean():>14.3f}{(n<2).mean():>14.3f}")
    print(f"{'median RMSD (A)':22}{np.median(o):>14.2f}{np.median(n):>14.2f}")
    print(f"{'mean RMSD (A)':22}{o.mean():>14.2f}{n.mean():>14.2f}")
    improved = int((n < o - 0.5).sum()); worse = int((n > o + 0.5).sum())
    print(f"\nimproved >0.5 A: {improved}   worse >0.5 A: {worse}   "
          f"unchanged: {len(ok)-improved-worse}")
    print("\nlargest improvements:")
    for i in np.argsort(n - o)[:8]:
        print(f"   {ok[i]['pdbid']}: {o[i]:.2f} -> {n[i]:.2f} A")
    json.dump(rows, open(OUTPUT_ROOT / "af3_mc_pilot_compare.json", "w"), indent=2)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="run only the N smallest")
    ap.add_argument("--variants", default="wt",
                    help="comma-separated: wt,rem,pack,inv (default wt)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--restore", action="store_true")
    ap.add_argument("--compare", action="store_true")
    args = ap.parse_args()
    variants = [v.strip() for v in args.variants.split(",") if v.strip()]

    sysrows = multichain_systems()
    runnable = [r for r in sysrows if r[2] <= MAX_TOTAL_LENGTH]
    if args.limit:
        runnable = runnable[: args.limit]

    if args.restore:
        n = 0
        if BACKUP.is_dir():
            for p in sorted(BACKUP.iterdir()):
                for v in variants:
                    n += restore_variant(p.name, v)
        print(f"restored {n} files from {BACKUP} (variants: {','.join(variants)})")
        return 0
    if args.compare:
        return compare()

    print(f"multi-chain systems: {len(sysrows)}  runnable (<= {MAX_TOTAL_LENGTH} aa): "
          f"{len([r for r in sysrows if r[2] <= MAX_TOTAL_LENGTH])}")
    skipped = [r for r in sysrows if r[2] > MAX_TOTAL_LENGTH]
    if skipped:
        print(f"  over length gate (skipped): {', '.join(f'{s}({l}aa)' for s,_,l in skipped)}")
    print(f"this run: {len(runnable)} systems\n")
    if args.dry_run:
        for s, nc, l in runnable:
            print(f"   {s}: {nc} chains, {l} aa")
        return 0

    results = []
    if LOG.exists():                       # resume: keep earlier variants' records
        try:
            results = json.loads(LOG.read_text())
        except Exception:
            results = []
    n_cells = len(runnable) * len(variants)
    k = 0
    for variant in variants:
        for (pdbid, nc, l) in runnable:
            k += 1
            nb = backup_variant(pdbid, variant)
            cif = OUTPUT_ROOT / pdbid / variant / f"af3msa_{pdbid}_{variant}_model_0.cif"
            if nb < 0:
                if cif.exists():
                    print(f"[{k}/{n_cells}] {pdbid}/{variant} — already re-run, skip",
                          flush=True)
                    results.append({"pdbid": pdbid, "variant": variant,
                                    "status": "skip_done"})
                    LOG.write_text(json.dumps(results, indent=2))
                    continue
                print(f"[{k}/{n_cells}] {pdbid}/{variant} — baseline backed up, "
                      f"re-running", flush=True)
            else:
                print(f"[{k}/{n_cells}] {pdbid}/{variant} ({nc} chains, {l} aa) "
                      f"— backed up {nb} files", flush=True)
            rec = run_one(pdbid, variant)
            results.append(rec)
            print(f"     -> {rec.get('status')} "
                  f"msa={rec.get('msa_s','-')}s af3={rec.get('af3_s','-')}s", flush=True)
            LOG.write_text(json.dumps(results, indent=2))
    ok = sum(1 for r in results if r.get("status") == "ok")
    tot_af3 = sum(r.get("af3_s", 0) or 0 for r in results)
    print(f"\ndone: {ok}/{len(results)} ok; AF3 wall {tot_af3/60:.1f} min "
          f"(mean {tot_af3/max(ok,1):.0f} s/system)")
    print(f"log: {LOG}")
    return 0 if ok == len(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
