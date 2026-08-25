"""Multi-seed confirmation for the corrected (multi-chain) AF3+MSA numbers.

The 2026-08-24 re-run moved the AF3+MSA memorization bars a lot (`inv`
0.078 -> 0.549 on the affected systems), but it is a SINGLE seed, and the WT arm
showed real run-to-run variance (5 systems regressed, one 1.77 -> 32.56 A). This
re-runs the same cells under additional AF3 seeds so the corrected rate can be
quoted with a spread instead of a point estimate.

Seed runs are written to a SEPARATE tree so they never collide with the
canonical CIFs:

    outputs/_af3_seeds/seed<N>/<pdbid>/<variant>/af3msa_<prefix>_model_0.cif

RMSD is computed by temporarily swapping the seed CIF into the canonical path
and calling the pipeline's own `analyze_prediction`, so the metric is identical
to every other number in the benchmark. The canonical file is always restored
(try/finally), and the swap refuses to run if the canonical file is missing.

Usage:
  # what would run
  python .../28_af3_multiseed_check.py --seeds 2,3 --variants inv --dry-run
  # run
  python .../28_af3_multiseed_check.py --seeds 2,3 --variants inv
  # summarise across seeds
  python .../28_af3_multiseed_check.py --seeds 2,3 --variants inv --summarize
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

SEED_ROOT = OUTPUT_ROOT / "_af3_seeds"
LOG = OUTPUT_ROOT / "af3_multiseed_log.json"
PILOT_LOG = OUTPUT_ROOT / "af3_mc_pilot_log.json"


def target_cells(variants: list[str]) -> list[tuple[str, str]]:
    """Cells the multi-chain re-run completed — the ones worth re-seeding."""
    if not PILOT_LOG.exists():
        return []
    log = json.loads(PILOT_LOG.read_text())
    return sorted({(r["pdbid"], r["variant"]) for r in log
                   if r.get("status") in ("ok", "skip_done") and r["variant"] in variants})


def run_cell(pdbid: str, variant: str, seed: int) -> dict:
    src = OUTPUT_ROOT / pdbid / variant
    dst = SEED_ROOT / f"seed{seed}" / pdbid / variant
    prefix = f"{pdbid}_{variant}"
    rec = {"pdbid": pdbid, "variant": variant, "seed": seed, "prefix": prefix}
    if (dst / f"af3msa_{prefix}_model_0.cif").exists():
        rec["status"] = "skip_done"
        return rec

    af3_json = src / "af3.json"
    wt_json = OUTPUT_ROOT / pdbid / "wt" / "af3.json"
    if not af3_json.exists() or not wt_json.exists():
        rec["status"] = "missing_input"
        return rec
    chains = _R._read_af3_protein_chains(af3_json)
    wt_chains = dict(_R._read_af3_protein_chains(wt_json))
    total_len = sum(len(s) for _, s in chains)
    rec.update(n_chains=len(chains), total_len=total_len)
    if total_len > _R.MAX_TOTAL_LENGTH:
        rec["status"] = "skip_too_long"
        return rec

    try:
        by_seq, chain_msas = {}, {}
        for cid, seq in chains:
            wt_seq = wt_chains.get(cid)
            if wt_seq is None:
                raise RuntimeError(f"chain {cid} not in WT")
            if wt_seq not in by_seq:
                by_seq[wt_seq] = _R.fetch_msa_via_boltz(wt_seq, cache_dir=_R.MSA_CACHE)
            a3m = _R.rewrite_a3m_query(by_seq[wt_seq], seq)
            chain_msas[cid] = _R.clean_a3m_for_af3(a3m, len(seq))
    except Exception as exc:
        rec.update(status="error", error=f"msa: {exc}")
        return rec

    dst.mkdir(parents=True, exist_ok=True)
    json_msa = dst / "af3_msa.json"
    _R.render_af3(name=prefix, chain_seqs=chains,
                  ligand_smiles=_R._read_ligand_smiles(af3_json),
                  out_path=json_msa, chain_msas=chain_msas, seed=seed)

    t0 = time.time()
    work = dst / "_af3_msa_work"
    if work.exists():
        shutil.rmtree(work)
    try:
        af3_sys = _R._run_af3(json_msa, work)
        _R._copy_all_samples(af3_sys, prefix, dst, "af3msa")
        rec.update(status="ok", af3_s=round(time.time() - t0, 1))
        shutil.rmtree(work, ignore_errors=True)
    except Exception as exc:
        rec.update(status="error", error=f"af3: {str(exc)[:200]}",
                   af3_s=round(time.time() - t0, 1))
    return rec


def rmsd_for(pdbid: str, variant: str, seed: int) -> float | None:
    """RMSD via the pipeline's own analyzer, by temporarily swapping the CIF in."""
    from casf_mutagenesis.analysis import analyze_prediction
    prefix = f"{pdbid}_{variant}"
    canon = OUTPUT_ROOT / pdbid / variant / f"af3msa_{prefix}_model_0.cif"
    seed_cif = SEED_ROOT / f"seed{seed}" / pdbid / variant / f"af3msa_{prefix}_model_0.cif"
    if not seed_cif.exists() or not canon.exists():
        return None
    stash = canon.with_suffix(".cif.__multiseed_stash__")
    try:
        shutil.move(str(canon), str(stash))
        shutil.copy(str(seed_cif), str(canon))
        rec = analyze_prediction(pdbid, variant, "AF3+MSA")
        return rec.ligand_rmsd_a if rec.status == "ok" else None
    except Exception:
        return None
    finally:
        if stash.exists():
            if canon.exists():
                canon.unlink()
            shutil.move(str(stash), str(canon))


def summarize(seeds: list[int], variants: list[str]) -> int:
    import numpy as np
    from casf_mutagenesis.analysis import analyze_prediction

    cells = target_cells(variants)
    per_seed: dict[int, dict[tuple[str, str], float]] = {}
    # seed 1 == the canonical corrected run already on disk
    base: dict[tuple[str, str], float] = {}
    for p, v in cells:
        try:
            r = analyze_prediction(p, v, "AF3+MSA")
            if r.status == "ok" and r.ligand_rmsd_a is not None:
                base[(p, v)] = r.ligand_rmsd_a
        except Exception:
            pass
    per_seed[1] = base
    for s in seeds:
        d = {}
        for p, v in cells:
            x = rmsd_for(p, v, s)
            if x is not None:
                d[(p, v)] = x
        per_seed[s] = d

    print(f"{'variant':<8}{'seed':>6}{'n':>5}{'rate<2A':>10}{'median':>9}")
    print("-" * 40)
    rates: dict[str, list[float]] = {}
    for v in variants:
        for s in [1] + seeds:
            keys = [k for k in per_seed[s] if k[1] == v]
            if not keys:
                continue
            arr = np.array([per_seed[s][k] for k in keys])
            rate = float((arr < 2).mean())
            rates.setdefault(v, []).append(rate)
            print(f"{v:<8}{s:>6}{len(keys):>5}{rate:>10.3f}{np.median(arr):>9.2f}")
        print()
    print("across-seed summary (corrected, multi-chain):")
    for v, rs in rates.items():
        a = np.array(rs)
        print(f"  {v:<6} mean {a.mean():.3f}  sd {a.std(ddof=1) if len(a)>1 else 0:.3f}  "
              f"range [{a.min():.3f}, {a.max():.3f}]  n_seeds={len(a)}")
    json.dump({f"seed{s}": {f"{p}|{v}": x for (p, v), x in d.items()}
               for s, d in per_seed.items()},
              open(OUTPUT_ROOT / "af3_multiseed_rmsd.json", "w"), indent=2)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="2,3")
    ap.add_argument("--variants", default="inv")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--summarize", action="store_true")
    args = ap.parse_args()
    seeds = [int(s) for s in args.seeds.split(",") if s.strip()]
    variants = [v.strip() for v in args.variants.split(",") if v.strip()]

    if args.summarize:
        return summarize(seeds, variants)

    cells = target_cells(variants)
    if args.limit:
        cells = cells[: args.limit]
    n = len(cells) * len(seeds)
    print(f"cells: {len(cells)} x seeds {seeds} = {n} AF3 jobs "
          f"(~{n*110/3600:.1f} h at 110 s/job)")
    if args.dry_run:
        for p, v in cells[:10]:
            print(f"   {p}/{v}")
        return 0

    results = json.loads(LOG.read_text()) if LOG.exists() else []
    k = 0
    for seed in seeds:
        for pdbid, variant in cells:
            k += 1
            rec = run_cell(pdbid, variant, seed)
            results.append(rec)
            print(f"[{k}/{n}] seed{seed} {pdbid}/{variant} -> {rec.get('status')} "
                  f"{rec.get('af3_s','-')}s", flush=True)
            LOG.write_text(json.dumps(results, indent=2))
    ok = sum(1 for r in results if r.get("status") == "ok")
    print(f"\ndone: {ok} ok / {len(results)} records")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
