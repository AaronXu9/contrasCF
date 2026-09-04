"""Pre-fetch every WT MSA a mutagenesis arm needs, into the shared a3m cache.

Phase 1 of `06_run_af3_msa_subset20.py` fetches MSAs by driving Boltz's
`--use_msa_server`, which needs OUTBOUND NETWORK. CARC GPU compute nodes are
not guaranteed to have any, so a job that fetches and folds in one go can die
in Phase 1 after waiting hours in the queue. Running this on a host that does
have network, then shipping `_msa_cache/` to the cluster, makes the GPU job a
pure cache read.

The cache is keyed by `sha1(sequence)[:16]`, so it is shared across arms: the
ligand arm perturbs the LIGAND and leaves the receptor untouched, so its
sequences are exactly the binding-site arm's and every a3m already fetched
there is a hit.

Run (ligand arm):
    source env/lab.sh
    CONTRASCF_OUTPUTS_ROOT=$CONTRASCF_ROOT/analysis/ligand_mutagenesis/outputs \\
    CONTRASCF_SCOPE=disk \\
        $CONTRASCF_PY analysis/casf_mutagenesis/scripts/36_prefetch_msas.py
"""
from __future__ import annotations
import hashlib
import json
import os
import sys
import time
from pathlib import Path

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
sys.path.insert(0, str(REPO_ROOT / "analysis"))

from casf_mutagenesis.config import OUTPUT_ROOT  # noqa: E402
from casf_mutagenesis.msa_via_boltz import fetch_msa_via_boltz  # noqa: E402

OUTPUTS_ROOT = Path(os.environ.get("CONTRASCF_OUTPUTS_ROOT", str(OUTPUT_ROOT)))
MSA_CACHE = Path(os.environ.get("CONTRASCF_MSA_CACHE",
                                str(OUTPUT_ROOT / "_msa_cache")))
# Must match 06_run_af3_msa_subset20.py, or we would fetch MSAs for systems
# that script will then skip.
MAX_TOTAL_LENGTH = int(os.environ.get("CONTRASCF_MAX_TOTAL_LENGTH", "800"))


def _protein_chains(af3_json: Path) -> list[tuple[str, str]]:
    d = json.loads(af3_json.read_text())
    out: list[tuple[str, str]] = []
    for s in d["sequences"]:
        pr = s.get("protein")
        if not pr:
            continue
        ids = pr["id"]
        for cid in (ids if isinstance(ids, list) else [ids]):
            out.append((str(cid), pr["sequence"]))
    return out


def main() -> int:
    MSA_CACHE.mkdir(parents=True, exist_ok=True)
    # sequence -> the systems that need it, for reporting
    need: dict[str, list[str]] = {}
    n_over = 0
    for wt in sorted(OUTPUTS_ROOT.glob("*/wt/af3.json")):
        pdbid = wt.parent.parent.name
        if pdbid.startswith("_"):
            continue
        chains = _protein_chains(wt)
        if sum(len(s) for _, s in chains) > MAX_TOTAL_LENGTH:
            n_over += 1
            continue
        for _, seq in chains:
            need.setdefault(seq, []).append(pdbid)

    todo = [s for s in need
            if not (MSA_CACHE / f"{hashlib.sha1(s.encode()).hexdigest()[:16]}.a3m").exists()]
    print(f"outputs_root : {OUTPUTS_ROOT}")
    print(f"msa_cache    : {MSA_CACHE}")
    print(f"systems over {MAX_TOTAL_LENGTH} aa (skipped): {n_over}")
    print(f"distinct sequences needed : {len(need)}")
    print(f"  already cached          : {len(need) - len(todo)}")
    print(f"  to fetch                : {len(todo)}\n", flush=True)

    n_ok = n_fail = 0
    failures: list[dict] = []
    for i, seq in enumerate(sorted(todo, key=len), 1):
        key = hashlib.sha1(seq.encode()).hexdigest()[:16]
        t0 = time.time()
        try:
            a3m = fetch_msa_via_boltz(seq, cache_dir=MSA_CACHE)
            n_seqs = sum(1 for ln in a3m.splitlines() if ln.startswith(">"))
            n_ok += 1
            print(f"  [{i}/{len(todo)}] {key} len={len(seq)} "
                  f"n_seqs={n_seqs} ({time.time()-t0:.1f}s) "
                  f"e.g. {need[seq][0]}", flush=True)
        except Exception as exc:
            n_fail += 1
            failures.append({"key": key, "len": len(seq),
                             "systems": need[seq][:5], "error": str(exc)[-400:]})
            print(f"  [{i}/{len(todo)}] {key} len={len(seq)} FAILED: "
                  f"{str(exc)[-200:]}", flush=True)

    log = OUTPUTS_ROOT / "msa_prefetch_log.json"
    log.write_text(json.dumps(
        {"cache": str(MSA_CACHE), "needed": len(need), "fetched": n_ok,
         "failed": n_fail, "failures": failures}, indent=2))
    print(f"\nfetched={n_ok} failed={n_fail}  cache now has "
          f"{len(list(MSA_CACHE.glob('*.a3m')))} a3m files")
    print(f"log: {log}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
