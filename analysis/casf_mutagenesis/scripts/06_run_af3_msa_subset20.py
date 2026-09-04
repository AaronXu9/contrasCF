"""Run AF3 on subset20 with MSAs (replaces no-MSA results).

Two-phase pipeline:
  1. For each unique WT sequence: fetch a uniref+bfd merged A3M via Boltz's
     `--use_msa_server` machinery (cached by sha1 in outputs/_msa_cache/).
     ~2 s + ~6 s Boltz minimum-inference per system; ~150 s total for 19.
  2. For each (pdbid, variant): take the WT A3M, rewrite the query row at
     mutation positions, render af3.json with the MSA injected, run AF3,
     copy outputs as `af3msa_<prefix>_*` (the no-MSA `af3_*` outputs from
     `04_run_af3_subset20.py` are preserved as the baseline).

Run:
    LD_LIBRARY_PATH=/home/aoxu/miniconda3/envs/rdkit_env/lib:$LD_LIBRARY_PATH \\
        /home/aoxu/miniconda3/envs/rdkit_env/bin/python \\
        analysis/casf_mutagenesis/scripts/06_run_af3_msa_subset20.py
"""
from __future__ import annotations
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

_repo_root_default = "/mnt/katritch_lab2/aoxu/contrasCF"
REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", _repo_root_default))
sys.path.insert(0, str(REPO_ROOT / "analysis"))

from casf_mutagenesis.config import (  # noqa: E402
    AF3_DIR, AF3_ENV, AF3_MODEL_DIR, CUDA_DEVICE,
    OUTPUT_ROOT, SPLIT_JSON, SUBSET20_JSON, VARIANTS,
)
from casf_mutagenesis.inputs_af3 import (  # noqa: E402
    clean_a3m_for_af3, render_af3, rewrite_a3m_query,
)
from casf_mutagenesis.msa_via_boltz import fetch_msa_via_boltz  # noqa: E402

NUM_DIFFUSION_SAMPLES = 5
NUM_RECYCLES = 10
TIMEOUT_S = 3600
MAX_TOTAL_LENGTH = 800

# Which outputs tree to walk. Default = the binding-site (casf) module, so
# protein-arm behaviour is unchanged when the env var is unset. Point it at
# analysis/ligand_mutagenesis/outputs to drive the ligand arm — the same
# contract 08_/11_/14_ already honour.
OUTPUTS_ROOT = Path(os.environ.get("CONTRASCF_OUTPUTS_ROOT", str(OUTPUT_ROOT)))

# The MSA cache stays anchored to the casf tree even when OUTPUTS_ROOT is
# overridden, and that is deliberate: the ligand arm perturbs the LIGAND, so
# its receptor sequences are identical to the binding-site arm's and every
# a3m already fetched there is reusable. Override with CONTRASCF_MSA_CACHE.
MSA_CACHE = Path(os.environ.get("CONTRASCF_MSA_CACHE",
                                str(OUTPUT_ROOT / "_msa_cache")))


def _build_af3_env() -> dict:
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = CUDA_DEVICE
    env["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
    # AF3 v3 explicitly requires this XLA flag on Volta (compute 7.x):
    #   ValueError: For devices with GPU compute capability 7.x ... the ENV
    #   XLA_FLAGS must include "--xla_disable_hlo_passes=custom-kernel-fusion-rewriter".
    # Setting unconditionally is fine on Ampere+ (compute 8+) — at worst a
    # redundant optimisation skip. See run-log of CARC job 8622798 for the
    # original failure.
    xla_flag = "--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"
    existing_xla = env.get("XLA_FLAGS", "")
    env["XLA_FLAGS"] = (
        f"{existing_xla} {xla_flag}".strip() if xla_flag not in existing_xla
        else existing_xla
    )
    nvidia_lib_dirs: list[str] = []
    for sp in AF3_ENV.glob("lib/python*/site-packages/nvidia"):
        if sp.is_dir():
            for child in sp.iterdir():
                if (child / "lib").is_dir():
                    nvidia_lib_dirs.append(str(child / "lib"))
    if nvidia_lib_dirs:
        existing = env.get("LD_LIBRARY_PATH", "")
        env["LD_LIBRARY_PATH"] = ":".join(
            nvidia_lib_dirs + ([existing] if existing else [])
        )
    return env


def _read_af3_protein_chains(af3_json_path: Path) -> list[tuple[str, str]]:
    """Every protein chain in an af3.json, as [(chain_id, sequence), ...].

    Replaces an earlier `_read_af3_protein_sequence` that returned only the
    FIRST protein and silently discarded the rest. Combined with a hardcoded
    `chain_seqs=[("A", seq)]` at the render call, that dropped every extra
    chain from the AF3+MSA job: 51 of 52 multi-chain systems were folded as a
    single chain, and for `1bcu` (chains L=26 then H=257) the kept chain was
    the 26-residue L while both mutations sat on the discarded H.

    An entry's `id` may be a scalar or a list (AF3 lets one entry cover several
    identical copies); both are expanded to one tuple per chain.
    """
    d = json.loads(af3_json_path.read_text())
    chains: list[tuple[str, str]] = []
    for s in d["sequences"]:
        pr = s.get("protein")
        if not pr:
            continue
        ids = pr["id"]
        for cid in (ids if isinstance(ids, list) else [ids]):
            chains.append((str(cid), pr["sequence"]))
    if not chains:
        raise RuntimeError(f"no protein in {af3_json_path}")
    return chains


def _read_ligand_smiles(af3_json_path: Path) -> str:
    d = json.loads(af3_json_path.read_text())
    for s in d["sequences"]:
        if "ligand" in s:
            return s["ligand"]["smiles"]
    raise RuntimeError(f"no ligand in {af3_json_path}")


def _run_af3(json_path: Path, work_dir: Path) -> Path:
    work_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(AF3_ENV / "bin" / "python"),
        str(AF3_DIR / "run_alphafold.py"),
        f"--json_path={json_path}",
        f"--output_dir={work_dir}",
        f"--model_dir={AF3_MODEL_DIR}",
        "--norun_data_pipeline",
        "--run_inference=true",
        f"--num_diffusion_samples={NUM_DIFFUSION_SAMPLES}",
        f"--num_recycles={NUM_RECYCLES}",
        # Required on Volta (compute 7.x); harmless on Ampere+ (compute 8+).
        # AF3 v3 raises if this isn't set on V100. Pair with the XLA env flag
        # in _build_af3_env above.
        "--flash_attention_implementation=xla",
    ]
    result = subprocess.run(
        cmd, capture_output=True, text=True,
        env=_build_af3_env(), timeout=TIMEOUT_S,
    )
    (work_dir / "af3_stdout.log").write_text(result.stdout)
    (work_dir / "af3_stderr.log").write_text(result.stderr)
    if result.returncode != 0:
        raise RuntimeError(
            f"AF3 failed (exit {result.returncode}); "
            f"see {work_dir}/af3_stderr.log"
        )
    name = json.loads(json_path.read_text())["name"]
    sys_dir = work_dir / name
    if not sys_dir.exists():
        raise FileNotFoundError(f"AF3 output dir missing: {sys_dir}")
    return sys_dir


def _copy_all_samples(af3_sys_dir: Path, name: str, dst_dir: Path,
                      prefix_str: str = "af3msa") -> dict:
    """Copy every AF3 sample, renamed `<prefix_str>_<name>_model_<rank>.cif`,
    ranked 0..N-1 by ranking_score (best first)."""
    import csv
    ranking = af3_sys_dir / f"{name}_ranking_scores.csv"
    if not ranking.exists():
        raise FileNotFoundError(f"missing ranking CSV: {ranking}")
    with ranking.open() as f:
        rows = list(csv.DictReader(f))
    rows.sort(key=lambda r: float(r["ranking_score"]), reverse=True)
    out: dict = {"cifs": [], "confs": [], "ranking": None}
    for rank, row in enumerate(rows):
        seed, sample = int(row["seed"]), int(row["sample"])
        sub = af3_sys_dir / f"seed-{seed}_sample-{sample}"
        src_cif = sub / f"{name}_seed-{seed}_sample-{sample}_model.cif"
        src_conf = sub / f"{name}_seed-{seed}_sample-{sample}_summary_confidences.json"
        if not src_cif.exists():
            raise FileNotFoundError(f"missing model cif: {src_cif}")
        dst_cif = dst_dir / f"{prefix_str}_{name}_model_{rank}.cif"
        shutil.copy(src_cif, dst_cif)
        out["cifs"].append(str(dst_cif))
        if src_conf.exists():
            dst_conf = dst_dir / f"{prefix_str}_summary_confidences_{name}_{rank}.json"
            shutil.copy(src_conf, dst_conf)
            out["confs"].append(str(dst_conf))
    dst_rank = dst_dir / f"{prefix_str}_ranking_scores_{name}.csv"
    shutil.copy(ranking, dst_rank)
    out["ranking"] = str(dst_rank)
    # Legacy compat: rank-0 confidence JSON without the `_0` suffix, so the
    # old single-pose analysis path keeps working.
    if out["confs"]:
        legacy_conf = dst_dir / f"{prefix_str}_summary_confidences_{name}.json"
        shutil.copy(out["confs"][0], legacy_conf)
    return out


def _resolve_ids() -> tuple[list[str], str]:
    scope = os.environ.get("CONTRASCF_SCOPE", "subset20")
    if scope == "full":
        ids = json.loads(SPLIT_JSON.read_text())["casf2016"]
    elif scope == "disk":
        # Discover systems from OUTPUTS_ROOT itself. Required for the ligand
        # arm, whose system list comes from its build manifest rather than a
        # CASF label file, and whose per-system variant set is dynamic.
        ids = sorted(
            q.parent.parent.name
            for q in OUTPUTS_ROOT.glob("*/wt/af3.json")
            if not q.parent.parent.name.startswith("_")
        )
    else:
        ids = json.loads(SUBSET20_JSON.read_text())["casf2016"]
    start = int(os.environ.get("CONTRASCF_START", "0"))
    end = int(os.environ.get("CONTRASCF_END", str(len(ids))))
    return ids[start:end], f"{scope}[{start}:{end}]"


def _variants_for(pdbid: str) -> list[str]:
    """Variant names for one system.

    The binding-site arm has a fixed 4-tuple (wt/rem/pack/inv). The ligand
    arm's set is dynamic per system — halo_*/meth_*/chrg_* exist only where
    the ligand carries the matching functional group — so under
    CONTRASCF_SCOPE=disk we read it off disk instead. `wt` is forced first so
    the MSA-bearing reference is folded before the variants that reuse it.
    """
    if os.environ.get("CONTRASCF_SCOPE") != "disk":
        return list(VARIANTS)
    found = sorted(q.parent.name for q in (OUTPUTS_ROOT / pdbid).glob("*/af3.json"))
    if "wt" in found:
        found = ["wt"] + [v for v in found if v != "wt"]
    return found


def main() -> int:
    ids, scope_label = _resolve_ids()
    print(f"AF3+MSA scope: {scope_label}, n_pdb={len(ids)}")

    # Phase 1: ensure a WT MSA is cached for EVERY protein chain of each system.
    # One fetch per DISTINCT sequence — a homodimer reuses a single MSA, a
    # hetero-complex (e.g. 1bcu L/H) gets one per chain.
    print("Phase 1: fetch MSAs for WT sequences\n")
    wt_msas: dict[str, dict[str, str]] = {}          # pdbid -> {chain_id: a3m}
    for i, pdbid in enumerate(ids, 1):
        wt_json = OUTPUTS_ROOT / pdbid / "wt" / "af3.json"
        if not wt_json.exists():
            print(f"  [{i}/{len(ids)}] {pdbid}: no af3.json — skip"); continue
        wt_chains = _read_af3_protein_chains(wt_json)
        total_len = sum(len(s) for _, s in wt_chains)
        if total_len > MAX_TOTAL_LENGTH:
            print(f"  [{i}/{len(ids)}] {pdbid}: total len={total_len} "
                  f"> {MAX_TOTAL_LENGTH} — skip")
            continue
        t0 = time.time()
        per_chain: dict[str, str] = {}
        by_seq: dict[str, str] = {}
        try:
            for cid, seq in wt_chains:
                if seq not in by_seq:
                    by_seq[seq] = fetch_msa_via_boltz(seq, cache_dir=MSA_CACHE)
                per_chain[cid] = by_seq[seq]
        except Exception as exc:
            print(f"  [{i}/{len(ids)}] {pdbid}: MSA fetch failed: {exc}")
            continue
        wt_msas[pdbid] = per_chain
        n_seqs = min(sum(1 for ln in a.splitlines() if ln.startswith(">"))
                     for a in by_seq.values())
        print(f"  [{i}/{len(ids)}] {pdbid}: MSA ok (chains={len(wt_chains)}, "
              f"distinct_seqs={len(by_seq)}, n_seqs>={n_seqs}, "
              f"total_len={total_len}, {time.time()-t0:.1f}s)")

    # Phase 2: AF3 with MSA for each (pdbid, variant)
    print("\nPhase 2: run AF3 with MSA for each (pdbid, variant)\n")
    variants_by_id = {pdbid: _variants_for(pdbid) for pdbid in ids}
    n_total = sum(len(v) for v in variants_by_id.values())
    n_done = 0
    n_skip = 0
    n_fail = 0
    log_path = OUTPUTS_ROOT / "af3_msa_run_log.json"
    runs: list[dict] = []

    for pdbid in ids:
        wt_chain_msas = wt_msas.get(pdbid)
        variants = variants_by_id[pdbid]
        if wt_chain_msas is None:
            for v in variants:
                runs.append({"pdbid": pdbid, "variant": v, "status": "skip_no_msa"})
                n_skip += 1
                n_done += 1
            continue
        for variant in variants:
            n_done += 1
            v_dir = OUTPUTS_ROOT / pdbid / variant
            af3_json = v_dir / "af3.json"
            prefix = f"{pdbid}_{variant}"
            af3msa_cif = v_dir / f"af3msa_{prefix}_model_0.cif"
            entry = {"pdbid": pdbid, "variant": variant, "prefix": prefix}
            if not af3_json.exists():
                entry["status"] = "missing_input"
                n_skip += 1
                runs.append(entry); continue
            if af3msa_cif.exists():
                entry["status"] = "skip_existing"
                n_skip += 1
                runs.append(entry); continue
            variant_chains = _read_af3_protein_chains(af3_json)
            smiles = _read_ligand_smiles(af3_json)
            entry["n_chains"] = len(variant_chains)
            try:
                # One MSA per chain, each rewritten onto that chain's own
                # (possibly mutated) sequence. Chains are matched to the WT by
                # id, so a mutation on any chain reaches the model.
                chain_msas: dict[str, str] = {}
                for cid, seq in variant_chains:
                    base = wt_chain_msas.get(cid)
                    if base is None:
                        raise RuntimeError(
                            f"chain {cid} present in variant af3.json but not in WT")
                    a3m = rewrite_a3m_query(base, seq)
                    chain_msas[cid] = clean_a3m_for_af3(a3m, len(seq))
            except Exception as exc:
                entry.update({"status": "error", "error": f"a3m rewrite: {exc}"})
                n_fail += 1
                runs.append(entry); continue

            # Render fresh JSON with MSA — ALL protein chains, not just the first
            json_msa = v_dir / "af3_msa.json"
            render_af3(
                name=prefix, chain_seqs=variant_chains,
                ligand_smiles=smiles, out_path=json_msa,
                chain_msas=chain_msas,
            )
            print(f"  [{n_done}/{n_total}] {prefix} ...", flush=True)
            t0 = time.time()
            try:
                work = v_dir / "_af3_msa_work"
                if work.exists():
                    shutil.rmtree(work)
                af3_sys = _run_af3(json_msa, work)
                outs = _copy_all_samples(af3_sys, prefix, v_dir, "af3msa")
                entry.update({
                    "status": "ok", "outputs": outs,
                    "wallclock_s": round(time.time() - t0, 1),
                })
                shutil.rmtree(work, ignore_errors=True)
                print(f"      done in {entry['wallclock_s']}s")
            except Exception as exc:
                entry.update({
                    "status": "error", "error": str(exc),
                    "wallclock_s": round(time.time() - t0, 1),
                })
                n_fail += 1
                print(f"      FAILED: {exc}")
            runs.append(entry)
            log_path.write_text(json.dumps({"runs": runs}, indent=2))

    print(f"\nTotal {n_total} | done={n_done - n_skip - n_fail} "
          f"skip={n_skip} fail={n_fail}")
    print(f"Run log: {log_path}")
    # Always return 0 — see note in 03_run_boltz2_subset20.py.
    return 0


if __name__ == "__main__":
    sys.exit(main())
