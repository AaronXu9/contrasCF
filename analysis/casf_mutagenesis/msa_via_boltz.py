"""Fetch MSAs by piggybacking on Boltz's `--use_msa_server` flag.

Boltz internally calls the ColabFold MMseqs2 API with the right pairing
strategy and writes the resulting A3M to:
    <out_dir>/boltz_results_<name>/msa/<name>_unpaired_tmp_env/uniref.a3m
                                                              /bfd.mgnify30.metaeuk30.smag30.a3m

This wrapper:
  1. Writes a minimal protein-only YAML.
  2. Runs `boltz predict --use_msa_server` with absolute-minimum settings
     (1 sample, 1 recycle, 5 sampling steps) so the inference time is
     trivial and the dominant cost is the MSA download (~2 sec/sequence).
  3. Reads the resulting A3M, strips null bytes (Boltz/MMseqs2 occasionally
     emits a `\\x00` in one row — AF3 then rejects it as "Unknown residues").
  4. Returns the cleaned A3M text.

Why not call the ColabFold HTTP API directly? Direct submission via
`POST /ticket/msa` was returning persistent `status: ERROR` for every
sequence. Boltz uses different endpoint parameters / pairing strategy
internally and works. This module reuses Boltz's known-working code path.
"""
from __future__ import annotations
import hashlib
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from .config import BOLTZ_BIN as _DEFAULT_BOLTZ_BIN, CUDA_DEVICE, assert_boltz2_binary

BOLTZ_BIN = str(_DEFAULT_BOLTZ_BIN)
_BOLTZ_VERSION_CHECKED = False  # lazily checked on first fetch_msa_via_boltz call


def _yaml_for_msa(seq: str) -> str:
    """Minimal Boltz YAML — protein-only, no ligand, no MSA spec.

    The A and ligand-id chains used by `inputs_boltz.render_boltz` aren't
    needed here; we only want the MSA side-effect.
    """
    return (
        "version: 1\n"
        "sequences:\n"
        "  - protein:\n"
        "      id: [A]\n"
        f"      sequence: {seq}\n"
    )


def _strip_null_bytes(a3m: str) -> str:
    return a3m.replace("\x00", "")


def fetch_msa_via_boltz(
    sequence: str,
    *,
    cache_dir: Path | str,
    boltz_bin: str = BOLTZ_BIN,
    diffusion_samples: int = 1,
    recycling_steps: int = 1,
    sampling_steps: int = 5,
) -> str:
    """Fetch and cache the uniref+bfd merged A3M for a single sequence.

    The returned A3M concatenates Boltz's `uniref.a3m` and `bfd.mgnify30
    .metaeuk30.smag30.a3m` (uniref first, query row deduplicated). Cached
    at `<cache_dir>/<sha1>.a3m`.
    """
    sequence = sequence.strip().upper()
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha1(sequence.encode()).hexdigest()[:16]
    cache_path = cache_dir / f"{key}.a3m"
    if cache_path.exists():
        return cache_path.read_text()

    # Defensive: Boltz-1 (v0.4.1) and Boltz-2 (v2.2.1) coexist as `boltz`
    # binaries on lab workstations. We rely on the v2 schema implicitly
    # (the YAML format and the --use_msa_server contract). Fail early if
    # CONTRASCF_BOLTZ_BIN points at a v1 install.
    global _BOLTZ_VERSION_CHECKED
    if not _BOLTZ_VERSION_CHECKED:
        assert_boltz2_binary(boltz_bin)
        _BOLTZ_VERSION_CHECKED = True

    with tempfile.TemporaryDirectory(prefix="boltz_msa_") as tmp:
        tmp_dir = Path(tmp)
        yaml_path = tmp_dir / f"{key}.yaml"
        yaml_path.write_text(_yaml_for_msa(sequence))
        out_dir = tmp_dir / "out"
        out_dir.mkdir()

        cmd = [
            boltz_bin, "predict", str(yaml_path),
            "--out_dir", str(out_dir),
            "--model", "boltz2",
            "--output_format", "mmcif",
            "--diffusion_samples", str(diffusion_samples),
            "--recycling_steps", str(recycling_steps),
            "--sampling_steps", str(sampling_steps),
            "--seed", "42",
            "--use_msa_server",
        ]
        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = CUDA_DEVICE
        result = subprocess.run(cmd, capture_output=True, text=True, env=env)
        if result.returncode != 0:
            raise RuntimeError(
                f"boltz MSA fetch failed (exit {result.returncode}):\n"
                f"{result.stderr[-2000:]}"
            )

        msa_dir = out_dir / f"boltz_results_{key}" / "msa" / f"{key}_unpaired_tmp_env"
        uniref = msa_dir / "uniref.a3m"
        bfd = msa_dir / "bfd.mgnify30.metaeuk30.smag30.a3m"
        if not uniref.exists():
            raise FileNotFoundError(f"Boltz did not produce {uniref}")

        a3m = _strip_null_bytes(uniref.read_text())
        if bfd.exists():
            # Append BFD/MGnify hits, skipping the duplicate query row.
            bfd_text = _strip_null_bytes(bfd.read_text())
            seen_query = False
            extra_lines = []
            for ln in bfd_text.splitlines():
                if ln.startswith(">"):
                    if not seen_query:
                        # first header is the query — skip it
                        seen_query = True
                        continue
                    extra_lines.append(ln)
                elif seen_query:
                    extra_lines.append(ln)
            if extra_lines:
                a3m = a3m.rstrip("\n") + "\n" + "\n".join(extra_lines) + "\n"

    cache_path.write_text(a3m)
    return a3m
