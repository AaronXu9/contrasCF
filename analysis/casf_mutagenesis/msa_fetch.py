"""Minimal ColabFold MSA-server client.

Fetches a uniref-merged A3M MSA for a single protein sequence using the
public MMseqs2 API at https://api.colabfold.com/. Polls until completion,
unpacks the tar.gz response, and returns the merged A3M text.

This is the lightweight alternative to installing AlphaFold3's full data
pipeline (which needs ~600 GB of BFD/UniRef30/MGnify databases). The
ColabFold API is the same one Boltz's `--use_msa_server` flag uses
internally.

Example:
    from casf_mutagenesis.msa_fetch import fetch_msa_a3m
    a3m = fetch_msa_a3m(seq)
    Path("uniref.a3m").write_text(a3m)
"""
from __future__ import annotations
import io
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

DEFAULT_API = "https://api.colabfold.com"
POLL_INTERVAL_S = 6
MAX_WAIT_S = 1800   # 30 min per query


def _post(url: str, data: dict[str, str]) -> dict:
    """POST form data, parse JSON response."""
    body = urllib.parse.urlencode(data).encode()
    req = urllib.request.Request(url, data=body, method="POST")
    req.add_header(
        "Content-Type", "application/x-www-form-urlencoded; charset=utf-8"
    )
    with urllib.request.urlopen(req, timeout=60) as r:
        import json as _json
        return _json.loads(r.read())


def _get_json(url: str) -> dict:
    with urllib.request.urlopen(url, timeout=60) as r:
        import json as _json
        return _json.loads(r.read())


def _get_bytes(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=300) as r:
        return r.read()


def fetch_msa_a3m(
    sequence: str,
    *,
    mode: str = "env",
    api: str = DEFAULT_API,
    cache_dir: Path | str | None = None,
    poll_interval_s: int = POLL_INTERVAL_S,
    max_wait_s: int = MAX_WAIT_S,
) -> str:
    """Submit a single sequence to ColabFold's MSA server and return A3M text.

    Args:
        sequence: query protein sequence (1-letter, no whitespace).
        mode: "env" (uniref + environmental DBs, default) or "all" or
              "uniref". "env" is what Boltz / AF2 use by default.
        cache_dir: if given, results are cached as f"{sha1(seq)}.a3m" so
            repeat calls with the same sequence skip the network round-trip.

    Returns:
        Merged A3M text (the ColabFold-merged uniref + env MSA), with the
        query sequence as the first record.
    """
    sequence = sequence.strip().upper()
    if not sequence.isalpha():
        raise ValueError(f"sequence must be alphabetic; got {sequence!r}")

    if cache_dir is not None:
        import hashlib
        cache_dir = Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
        key = hashlib.sha1(f"{mode}:{sequence}".encode()).hexdigest()[:16]
        cache_path = cache_dir / f"{key}.a3m"
        if cache_path.exists():
            return cache_path.read_text()

    submit = _post(f"{api}/ticket/msa", {"q": sequence, "mode": mode})
    ticket_id = submit["id"]
    status = submit.get("status", "PENDING")

    deadline = time.monotonic() + max_wait_s
    while status in ("PENDING", "RUNNING"):
        if time.monotonic() > deadline:
            raise TimeoutError(f"MSA ticket {ticket_id} did not finish in {max_wait_s}s")
        time.sleep(poll_interval_s)
        status_resp = _get_json(f"{api}/ticket/{ticket_id}")
        status = status_resp.get("status", "ERROR")
    if status != "COMPLETE":
        raise RuntimeError(f"MSA ticket {ticket_id} ended with status {status}")

    raw = _get_bytes(f"{api}/result/download/{ticket_id}")
    a3m = _extract_merged_a3m(raw)

    if cache_dir is not None:
        cache_path.write_text(a3m)
    return a3m


def _extract_merged_a3m(raw_tar_bytes: bytes) -> str:
    """ColabFold returns a tar.gz with several A3M files (uniref30, bfd
    mgnify, etc.). We concatenate them — common ColabFold pattern — and
    keep the query as the first record once.
    """
    with tarfile.open(fileobj=io.BytesIO(raw_tar_bytes), mode="r:*") as tar:
        a3m_members = [
            m for m in tar.getmembers()
            if m.isfile() and m.name.endswith(".a3m")
        ]
        if not a3m_members:
            raise RuntimeError("no .a3m in MSA result tar")
        # Sort by typical priority (uniref first, then env DBs)
        order = ["uniref", "bfd.mgnify30.metaeuk30.smag30", "pdb70"]
        def rank(m):
            for i, p in enumerate(order):
                if p in m.name:
                    return i
            return 99
        a3m_members.sort(key=rank)

        merged_lines: list[str] = []
        seen_query = False
        for m in a3m_members:
            data = tar.extractfile(m).read().decode()
            for hdr_seq in _iter_a3m(data):
                hdr, seq = hdr_seq
                if hdr.startswith(">101") or hdr.startswith(">query") or (
                    not seen_query and hdr.startswith(">")
                ):
                    if seen_query:
                        # skip subsequent query rows (one per .a3m file)
                        continue
                    seen_query = True
                    merged_lines.append(">query")
                    merged_lines.append(seq)
                    continue
                merged_lines.append(hdr)
                merged_lines.append(seq)
        return "\n".join(merged_lines) + "\n"


def _iter_a3m(text: str):
    """Yield (header, sequence) tuples from an A3M string. Wrapped sequences
    are concatenated; non-AA characters are preserved (a3m keeps lowercase
    insertions)."""
    hdr = None
    seq_parts: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if hdr is not None:
                yield hdr, "".join(seq_parts)
            hdr = line
            seq_parts = []
        else:
            if line:
                seq_parts.append(line)
    if hdr is not None:
        yield hdr, "".join(seq_parts)
