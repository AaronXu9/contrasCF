"""Repair cached A3M files written before the 2026-11-02 msa_via_boltz merge fix.

The old merge skipped the BFD file's query HEADER but not its sequence line, so
every cached A3M has the query sequence glued onto the record that ended
uniref.a3m:
  * usually the last UniRef hit -> that hit gets dropped by clean_a3m_for_af3
    (wrong length), costing one row of thousands;
  * when UniRef returned no hits, the query record itself -> the query row is
    doubled and AF3 refuses the job (1o5b, 2c3i: 8-residue peptide chains).

Repair: un-double a doubled query row; drop a trailing sequence line equal to the
query from any multi-line hit record. Originals are copied to
<cache>/../_msa_cache_pre_repair/ first. Idempotent: a repaired file is unchanged
on a second run.

Run:  $CONTRASCF_PY analysis/casf_mutagenesis/scripts/build/10_repair_msa_cache.py [--cache DIR] [--dry-run]
"""
from __future__ import annotations
import argparse, os, shutil, sys
from pathlib import Path

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))


def _records(text: str) -> list[tuple[str, list[str]]]:
    out, hdr, seq = [], None, []
    for ln in text.splitlines():
        if ln.startswith(">"):
            if hdr is not None:
                out.append((hdr, seq))
            hdr, seq = ln, []
        elif ln:
            seq.append(ln)
    if hdr is not None:
        out.append((hdr, seq))
    return out


def repair(text: str) -> tuple[str, dict]:
    recs = _records(text)
    if not recs:
        return text, {"doubled": 0, "unglued": 0}
    q = "".join(recs[0][1]); doubled = 0
    half = len(q) // 2
    if len(q) % 2 == 0 and half and q[:half] == q[half:]:
        q = q[:half]; recs[0] = (recs[0][0], [q]); doubled = 1
    unglued = 0
    for i in range(1, len(recs)):
        h, s = recs[i]
        if len(s) > 1 and s[-1] == q:
            recs[i] = (h, s[:-1]); unglued += 1
    lines = []
    for h, s in recs:
        lines.append(h); lines.extend(s)
    return "\n".join(lines) + "\n", {"doubled": doubled, "unglued": unglued}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", default=str(REPO_ROOT / "analysis/casf_mutagenesis/outputs/_msa_cache"))
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    cache = Path(a.cache); backup = cache.parent / "_msa_cache_pre_repair"
    tot = {"files": 0, "changed": 0, "doubled": 0, "unglued": 0}
    for f in sorted(cache.glob("*.a3m")):
        old = f.read_text(); new, st = repair(old)
        tot["files"] += 1; tot["doubled"] += st["doubled"]; tot["unglued"] += st["unglued"]
        if new != old:
            tot["changed"] += 1
            if not a.dry_run:
                backup.mkdir(exist_ok=True)
                if not (backup / f.name).exists():
                    shutil.copy2(f, backup / f.name)
                tmp = f.with_suffix(".a3m.tmp"); tmp.write_text(new); tmp.replace(f)
    print(("DRY RUN " if a.dry_run else "") + str(tot))
    return 0


if __name__ == "__main__":
    sys.exit(main())
