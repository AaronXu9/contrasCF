#!/usr/bin/env python3
"""Fill AF3+MSA confidence (iptm/ptm/ranking_score) into results_full.csv.

The full-CASF AF3+MSA run on CARC wrote per-cell `af3msa_summary_confidences_*`
JSONs, but only the model cifs were rsync'd back, so `05_analyze` left the
AF3+MSA confidence columns blank for all but the 19 subset20 systems. This
patches results_full.csv in place (with a .bak) from the rsync'd JSONs so the
confidence analyses (16/17) see AF3+MSA at full CASF without a full rdkit_env
re-run of 05_analyze. The canonical fix is to re-run 05_analyze on the lab box;
this is the lightweight confidence-columns-only equivalent. Stdlib only.
"""
import csv, json, glob, shutil
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "outputs"
RES = OUT / "results_full.csv"

conf = {}
for f in glob.glob(str(OUT / "*" / "*" / "af3msa_summary_confidences_*.json")):
    p = Path(f).parts
    pdbid, variant = p[-3], p[-2]
    try:
        conf[(pdbid, variant)] = json.load(open(f))
    except Exception:
        pass
print(f"loaded af3msa confidence for {len(conf)} cells "
      f"({len({k[0] for k in conf})} systems)")

with open(RES) as fh:
    reader = csv.DictReader(fh)
    fields = reader.fieldnames
    rows = list(reader)

shutil.copy(RES, str(RES) + ".bak")
filled = 0
for r in rows:
    if r["model"] != "AF3+MSA":
        continue
    d = conf.get((r["pdbid"], r["variant"]))
    if d is None or r.get("iptm", "") != "":   # don't clobber the 19 already there
        continue
    r["iptm"] = d.get("iptm", "")
    r["ptm"] = d.get("ptm", "")
    r["ranking_score"] = d.get("ranking_score", "")
    filled += 1

with open(RES, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields)
    w.writeheader()
    w.writerows(rows)
print(f"filled {filled} AF3+MSA rows  (backup: {RES.name}.bak)")
