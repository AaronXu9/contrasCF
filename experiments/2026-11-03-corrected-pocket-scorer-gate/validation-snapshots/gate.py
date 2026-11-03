"""Gate: the corrected six-method scorer reproduces the 2026-10-21 audit on Set B (96 systems, pack/inv).
References: GNINA from the corrected-replay table (merged_gnina192.csv), UniDock2/ICM from
pose_diagnostics.csv (exact_local5A_rmsd_A), SurfDock from surfdock_pose_diagnostics.csv.
Pass: per-row retained/not agreement >= 95% for each method, and every disagreement is one
where our symmetry-aware RMSD is LOWER than the audit's fixed historical mapping."""
import csv, json, sys
from pathlib import Path
R = Path("/mnt/katritch_lab2/aoxu/contrasCF"); A = R / "experiments/2026-10-21-pipeline-audit/pocket"
setb = set(json.load(open(R / "experiments/2026-10-20-setb-local-rmsd/cohort.json"))["set_b"])
mine = {(r["system"], r["variant"], r["method"]): r for r in csv.DictReader(open(R / "analysis/casf_mutagenesis/outputs/corrected_pocket_scores.csv"))}
ref = {}
for r in csv.DictReader(open(A / "gnina_corrected_analysis/merged_gnina192.csv")):
    if r.get("rmsd_A"): ref[(r["system"], r["variant"], "GNINA")] = float(r["rmsd_A"])
for r in csv.DictReader(open(A / "pose_diagnostics.csv")):
    if r["method"] in ("UniDock2", "ICM") and r.get("exact_local5A_rmsd_A"):
        ref[(r["system"], r["variant"], r["method"])] = float(r["exact_local5A_rmsd_A"])
sd = A / "surfdock_pose_diagnostics.csv"
if sd.exists():
    rows = list(csv.DictReader(open(sd))); col = "exact_local5A_rmsd_A"
    for r in rows:
        if r.get(col): ref[(r["system"], r["variant"], "SurfDock")] = float(r[col])
ok = True; report = {}
for meth in ("GNINA", "UniDock2", "ICM", "SurfDock"):
    keys = [k for k in ref if k[2] == meth and k[0] in setb]
    agree = 0; bad = []; unexplained = []; counts = {"pack": [0, 0], "inv": [0, 0]}
    for k in keys:
        m = mine.get(k)
        if not m or m.get("status") != "ok": bad.append((k, (m or {}).get("error", "missing"))); continue
        a, b = ref[k] < 2, float(m["ligand_rmsd_a"]) < 2
        counts[k[1]][0] += a; counts[k[1]][1] += b
        if a == b: agree += 1
        elif not (b and float(m["ligand_rmsd_a"]) < ref[k]): unexplained.append((k, ref[k], m["ligand_rmsd_a"]))
    n = len(keys) - len(bad); frac = agree / n if n else 0
    passed = frac >= 0.95 and not unexplained
    ok &= passed
    report[meth] = dict(rows=len(keys), scored=n, agreement=round(frac, 4), not_scored=bad[:5], unexplained=unexplained,
                        retained_audit_vs_ours={v: f"{c[0]} vs {c[1]}" for v, c in counts.items()}, passed=passed)
    print(meth, json.dumps(report[meth]))
print("GATE", "PASS" if ok else "FAIL")
json.dump(report, open(Path(__file__).with_name("gate_result.json"), "w"), indent=2)
