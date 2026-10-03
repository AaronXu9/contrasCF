"""Compare SurfDock input surfaces built from a PROTONATED vs an H-STRIPPED
receptor, on the identical pocket.

Motivation: the CASF mutagenesis WT arm docks into a protonated crystal receptor
(HiQBind ships ~4k hydrogens) while every AF3-predicted mutant receptor has none.
SurfDock consumes an MSMS molecular surface plus MaSIF-style chemical channels,
both of which read atomic radii and H-bond donors -- so unlike the Vina family
(bit-identical scores, verified) it may be sensitive to protonation. If it is,
part of SurfDock's WT->mutant gap is a preparation artifact, not a mutation
response.
"""
import sys, glob, os
import numpy as np
from plyfile import PlyData
from scipy.spatial import cKDTree

T = "/home/aoxu/tmp/h_surface_test"

def load(path):
    p = PlyData.read(path)
    v = p["vertex"]
    xyz = np.stack([v["x"], v["y"], v["z"]], axis=1)
    ch = {n: np.asarray(v[n]) for n in v.data.dtype.names if n not in ("x", "y", "z")}
    nf = len(p["face"].data) if "face" in p else 0
    return xyz, ch, nf

def find(mode, sysid):
    hits = glob.glob(f"{T}/out_{mode}/{sysid}/*.ply")
    return hits[0] if hits else None

print(f"{'system':<7s} {'mode':<6s} {'verts':>7s} {'faces':>7s}  "
      f"{'channel stats (mean +- sd)'}")
rows = {}
for sysid in ("1e66", "3u5j", "2brb"):
    for mode in ("withH", "noH"):
        f = find(mode, sysid)
        if not f:
            print(f"{sysid:<7s} {mode:<6s} {'MISSING':>7s}")
            continue
        xyz, ch, nf = load(f)
        rows[(sysid, mode)] = (xyz, ch)
        desc = "  ".join(f"{k}={ch[k].mean():+.3f}+-{ch[k].std():.3f}"
                         for k in sorted(ch) if ch[k].dtype.kind == "f")
        print(f"{sysid:<7s} {mode:<6s} {len(xyz):>7d} {nf:>7d}  {desc}")

print("\n--- protonated vs stripped, same pocket ---")
print(f"{'system':<7s} {'d_verts':>9s} {'%':>7s} {'median NN dist':>15s} {'p95 NN':>8s} {'hbond change':>14s}")
for sysid in ("1e66", "3u5j", "2brb"):
    a = rows.get((sysid, "withH")); b = rows.get((sysid, "noH"))
    if not a or not b:
        continue
    (xa, ca), (xb, cb) = a, b
    d, _ = cKDTree(xa).query(xb)           # each stripped vertex -> nearest protonated
    dv = len(xb) - len(xa)
    hb = ""
    if "hbond" in ca and "hbond" in cb:
        hb = f"{ca['hbond'].mean():+.3f} -> {cb['hbond'].mean():+.3f}"
    print(f"{sysid:<7s} {dv:>+9d} {100*dv/len(xa):>+6.1f}% "
          f"{np.median(d):>14.3f}A {np.percentile(d,95):>7.3f}A {hb:>14s}")
