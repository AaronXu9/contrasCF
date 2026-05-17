"""Compute the pbcatom (1-based atom index) for a named index-file group.

When a receptor or ligand group spans more than 0.25 × box-side, GROMACS-RAMD
requires an explicit `ramd-group1-{receptor,ligand}-pbcatom` — a single atom
near the group's geometric centre, used as the PBC reference for COM
calculations.

This helper picks the atom in the group that is closest to the group's
geometric centroid (no mass weighting; the pbcatom doesn't need mass-COM,
just a "central" reference). Returns a 1-based index suitable for direct
substitution into the production mdp.
"""
from __future__ import annotations
from pathlib import Path

import numpy as np


def parse_gro_coords(gro_path: Path) -> np.ndarray:
    """Return (N, 3) coords array. Atom order is file order = 1-based gro index."""
    lines = Path(gro_path).read_text().splitlines()
    natoms = int(lines[1].strip())
    coords = np.empty((natoms, 3), dtype=float)
    for i, line in enumerate(lines[2:2 + natoms]):
        coords[i] = (float(line[20:28]), float(line[28:36]), float(line[36:44]))
    return coords


def parse_ndx_group(ndx_path: Path, group_name: str) -> list[int]:
    """Return 1-based atom indices for the FIRST occurrence of the named group."""
    text = Path(ndx_path).read_text()
    in_group = False
    indices: list[int] = []
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("[") and s.endswith("]"):
            current = s.strip("[ ]").strip()
            if in_group:
                break  # next group hit — first match captured, stop
            if current == group_name:
                in_group = True
            continue
        if in_group and s:
            indices.extend(int(x) for x in s.split())
    if not indices:
        raise KeyError(f"group {group_name!r} not found or empty in {ndx_path}")
    return indices


def compute_pbcatom(gro_path: Path, ndx_path: Path, group_name: str) -> int:
    """1-based atom index closest to the geometric centre of `group_name`."""
    coords = parse_gro_coords(gro_path)
    group_idx = parse_ndx_group(ndx_path, group_name)
    # 1-based → 0-based array indices
    rows = np.array(group_idx, dtype=int) - 1
    group_coords = coords[rows]
    centre = group_coords.mean(axis=0)
    dists = np.linalg.norm(group_coords - centre, axis=1)
    return int(group_idx[int(np.argmin(dists))])


def compute_ramd_pbcatoms(
    system_dir: Path, receptor: str = "Protein", ligand: str = "ATP",
) -> dict[str, int]:
    """Convenience: returns {receptor: idx, ligand: idx} for a built system_dir."""
    gro = Path(system_dir) / "solv_ions.gro"
    ndx = Path(system_dir) / "index.ndx"
    return {
        receptor: compute_pbcatom(gro, ndx, receptor),
        ligand:   compute_pbcatom(gro, ndx, ligand),
    }


if __name__ == "__main__":
    import argparse, json, sys
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--system-dir", type=Path, required=True)
    p.add_argument("--receptor", default="Protein")
    p.add_argument("--ligand", default="ATP")
    args = p.parse_args()
    out = compute_ramd_pbcatoms(args.system_dir, args.receptor, args.ligand)
    print(json.dumps(out, indent=2))
    sys.exit(0)
