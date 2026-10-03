#!/usr/bin/env python3
"""Batch structural views of the pocket-mutagenesis contrast (SurfDock + Boltz-2).

Everything is exported in the CRYSTAL frame, using the SAME transforms the RMSD
analysis applies, so what you see is the pose whose number we report:

  SurfDock  rem/pack/inv docked into an AF3-predicted mutant receptor, so the
            on-disk SDF is in the AF3 frame. We apply `gnina_analysis._ca_transform`
            (receptor -> crystal) to BOTH the pose and the mutant receptor.
  SurfDock  wt docked into the crystal receptor -> identity, no transform.
  Boltz-2   folds its own protein, so the whole predicted complex is superposed
            onto the crystal by ligand-near-chain Ca (`superpose_by_index`),
            exactly as `analysis._analyze_single_pose` does.

Skipping either transform puts mutant poses ~80-90 A from the crystal ligand.

Writes per system, under --out/<system>/:
  crystal_protein.pdb              crystal receptor (the common frame)
  crystal_ligand.sdf               the RMSD reference for every cell
  surfdock_<variant>_receptor.pdb  mutant receptor, aligned  (wt: the crystal one)
  surfdock_<variant>_pose.sdf      rank-1 pose, aligned
  boltz2_<variant>_complex.pdb     predicted complex, aligned
  boltz2_<variant>_ligand.sdf      predicted ligand, aligned
  view.pml                         PyMOL session: one scene per (method, variant)
  manifest.json                    per-cell RMSD + provenance

Run (rdkit_env):
  .../rdkit_env/bin/python plot/07_render_mutation_views.py --pick memorized --limit 12 --render
  .../rdkit_env/bin/python plot/07_render_mutation_views.py --systems 1e66,1gpk --render --contact-sheet
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
sys.path.insert(0, str(REPO_ROOT / "analysis"))

import gemmi  # noqa: E402
from rdkit import Chem  # noqa: E402

from casf_mutagenesis.config import (  # noqa: E402
    CASF_LIGANDS, CASF_RAW, OUTPUT_ROOT, VARIANTS,
)
from casf_mutagenesis.analysis import (  # noqa: E402
    MODEL_FILES, _crystal_ligand_mol, _heavy_coords, _heavy_indices,
    _predicted_ligand, extract_protein_ca_near, read_structure,
    superpose_by_index,
)
from casf_mutagenesis.gnina_analysis import _ca_transform, _read_top_pose  # noqa: E402

_SUP_CACHE: dict = {}

PYMOL_BIN = Path(os.environ.get(
    "CONTRASCF_PYMOL_BIN", "/home/aoxu/miniconda3/envs/PyMOL-PoseBench/bin/pymol"))

# Two docking engines then two co-folding models, so the grid reads
# physics-on-a-given-pocket (top) vs builds-its-own-pocket (bottom).
METHODS = ("icm", "surfdock", "boltz2", "af3msa")
PRETTY = {"icm": "ICM", "surfdock": "SurfDock",
          "boltz2": "Boltz-2", "af3msa": "AF3+MSA"}
COFOLD = {"boltz2": "{s}_{v}_model_*.cif", "af3msa": "af3msa_{s}_{v}_model_*.cif"}
ICM_ROOT = OUTPUT_ROOT / "_icm_poses"

# Okabe-Ito. Crystal reference is yellow so it never collides with a variant.
COLORS = {
    "crystal": "#F0E442",
    "wt":      "#0072B2",
    "rem":     "#E69F00",
    "pack":    "#009E73",
    "inv":     "#CC79A7",
}


# --------------------------------------------------------------------------
# geometry helpers
# --------------------------------------------------------------------------

def _apply_to_structure(st: gemmi.Structure, sup) -> None:
    """Transform every atom of `st` in place with a Superposition."""
    for model in st:
        for chain in model:
            for res in chain:
                for atom in res:
                    p = sup.apply(np.array([[atom.pos.x, atom.pos.y, atom.pos.z]]))[0]
                    atom.pos = gemmi.Position(float(p[0]), float(p[1]), float(p[2]))


def _apply_to_mol(mol: Chem.Mol, sup) -> Chem.Mol:
    """Return a copy of `mol` with every conformer atom transformed."""
    out = Chem.Mol(mol)
    conf = out.GetConformer()
    xyz = np.array([list(conf.GetAtomPosition(i)) for i in range(out.GetNumAtoms())])
    xyz = sup.apply(xyz)
    for i, p in enumerate(xyz):
        conf.SetAtomPosition(i, (float(p[0]), float(p[1]), float(p[2])))
    return out


def _write_sdf(mol: Chem.Mol, path: Path, name: str) -> None:
    mol = Chem.Mol(mol)
    try:
        mol = Chem.RemoveHs(mol, sanitize=False)   # heavy atoms only, see _write_pdb
    except Exception:
        pass
    mol.SetProp("_Name", name)
    w = Chem.SDWriter(str(path))
    # ICM poses are parsed leniently (valences RDKit rejects), so they carry
    # aromatic flags without a valid kekule structure and SDWriter's default
    # kekulization raises. Aromatic bond records are fine for viewing.
    w.SetKekulize(False)
    try:
        w.write(mol)
    finally:
        w.close()


def _write_pdb(st: gemmi.Structure, path: Path) -> None:
    """Write heavy atoms only.

    Protonation is NOT uniform across the arms: the wt docking receptor is the
    PROTONATED crystal structure (1059 H of 2121 atoms on 3u5j) while every
    AF3-predicted mutant receptor and every co-folding complex has none. Left
    as-is, the wt panel of each docking row renders roughly twice the sticks of
    its own rem/pack/inv panels, so the eye reads a difference in pocket density
    that is really a difference in file provenance. RMSD is heavy-atom
    throughout, so dropping H changes nothing quantitative.
    """
    st.setup_entities()
    st.remove_hydrogens()
    st.write_pdb(str(path))


# --------------------------------------------------------------------------
# per-method exporters
# --------------------------------------------------------------------------

def docking_receptor(system: str, variant: str, outdir: Path):
    """Write (once) the docking receptor for this cell, in the CRYSTAL frame.

    Both docking engines were handed the SAME receptor for a given cell, so the
    aligned copy is shared rather than duplicated per engine.
    Returns (path, Superposition|None) or (None, None).
    """
    out = outdir / f"mutant_{variant}_receptor.pdb"
    receptor = OUTPUT_ROOT / system / variant / "docking" / "receptor.pdb"
    if not receptor.exists():
        return None, None
    if out.exists():
        return out, _SUP_CACHE.get((system, variant), None)

    st_rec = read_structure(receptor)
    sup = None
    if variant != "wt":
        # Same condition and same transform as gnina_analysis.analyze_gnina.
        crystal_pdb = CASF_RAW / system / f"{system}_protein.pdb"
        crystal = _read_top_pose(CASF_LIGANDS / f"{system}_ligand.sdf")
        if not crystal_pdb.exists() or crystal is None:
            return None, None
        try:
            sup = _ca_transform(receptor, crystal_pdb, _heavy_coords(crystal))
        except Exception as exc:
            print(f"    {variant}: receptor align failed "
                  f"({type(exc).__name__}: {exc})")
            return None, None
        _apply_to_structure(st_rec, sup)
    _SUP_CACHE[(system, variant)] = sup
    _write_pdb(st_rec, out)
    return out, sup


def export_surfdock(system: str, variant: str, outdir: Path) -> dict | None:
    """SurfDock rank-1 pose, moved into the crystal frame."""
    pose_sdf = OUTPUT_ROOT / system / variant / "surfdock" / "poses.sdf"
    if not pose_sdf.exists():
        return None
    rec_path, sup = docking_receptor(system, variant, outdir)
    if rec_path is None:
        return None
    pose = _read_top_pose(pose_sdf)
    if pose is None:
        return None
    if sup is not None:
        pose = _apply_to_mol(pose, sup)
    _write_sdf(pose, outdir / f"surfdock_{variant}_pose.sdf",
               f"{system}_{variant}_surfdock")
    return {
        "aligned": sup is not None,
        "ca_rmsd_a": round(sup.rmsd, 3) if sup is not None else 0.0,
        "n_ca_paired": sup.n_paired if sup is not None else None,
        "source": str(pose_sdf.relative_to(REPO_ROOT)),
    }


def _read_icm_pose(path: Path):
    """ICM writes valences RDKit rejects (phosphates); fall back to a lenient
    parse exactly as analyze/04_analyze_icm.py does. Only coordinates are used."""
    m = next(iter(Chem.SDMolSupplier(str(path), sanitize=True)), None)
    if m is not None:
        return m
    m = next(iter(Chem.SDMolSupplier(str(path), sanitize=False)), None)
    if m is None:
        return None
    try:
        m = Chem.Mol(m)
        Chem.SanitizeMol(m, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
                         ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        return m
    except Exception:
        return None


def export_icm(system: str, variant: str, outdir: Path) -> dict | None:
    """ICM rank-1 pose. ICM docked into a PRE-ALIGNED receptor, so its poses are
    already in the crystal frame and must NOT be transformed again."""
    pose_sdf = ICM_ROOT / system / variant / "ICM" / f"D_{system}_{variant}_pose1.sdf"
    if not pose_sdf.exists():
        return None
    rec_path, sup = docking_receptor(system, variant, outdir)
    if rec_path is None:
        return None
    pose = _read_icm_pose(pose_sdf)
    if pose is None or not pose.GetNumConformers():
        return None
    _write_sdf(pose, outdir / f"icm_{variant}_pose.sdf", f"{system}_{variant}_icm")
    return {
        "aligned": False,          # already in the crystal frame by construction
        "ca_rmsd_a": round(sup.rmsd, 3) if sup is not None else 0.0,
        "n_ca_paired": sup.n_paired if sup is not None else None,
        "source": str(pose_sdf.relative_to(REPO_ROOT)),
    }


def export_cofold(model: str, system: str, variant: str, outdir: Path) -> dict | None:
    """Export a co-folding model's rank-0 complex + ligand in the crystal frame."""
    cell = OUTPUT_ROOT / system / variant
    cifs = sorted(cell.glob(COFOLD[model].format(s=system, v=variant)))
    if not cifs:
        return None
    cif = cifs[0]                     # rank 0 = top-ranked by model confidence

    crystal_pdb = CASF_RAW / system / f"{system}_protein.pdb"
    crystal_sdf = CASF_LIGANDS / f"{system}_ligand.sdf"
    if not crystal_pdb.exists() or not crystal_sdf.exists():
        return None
    st_native = read_structure(crystal_pdb)
    crystal_mol = _crystal_ligand_mol(crystal_sdf)
    crystal_lig_xyz = _heavy_coords(crystal_mol)

    st_pred = read_structure(cif)
    lig_block, pred_mol = _predicted_ligand(st_pred, Chem.MolToSmiles(crystal_mol))
    if lig_block is None:
        return None

    try:
        sup = superpose_by_index(
            extract_protein_ca_near(st_pred, lig_block.heavy_coords),
            extract_protein_ca_near(st_native, crystal_lig_xyz),
        )
    except Exception as exc:
        print(f"    {model}/{variant}: align failed ({type(exc).__name__}: {exc})")
        return None

    _apply_to_structure(st_pred, sup)
    _write_pdb(st_pred, outdir / f"{model}_{variant}_complex.pdb")
    if pred_mol is not None and pred_mol.GetNumConformers():
        _write_sdf(_apply_to_mol(pred_mol, sup),
                   outdir / f"{model}_{variant}_ligand.sdf",
                   f"{system}_{variant}_{model}")
        has_sdf = True
    else:
        has_sdf = False       # complex PDB still carries the ligand as HETATM
    return {
        "aligned": True,
        "ca_rmsd_a": round(sup.rmsd, 3),
        "n_ca_paired": sup.n_paired,
        "ligand_sdf": has_sdf,
        "source": str(cif.relative_to(REPO_ROOT)),
    }


# --------------------------------------------------------------------------
# crystal reference
# --------------------------------------------------------------------------

def export_crystal(system: str, outdir: Path) -> bool:
    crystal_pdb = CASF_RAW / system / f"{system}_protein.pdb"
    crystal_sdf = CASF_LIGANDS / f"{system}_ligand.sdf"
    if not crystal_pdb.exists() or not crystal_sdf.exists():
        return False
    _write_pdb(read_structure(crystal_pdb), outdir / "crystal_protein.pdb")
    mol = _read_top_pose(crystal_sdf)
    if mol is None:
        return False
    _write_sdf(mol, outdir / "crystal_ligand.sdf", f"{system}_crystal")
    return True


# --------------------------------------------------------------------------
# published RMSDs, so the views are labelled with the numbers we report
# --------------------------------------------------------------------------

def load_rmsds() -> dict[tuple[str, str, str], float]:
    """(method, system, variant) -> reported top-1 RMSD."""
    out: dict[tuple[str, str, str], float] = {}
    dock = OUTPUT_ROOT / "docking_results.csv"
    if dock.exists():
        for r in csv.DictReader(dock.open()):
            if r["module"] == "casf" and r["engine"] == "surfdock" \
                    and r["status"] == "ok" and r["rmsd_a"]:
                out[("surfdock", r["system"], r["variant"])] = float(r["rmsd_a"])
    icm = OUTPUT_ROOT / "icm_results.csv"
    if icm.exists():
        for r in csv.DictReader(icm.open()):
            if r["status"] == "ok" and r["rmsd_a"]:
                out[("icm", r["system"], r["variant"])] = float(r["rmsd_a"])
    cof = OUTPUT_ROOT / "results_full.csv"
    if cof.exists():
        model_key = {"Boltz2": "boltz2", "AF3+MSA": "af3msa"}
        for r in csv.DictReader(cof.open()):
            k = model_key.get(r["model"])
            if k and r.get("pose_idx") == "0" and r["status"] == "ok" \
                    and r.get("ligand_rmsd_a"):
                out[(k, r["pdbid"], r["variant"])] = float(r["ligand_rmsd_a"])
    return out


def load_mutation_counts() -> dict[str, int]:
    """{system: n_mutations} for the rem variant, from the build manifest."""
    man = OUTPUT_ROOT / "manifest_full_casf.json"
    if not man.exists():
        return {}
    out = {}
    for s in json.loads(man.read_text()).get("systems", []):
        if s.get("status") != "ok":
            continue
        out[s["pdbid"]] = s.get("variants", {}).get("rem", {}).get("n_mutations", 0)
    return out


def pick_systems(mode: str, rmsds: dict, limit: int | None,
                 min_mutations: int = 5) -> list[str]:
    """Select systems from the result CSVs rather than a hand-typed list.

    memorized -- WT solved AND at least one adversarial variant still < 2 A.
                 The cases that carry the memorization argument.
    collapse  -- WT solved AND every adversarial variant > 8 A. Clean physics.
    all       -- every system with a WT cell for either method.

    CRITICAL: pocket size varies a lot (median 5 mutations, range 0-13, and two
    systems have ZERO -- the known generation no-op). Retention is strongly
    dose-dependent, so an unfiltered `memorized` pick surfaces exactly the
    weakest cases: a system with 1 mutation trivially keeps its pose and proves
    nothing. `--min-mutations` (default 5) guards against that, and systems are
    ranked most-mutated first so the strongest exhibits come out on top.
    """
    nmut = load_mutation_counts()
    systems = sorted(s for s in {s for (_, s, _) in rmsds}
                     if nmut.get(s, 0) >= min_mutations
                     or mode in ("low-mutation", "extremes"))
    picked: list[tuple[float, float, str]] = []
    for s in systems:
        for method in METHODS:
            wt = rmsds.get((method, s, "wt"))
            adv = [rmsds.get((method, s, v)) for v in ("rem", "pack", "inv")]
            adv = [a for a in adv if a is not None]
            if wt is None or len(adv) < 3 or wt >= 2.0:
                continue
            # rank most-mutated first: the strongest evidence, not the weakest
            if mode == "memorized" and min(adv) < 2.0:
                picked.append((-nmut.get(s, 0), min(adv), s))
                break
            if mode == "collapse" and min(adv) > 8.0:
                picked.append((-nmut.get(s, 0), -min(adv), s))
                break
    if mode in ("low-mutation", "high-mutation", "extremes"):
        # Representative selection by perturbation size. Only systems with a
        # COMPLETE 4x4 grid, trustworthy frames and every method solving WT
        # qualify, so the panel differences are the mutation and nothing else.
        ok = []
        for sy in {s for (_, s, _) in rmsds}:
            cells = {m: {v for (mm, ss, v) in rmsds if mm == m and ss == sy}
                     for m in METHODS}
            if not all(len(cells[m]) == len(VARIANTS) for m in METHODS):
                continue
            if not all(rmsds.get((m, sy, "wt"), 99) < 2.0 for m in METHODS):
                continue
            ok.append((nmut.get(sy, 0), sy))
        ok.sort()
        k = max(1, (limit or 6) // 2)
        if mode == "low-mutation":
            return [s for _, s in ok[:limit or 3]]
        if mode == "high-mutation":
            return [s for _, s in ok[-(limit or 3):]]
        return [s for _, s in ok[:k]] + [s for _, s in ok[-k:]]

    if mode == "all":
        out = systems
    else:
        out = [s for _, _, s in sorted(picked)]
    return out[:limit] if limit else out


# --------------------------------------------------------------------------
# PyMOL session
# --------------------------------------------------------------------------

PML_HEADER = """# {system} -- pocket-mutagenesis structural contrast
# Generated by plot/07_render_mutation_views.py. Every object is in the CRYSTAL frame.
reinitialize
bg_color white
set ray_opaque_background, 1
set cartoon_transparency, 0.55
set stick_radius, 0.15
set label_size, 16
"""


def write_pml(system: str, outdir: Path, cells: dict, rmsds: dict) -> Path:
    L: list[str] = [PML_HEADER.format(system=system)]
    for name, hexcol in COLORS.items():
        r, g, b = (int(hexcol[i:i + 2], 16) / 255 for i in (1, 3, 5))
        L.append(f"set_color c_{name}, [{r:.3f}, {g:.3f}, {b:.3f}]")

    L += [
        "",
        "load crystal_protein.pdb, crystal",
        "load crystal_ligand.sdf, ref_ligand",
        "hide everything",
        "show cartoon, crystal and polymer",
        "color grey70, crystal",
        "show sticks, ref_ligand",
        "color c_crystal, ref_ligand",
        "util.cnc('ref_ligand')",
        "set_bond stick_radius, 0.22, ref_ligand",
        "",
    ]

    scenes: list[str] = []
    for method in METHODS:
        for variant in VARIANTS:
            if variant not in cells.get(method, {}):
                continue
            obj = f"{method}_{variant}"
            if method in ("surfdock", "icm"):
                L.append(f"load mutant_{variant}_receptor.pdb, {obj}_rec")
                L.append(f"load {method}_{variant}_pose.sdf, {obj}_lig")
            else:
                L.append(f"load {method}_{variant}_complex.pdb, {obj}_rec")
                if cells[method][variant].get("ligand_sdf"):
                    L.append(f"load {method}_{variant}_ligand.sdf, {obj}_lig")
                else:
                    L.append(f"create {obj}_lig, {obj}_rec and not polymer and not solvent")
            L += [
                f"hide everything, {obj}_rec or {obj}_lig",
                f"show cartoon, {obj}_rec and polymer",
                f"color c_{variant}, {obj}_rec",
                f"show sticks, {obj}_lig",
                f"color c_{variant}, {obj}_lig",
                f"util.cnc('{obj}_lig')",
                f"set_bond stick_radius, 0.22, {obj}_lig",
                # the pocket that is supposed to be destroyed: thinner sticks
                # so the pose stays visually dominant
                f"select {obj}_pocket, byres ({obj}_rec and polymer within 5 of ref_ligand)",
                f"show sticks, {obj}_pocket and not (name C+N+O)",
                f"color c_{variant}, {obj}_pocket and elem C",
                f"set_bond stick_radius, 0.10, {obj}_pocket",
                "",
            ]
            scenes.append(obj)

    L += [
        "orient ref_ligand",
        "zoom ref_ligand, 8",
        "deselect",
        "",
        "# one scene per cell: crystal reference + that cell's receptor and pose",
    ]
    for obj in scenes:
        meth, _, var = obj.rpartition("_")
        rmsd = rmsds.get((meth, system, var))
        tag = f"{rmsd:.2f} A" if rmsd is not None else "n/a"
        L += [
            "disable *",
            "enable crystal",
            "enable ref_ligand",
            f"enable {obj}_rec",
            f"enable {obj}_lig",
            # include the pose itself: an ejected pose (>8 A) falls outside a
            # crystal-ligand-only zoom and would render as an empty pocket
            f"zoom (ref_ligand or {obj}_lig), 4",
            f"scene {obj}, store, message={system} {obj} rank1 RMSD {tag}",
        ]
    L += ["", "enable *", "scene all, store", f"scene {scenes[0]}, recall" if scenes else ""]

    pml = outdir / "view.pml"
    pml.write_text("\n".join(L) + "\n")
    return pml


RENDER_TAIL = """
python
from pymol import cmd
import os
w, h = {w}, {h}
for name in cmd.get_scene_list():
    if name == 'all':
        continue
    cmd.scene(name, 'recall')
    cmd.ray(w, h)
    cmd.png(os.path.join(r'{outdir}', name + '.png'), dpi=150)
python end
"""


def render(outdir: Path, width: int, height: int) -> list[Path]:
    script = outdir / "_render.pml"
    script.write_text(
        (outdir / "view.pml").read_text()
        + RENDER_TAIL.format(w=width, h=height, outdir=str(outdir))
    )
    res = subprocess.run(
        [str(PYMOL_BIN), "-cq", str(script)],
        cwd=str(outdir), capture_output=True, text=True, timeout=900,
    )
    if res.returncode != 0:
        print(f"    PyMOL failed (rc={res.returncode}): {res.stderr.strip()[:400]}")
    return sorted(outdir.glob("*_*.png"))


def contact_sheet(system: str, outdir: Path, cells: dict, rmsds: dict) -> Path | None:
    """Assemble the rendered scenes into one labelled methods x variants grid."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = [m for m in METHODS if cells.get(m)]
    if not rows:
        return None
    fig, axes = plt.subplots(len(rows), len(VARIANTS),
                             figsize=(3.6 * len(VARIANTS), 3.2 * len(rows)),
                             squeeze=False, constrained_layout=True)
    for i, method in enumerate(rows):
        for j, variant in enumerate(VARIANTS):
            ax = axes[i][j]
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(False)
            png = outdir / f"{method}_{variant}.png"
            if png.exists():
                ax.imshow(plt.imread(png))
            else:
                ax.text(0.5, 0.5, "no cell", ha="center", va="center",
                        transform=ax.transAxes, color="grey", fontsize=11)
            rmsd = rmsds.get((method, system, variant))
            tag = f"{rmsd:.2f} Å" if rmsd is not None else "—"
            ax.set_title(f"{variant} · {tag}", fontsize=11, color=COLORS[variant])
            if j == 0:
                ax.set_ylabel(PRETTY.get(method, method), fontsize=12, labelpad=8)
    nmut = load_mutation_counts().get(system)
    nm_txt = (f" · {nmut} pocket mutation{'' if nmut == 1 else 's'}"
              if nmut is not None else "")
    fig.suptitle(
        f"{system}{nm_txt} — top-1 pose vs crystal ligand (yellow), "
        f"all in the crystal frame", fontsize=14)
    out = outdir / f"{system}_contact.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


# --------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--systems", help="comma-separated PDB ids")
    g.add_argument("--pick",
                   choices=("memorized", "collapse", "all",
                            "low-mutation", "high-mutation", "extremes"),
                   help="select systems from the result CSVs. `extremes` returns "
                        "an equal number of minimally- and maximally-perturbed "
                        "systems -- the representative pair for a figure")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--min-mutations", type=int, default=5,
                    help="skip systems whose pocket has fewer than this many "
                         "mutations (default 5 = the dataset median); retention "
                         "in a 1-2 residue pocket is not evidence of anything")
    ap.add_argument("--out", default=str(REPO_ROOT / "analysis" / "casf_mutagenesis"
                                         / "figures" / "structures"))
    ap.add_argument("--render", action="store_true", help="run headless PyMOL -> PNG")
    ap.add_argument("--contact-sheet", action="store_true",
                    help="assemble rendered scenes into one labelled grid (implies --render)")
    ap.add_argument("--width", type=int, default=1200)
    ap.add_argument("--height", type=int, default=900)
    args = ap.parse_args()

    rmsds = load_rmsds()
    if args.systems:
        systems = [s.strip().lower() for s in args.systems.split(",") if s.strip()]
    else:
        systems = pick_systems(args.pick, rmsds, args.limit, args.min_mutations)
    if not systems:
        print("No systems selected."); return 1
    # Deliberately does NOT imply --render: rebuilding a sheet from PNGs that
    # already exist is seconds, re-ray-tracing 16 scenes x 17 systems is ~15 min.
    # Pass --render explicitly when the scenes themselves need regenerating.

    root = Path(args.out)
    root.mkdir(parents=True, exist_ok=True)
    print(f"{len(systems)} system(s) -> {root}")

    n_ok = 0
    for system in systems:
        outdir = root / system
        outdir.mkdir(parents=True, exist_ok=True)
        print(f"\n[{system}]")
        if not export_crystal(system, outdir):
            print("    no crystal reference — skipped"); continue

        cells: dict[str, dict] = {m: {} for m in METHODS}
        exporters = {
            "icm": export_icm,
            "surfdock": export_surfdock,
            "boltz2": lambda sy, v, o: export_cofold("boltz2", sy, v, o),
            "af3msa": lambda sy, v, o: export_cofold("af3msa", sy, v, o),
        }
        for variant in VARIANTS:
            for method in METHODS:
                try:
                    info = exporters[method](system, variant, outdir)
                except Exception as exc:
                    print(f"    {method:<9s} {variant:<5s} export failed: "
                          f"{type(exc).__name__}: {exc}")
                    continue
                if info is None:
                    continue
                cells[method][variant] = info
                rmsd = rmsds.get((method, system, variant))
                tag = f"{rmsd:6.2f} Å" if rmsd is not None else "     —"
                print(f"    {method:<9s} {variant:<5s} rmsd={tag}"
                      f"  Cα-fit={info['ca_rmsd_a']:.2f} Å")
        cells = {m: c for m, c in cells.items() if c}
        if not cells:
            print("    no exportable cells — skipped"); continue

        pml = write_pml(system, outdir, cells, rmsds)
        (outdir / "manifest.json").write_text(json.dumps({
            "system": system,
            "frame": "crystal",
            "reference_ligand": str((CASF_LIGANDS / f"{system}_ligand.sdf")
                                    .relative_to(REPO_ROOT)),
            "cells": {m: {v: {**i, "reported_rmsd_a": rmsds.get((m, system, v))}
                          for v, i in c.items()} for m, c in cells.items()},
        }, indent=2) + "\n")
        print(f"    session: {pml}")

        if args.render:
            pngs = render(outdir, args.width, args.height)
            print(f"    rendered {len(pngs)} scene(s)")
        if args.contact_sheet:
            sheet = contact_sheet(system, outdir, cells, rmsds)
            if sheet:
                print(f"    contact sheet: {sheet}")
        n_ok += 1

    print(f"\nDone: {n_ok}/{len(systems)} system(s).")
    return 0 if n_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
