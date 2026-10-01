#!/usr/bin/env python3
"""Is the retained pose still PHYSICALLY favourable in the mutated pocket?

`34_validate_cofold_pockets.py` showed the surprising result that when Boltz-2
keeps the crystal pose in a mutated pocket, the resulting complex is as clean as
its own WT prediction (88.0% PoseBusters pass vs 87.5% for WT). Two readings:

  (a) the model imposes a remembered pose and the mutation simply did not
      destroy the site -- retention is defensible, or
  (b) the site IS degraded and the model is insensitive to it.

Those are distinguished by an INDEPENDENT physics score. This script holds the
CRYSTAL ligand pose fixed and rescores it (GNINA `--score_only`) in each
receptor:

    wt          -> the crystal receptor         (baseline)
    rem/pack/inv-> the AF3-predicted mutant receptor that docking was handed

The ligand is never redocked and never moved relative to the protein: the
crystal pose is mapped into the mutant receptor's frame with the INVERSE of the
same Ca transform the RMSD analysis uses, so the only thing that changes between
rows is the pocket. A large drop in CNNaffinity/Vina score means the site really
was degraded; a small drop means it was not.

Outputs `outputs/crystal_pose_rescore.csv` with, per (system, variant):
  vina_affinity, cnn_score, cnn_affinity, ca_fit_rmsd_a

Run (any env with rdkit + gemmi; GNINA is a subprocess):
  .../rdkit_env/bin/python 35_rescore_crystal_pose.py [--limit N]
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
sys.path.insert(0, str(REPO_ROOT / "analysis"))

from rdkit import Chem, RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")

import gemmi  # noqa: E402

from casf_mutagenesis.config import CASF_LIGANDS, CASF_RAW, OUTPUT_ROOT, VARIANTS  # noqa: E402

AA3 = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR "
          "TRP TYR VAL MSE".split())
BACKBONE = {"N", "CA", "C", "O", "OXT"}
MODEL_CIF = {"boltz2": "{s}_{v}_model_*.cif", "af3msa": "af3msa_{s}_{v}_model_*.cif"}
from casf_mutagenesis.gnina_analysis import _ca_transform, _heavy_coords, _read_top_pose  # noqa: E402

GNINA_BIN = os.environ.get(
    "CONTRASCF_GNINA_BIN", "/home/aoxu/projects/PoseBench/forks/GNINA/gnina")
TIMEOUT_S = 300

PATTERNS = {
    "vina_affinity": re.compile(r"^Affinity:\s+(-?[\d.]+)", re.M),
    "cnn_score": re.compile(r"^CNNscore:\s+([\d.]+)", re.M),
    "cnn_affinity": re.compile(r"^CNNaffinity:\s+(-?[\d.]+)", re.M),
}


def score_only(receptor: Path, ligand: Path, cleanup: tuple[Path, ...] = ()) -> dict:
    """Score one complex, then delete any temp inputs we created for it.

    `cleanup` exists because an earlier run wrote 2 files per cell x ~3000 cells
    into TMPDIR and never removed them, filling the root partition and dying at
    240/251 with ENOSPC. Temp files are per-cell; they must not outlive the cell.
    """
    try:
        return _score_only(receptor, ligand)
    finally:
        for f in cleanup:
            try:
                f.unlink()
            except OSError:
                pass


def _score_only(receptor: Path, ligand: Path) -> dict:
    cmd = [GNINA_BIN, "--score_only", "-r", str(receptor), "-l", str(ligand),
           "--cnn_scoring", "rescore", "--seed", "42"]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT_S)
    except subprocess.TimeoutExpired:
        return {"status": "timeout"}
    if res.returncode != 0:
        return {"status": f"gnina_exit_{res.returncode}"}
    out = {"status": "ok"}
    for key, pat in PATTERNS.items():
        m = pat.search(res.stdout)
        out[key] = float(m.group(1)) if m else None
    if out.get("cnn_affinity") is None:
        return {"status": "unparsed"}
    return out


def self_complex(model: str, system: str, variant: str, tmp: Path):
    """(protein pdb, ligand sdf) for a model's OWN rank-0 prediction.

    Protein and ligand are mutually consistent by construction (the model built
    them together), so this asks a different question from the crystal-pose test:
    not "is the native pose viable in a mutant pocket" but "is the arrangement
    Boltz-2 actually produced physically sound".
    """
    cifs = sorted((OUTPUT_ROOT / system / variant).glob(
        MODEL_CIF[model].format(s=system, v=variant)))
    if not cifs:
        return None, None
    st = gemmi.read_structure(str(cifs[0]))
    st.setup_entities()
    st.remove_hydrogens()
    lig, lig_n = None, 0
    for model in st:
        for ch in model:
            for res in ch:
                if res.name in AA3 or res.name in ("HOH", "WAT"):
                    continue
                n = sum(1 for a in res if a.element.name != "H")
                if n > lig_n:
                    lig, lig_n = res, n
    if lig is None:
        return None, None
    lp = tmp / f"{system}_{variant}_{model}_lig.sdf"
    blk = ["HETATM%5d %-4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s" %
           (i + 1, a.name[:4], "LIG", 1, a.pos.x, a.pos.y, a.pos.z,
            a.element.name.upper())
           for i, a in enumerate(lig) if a.element.name != "H"]
    m = Chem.MolFromPDBBlock("\n".join(blk) + "\nEND\n", removeHs=True, sanitize=False)
    if m is None:
        return None, None
    w = Chem.SDWriter(str(lp)); w.write(m); w.close()
    pp = tmp / f"{system}_{variant}_{model}_prot.pdb"
    s2 = st.clone(); s2.setup_entities()
    s2.remove_ligands_and_waters(); s2.remove_empty_chains()
    s2.write_pdb(str(pp))
    return pp, lp


def crystal_rem_receptor(system: str, tmp: Path) -> Path | None:
    """The CRYSTAL receptor with every pocket side chain deleted (= all-Gly).

    This is the only variant we can build without a model and without a
    superposition: `rem` is a pure deletion, so truncating the crystal residues
    to backbone reproduces it exactly, in the crystal frame. Scoring the crystal
    ligand here answers "what did removing these side chains actually cost?"
    with no prediction and no alignment error anywhere in the path.
    (`pack`/`inv` would need side-chain packing, so they get no such control.)
    """
    man = OUTPUT_ROOT / "manifest_full_casf.json"
    if not man.exists():
        return None
    entry = next((e for e in json.loads(man.read_text())["systems"]
                  if e.get("pdbid") == system and e.get("status") == "ok"), None)
    if not entry or not entry.get("pocket"):
        return None
    targets = {(p["chain"], p["resnum"], p.get("ins", "").strip())
               for p in entry["pocket"]}
    src = CASF_RAW / system / f"{system}_protein.pdb"
    if not src.exists():
        return None
    out = tmp / f"{system}_rem_truncated.pdb"
    kept = []
    for line in src.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            key = (line[21], int(line[22:26]), line[26].strip())
        except ValueError:
            continue
        if key in targets and line[12:16].strip() not in BACKBONE:
            continue                       # delete this side-chain atom
        kept.append(line)
    if not kept:
        return None
    out.write_text("\n".join(kept) + "\nEND\n")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--systems", default=None)
    ap.add_argument("--out", default=str(OUTPUT_ROOT / "crystal_pose_rescore.csv"))
    args = ap.parse_args()

    if args.systems:
        systems = [s.strip().lower() for s in args.systems.split(",") if s.strip()]
    else:
        systems = sorted(d.name for d in OUTPUT_ROOT.iterdir()
                         if d.is_dir() and len(d.name) == 4)
    if args.limit:
        systems = systems[:args.limit]

    tmp = Path(tempfile.mkdtemp(prefix="rescore_"))
    rows: list[dict] = []
    keys = ["system", "variant", "source", "status", "vina_affinity", "cnn_score",
            "cnn_affinity", "ca_fit_rmsd_a"]
    fh = open(args.out, "w", newline="")
    writer = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
    writer.writeheader(); fh.flush()

    def emit(row: dict) -> None:
        rows.append(row)
        writer.writerow(row)
        fh.flush()          # a multi-hour run must be readable while running

    t0 = time.time()
    for n, system in enumerate(systems, 1):
        crystal_sdf = CASF_LIGANDS / f"{system}_ligand.sdf"
        crystal_pdb = CASF_RAW / system / f"{system}_protein.pdb"
        if not crystal_sdf.exists() or not crystal_pdb.exists():
            continue
        cm = _read_top_pose(crystal_sdf)
        if cm is None:
            continue
        crys_xyz = _heavy_coords(cm, list(range(cm.GetNumHeavyAtoms())))

        for variant in VARIANTS:
            receptor = OUTPUT_ROOT / system / variant / "docking" / "receptor.pdb"
            if not receptor.exists():
                continue
            rec = {"system": system, "variant": variant, "ca_fit_rmsd_a": None}

            if variant == "wt":
                # receptor IS the crystal protein -> pose already in its frame
                lig_path = crystal_sdf
                rec["ca_fit_rmsd_a"] = 0.0
            else:
                # map the crystal pose INTO the mutant receptor's frame using the
                # inverse of the receptor->crystal transform (x_rec = (x_cry - t) R)
                try:
                    sup = _ca_transform(receptor, crystal_pdb, crys_xyz)
                except Exception as exc:
                    emit({**rec, "source": "crystal_in_docking_receptor",
                          "status": f"align_failed:{type(exc).__name__}"})
                    continue
                rec["ca_fit_rmsd_a"] = round(sup.rmsd, 3)
                mol = Chem.Mol(cm)
                conf = mol.GetConformer()
                xyz = np.array([list(conf.GetAtomPosition(i))
                                for i in range(mol.GetNumAtoms())])
                xyz = (xyz - sup.t) @ sup.R
                for i, p in enumerate(xyz):
                    conf.SetAtomPosition(i, (float(p[0]), float(p[1]), float(p[2])))
                lig_path = tmp / f"{system}_{variant}.sdf"
                w = Chem.SDWriter(str(lig_path)); w.write(mol); w.close()

            emit({**rec, "source": "crystal_in_docking_receptor",
                  **score_only(receptor, lig_path,
                               cleanup=() if variant == "wt" else (lig_path,))})

            # Alignment-free arm: each model's OWN complex (protein and ligand
            # mutually consistent by construction, so no superposition is used).
            for model in ("boltz2", "af3msa"):
                pp, lp = self_complex(model, system, variant, tmp)
                if pp is not None:
                    emit({"system": system, "variant": variant,
                          "ca_fit_rmsd_a": None, "source": f"{model}_self",
                          **score_only(pp, lp, cleanup=(pp, lp))})

            # Model-free control, `rem` only: crystal receptor truncated to Gly.
            if variant == "rem":
                tr = crystal_rem_receptor(system, tmp)
                if tr is not None:
                    emit({"system": system, "variant": "rem",
                          "ca_fit_rmsd_a": 0.0,
                          "source": "crystal_truncated_in_place",
                          **score_only(tr, crystal_sdf, cleanup=(tr,))})

        if n % 10 == 0:
            el = time.time() - t0
            print(f"  [{n}/{len(systems)}] {len(rows)} cells  "
                  f"{el/60:.1f} min  eta {el/n*(len(systems)-n)/60:.0f} min", flush=True)

    fh.close()
    ok = sum(1 for r in rows if r.get("status") == "ok")
    print(f"\n{ok}/{len(rows)} cells scored -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
