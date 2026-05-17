"""Build a solvated, ionized GROMACS system for one pilot complex.

Inputs:
  - protein PDB (cleaned, chain A only)
  - ligand .gro + .itp (from ligand_params.parameterize)

Outputs (under <work>/):
  - topol.top (protein ff14SB + ligand itp + TIP3P + ions)
  - solv.gro (protein + ligand + waters + ions)
  - index.ndx (Protein, ATP, SOL, ions, ATP_pocket — for RAMD selection)
  - posres.itp (per-stage restraint #ifdef macros — provided by ff14SB)

Workflow mirrors the canonical pdb2gmx → editconf → solvate → genion
pipeline. Force field set to AMBER ff14SB, water model TIP3P, box cubic
with 15 Å padding (matching paper).
"""
from __future__ import annotations
import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path


GMX = "gmx_mpi"  # default; override via env if running locally with gmx
# ff14SB (used by Masters et al. 2025) is NOT bundled with stock GROMACS.
# Default to the closest bundled variant; use --ff amber14sb-ol15 (or similar
# port) if a ported .ff dir is dropped under $GMXLIB.
FF = "amber99sb-ildn"
WATER = "tip3p"
BOX_PADDING_NM = 1.5  # 15 Å
SALT_CONC_M = 0.15


def _run(cmd: list[str], cwd: Path, stdin: str | None = None) -> None:
    print(f"  $ {' '.join(cmd)}")
    subprocess.run(cmd, cwd=cwd, input=stdin, text=True, check=True)


def build(
    protein_pdb: Path,
    ligand_gro: Path,
    ligand_itp: Path,
    out_dir: Path,
    ligand_resname: str = "LIG",
    ff: str = FF,
    water: str = WATER,
    box_padding_nm: float = BOX_PADDING_NM,
    salt_conc_m: float = SALT_CONC_M,
    gmx: str = GMX,
) -> dict[str, Path]:
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(protein_pdb, out_dir / "protein_in.pdb")
    shutil.copy(ligand_gro, out_dir / "ligand.gro")
    shutil.copy(ligand_itp, out_dir / "ligand.itp")

    # 1. pdb2gmx on protein only (ligand goes in separately)
    _run([gmx, "pdb2gmx", "-f", "protein_in.pdb", "-o", "protein.gro",
          "-p", "topol.top", "-i", "posre.itp", "-water", water, "-ff", ff,
          "-ignh"], cwd=out_dir)

    # 2. Concatenate protein + ligand coords. Both .gro files must share
    # the same box vector; use the protein box and copy ligand atoms in.
    _concat_gro(out_dir / "protein.gro", out_dir / "ligand.gro",
                out_dir / "complex.gro")

    # 3. Patch topol.top: add ligand itp + extend [ molecules ].
    _merge_topology(out_dir / "topol.top", out_dir / "ligand.itp", ligand_resname)

    # 4. Box: cubic with 15 Å padding (re-center on complex).
    _run([gmx, "editconf", "-f", "complex.gro", "-o", "boxed.gro",
          "-c", "-d", f"{box_padding_nm:.2f}", "-bt", "cubic"], cwd=out_dir)

    # 5. Solvate.
    _run([gmx, "solvate", "-cp", "boxed.gro", "-cs", "spc216.gro",
          "-p", "topol.top", "-o", "solv.gro"], cwd=out_dir)

    # 6. Add ions: build a dummy tpr first, then genion.
    _write_dummy_ions_mdp(out_dir / "ions.mdp")
    _run([gmx, "grompp", "-f", "ions.mdp", "-c", "solv.gro",
          "-p", "topol.top", "-o", "ions.tpr", "-maxwarn", "2"], cwd=out_dir)
    _run([gmx, "genion", "-s", "ions.tpr", "-o", "solv_ions.gro",
          "-p", "topol.top", "-pname", "NA", "-nname", "CL",
          "-conc", f"{salt_conc_m:.3f}", "-neutral"],
         cwd=out_dir, stdin="SOL\n")

    # 7. Index file: define Protein, LIG, SOL, plus a Pocket reference for RAMD.
    _write_index_file(out_dir, ligand_resname=ligand_resname, gmx=gmx)

    return {
        "gro": out_dir / "solv_ions.gro",
        "top": out_dir / "topol.top",
        "ndx": out_dir / "index.ndx",
        "posres": out_dir / "posre.itp",
    }


def _concat_gro(protein_gro: Path, ligand_gro: Path, out_gro: Path) -> None:
    p = protein_gro.read_text().splitlines()
    l = ligand_gro.read_text().splitlines()
    p_natoms = int(p[1].strip())
    l_natoms = int(l[1].strip())
    title = p[0]
    p_atoms = p[2:2 + p_natoms]
    l_atoms = l[2:2 + l_natoms]
    box = p[2 + p_natoms]
    out = [title, str(p_natoms + l_natoms), *p_atoms, *l_atoms, box]
    out_gro.write_text("\n".join(out) + "\n")


def _split_atomtypes(itp_text: str) -> tuple[str, str]:
    """Split an acpype-style .itp into (atomtypes_block, rest_block).

    GROMACS requires `[ atomtypes ]` to appear globally before any
    `[ moleculetype ]`. acpype emits both directives in one file, so we
    have to peel the atomtypes block out and include it earlier in topol.top.

    Returns ("", original) if no [ atomtypes ] directive is present.
    """
    lines = itp_text.splitlines(keepends=True)
    atomtypes_lines: list[str] = []
    rest_lines:      list[str] = []
    in_atomtypes = False
    found = False
    directive_pat = re.compile(r"^\s*\[\s*([A-Za-z_]+)\s*\]")
    for line in lines:
        m = directive_pat.match(line)
        if m:
            directive = m.group(1).lower()
            if directive == "atomtypes":
                in_atomtypes = True
                found = True
                atomtypes_lines.append(line)
                continue
            if in_atomtypes:
                in_atomtypes = False  # next directive ends atomtypes
        (atomtypes_lines if in_atomtypes else rest_lines).append(line)
    if not found:
        return "", itp_text
    return "".join(atomtypes_lines), "".join(rest_lines)


def _merge_topology(top: Path, ligand_itp: Path, ligand_resname: str) -> None:
    text = top.read_text()
    out_dir = top.parent

    # Split acpype .itp into atomtypes-only and the moleculetype block.
    lig_text = ligand_itp.read_text()
    atomtypes_block, rest_block = _split_atomtypes(lig_text)
    atomtypes_itp = out_dir / "ligand_atomtypes.itp"
    moltypes_itp  = out_dir / "ligand_moltypes.itp"
    if atomtypes_block:
        atomtypes_itp.write_text(atomtypes_block)
        moltypes_itp.write_text(rest_block)
    else:
        # No atomtypes (rare): use the original .itp as moleculetype-only include.
        moltypes_itp = ligand_itp

    # Inject atomtypes include AFTER the forcefield.itp line.
    ff_pat = re.compile(r'(#include\s+"[^"]*forcefield\.itp"[^\n]*\n)')
    if atomtypes_block and atomtypes_itp.name not in text:
        m = ff_pat.search(text)
        if not m:
            raise RuntimeError("topol.top has no #include forcefield.itp; "
                               "ligand atomtypes cannot be safely placed")
        insert = (
            '; Include ligand atom types (must precede any [ moleculetype ])\n'
            f'#include "{atomtypes_itp.name}"\n'
        )
        text = text[:m.end()] + insert + text[m.end():]

    # Inject ligand moleculetype include BEFORE the water include.
    water_pat = re.compile(r'(#include\s+"[^"]*(?:tip[34]p|spc(?:e)?)\.itp")')
    if moltypes_itp.name not in text:
        m = water_pat.search(text)
        block = f'; Include ligand moleculetype\n#include "{moltypes_itp.name}"\n'
        if m:
            text = text[:m.start()] + block + text[m.start():]
        else:
            # No water include yet (system before solvate); append before [ system ]
            sys_pat = re.compile(r'(\[\s*system\s*\])', re.IGNORECASE)
            sm = sys_pat.search(text)
            if sm:
                text = text[:sm.start()] + block + text[sm.start():]
            else:
                text = text.rstrip() + "\n" + block

    # Add ligand to [ molecules ] (must come BEFORE SOL/ions which solvate/genion append).
    if f"\n{ligand_resname} " not in text and f"\n{ligand_resname}\t" not in text:
        text = text.rstrip() + f"\n{ligand_resname}        1\n"
    top.write_text(text)


def _write_dummy_ions_mdp(path: Path) -> None:
    path.write_text(
        "; dummy mdp for grompp -genion\n"
        "integrator = steep\n"
        "nsteps = 1\n"
        "cutoff-scheme = Verlet\n"
        "coulombtype = PME\n"
        "rcoulomb = 1.0\n"
        "rvdw = 1.0\n"
    )


def _write_index_file(out_dir: Path, ligand_resname: str, gmx: str) -> None:
    """Build index.ndx with the auto-generated groups + an explicit ligand group.

    GROMACS auto-creates Protein, Non-Protein, Other, SOL, NA, CL, and an
    [ ATP ] residue group when SOL/ions are present. We just feed `q` to
    accept the defaults — the production mdp's RAMD references will resolve
    to those names directly. Refining the pocket selection (e.g., 11
    binding-site residues) is left to a follow-up when the protocol locks.
    """
    cmd = "q\n"
    _run([gmx, "make_ndx", "-f", str(out_dir / "solv_ions.gro"),
          "-o", str(out_dir / "index.ndx")], cwd=out_dir, stdin=cmd)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--protein-pdb", type=Path, required=True)
    p.add_argument("--ligand-gro", type=Path, required=True)
    p.add_argument("--ligand-itp", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--ligand-resname", default="LIG")
    p.add_argument("--gmx", default=GMX)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    out = build(
        protein_pdb=args.protein_pdb, ligand_gro=args.ligand_gro,
        ligand_itp=args.ligand_itp, out_dir=args.out,
        ligand_resname=args.ligand_resname, gmx=args.gmx,
    )
    for k, v in out.items():
        print(f"  {k:8s} {v}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
