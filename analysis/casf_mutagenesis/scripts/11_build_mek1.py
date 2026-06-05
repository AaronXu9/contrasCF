"""Build cofolding inputs for the MEK1 / 7XLP binding-site mutagenesis
experiment (paper Fig. 2). Uses the paper's explicit 7-residue list rather
than the strict 3.5-Å side-chain rule because that rule is documented to
miss 5 of the 7 paper residues on this system (see
scripts/00_verify_reference_systems.py for the +36 PDB offset + side-chain
distance analysis).

Outputs:
  - analysis/casf_mutagenesis/outputs/7xlp_mek1/<variant>/af3.json
                                              /boltz.yaml
                                              /docking/{ligand.sdf, box.json,
                                                        receptor_crystal.pdb}
  - contrasCF/data/AF3/mek1_<variant>/<prefix>_af3.json    (input only — runner uses this)
  - contrasCF/data/Boltz/mek1_<variant>/<prefix>_boltz.yaml
  - contrasCF/data/Boltz2/mek1_<variant>/<prefix>_boltz.yaml

The protein chains are constant across variants; only the 7 pocket residues
flip per the mutation table. AF3 / Boltz prediction outputs land later via
their respective runner scripts.

Run:
    LD_LIBRARY_PATH=/home/aoxu/miniconda3/envs/rdkit_env/lib:$LD_LIBRARY_PATH \\
        /home/aoxu/miniconda3/envs/rdkit_env/bin/python \\
        analysis/casf_mutagenesis/scripts/11_build_mek1.py
"""
from __future__ import annotations
import json
import os
import sys
from pathlib import Path

REPO_ROOT = Path(os.environ.get("CONTRASCF_ROOT", "/mnt/katritch_lab2/aoxu/contrasCF"))
sys.path.insert(0, str(REPO_ROOT / "analysis"))

from casf_mutagenesis.config import OUTPUT_ROOT, NATIVE_DIR  # noqa: E402
from casf_mutagenesis.inputs_af3 import render_af3  # noqa: E402
from casf_mutagenesis.inputs_boltz import render_boltz  # noqa: E402
from casf_mutagenesis.mutate import apply_mutations  # noqa: E402
from casf_mutagenesis.pocket import PocketResidue  # noqa: E402
from casf_mutagenesis.sequence import extract_chain_sequences  # noqa: E402

DATA_ROOT = REPO_ROOT / "contrasCF" / "data"
VARIANTS = ("wt", "rem", "pack", "inv")

# Paper's MEK1 pocket residues (paper Methods §"Removing interacting residues"):
#   A40, A59, I105, E108, M110, S158, F173.
# 7XLP auth_seq numbering is +36 relative to paper:
#   paper R → 7XLP R+36; AA identity matches at all 7 positions.
PAPER_RESIDUES = [
    ("A", 40, "ALA"),
    ("A", 59, "ALA"),
    ("I", 105, "ILE"),
    ("E", 108, "GLU"),
    ("M", 110, "MET"),
    ("S", 158, "SER"),
    ("F", 173, "PHE"),
]
PDB_OFFSET = 36   # 7XLP auth_seq = paper + 36


# FZC SMILES (CCD canonical, 33 heavy atoms). Defined here AND in
# analysis/src/config.py so the input generation and the analysis pipeline
# stay consistent. If you change one, change the other.
FZC_SMILES = "Cc1ccnc(Oc2ccc(c(Cl)c2)c3cc4[nH]nc(C)c4c(O[C@H]5CCC[C@@H](N)C5)c3)n1"


def _ensure_7xlp_cif() -> Path:
    path = NATIVE_DIR / "7XLP.cif"
    if path.exists():
        return path
    import urllib.request
    NATIVE_DIR.mkdir(parents=True, exist_ok=True)
    print(f"  downloading https://files.rcsb.org/download/7XLP.cif → {path}")
    urllib.request.urlretrieve("https://files.rcsb.org/download/7XLP.cif", path)
    return path


def _extract_crystal_protein_pdb(cif_path: Path) -> str:
    """7XLP as a plain PDB string with standard-AA atoms only.

    The casf_mutagenesis.sequence loader can't read 7XLP.cif because it's a CIF
    file format; for sequence extraction we use a gemmi structure. But we ALSO
    need a PDB file to drive `extract_chain_sequences` because that function
    accepts paths only. Quick fix: convert via gemmi.
    """
    import gemmi
    st = gemmi.read_structure(str(cif_path))
    st.setup_entities()
    # Round-trip through PDB via gemmi.
    return st.make_pdb_string()


def _build_pocket(chain_id: str = "A") -> list[PocketResidue]:
    """Construct paper-explicit pocket-residue list for 7XLP."""
    return [
        PocketResidue(
            chain=chain_id,
            resnum=paper_n + PDB_OFFSET,
            ins_code="",
            aa3=aa3,
        )
        for _, paper_n, aa3 in PAPER_RESIDUES
    ]


def _make_crystal_fzc_sdf(out_path: Path) -> None:
    """Build a 3D SDF for FZC from the 7XLP crystal geometry.

    We pull the FZC residue's heavy atoms straight from the CIF and emit an SDF
    with elements + xyz; bond orders are then perceived by RDKit. We could also
    template-match against the SMILES, but for box placement the geometry alone
    (centroid) is enough; the SDF here is mostly for docking input prep later
    (which embeds from SMILES anyway, so this file is informational).
    """
    import gemmi
    from rdkit import Chem
    from rdkit.Chem import AllChem
    cif = _ensure_7xlp_cif()
    st = gemmi.read_structure(str(cif))
    atoms = []
    for model in st:
        for chain in model:
            for res in chain:
                if res.name.strip() == "FZC":
                    for a in res:
                        if a.element.name == "H":
                            continue
                        atoms.append((a.name, a.element.name,
                                      a.pos.x, a.pos.y, a.pos.z))
                    break
    if not atoms:
        raise RuntimeError("FZC not found in 7XLP")
    # Write minimal SDF v2000
    n = len(atoms)
    lines = [
        "FZC_crystal_7XLP",
        "  contrasCF generated from 7XLP",
        "",
        f"{n:>3d}  0  0  0  0  0  0  0  0  0999 V2000",
    ]
    for _, el, x, y, z in atoms:
        sym = el[0].upper() + el[1:].lower() if len(el) > 1 else el.upper()
        lines.append(f"{x:>10.4f}{y:>10.4f}{z:>10.4f} {sym:<3s} 0  0  0  0  0  0  0  0  0  0  0  0")
    lines.append("M  END")
    lines.append("$$$$")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n")


def main() -> int:
    print("=== MEK1 / 7XLP cofolding input build ===\n")

    cif_path = _ensure_7xlp_cif()
    pdb_text = _extract_crystal_protein_pdb(cif_path)
    pdb_path = NATIVE_DIR / "7XLP.pdb"
    pdb_path.write_text(pdb_text)

    chains = extract_chain_sequences(pdb_path)
    if not chains:
        raise RuntimeError("no protein chains parsed from 7XLP")
    print(f"chains: {[(cid, len(recs)) for cid, recs in chains]}")

    # 7XLP is a single chain A (315 standard AAs); double-check.
    chain_a = next((cid, recs) for cid, recs in chains if cid == "A")
    print(f"chain A length: {len(chain_a[1])} residues")
    print(f"resnum range:   {chain_a[1][0].resnum}..{chain_a[1][-1].resnum}")

    pocket = _build_pocket("A")
    print(f"pocket residues (paper-explicit, +36 offset applied):")
    for p in pocket:
        print(f"  A {p.aa3} {p.resnum}")

    # Sanity check: the AAs in the chain at each pocket residue's PDB num
    # match the paper's expected identity.
    by_num = {(r.resnum, r.ins_code): r for r in chain_a[1]}
    for p in pocket:
        rec = by_num.get((p.resnum, p.ins_code))
        if rec is None:
            raise RuntimeError(
                f"residue {p.resnum} missing in chain A — bad numbering"
            )
        if rec.aa3 != p.aa3:
            raise RuntimeError(
                f"residue {p.resnum}: paper expects {p.aa3}, "
                f"got {rec.aa3} from 7XLP"
            )
    print("✓ pocket-residue AA identities verified against 7XLP\n")

    # Per-variant: render AF3 + Boltz inputs in two locations
    # (1) casf_mutagenesis outputs/ for the module's own pipeline
    # (2) contrasCF/data/<model>/mek1_<variant>/ for the original 16-case
    #     analysis pipeline to pick up.
    mod_out = OUTPUT_ROOT / "7xlp_mek1"

    manifest: dict = {
        "system": "7xlp_mek1",
        "pocket_residues_7xlp_numbering": [(p.chain, p.resnum, p.aa3) for p in pocket],
        "pocket_residues_paper_numbering": [(c, n, aa) for c, n, aa in PAPER_RESIDUES],
        "pdb_offset": PDB_OFFSET,
        "ligand_smiles": FZC_SMILES,
        "ligand_source": "CCD canonical (7XLP HET FZC)",
        "variants": {},
    }

    for variant in VARIANTS:
        seqs, applied = apply_mutations([chain_a], pocket, variant)
        v_dir = mod_out / variant
        v_dir.mkdir(parents=True, exist_ok=True)

        # AF3 — single-sequence placeholder (no MSA yet). Runner will overwrite
        # with the MSA-injected JSON via render_af3 + chain_msas.
        af3_json_path = v_dir / "af3.json"
        render_af3(
            name=f"mek1_{variant}",
            chain_seqs=seqs,
            ligand_smiles=FZC_SMILES,
            out_path=af3_json_path,
        )

        # Boltz — single-sequence (msa: empty).
        boltz_yaml_path = v_dir / "boltz.yaml"
        render_boltz(
            chain_seqs=seqs,
            ligand_smiles=FZC_SMILES,
            out_path=boltz_yaml_path,
        )

        # Drop copies into the 16-case layout for the runners to pick up.
        for model_dir in ("AF3", "Boltz", "Boltz2"):
            case_dir = DATA_ROOT / model_dir / f"mek1_{variant}"
            case_dir.mkdir(parents=True, exist_ok=True)
        # AF3 input goes only in the AF3 case dir as `af3_input.json`
        (DATA_ROOT / "AF3" / f"mek1_{variant}" / "af3_input.json").write_text(
            af3_json_path.read_text()
        )
        # Boltz inputs go in both Boltz and Boltz2 case dirs as
        # `mek1_<variant>.yaml` (the runner uses the YAML stem as the
        # output prefix).
        for model_dir in ("Boltz", "Boltz2"):
            (DATA_ROOT / model_dir / f"mek1_{variant}" / f"mek1_{variant}.yaml").write_text(
                boltz_yaml_path.read_text()
            )

        manifest["variants"][variant] = {
            "n_mutations": len(applied),
            "applied": [m.code() for m in applied],
            "af3_input": str(DATA_ROOT / "AF3" / f"mek1_{variant}" / "af3_input.json"),
            "boltz1_input": str(DATA_ROOT / "Boltz" / f"mek1_{variant}" / f"mek1_{variant}.yaml"),
            "boltz2_input": str(DATA_ROOT / "Boltz2" / f"mek1_{variant}" / f"mek1_{variant}.yaml"),
            "mod_dir": str(v_dir),
        }
        print(f"  {variant}: {len(applied)} mutations — {[m.code() for m in applied]}")

    # Save a crystal-FZC SDF in the casf_mutagenesis output (informational; the
    # actual docking ligand.sdf is RDKit-embedded by the docking-prep step).
    sdf_path = mod_out / "wt" / "fzc_crystal.sdf"
    _make_crystal_fzc_sdf(sdf_path)
    manifest["crystal_fzc_sdf"] = str(sdf_path)

    manifest_path = mod_out / "manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print(f"\nmanifest: {manifest_path}")
    print("done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
