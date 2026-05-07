"""Single source of truth for cases, models, and native references.

All paths in this module are absolute so scripts can be run from anywhere.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path("/mnt/katritch_lab2/aoxu/contrasCF")
DATA_ROOT = REPO_ROOT / "contrasCF" / "data"          # shipped predictions
ANALYSIS_ROOT = REPO_ROOT / "analysis"                # our work
NATIVE_DIR = ANALYSIS_ROOT / "native"
RESULTS_DIR = ANALYSIS_ROOT / "results"
FIGURES_DIR = RESULTS_DIR / "figures"


# -- Native crystal references ------------------------------------------------
# Resname and chain from PDB convention, not from predicted outputs.

@dataclass(frozen=True)
class NativeRef:
    pdb_id: str
    ligand_resname: str     # 3-letter code of target ligand in the crystal
    ligand_chain: str | None = None   # None = pick any chain matching resname


NATIVES = {
    "cdk2_atp": NativeRef(pdb_id="1B38", ligand_resname="ATP"),   # CDK2-ATP
    "gdh_glc":  NativeRef(pdb_id="2VWH", ligand_resname="BGC"),   # GDH-β-D-glucose (CCD BGC)
}


# -- Cases --------------------------------------------------------------------
# `target_smiles` is the ligand SMILES used as input to the cofolding models.
# `common_smarts_key` selects a substructure for the RMSD_common metric.

@dataclass(frozen=True)
class CaseSpec:
    name: str
    family: str                       # "bindingsite" | "atp_charge" | "glucose"
    native_key: str                   # key into NATIVES
    target_smiles: str                # SMILES of the target (modified) ligand
    common_smarts_key: str            # key into COMMON_SUBSETS
    has_mg_input: bool = False        # Mg2+ in AF3 job request
    has_zn_input: bool = False        # Zn2+ (GDH crystal contains Zn)
    has_nadp_input: bool = False      # NADP cofactor (glucose family)


# ATP ligand SMILES (as used in 1B38 crystal; all-atom form with formal charges).
# The job requests use CCD_ATP for binding-site cases, so we embed the SMILES
# here for consistent downstream handling.
ATP_SMILES = (
    "c1nc(N)c2ncn([C@@H]3O[C@H](CO[P@@](=O)([O-])O[P@@](=O)([O-])OP(=O)([O-])[O-])"
    "[C@@H](O)[C@H]3O)c2n1"
)

# Native glucose in GDH (2VWH): alpha-D-glucose CCD "GLC".
GLC_SMILES = "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"

# ATP_charge variants (triphosphate replaced by neutral alkyls / quaternary amines).
ATP_CHARGE_SMILES = {
    "atp_charge_methyl": "CCCC(C)(C)COC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(N)ncnc12",
    "atp_charge_ethyl":  "CC(C)(C)CC(C)(C)COC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(N)ncnc12",
    "atp_charge_propyl": "CC(C)(C)CC(C)(C)CC(C)(C)COC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(N)ncnc12",
    "atp_charge_1": "C[N+](C)(C)COC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(N)ncnc12",
    "atp_charge_2": "C[N+](C)(C)C[N+](C)(C)COC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(N)ncnc12",
    "atp_charge_3": "C[N+](C)(C)C[N+](C)(C)C[N+](C)(C)COC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(N)ncnc12",
}

# Glucose methylation variants (SMILES from AF3 *_data.json files).
GLUCOSE_SMILES = {
    "glucose_0": "C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O",
    "glucose_1": "CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O",
    "glucose_2": "CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1OC",
    "glucose_3": "CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](OC)[C@H]1OC",
    "glucose_4": "CO[C@@H]1O[C@H](CO)[C@@H](OC)[C@H](OC)[C@H]1OC",
    "glucose_5": "COC[C@H]1O[C@@H](OC)[C@H](OC)[C@@H](OC)[C@@H]1OC",
}


CASES: dict[str, CaseSpec] = {}

# bindingsite_* — native ATP with protein mutations.
for n, has_mg in [
    ("bindingsite_wt", True),
    ("bindingsite_rem", True),
    ("bindingsite_pack", False),        # no Mg2+ in job request (data quirk)
    ("bindingsite_inv", True),
]:
    CASES[n] = CaseSpec(
        name=n, family="bindingsite", native_key="cdk2_atp",
        target_smiles=ATP_SMILES, common_smarts_key="ATP_FULL",
        has_mg_input=has_mg,
    )

# atp_charge_* — modified ATP with CDK2 unchanged.
for n, smi in ATP_CHARGE_SMILES.items():
    CASES[n] = CaseSpec(
        name=n, family="atp_charge", native_key="cdk2_atp",
        target_smiles=smi, common_smarts_key="ADENOSINE",
    )

# glucose_* — methylated glucose, GDH + NADP + Zn in input.
for n, smi in GLUCOSE_SMILES.items():
    CASES[n] = CaseSpec(
        name=n, family="glucose", native_key="gdh_glc",
        target_smiles=smi, common_smarts_key="GLC_CORE",
        has_zn_input=True, has_nadp_input=True,
    )


# -- Common subsets (SMARTS) --------------------------------------------------
# For cross-case RMSD we want the SAME set of heavy atoms paired across
# native <-> prediction, even when the ligand has extra or altered atoms.
#
# ATP_FULL     : all heavy atoms of ATP (graph is unchanged in bindingsite_*).
# ADENOSINE    : adenine + ribose + 5' CH2-O (excludes tri/alkyl/amine tail).
# GLC_CORE     : pyranose ring skeleton + ring-oxygen substituent positions
#                (excludes O-methyl carbons).
COMMON_SUBSETS = {
    "ATP_FULL":  ATP_SMILES,                               # full match (bindingsite_*)
    # Adenosine = adenine + ribose + 5'-CH2-OH (19 heavy atoms). Conserved
    # across native ATP and all atp_charge_* variants (native has α-phosphate
    # in place of the terminal OH; pairing is done on the 19 shared atoms).
    "ADENOSINE": "Nc1ncnc2n(cnc12)C1OC(CO)C(O)C1O",
    # Pyranose skeleton only (6 atoms). Conserved across native BGC and
    # glucose_0..5 methylated variants.
    "GLC_CORE":  "[C;R]1[O;R][C;R][C;R][C;R][C;R]1",
}


# -- Models -------------------------------------------------------------------

@dataclass(frozen=True)
class ModelSpec:
    name: str
    structure_globs: tuple[str, ...]   # tried in order to find the canonical file
    confidence_globs: tuple[str, ...]  # paired with structure; optional (empty => none)


MODELS: dict[str, ModelSpec] = {
    "AF3": ModelSpec(
        name="AF3",
        # bindingsite_wt/rem: `*_model_0.cif`; atp_charge/glucose: `*_model.cif`.
        # bindingsite_pack: `*_model_0_6.cif`; bindingsite_inv: `*_model_0_2.cif`.
        structure_globs=("*_model_0.cif", "*_model_0_*.cif", "*_model.cif"),
        confidence_globs=(
            "*_summary_confidences_0.json",
            "*_summary_confidences_0_*.json",
            "*_summary_confidences.json",
        ),
    ),
    "Boltz": ModelSpec(
        name="Boltz",
        structure_globs=("*_model_0.cif",),
        confidence_globs=("confidence_*_model_0.json",),
    ),
    "Boltz2": ModelSpec(
        name="Boltz2",
        structure_globs=("*_model_0.cif",),
        confidence_globs=("confidence_*_model_0.json",),
    ),
    "Chai": ModelSpec(
        name="Chai",
        structure_globs=("pred.rank_0.cif",),
        confidence_globs=("scores.rank_0.json",),
    ),
    "RFAA": ModelSpec(
        name="RFAA",
        structure_globs=("*.pdb",),
        confidence_globs=(),  # no JSON; B-factor only
    ),
    # Physics-based docking outputs. `combined_top1.pdb` is the AF3 receptor
    # with the top-1 docked pose attached as a HETATM block (written by
    # analysis/scripts/11_run_docking.py). Scores are read from `poses.sdf` tags.
    "UniDock2": ModelSpec(
        name="UniDock2",
        structure_globs=("combined_top1.pdb",),
        confidence_globs=("poses.sdf",),
    ),
    "GNINA": ModelSpec(
        name="GNINA",
        structure_globs=("combined_top1.pdb",),
        confidence_globs=("poses.sdf",),
    ),
    # SurfDock is a structure-based pose predictor; it docks into the AF3
    # receptor, so the output layout mirrors UniDock2/GNINA.
    "SurfDock": ModelSpec(
        name="SurfDock",
        structure_globs=("combined_top1.pdb",),
        confidence_globs=("poses.sdf",),
    ),
}


# Binding-site residues of CDK2 (PDB 1B38), as used by the paper for the
# binding-site mutagenesis challenges and for paper-style rendering.
CDK2_BS_RESIDUES_LIST = [10, 14, 18, 31, 33, 86, 129, 131, 132, 134, 145]
CDK2_BS_RESIDUES = "+".join(str(r) for r in CDK2_BS_RESIDUES_LIST)


def case_dir(model: str, case: str) -> Path:
    return DATA_ROOT / model / case


def find_first(d: Path, globs: tuple[str, ...]) -> Path | None:
    for pat in globs:
        hits = sorted(d.glob(pat))
        if hits:
            return hits[0]
    return None
