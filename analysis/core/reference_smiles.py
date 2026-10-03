"""Reference SMILES from Masters et al. 2025 — shared by the paper-reproduction arm
(analysis/paper_repro/lib/config.py re-exports them) and the ligand-mutagenesis
arm's verify gate. Moved verbatim from analysis/paper_repro/lib/config.py in S3 (2026-10)."""
from __future__ import annotations

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


# FZC inhibitor (allosteric MEK1 binder, paper Fig. 2). CCD canonical SMILES;
# 33 heavy atoms. Constant across the 4 mek1_* variants (only the protein
# changes — binding-site mutagenesis).
FZC_SMILES = "Cc1ccnc(Oc2ccc(c(Cl)c2)c3cc4[nH]nc(C)c4c(O[C@H]5CCC[C@@H](N)C5)c3)n1"


# Glucose methylation variants (SMILES from AF3 *_data.json files).
GLUCOSE_SMILES = {
    "glucose_0": "C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O",
    "glucose_1": "CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O",
    "glucose_2": "CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1OC",
    "glucose_3": "CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](OC)[C@H]1OC",
    "glucose_4": "CO[C@@H]1O[C@H](CO)[C@@H](OC)[C@H](OC)[C@H]1OC",
    "glucose_5": "COC[C@H]1O[C@@H](OC)[C@H](OC)[C@@H](OC)[C@@H]1OC",
}
