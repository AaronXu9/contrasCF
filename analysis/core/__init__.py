"""Shared core used by every contrasCF arm.

  ligand_rmsd      atom correspondence + in-place / best-fit ligand RMSD
  loaders          structure / ligand loading (gemmi, RDKit)
  surfdock_engine  SurfDock 4-step pipeline helpers (surface -> CSV -> ESM -> diffusion)
  reference_smiles Masters et al. 2025 reference SMILES (ATP, glucose, FZC, ...)

Created in S3 of plans/2026-10-01-organize-codes-results.md so that no arm reaches
into another arm's scripts/ by file path. Import as `core.<module>` with
`analysis/` on sys.path.
"""
