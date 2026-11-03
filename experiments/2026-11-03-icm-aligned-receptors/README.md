# ICM aligned receptors (copied 2026-11-03)

717 `<pdbid>/<variant>/docking/receptor_aligned.pdb` files: the receptors ICM docked
into on CARC. ICM poses are in this frame. Source:
`aoxu@discovery2.usc.edu:/project2/katritch_223/aoxu/contrasCF/analysis/casf_mutagenesis/outputs/`
(rem/pack/inv 239 each). Used by `analysis/casf_mutagenesis/scripts/analyze/15_corrected_pocket_scores.py`.
The PDBs are git-ignored (84 MB); re-fetch with
`rsync -a --prune-empty-dirs --include='*/' --include='receptor_aligned.pdb' --exclude='*' <source> <here>/`.
