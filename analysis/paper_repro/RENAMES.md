# Paper-reproduction arm moved (S3, 2026-10)

| old | new |
|---|---|
| `analysis/src/` | `analysis/paper_repro/lib/` |
| `analysis/scripts/` | `analysis/paper_repro/scripts/` |
| `analysis/src/loaders.py` | `analysis/core/loaders.py` (shim left in `lib/`) |
| reference SMILES in `analysis/src/config.py` | `analysis/core/reference_smiles.py` (re-exported by `lib/config.py`) |
| SurfDock helpers in `analysis/scripts/13_run_surfdock.py` | `analysis/core/surfdock_engine.py` |

`analysis/native/` and `analysis/results/` did not move.
