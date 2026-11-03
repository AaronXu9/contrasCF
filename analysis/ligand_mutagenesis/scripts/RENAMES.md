# Script renames (S3, 2026-10)

Scripts were grouped by stage (`build/`, `run/`, `analyze/`, `plot/`, `export/`) and renumbered inside each stage. Journal entries and older plans cite the old names; this table resolves them. See `plans/2026-10-01-organize-codes-results.md`.

| old | new |
|---|---|
| `00_verify_reference_systems.py` | `build/00_verify_reference_systems.py` |
| `01_build_subset20.py` | `build/01_build_subset20.py` |
| `02_build_full_casf.py` | `build/02_build_full_casf.py` |
| `03_run_boltz2.py` | `run/01_run_boltz2.py` |
| `05_analyze.py` | `analyze/01_analyze_cofold.py` |
| `06_plot_affinity.py` | `plot/01_plot_affinity.py` |
| `07_confidence_response.py` | `analyze/02_confidence_response.py` |
| `08_retention_table.py` | `analyze/03_retention_table.py` |
