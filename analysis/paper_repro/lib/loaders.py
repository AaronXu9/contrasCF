"""Compatibility shim — implementation moved to analysis/core/loaders.py (S3, 2026-10).

Re-exports every name, private ones included, so `from loaders import X` keeps
working for the paper-reproduction arm."""
import sys as _sys
from pathlib import Path as _Path

_ANALYSIS = _Path(__file__).resolve().parents[2]   # analysis/paper_repro/lib -> analysis
if str(_ANALYSIS) not in _sys.path:
    _sys.path.insert(0, str(_ANALYSIS))
from core import loaders as _impl  # noqa: E402

globals().update({k: v for k, v in vars(_impl).items() if not k.startswith("__")})
