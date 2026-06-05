#!/usr/bin/env python3
"""Apply the pocket_center shape fix to SurfDock's inference_accelerate.py

The bug: line 208 in /home/aoxu/projects/SurfDock/inference_accelerate.py creates
a numpy array with shape (1,3) instead of (3,), causing dimensional mismatch when
tensor operations attempt to subtract original_center.

This fix rewrites that line to create the correct shape.
"""
import re
import sys
from pathlib import Path

SURFDOCK_INFERENCE = Path("/home/aoxu/projects/SurfDock/inference_accelerate.py")

def apply_fix():
    """Apply the shape fix to inference_accelerate.py"""
    if not SURFDOCK_INFERENCE.exists():
        print(f"ERROR: {SURFDOCK_INFERENCE} not found", file=sys.stderr)
        return False

    content = SURFDOCK_INFERENCE.read_text()

    # Check if fix is already applied (look for the correct version)
    if "np.array([float(x),float(y),float(z)])" in content:
        print(f"[apply_surfdock_fix] Fix already applied to {SURFDOCK_INFERENCE}")
        return True

    # Apply the fix: change np.array([(float(x),float(y),float(z))])
    # to np.array([float(x),float(y),float(z)])
    pattern = r'new_pocket_centers\.append\(np\.array\(\[\(float\(x\),float\(y\),float\(z\)\)\]\)\)'
    replacement = 'new_pocket_centers.append(np.array([float(x),float(y),float(z)]))'

    new_content = re.sub(pattern, replacement, content)

    if new_content == content:
        print(f"ERROR: Could not find pattern to fix in {SURFDOCK_INFERENCE}", file=sys.stderr)
        print("Pattern:", pattern, file=sys.stderr)
        return False

    SURFDOCK_INFERENCE.write_text(new_content)
    print(f"[apply_surfdock_fix] Successfully patched {SURFDOCK_INFERENCE}")
    print("[apply_surfdock_fix]   Changed: np.array([(float(x),float(y),float(z))])")
    print("[apply_surfdock_fix]   To:      np.array([float(x),float(y),float(z)])")
    return True

if __name__ == "__main__":
    success = apply_fix()
    sys.exit(0 if success else 1)
