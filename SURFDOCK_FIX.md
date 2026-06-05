# SurfDock Rank-1 Positioning Fix

## Problem
SurfDock rank-1 poses were placing ~30 Å from the binding pocket instead of anchoring to the pocket center.

## Root Cause
**Bug in `/home/aoxu/projects/SurfDock/inference_accelerate.py` line 208**

The pocket center coordinates from the CSV were being parsed incorrectly:

```python
# BUGGY (creates shape (1,3) instead of (3,))
new_pocket_centers.append(np.array([(float(x),float(y),float(z))]))
```

The extra parentheses `[(...)]` create a 2D array with shape (1,3):
```
array([[-9.974, -4.702, -1.542]])  # shape (1,3)
```

When this was passed to `ScreenDataset`, it was supposed to be flattened to shape (3,):
```
array([-9.974, -4.702, -1.542])  # shape (3,) — CORRECT
```

This shape mismatch caused dimensional errors in the tensor operations in `score_dataset.py`:
```python
complex_graph['receptor'].pocket_center = torch.tensor(self.pocket_center).float() - complex_graph.original_center
```

The `--ligand_to_pocket_center` flag in SurfDock's sampling code was never able to apply the constraint because `pocket_center` was malformed, so ligands were sampled with random translations (Normal(0, tr_sigma_max)) that could drift far from the true binding site.

## Solution
Fix the array creation to remove the extra parentheses:

```python
# FIXED (creates shape (3,))
new_pocket_centers.append(np.array([float(x),float(y),float(z)]))
```

This produces:
```
array([-9.974, -4.702, -1.542])  # shape (3,) — CORRECT
```

## How to Apply
Run the fix script **before running SurfDock**:
```bash
python3 analysis/scripts/apply_surfdock_fix.py
```

This script:
1. Checks if the fix is already applied (safe to run multiple times)
2. Uses regex to find and replace the buggy line
3. Verifies the fix was successful

The fix is idempotent and will report "already applied" if run again.

## Verification
After running the fix, check that the file was patched:
```bash
grep "np.array\(\[float(x),float(y),float(z)\]\)" /home/aoxu/projects/SurfDock/inference_accelerate.py
```

You should see one match on line 208.

## Re-running SurfDock
After applying the fix, delete previous SurfDock output and re-run:
```bash
rm -rf contrasCF/data/SurfDock
python3 analysis/scripts/13_run_surfdock.py
```

The rank-1 poses should now be correctly positioned at the binding pocket center instead of 30+ Å away.

## Technical Details
- **Affected file**: `/home/aoxu/projects/SurfDock/inference_accelerate.py` (external tool, not contrasCF repo)
- **Affected line**: 208
- **Root cause**: Numpy array shape mismatch (1,3) vs (3,)
- **Impact**: All SurfDock runs with `--ligand_to_pocket_center` flag had corrupted pocket center data
- **Workaround**: Only the ligand position hint (via `_translate_ligand_to_pocket`) was being used, which provides some constraint but isn't as reliable as explicit pocket center anchoring

## Files Modified
- `/home/aoxu/projects/SurfDock/inference_accelerate.py` (patched via `apply_surfdock_fix.py`)

## Scripts Involved
- `/mnt/katritch_lab2/aoxu/contrasCF/analysis/scripts/13_run_surfdock.py` - Injects pocket_center into CSV (works correctly once SurfDock is patched)
- `/mnt/katritch_lab2/aoxu/contrasCF/analysis/scripts/apply_surfdock_fix.py` - Applies the SurfDock fix
- `/mnt/katritch_lab2/aoxu/contrasCF/analysis/casf_mutagenesis/scripts/14_run_surfdock_variants.py` - Similar runner for variants (same bug applies)
