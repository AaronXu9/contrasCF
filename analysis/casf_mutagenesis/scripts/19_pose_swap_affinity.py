#!/usr/bin/env python3
"""Pose-swap test for Boltz-2's affinity head (hook B from the design spec
docs/superpowers/specs/2026-06-06-pose-swap-test-design.md).

Question: is the affinity head functionally pose-insensitive, or pose-sensitive
but the CASF mutations just don't move the ligand enough? We hold the trunk fixed
and vary ONLY the coordinates the head distograms over (x_pred): translate the
ligand out of the pocket by 0/2/5/10/20 A and read affinity_pred_value.

Mechanism: the affinity head reads pose only via cdist of x_pred over
protein-ligand cross-pairs (affinity.py:104-108); s_inputs/z are pose-independent
trunk outputs. We monkeypatch AffinityModule.forward so that, on the production
call (native x_pred), we ALSO score a ladder of decoy x_pred with the same module
and inputs, then return the native result unchanged. The native row reproducing
production affinity is the identity check (automatic by construction).

Run (boltzina_env, GPU):
    CUDA_VISIBLE_DEVICES=0 /mnt/katritch_lab2/aoxu/envs/boltzina_env/bin/python \
        analysis/casf_mutagenesis/scripts/19_pose_swap_affinity.py \
        --yaml <affinity-enabled boltz.yaml> --out /tmp/poseswap_<id> [--tag <id>]
"""
from __future__ import annotations
import argparse
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("CUDA_VISIBLE_DEVICES", "0")

import numpy as np
import torch

# MW-correction constants (boltz2.py:687-697)
MW_MODEL_COEF, MW_COEF, MW_BIAS = 1.03525938, -0.59992683, 2.83288489
EXIT_LADDER = [0.0, 2.0, 5.0, 10.0, 20.0]   # A along pocket-exit vector
RESULTS: list[dict] = []                    # filled inside the patched forward
_DECOYS_CACHE: dict = {}


def _rep_atom_indices(feats, tok_mask):
    """Atom-row indices in x_pred that each masked token gathers (argmax of the
    token_to_rep_atom gather matrix)."""
    t2a = feats["token_to_rep_atom"][0].float()          # [N_tok, N_atom]
    idx = t2a.argmax(dim=-1)                              # [N_tok]
    return idx[tok_mask[0].bool()].tolist()


def _build_decoys(x_pred, feats):
    """Return list of (pose_id, disp_A, axis, x_pred_decoy). Translates only the
    ligand rep-atom rows; receptor rows fixed. x_pred: [B,mult,N,3] or [BM,N,3]."""
    g = torch.Generator(device="cpu").manual_seed(0)
    pad = feats["token_pad_mask"][0].bool()
    lig_tok = feats["affinity_token_mask"][0].bool() & pad
    rec_tok = (feats["mol_type"][0] == 0) & pad
    lig_atoms = _rep_atom_indices(feats, lig_tok.unsqueeze(0))
    rec_atoms = _rep_atom_indices(feats, rec_tok.unsqueeze(0))

    x = x_pred.detach().float().cpu()
    flat = x.reshape(-1, x.shape[-2], 3)[0]               # [N_atom,3], first sample
    lig_xyz = flat[lig_atoms]
    rec_xyz = flat[rec_atoms]
    lig_com = lig_xyz.mean(0)
    centroid = rec_xyz.mean(0)
    Rprot = float((rec_xyz - centroid).norm(dim=1).max())        # protein bounding radius
    outward = lig_com - centroid
    outward = outward / outward.norm() if float(outward.norm()) > 2.0 else torch.tensor([0., 0., 1.])

    def translate_lig(vec):
        xd = x.clone(); view = xd.reshape(-1, xd.shape[-2], 3)
        for a in lig_atoms:
            view[:, a, :] += vec
        return xd.to(x_pred.dtype).to(x_pred.device)

    def mindist(xp):
        f = xp.detach().float().cpu().reshape(-1, xp.shape[-2], 3)[0]
        return float(torch.cdist(f[lig_atoms], f[rec_atoms]).min())

    decoys = [("native", 0.0, "none", x_pred)]
    # EJECTION ladder: place the ligand COM BEYOND the protein bounding sphere so
    # it genuinely clears all protein contacts (min-dist ~ clearance), regardless
    # of how buried the pocket is. The achieved min ligand-protein distance is
    # recorded per pose and is the real x-axis (raw translation is not).
    for c in (5.0, 15.0, 30.0):
        target = centroid + outward * (Rprot + c)
        decoys.append((f"eject{int(c)}", c, "radial", translate_lig(target - lig_com)))
    rd = torch.randn(3, generator=g); rd = rd / rd.norm()         # direction control
    decoys.append(("eject_rand", 15.0, "rand", translate_lig((centroid + rd * (Rprot + 15.0)) - lig_com)))
    # in-pocket rotation (orientation scramble, no translation)
    th = float(torch.rand(1, generator=g) * 2 * np.pi)
    ax = torch.randn(3, generator=g); ax = ax / ax.norm()
    K = torch.tensor([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
    Rot = torch.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)
    xr = x.clone(); vr = xr.reshape(-1, xr.shape[-2], 3)
    for a in lig_atoms:
        vr[:, a, :] = (Rot @ (vr[:, a, :].T - lig_com[:, None])).T + lig_com
    decoys.append(("rot_inpocket", 0.0, "rot", xr.to(x_pred.dtype).to(x_pred.device)))
    # SANITY: translate the WHOLE complex far. cdist is translation-invariant, so
    # this MUST give Delta-aff ~ 0 and min-dist ~ native (proves the harness).
    xw = x.clone(); xw.reshape(-1, xw.shape[-2], 3)[:, :, :] += outward * (Rprot + 30.0)
    decoys.append(("whole_far", 0.0, "sanity", xw.to(x_pred.dtype).to(x_pred.device)))
    decoys = [(pid, disp, axis, mindist(xp), xp) for pid, disp, axis, xp in decoys]
    return decoys, len(lig_atoms), len(rec_atoms)


def _install_patch():
    import boltz.model.modules.affinity as am
    orig = am.AffinityModule.forward
    state = {"call": 0}

    def patched(self, s_inputs, z, x_pred, feats, multiplicity=1, use_kernels=False):
        native_out = orig(self, s_inputs, z, x_pred, feats, multiplicity, use_kernels)
        call = state["call"]; state["call"] += 1
        key = int(x_pred.shape[-2])
        if key not in _DECOYS_CACHE:
            _DECOYS_CACHE[key] = _build_decoys(x_pred, feats)
        decoys, n_lig, n_rec = _DECOYS_CACHE[key]
        mw = float(feats["affinity_mw"][0]) if "affinity_mw" in feats else float("nan")
        with torch.no_grad():
            for pose_id, disp, axis, min_dist, xp in decoys:
                out = orig(self, s_inputs, z, xp, feats, multiplicity, use_kernels)
                RESULTS.append(dict(
                    module=call, pose_id=pose_id, disp_A=disp, axis=axis, min_lig_prot=min_dist,
                    aff_value=float(out["affinity_pred_value"].reshape(-1)[0]),
                    logits=float(out["affinity_logits_binary"].reshape(-1)[0]),
                    mw=mw, n_lig_atoms=n_lig, n_rec_atoms=n_rec,
                ))
        return native_out

    am.AffinityModule.forward = patched


def _aggregate(tag, out_dir):
    import pandas as pd
    df = pd.DataFrame(RESULTS)
    if df.empty:
        print("[poseswap] NO affinity calls captured — did the yaml enable affinity?")
        return
    # ensemble over modules (module 0/1), then MW correction + probability
    def sig(x): return 1 / (1 + np.exp(-x))
    rows = []
    for pose_id, g in df.groupby("pose_id", sort=False):
        ens = g.aff_value.mean()                          # raw ensemble = production value
        mw = g.mw.iloc[0]
        aff_mw = MW_MODEL_COEF * ens + MW_COEF * (mw ** 0.3) + MW_BIAS   # constant offset/ladder
        rows.append(dict(tag=tag, pose_id=pose_id, disp_A=g.disp_A.iloc[0], axis=g.axis.iloc[0],
                         min_lig_prot=float(g.min_lig_prot.iloc[0]) if "min_lig_prot" in g else float("nan"),
                         boltz_aff=ens, boltz_aff_mw=aff_mw, boltz_prob=float(sig(g.logits).mean()),
                         mw=mw, n_lig=int(g.n_lig_atoms.iloc[0]),
                         aff_m0=g[g.module == 0].aff_value.mean(),
                         aff_m1=g[g.module == 1].aff_value.mean()))
    summ = pd.DataFrame(rows)
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    summ.to_csv(Path(out_dir) / f"poseswap_{tag}.csv", index=False)
    nat = summ[summ.pose_id == "native"].boltz_aff.iloc[0]
    print(f"\n=== POSE-SWAP affinity (system={tag}; lower=tighter; native={nat:.3f}) ===")
    for _, r in summ.iterrows():
        print(f"  {r.pose_id:13s} d={r.disp_A:4.0f}A {r.axis:5s} minD={r.min_lig_prot:5.1f}A  "
              f"aff={r.boltz_aff:+.3f}  Δvs_native={r.boltz_aff-nat:+.3f}  P(bind)={r.boltz_prob:.3f}")
    far = summ[summ.pose_id == "eject30"]
    if len(far):
        gap = far.boltz_aff.iloc[0] - nat
        print(f"\n  native→eject30 gap = {gap:+.3f} log units  (ligand now {far.min_lig_prot.iloc[0]:.0f} A "
              f"from protein; physics wants strongly +; ~0 ⇒ pose-insensitive)")
    print(f"[wrote] {Path(out_dir)/f'poseswap_{tag}.csv'}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--yaml", required=True, help="affinity-enabled boltz.yaml")
    ap.add_argument("--out", required=True, help="output dir (boltz work + results)")
    ap.add_argument("--tag", default=None)
    ap.add_argument("--diffusion_samples", type=int, default=5)
    ap.add_argument("--recycling_steps", type=int, default=5)
    ap.add_argument("--sampling_steps", type=int, default=200)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    tag = args.tag or Path(args.yaml).resolve().parent.name

    work = Path(args.out); work.mkdir(parents=True, exist_ok=True)
    yaml_renamed = work / f"{tag}.yaml"
    import shutil; shutil.copy(args.yaml, yaml_renamed)

    _install_patch()
    from boltz.main import cli
    argv = ["predict", str(yaml_renamed), "--out_dir", str(work), "--model", "boltz2",
            "--output_format", "mmcif", "--diffusion_samples", str(args.diffusion_samples),
            "--recycling_steps", str(args.recycling_steps), "--sampling_steps", str(args.sampling_steps),
            "--seed", str(args.seed)]
    print(f"[poseswap] boltz predict {' '.join(argv[1:])}")
    try:
        cli(argv, standalone_mode=False)
    except SystemExit:
        pass
    _aggregate(tag, work)


if __name__ == "__main__":
    main()
