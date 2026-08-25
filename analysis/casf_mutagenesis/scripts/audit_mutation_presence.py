"""4b deep-dive v2 — does each mutant receptor.pdb actually CONTAIN its mutations?

v1 was WRONG: it required the full spec chain to sit inside the receptor
sequence, so AF3's routine 2-3 residue terminal trimming was misread as
"chain missing" (false NOCHAIN on 1q8t, 4ivb, 1ydr, ...).

v2 aligns the RECEPTOR (shorter, actual model) onto the SPEC sequence with
difflib opcodes, giving a spec_index -> receptor_index map that tolerates
indels. Each specified mutation is then classified:

  present      position modeled AND carries the mutant residue   (correct)
  wrong_aa     position modeled but carries some other residue
  not_modeled  position absent from the receptor                 (silently lost)
  no_chain     the whole mutated chain is absent from the receptor
"""
from __future__ import annotations
import glob, json, os
from collections import defaultdict
from difflib import SequenceMatcher

import gemmi
import yaml

ROOT = "analysis/casf_mutagenesis/outputs"
VARIANTS = ("rem", "pack", "inv")
AA3to1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
    "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
    "THR":"T","TRP":"W","TYR":"Y","VAL":"V","MSE":"M","SEC":"U","PYL":"O",
}


def yaml_seqs(path):
    try:
        d = yaml.safe_load(open(path))
    except Exception:
        return {}
    out = {}
    for s in d.get("sequences", []) or []:
        for k, v in (s or {}).items():
            if k.lower() != "protein" or not isinstance(v, dict) or not v.get("sequence"):
                continue
            ids = v.get("id")
            ids = ids if isinstance(ids, list) else [ids]
            for cid in ids:
                out[str(cid)] = v["sequence"]
    return out


def pdb_chain_seqs(path):
    try:
        st = gemmi.read_structure(path)
    except Exception:
        return {}
    st.remove_ligands_and_waters()
    out = {}
    for ch in st[0]:
        s = "".join(AA3to1[r.name] for r in ch if r.name in AA3to1)
        if s:
            out[ch.name] = s
    return out


def spec_to_rec_map(spec: str, rec: str):
    """{spec_index: rec_index} via difflib matching blocks (tolerates indels)."""
    m = {}
    for a, b, size in SequenceMatcher(None, spec, rec, autojunk=False).get_matching_blocks():
        for k in range(size):
            m[a + k] = b + k
    return m


def main():
    rows = []
    systems = sorted({p.split("/")[-4] for p in glob.glob(f"{ROOT}/*/wt/docking/receptor.pdb")})
    for sysid in systems:
        wt_seqs = yaml_seqs(f"{ROOT}/{sysid}/wt/boltz.yaml")
        if not wt_seqs:
            continue
        for var in VARIANTS:
            rec_path = f"{ROOT}/{sysid}/{var}/docking/receptor.pdb"
            if not os.path.exists(rec_path):
                continue
            mut_seqs = yaml_seqs(f"{ROOT}/{sysid}/{var}/boltz.yaml")
            muts = defaultdict(list)
            for cid, ws in wt_seqs.items():
                ms = mut_seqs.get(cid)
                if not ms or len(ms) != len(ws):
                    continue
                for i, (a, b) in enumerate(zip(ws, ms)):
                    if a != b:
                        muts[cid].append((i, a, b))
            n_mut = sum(len(v) for v in muts.values())
            if n_mut == 0:
                rows.append(dict(pdbid=sysid, variant=var, n_mut=0, present=0,
                                 wrong_aa=0, not_modeled=0, no_chain=0, verdict="NOMUT"))
                continue

            rec = pdb_chain_seqs(rec_path)
            present = wrong = notmod = nochain = 0
            for cid, mlist in muts.items():
                spec = mut_seqs[cid]
                # best receptor chain for this spec chain
                best, best_r = None, 0.0
                for rname, rseq in rec.items():
                    r = SequenceMatcher(None, spec, rseq, autojunk=False).ratio()
                    if r > best_r:
                        best, best_r = rname, r
                if best is None or best_r < 0.50:
                    nochain += len(mlist)
                    continue
                m = spec_to_rec_map(spec, rec[best])
                rseq = rec[best]
                for (i, a, b) in mlist:
                    j = m.get(i)
                    if j is None:
                        notmod += 1
                    elif rseq[j] == b:
                        present += 1
                    else:
                        wrong += 1
            lost = wrong + notmod + nochain
            verdict = ("OK" if lost == 0 else
                       "NOCHAIN" if nochain == n_mut else
                       "ABSENT" if present == 0 else "PARTIAL")
            rows.append(dict(pdbid=sysid, variant=var, n_mut=n_mut, present=present,
                             wrong_aa=wrong, not_modeled=notmod, no_chain=nochain,
                             verdict=verdict))

    out = "analysis/casf_mutagenesis/mutation_presence_audit.json"
    json.dump(rows, open(out, "w"), indent=2)

    tot = len(rows)
    by = defaultdict(int)
    for r in rows:
        by[r["verdict"]] += 1
    print(f"cells audited: {tot}")
    for k in ("OK", "PARTIAL", "ABSENT", "NOCHAIN", "NOMUT"):
        if by[k]:
            print(f"  {k:8s} {by[k]:4d}  ({100*by[k]/tot:.1f}%)")
    tm = sum(r["n_mut"] for r in rows)
    tp = sum(r["present"] for r in rows)
    print(f"\nmutations specified {tm} | present {tp} ({100*tp/max(tm,1):.1f}%) | "
          f"lost {tm-tp}  [no_chain {sum(r['no_chain'] for r in rows)}, "
          f"not_modeled {sum(r['not_modeled'] for r in rows)}, "
          f"wrong_aa {sum(r['wrong_aa'] for r in rows)}]")

    dead = sorted({r["pdbid"] for r in rows if r["verdict"] in ("ABSENT", "NOCHAIN")})
    part = sorted({r["pdbid"] for r in rows if r["verdict"] == "PARTIAL"})
    print(f"\nsystems with >=1 cell where NO mutation survived ({len(dead)}):")
    print("  " + " ".join(dead))
    print(f"\nsystems with PARTIAL loss ({len(part)}):")
    print("  " + " ".join(part))
    print(f"\nCLEAN systems (all 3 variants OK): "
          f"{len({r['pdbid'] for r in rows} - set(dead) - set(part))}")


if __name__ == "__main__":
    main()
