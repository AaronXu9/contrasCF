#!/usr/bin/env bash
# Verify every factual claim in docs/data_store_map.md (2026-08-20 revision).
# Exit 1 if any claim fails. Prints OK/FAIL per claim with expected vs actual.
set -uo pipefail
FAIL=0
PY=/home/aoxu/miniconda3/envs/rdkit_env/bin/python

LAB_REPO=/mnt/katritch_lab2/aoxu/contrasCF
LAB_CASF=/home/aoxu/projects/VLS-Benchmark-Dataset/data/pdbbind_cleansplit
LAB_HIQBIND=/mnt/katritch_lab2/aoxu/data/hiqbind
LAB_CACHE=$LAB_REPO/.claude/worktrees/counterfold/counterfold/_cache
OUT=$LAB_REPO/analysis/casf_mutagenesis/outputs
export OUT   # the truncated-receptor check reads it from the environment
CARC=aoxu@discovery.usc.edu
CARC_REPO=/project2/katritch_223/aoxu/contrasCF
CARC_CASF=/project2/katritch_223/aoxu/projects/VLS-Benchmark-Dataset/data/pdbbind_cleansplit

eq() { # eq <label> <expected> <actual>
  if [ "$2" = "$3" ]; then printf "  OK    %-52s %s\n" "$1" "$3"
  else printf "  FAIL  %-52s expected=%s actual=%s\n" "$1" "$2" "$3"; FAIL=1; fi
}
tf() { # tf <label> <0-or-1 from test>
  if [ "$2" = "0" ]; then printf "  OK    %s\n" "$1"
  else printf "  FAIL  %s\n" "$1"; FAIL=1; fi
}

echo "################ KATLAB ################"
eq "raw/ dirs = 281 (251 HiQBind + 30 recovered)" 281 "$(ls -1 $LAB_CASF/raw 2>/dev/null | wc -l)"
eq "crystal_ligands entries=14662" 14662 "$(ls -1 $LAB_CASF/crystal_ligands 2>/dev/null | wc -l)"
eq "crystal_ligands *.sdf = 14661"  14661 "$(ls -1 $LAB_CASF/crystal_ligands/*_ligand.sdf 2>/dev/null | wc -l)"

# crystal_ligands: 0 symlinks / 14661 real  (full scan)
read -r CLSL CLRF < <($PY -c "
import os
d='$LAB_CASF/crystal_ligands'; sl=rf=0
for f in os.scandir(d):
    if f.name.endswith('_ligand.sdf'):
        sl+=1 if f.is_symlink() else 0; rf+=0 if f.is_symlink() else 1
print(sl,rf)")
eq "crystal_ligands symlinks = 0"   0     "$CLSL"
eq "crystal_ligands real files"     14661 "$CLRF"

# raw/<id>/ now carries TWO provenances (docs/data_prep_todo.md item 1):
#   251 HiQBind-backed -- <id>_protein.pdb + <id>_ligand.sdf, both symlinks into $HIQBIND
#    30 RCSB-recovered -- <id>_protein.pdb only, a REAL file, unrefined deposited
#                         coordinates written by 25_recover_missing_receptors.py
REC=$LAB_REPO/analysis/casf_mutagenesis/recovered_receptors.json
recids=$($PY -c "
import json
print(' '.join(r['pdbid'] for r in json.load(open('$REC'))['records']))" 2>/dev/null)
sl=0; res=0; tot=0; rec=0; recreal=0
for d in $LAB_CASF/raw/*/; do
  id=$(basename "$d")
  case " $recids " in
    *" $id "*)
      rec=$((rec+1))
      f="$d/${id}_protein.pdb"
      { [ -f "$f" ] && [ ! -L "$f" ]; } && recreal=$((recreal+1))
      continue ;;
  esac
  for f in "$d/${id}_protein.pdb" "$d/${id}_ligand.sdf"; do
    tot=$((tot+1)); [ -L "$f" ] && sl=$((sl+1)); [ -e "$f" ] && res=$((res+1))
  done
done
eq "raw/ HiQBind entries (251x2)"       502 "$tot"
eq "raw/ HiQBind entries symlinked"     502 "$sl"
eq "raw/ HiQBind symlinks RESOLVE"      502 "$res"
eq "raw/ RCSB-recovered dirs"            30 "$rec"
eq "raw/ RCSB-recovered are REAL files"  30 "$recreal"

[ -L "$LAB_REPO/data/casf2016" ]; tf "data/casf2016 is a symlink (lab)" $?

# docking cells
# NOTE: -mindepth/-maxdepth 4 restricts to outputs/<id>/<variant>/docking/, so
# backup trees such as outputs/_docking_prefix_backup/<id>/... are not counted.
eq "docking wt cells"   251 "$(find $OUT -mindepth 4 -maxdepth 4 -path '*/wt/docking/receptor.pdb' 2>/dev/null | wc -l)"
for v in rem pack inv; do
  eq "docking $v cells"  239 "$(find $OUT -mindepth 4 -maxdepth 4 -path "*/$v/docking/receptor.pdb" 2>/dev/null | wc -l)"
done
eq "docking cells total = 968" 968 "$(find $OUT -mindepth 4 -maxdepth 4 -path '*/docking/receptor.pdb' 2>/dev/null | wc -l)"
for e in gnina unidock2 surfdock; do
  eq "$e output dirs = 968" 968 "$(find $OUT -mindepth 3 -maxdepth 3 -type d -name $e 2>/dev/null | wc -l)"
done
eq "outputs/ system dirs = 251" 251 "$(ls -1 $OUT 2>/dev/null | grep -Ec '^[0-9a-z]{4}$')"

# truncated-receptor list
# Was 16 systems; 0169b42 (AF3+MSA all-chains) + the docking rebuild fixed 15.
# 2vw5 remains: the prediction drops 3 of its 4 homotetramer chains.
EXPECT="2vw5"
ACTUAL=$($PY - <<'EOF'
import glob,os
root=os.environ.get("OUT","")
def nres(p):
    try:
        return sum(1 for l in open(p) if l.startswith("ATOM") and l[12:16].strip()=="CA")
    except Exception: return None
bad=[]
for wt in glob.glob(f"{root}/*/wt/docking/receptor.pdb"):
    s=wt.split("/")[-4]; a=nres(wt)
    ms=[nres(f"{root}/{s}/{v}/docking/receptor.pdb") for v in ("rem","pack","inv")
        if os.path.exists(f"{root}/{s}/{v}/docking/receptor.pdb")]
    if a and ms and min(ms) < 0.5*a: bad.append(s)
print(" ".join(sorted(bad)))
EOF
)
eq "truncated (<50%) systems" "$EXPECT" "$ACTUAL"

# scope ladder + clusters typo
$PY - <<'EOF'
import json,os
R="/home/aoxu/projects/VLS-Benchmark-Dataset/data/pdbbind_cleansplit"
def ids(p):
    s=set()
    def fl(x):
        if isinstance(x,str) and len(x)==4: s.add(x.lower())
        elif isinstance(x,(list,tuple)):
            for y in x: fl(y)
        elif isinstance(x,dict):
            for k,v in x.items(): fl(k); fl(v)
    fl(json.load(open(p))); return s
orig=set(); clean=set()
for f in os.listdir(f"{R}/labels"):
    if "original_train_val_split" in f: orig|=ids(f"{R}/labels/{f}")
    if "cleansplit_train_val_split" in f: clean|=ids(f"{R}/labels/{f}")
d=json.load(open(f"{R}/labels/clusters_casf2016.json"))
core=[r[0].lower() for c in d.values() for r in c]
print(f"  {'OK' if len(orig)==18623 else 'FAIL'}    scope: original splits = {len(orig)} (expect 18623)")
print(f"  {'OK' if len(clean)==16491 else 'FAIL'}    scope: cleansplit splits = {len(clean)} (expect 16491)")
print(f"  {'OK' if len(core)==285 else 'FAIL'}    scope: CASF core = {len(core)} in {len(d)} clusters (expect 285/57)")
typo=sorted(set(core)&{'105b','10wh'})
print(f"  {'OK' if typo==['105b','10wh'] else 'FAIL'}    clusters json contains corrupt ids {typo}")
real=[i for i in ('1o5b','1owh') if os.path.isdir(f"{R}/raw/{i}")]
print(f"  {'OK' if real==['1o5b','1owh'] else 'FAIL'}    correct ids exist on disk: {real}")
EOF

for d in plip_labels plip flowr affinity; do
  [ -d "$LAB_CACHE/$d" ]; tf "LAB_CACHE/$d exists" $?
done

echo ""
echo "################ CARC ################"
CO=$(timeout 150 ssh -o BatchMode=yes -o ConnectTimeout=15 "$CARC" "
R=$CARC_REPO/analysis/casf_mutagenesis/outputs
C=$CARC_CASF
echo \"rawdirs \$(ls -1 \$C/raw 2>/dev/null | wc -l)\"
e=0; for d in \$C/raw/*/; do [ -z \"\$(ls -A \$d 2>/dev/null)\" ] && e=\$((e+1)); done
echo \"rawempty \$e\"
echo \"wtdock \$(find \$R -path '*/wt/docking/receptor.pdb' 2>/dev/null | wc -l)\"
echo \"mutdock \$(find \$R -path '*/rem/docking/receptor.pdb' -o -path '*/pack/docking/receptor.pdb' -o -path '*/inv/docking/receptor.pdb' 2>/dev/null | wc -l)\"
echo \"gninadirs \$(find \$R -maxdepth 3 -type d -name gnina 2>/dev/null | wc -l)\"
ne=0; for d in \$(find \$R -maxdepth 3 -type d -name gnina 2>/dev/null); do [ -n \"\$(ls -A \$d 2>/dev/null)\" ] && ne=\$((ne+1)); done
echo \"gninanonempty \$ne\"
echo \"af3cif \$(find \$R -name 'af3msa_*_rem_model_0.cif' 2>/dev/null | wc -l)\"
echo \"outents \$(ls -1 \$R 2>/dev/null | wc -l)\"
[ -f \$C/crystal_ligands/1bcu_ligand.sdf ] && [ ! -L \$C/crystal_ligands/1bcu_ligand.sdf ] && echo 'clreal 1' || echo 'clreal 0'
{ [ -d $CARC_REPO/data/casf2016 ] && [ ! -L $CARC_REPO/data/casf2016 ]; } && echo 'realdir 1' || echo 'realdir 0'
" 2>&1)
g() { printf '%s\n' "$CO" | awk -v k="$1" '$1==k{print $2; exit}'; }
eq "CARC raw/ dirs = 285"          285 "$(g rawdirs)"
eq "CARC raw/ EMPTY dirs = 34"      34 "$(g rawempty)"
eq "CARC docking wt cells = 251"   251 "$(g wtdock)"
eq "CARC docking mutant cells = 0"   0 "$(g mutdock)"
eq "CARC gnina dirs = 968"         968 "$(g gninadirs)"
eq "CARC gnina NON-EMPTY = 0"        0 "$(g gninanonempty)"
eq "CARC af3msa rem CIFs = 239"    239 "$(g af3cif)"
eq "CARC outputs/ entries = 273"   273 "$(g outents)"
eq "CARC crystal ligand is real file" 1 "$(g clreal)"
eq "CARC data/casf2016 is real dir"   1 "$(g realdir)"

echo ""
[ "$FAIL" = 0 ] && echo "RESULT: PASS — every claim in data_store_map.md verified." \
                || echo "RESULT: FAIL — see FAIL lines above."
exit $FAIL
