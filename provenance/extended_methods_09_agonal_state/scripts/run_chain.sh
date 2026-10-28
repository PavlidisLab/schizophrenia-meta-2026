#!/bin/bash
# run_chain.sh NAME SETNAME : the correlation-corrected chain for Extended Fig. 8 a, b, c on one tree (one worker; threads = 1; nice; load gate 28; checkpointed)
NAME=$1; SET=$2
R=$HOME/odc375_2026-10-15
T=$R/trees/$NAME/repo/analyses
W=$R/work/$NAME; mkdir -p $W/od $W/scratch $W/tmp $W/home $R/inputs/ogdata
[ -e $R/inputs/ogdata/GSE12649.txt.gz ] || ln -s Ogdata/expression/GSE12649.txt.gz $R/inputs/ogdata/GSE12649.txt.gz
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
export HOME=$W/home TMPDIR=$W/tmp
export OD_SET=study_brain OD_SETNAME=$SET ITEM_E_SETNAME=$SET OD_OUT=$W/od OD_RESULTS=$W/od
export ITEM_E_SCRATCH=$W/scratch ITEM_E_CACHE=$R/inputs/cache_asd ITEM_E_OGDATA=$R/inputs/ogdata
export ITEM_E_PANEL_A=$W/od/enrichment_209_by_disorder.tsv ITEM_E_PANEL_B=$W/od/direction_209_by_disorder.tsv
PY=$R/venv/bin/python
IE=$T/stat012_vif_2026-10-12/addendum4_2026-10-12/axis_prg_proteomics_rodent/item_e_de_restricted
A3=$T/stat012_vif_2026-10-12/addendum3_2026-10-12/other_disorders_and_scan
cd $W
gate() { while [ "$(cut -d' ' -f1 /proc/loadavg | cut -d. -f1)" -ge 28 ]; do echo "$(date +%T) load gate" >> $W/gate.log; sleep 120; done; }
step() { # step NAME cmd...
  local n=$1; shift
  [ -e $W/done_$n ] && { echo "skip $n"; return 0; }
  gate; echo "$(date +%F_%T) START $n" >> $W/steps.log; local t0=$(date +%s)
  nice -n 19 "$@" > $W/log_$n.txt 2>&1 || { echo "$(date +%F_%T) FAIL $n" >> $W/steps.log; exit 1; }
  echo "$(date +%F_%T) DONE $n $(( $(date +%s) - t0 ))s" >> $W/steps.log; touch $W/done_$n
}
step od      $PY $T/other_disorders_176_2026-10-09/other_disorders_176.py
step tables  $PY $IE/item_e_tables.py
# the stored selfcontained table has the GSE12649 rows from a separate invocation (ONLY=GSE12649: seed index 0 and 1); the main invocation holds the other 44 contrasts
mkdir -p $W/emptyog
step perms   env ITEM_E_OGDATA=$W/emptyog $PY $IE/item_e_series_permutations.py
step perms_og env ONLY=GSE12649 $PY $IE/item_e_series_permutations.py
step rows    $PY $IE/item_e_rows.py
step cluster $PY $IE/item_e_cluster_jackknife.py
step dir     $PY $IE/item_e_direction.py
step assemble $PY $IE/item_e_assemble.py
step axisz   $PY $A3/axis_z_vif.py
echo "$(date +%F_%T) ALL DONE" >> $W/steps.log
