#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
PY=rerun-envs/py/bin/python
B=numbers_batch_B_2026-10-15/item10; D=core80_work_2026-10-16/P4/L09/item10; cd $D
t0=$(date +%s)
case $1 in
 v1) nice $PY scripts_edit.py snap_sz/pb $B/donors_old $D/axis_v1 $D/run_v1 1000 > ../logs_i10_v1.log 2>&1;;
 v2) nice $PY scripts_edit.py snap_sz/pb $B/donors_old $D/axis_v2 $D/run_v2 1000 > ../logs_i10_v2.log 2>&1;;
 val) nice $PY scripts_edit.py snap_sz/pb $B/donors_old $B/axis_v1 $D/val_v1 1000 > ../logs_i10_val.log 2>&1;;
esac
echo -e "L09\ti10_$1\t$(( $(date +%s)-t0 ))\texit $?" >> ../walls.tsv
