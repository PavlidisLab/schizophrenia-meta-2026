#!/bin/bash
# A2 item 'The 383 recurrent genes ... 11 genes remain without Voisey 2017, 6 without Xiao 2014'.
# Run on remote 2026-10-19. Runs the core-80 recurrence code (recurrence_fork.py, unchanged copy of
# base78_2026-10-16/B80/ws) three times with the same settings as base78 run.sh (NK=membership,
# EMPTY_LISTS=drop, UNIVERSE_LEVEL=current, include file = 80 primary PMIDs):
#   none (validation: must reproduce stored recurrence_genes.tsv), Voisey 2017 (27940260), Xiao 2014 (24776767) left out.
# No randomness in recurrence_fork.py (exact Poisson-binomial). Base: core-80 (80 studies, 72 lists; Shao out).
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
R=$HOME/factcheck_A_2026-10-19
[ -d $R/ws ] || cp -a $HOME/base78_2026-10-16/B80/ws $R/ws
cd $R/ws/analyses
for tag in none:"" voisey:27940260 xiao:24776767; do
  name=${tag%%:*}; pm=${tag#*:}; OUT=$R/out_$name; mkdir -p $OUT
  RECURRENCE_OUT=$OUT RECURRENCE_NK=membership RECURRENCE_INCLUDE_FILE=$R/ws/pipeline_inputs/primary_pmids.txt \
    RECURRENCE_UNIVERSE_LEVEL=current RECURRENCE_EXCLUDE=$pm RECURRENCE_EMPTY_LISTS=drop \
    nice rerun-envs/py/bin/python recurrence_fork.py > $OUT/recurrence.log 2>&1 &
done
wait
echo ALLDONE > $R/done.txt
