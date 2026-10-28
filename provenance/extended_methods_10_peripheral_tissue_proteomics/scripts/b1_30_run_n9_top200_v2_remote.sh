#!/bin/bash
# B1 item 'Blood and patient-derived neural cells' (Ext. Methods 10): rerun N9 (blood / iPSC lists vs the pH top-200 sets) with the top 200 built from pH axis v2.
# Copy of core80_work_2026-10-16/P4/L14/run_N9.sh (remote: core80_work_2026-10-16/P4/L14/run_N9.sh); only PH and OUT dirs differ. date 2026-10-19; no random draws added here (N9 is exact Poisson-binomial; seed n/a).
# validation: tag v1 uses the stored v1 table (combined_controls_pH.tsv sha256 ebab86afedd3) and must reproduce core80 N9 out_new/results/set_overlap.tsv.
# base: core-80 (72 brain lists, 25 blood, 18 iPSC-neural modelled; Shao 2008 out via RECURRENCE_EXCLUDE=18191109).
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 RECURRENCE_NK=membership
PY=rerun-envs/py/bin/python; REPO=sztr-rebuild-infra/repo; P=core80_work_2026-10-16/P4; D=$P/L14/N9; W=factcheck_B_2026-10-19; mkdir -p $W; cd $W
PHV1=$REPO/analyses/confounds/results/combined_controls_pH.tsv
PHV2=fdv2_2026-10-15/work/axis_v2_tables/combined_controls_pH.tsv
one(){ tag=$1; PH=$2; mkdir -p $W/out_$tag/results
 SETSFILE=$P/NEW/analyses/confounds/recurrent_sets.tsv PH=$PH BASE=$D/in/periphery_209_sets.tsv SLOT133=study_brain SLOT5=study_brain_ge5 EXTRA=study_brain_ge10 TOPPH=1 OUT=$W/out_$tag/sets.tsv nice $PY $D/mod/n9_build_sets.py > $W/out_$tag/build.log 2>&1
 env RECURRENCE_EXCLUDE=18191109 SZTR_PARAMS=base78_2026-10-16/B80/params.tsv OUTDIR=$W/out_$tag ANALYSES=$P/NEW/analyses SETSTSV=$W/out_$tag/sets.tsv nice $PY $D/mod/overlap_375.py > $W/out_$tag/run.log 2>&1; echo "$tag exit $?"; }
one v1 $PHV1
one v2 $PHV2
sha256sum $PHV1 $PHV2 > $W/input_sha256.txt
