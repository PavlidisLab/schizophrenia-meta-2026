#!/usr/bin/env python3
"""Extended Fig. 1f on pH axis v2: panel_c.tsv (pH-axis rows replaced by v2 between-collection reliability) and split_half_ph_axis.tsv (v2 one-collection split halves). Source: ph_axis_v2_2026-10-14/reliability_summary.tsv (row 'v2 bulk (primary assignment)'), reliability_single_dataset_halves.tsv (same axis label). Non-axis rows of panel_c.tsv are copied unchanged from the manuscript repo file."""
import pandas as pd, hashlib, os, re
H=os.path.expanduser('~/analysis-delivery/ph_axis_v2_2026-10-14'); M=os.path.expanduser('~/Documents/manuscripts/schizophrenia-transcriptome-reproducibility-manuscript/figures/data/figures_data_quality')
sha=lambda p:hashlib.sha256(open(p,'rb').read()).hexdigest()
AX='v2 bulk (primary assignment)'
S=pd.read_csv(f'{H}/reliability_summary.tsv',sep='\t',comment='#'); S=S[S.axis==AX].iloc[0]
Hh=pd.read_csv(f'{H}/reliability_single_dataset_halves.tsv',sep='\t',comment='#'); Hh=Hh[Hh.axis==AX]
Bp=pd.read_csv(f'{H}/reliability_between_collections_point.tsv',sep='\t',comment='#'); Bp=Bp[Bp.axis=='v2_bulk'].iloc[0]
raw=open(f'{M}/panel_c.tsv').read().split('\n'); hdr=[l for l in raw if l.startswith('#')]; body=[l for l in raw if not l.startswith('#') and l.strip()]
cols=body[0].split('\t'); rows=[l.split('\t') for l in body[1:]]
new=[]
for r in rows:
    m=r[0]
    if m.startswith('Control-pH axis, 5 collections'): r=[m,f'{S.between_SB_point:.4f}',f'{S.between_SB_lo:.4f}',f'{S.between_SB_hi:.4f}','pH axis v2; Spearman-Brown of the mean pairwise Spearman between the 5 collections (Stanley, Pitt, BA10, Pritzker, NSW), one Stouffer z per collection; lo/hi = point +- 1.96 bootstrap SE (200 donor bootstraps, Fisher z scale); ph_axis_v2_2026-10-14/reliability_summary.tsv row "'+AX+'"']
    elif m.startswith('4 collections combined'): r=[m,f'{S.between_SB4_point:.4f}',f'{S.between_SB4_lo:.4f}',f'{S.between_SB4_hi:.4f}','pH axis v2; Spearman-Brown to 4 collections of the same mean pairwise Spearman; interval as above']
    elif m.startswith('3 collections combined'): r=[m,f'{S.between_SB3_point:.4f}',f'{S.between_SB3_lo:.4f}',f'{S.between_SB3_hi:.4f}','pH axis v2; Spearman-Brown to 3 collections of the same mean pairwise Spearman; interval as above']
    elif m.startswith('1 collection alone'): r=[m,f'{S.between_mean_pairwise_rho_point:.4f}',f'{S.between_rho_lo:.4f}',f'{S.between_rho_hi:.4f}','pH axis v2; mean pairwise Spearman between the five collections; interval as above']
    new.append(r)
stamp=['# STATUS 2026-10-18: PH_AXIS_V2 - the four pH-axis rows (5, 4, 3 collections combined; 1 collection alone) hold the v2 between-collection values (v1 held 0.5457, 0.4901, 0.4189, 0.1937); all other rows are the rows of the stored file (GTEx rows and the control pH-effect rows on the independent v1 five datasets are unchanged: the latter measure the pH slope reliability of the five independent datasets, not the axis)',
 f'# input: ph_axis_v2_2026-10-14/reliability_summary.tsv sha256 {sha(H+"/reliability_summary.tsv")[:12]}; reliability_between_collections_point.tsv sha256 {sha(H+"/reliability_between_collections_point.tsv")[:12]}; pH axis v2 = control z over ten dataset-regions (ph_axis_v2_2026-10-14/axis_v2_controls.tsv.gz)',
 f'# input: figures/data/figures_data_quality/panel_c.tsv of the manuscript repo sha256 {sha(M+"/panel_c.tsv")[:12]} (non-axis rows copied)']
with open('figures_data_quality/panel_c.tsv','w') as f:
    f.write('\n'.join(stamp+hdr)+'\n'+'\t'.join(cols)+'\n'+'\n'.join('\t'.join(r) for r in new)+'\n')
rows=[]
for _,r in Hh.sort_values('dataset').iterrows():
    n=int(r.n_ctrl); h1=n//2
    rows.append([f'Control-pH axis, two halves of the {r.dataset} controls (n={n})',f'{r.SB_mean}',f'{r.SB_lo}',f'{r.SB_hi}',f'pH axis v2 within one collection ({r.dataset}, {n} control donors with pH, age, sex); per-half axis = that collection\'s own per-gene control pH t (limma, expr ~ pH + age + sex), no Stouffer across collections; Spearman-Brown reliability from 200 random donor splits (mean; interval as stored in reliability_single_dataset_halves.tsv: 2.5% to 97.5% of splits); raw split-half rho {r.raw_rho_mean}'])
st=['# STATUS 2026-10-18: PH_AXIS_V2 - v2 controls (persons counted once across datasets, so n differs from the v1 file: GSE35977 48 against 50, GSE92538 GPL96 33 against 45); datasets with at least 30 controls',
    f'# input: ph_axis_v2_2026-10-14/reliability_single_dataset_halves.tsv sha256 {sha(H+"/reliability_single_dataset_halves.tsv")[:12]}, rows axis "{AX}"; script scripts/s5_compile.py of that folder (from s3b_boot.py)',
    f'# replaces figures/data/figures_data_quality/split_half_ph_axis.tsv sha256 {sha(M+"/split_half_ph_axis.tsv")[:12]} (v1: GSE35977 0.73, GSE92538 GPL96 0.55)']
with open('figures_data_quality/split_half_ph_axis.tsv','w') as f:
    f.write('\n'.join(st)+'\nmeasure\treliability\tlo\thi\twhat\n'+'\n'.join('\t'.join(map(str,r)) for r in rows)+'\n')
print(open('figures_data_quality/panel_c.tsv').read().split('measure\t')[1][:2500]); print(open('figures_data_quality/split_half_ph_axis.tsv').read()[-1800:])
