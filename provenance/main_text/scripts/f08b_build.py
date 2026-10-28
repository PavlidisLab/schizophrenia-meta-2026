import sys, os, shutil, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ['RECURRENCE_NK'] = 'membership'; os.environ['EXPORT_RUN_DATE'] = '2026-10-15'
from fdlib import *
W = f'{RUN}/work'; D = f'{RUN}/in/fd375/fig08'
sys.path.insert(0, f'{AN}/figure_exports'); sys.path.insert(0, AN)
import export_fig08_fig09_fig10 as E
def stage(tag):
    A = f'{W}/axroot_{tag}/analyses/axes_pritzker_agonal_2026-10-09'; C = f'{W}/axroot_{tag}/analyses/claims_vs_ph_axis_2026-10-08'
    S = f'{W}/b_{tag}'; shutil.rmtree(S, ignore_errors=True); os.makedirs(f'{S}/axes/results_remote'); os.makedirs(f'{S}/claims'); os.makedirs(f'{S}/out')
    for f in ('alignment_by_dataset.tsv', 'pairs_partial.tsv', 'summary.txt'): shutil.copy(f'{A}/out_{tag}/{f}', f'{S}/axes/results_remote/{f}')
    shutil.copy(f'{A}/single_axis_same_genes_remote.tsv', f'{S}/axes/single_axis_same_genes_remote.tsv')
    for f in ('alignment_two_axes.tsv', 'alignment_by_dataset.tsv', 'pairs_two_axes.tsv'): shutil.copy(f'{C}/{f}', f'{S}/claims/{f}')
    E.P_AX_REMOTE = f'{S}/axes'; E.P_CLAIMS_AX = f'{S}/claims'; E.OUT = f'{S}/out'
    F = E.Figure('fig08'); E.fig8_b(F)
    return S
for tag in ('v1', 'v2'): stage(tag)
FILES = {'fig08_b_dataset_axis_alignment.tsv': (['dataset'], None), 'fig08_b_pair_agreement_gtex_axis.tsv': (['axis_version', 'pair_class', 'level', 'axis_removed'], None),
         'fig08_b_pair_agreement_single_axis_same_genes.tsv': (['axis_version', 'set', 'axis_removed', 'level'], None), 'fig08_b_pairs_remote.tsv': (['dataset_a', 'dataset_b'], None),
         'fig08_b_axis_correlations.tsv': (['pair_of_axes'], None), 'fig08_b_alignment_medians_by_target_set.tsv': (['target_set', 'axis', 'level'], None),
         'fig08_b_pair_medians_by_target_set.tsv': (['target_set', 'pair_class', 'level', 'axes_removed'], None)}
def rdb(p): return pd.read_csv(p, sep='\t', comment='#', dtype=str, keep_default_na=False)
def cmts(p): return [l.rstrip('\n') for l in open(p) if l.startswith('#')]
VAL = []
AXH = '# axis: unsuffixed values = pH axis v2 (primary; inputs_ph_axis_v2/axis_v2_controls.tsv.gz sha256 %s; leave-collection-out construction of the decomposition scripts with the v2 axis from the stored per-dataset limma fits, scripts v2axis.py); columns ending _axis_v1 = v1 fork (five-collection control pH axis; stored 375 values reproduced by unchanged code)' % sha(P_V2CTL)
for nm, (key, _) in FILES.items():
    st = rdb(f'{D}/{nm}'); v1 = rdb(f'{W}/b_v1/out/fig08/{nm}'); v2 = rdb(f'{W}/b_v2/out/fig08/{nm}')
    same = st.shape == v1.shape and bool((st.values == v1.values).all())
    if not same:   # tolerant numeric comparison
        try:
            num_eq = all((a == b) or (a != '' and b != '' and abs(float(a) - float(b)) < 1e-9) for a, b in zip(st.values.ravel(), v1.values.ravel())) and st.shape == v1.shape
        except ValueError: num_eq = False
        same = num_eq
    VAL.append((f'fig08/{nm}', 'unchanged export_fig08_fig09_fig10.py fig8_b on the v1 decomposition outputs (rerun of the unchanged scripts; sha of outputs equal to the stored ones) reproduces the stored 375 file: shape; all cells equal', f'{v1.shape[0]}x{v1.shape[1]}; {same}', f'figure_data_375_2026-10-15/fig08/{nm} {st.shape[0]}x{st.shape[1]}', 'match' if same else 'DIFF', ''))
    vals = [c for c in v2.columns if c not in key]
    # keep only value columns that can differ; add fork columns for all non-key, non-descriptive ones
    desc = {'axis_collections', 'in_ALL', 'in_NOPRITZ', 'in_NODAI', 'set_description', 'n_pairs', 'n_datasets', 'n_expected_sign', 'n_of', 'both_bulk'}
    fork = [c for c in vals if c not in desc]
    m = v2.merge(v1[key + fork].rename(columns={c: c + '_axis_v1' for c in fork}), on=key, how='left')
    nch = {c: int((m[c] != m[c + '_axis_v1']).sum()) for c in fork}
    cm = cmts(f'{W}/b_v2/out/fig08/{nm}')
    p = f'{OUT}/fig08/{nm}'; os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, 'w') as fh:
        fh.write('# STATUS 2026-10-15: REGENERATED_AXIS_V2 - pH-axis quantities on pH axis v2 (primary), v1 as *_axis_v1 fork columns\n')
        for l in cm:
            if l.startswith('# script:'): l = '# script: fdv2_2026-10-15/work/f08b_build.py running analyses/figure_exports/export_fig08_fig09_fig10.py fig8_b unchanged on decomposition outputs of axes_pritzker_agonal_2026-10-09/decomposition_axes.py, single_axis_same_genes.py and claims_vs_ph_axis_2026-10-08/decomposition.py (copies; only the pH-axis function edited), run 2026-10-15'
            fh.write(l + '\n')
        fh.write(AXH + '\n')
        m.to_csv(fh, sep='\t', index=False, lineterminator='\n')
    print(nm, m.shape, {k: v for k, v in nch.items() if v})
pd.DataFrame(VAL, columns=['file','check','value_recomputed','reference','result','note']).to_csv(f'{OUT}/VALIDATION_fig08b.tsv', sep='\t', index=False)
print(pd.DataFrame(VAL)[[0, 2, 4]].to_string())
