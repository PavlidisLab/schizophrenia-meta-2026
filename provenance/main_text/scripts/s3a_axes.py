"""S3a: v2 axes (controls, cases, both), leave-collection-out axes, and point-estimate comparisons with v1. Run on remote: rerun-envs/vpy/bin/python scripts/s3a_axes.py"""
import sys, os, json, datetime
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from phlib import *
from scipy.stats import norm
OUT = f'{B}/out/tables'; os.makedirs(OUT, exist_ok=True)
fits = {m: load_fits(m) for m in ['primary', 'v1first', 'v1repro']}
reg = registry('primary'); reg1 = registry('v1first'); fold = reg.fold.to_dict(); tier = reg.tier.to_dict()
REF = f'{B}/ref'
v1c = pd.read_csv(f'{REF}/combined_controls_pH.tsv', sep='\t', comment='#'); v1c = v1c[v1c.factor == 'pH'].set_index('gene')
rg = pd.read_csv(f'{REF}/recurrent_genes_pH.tsv', sep='\t')
agon = pd.read_csv(f'{B}/../agonal_ph_axis_tool/input/pooled_axes.tsv.gz', sep='\t', index_col=0)['w_agonal'].dropna()
def variant_ds(var, model):
    mode, tiers = VARIANTS[var]; return [d for d in fits[mode][model] if reg.loc[d, 'tier'] in tiers]
def build_axis(fm, dsl, kfloor=4):
    genes, T, Z, Y, V, N = stack(fm, dsl); K = len(dsl); km = kmin_rule(K, kfloor)
    z, k = stouffer(Z, N, km); mu, se, tau2, I2, Q = dl(Y, V, km)
    W = np.where(~np.isnan(Z), N[None, :], 0).sum(1)
    df = pd.DataFrame({'z_stouffer': z, 'n_datasets': k, 'n_persons': W, 'beta_re_sdunits': mu, 'se_re': se, 'z_re': mu / se, 'tau2': tau2, 'I2': I2}, index=genes)
    df['p'] = 2 * norm.sf(np.abs(df.z_stouffer)); df['q_bh'] = bh(df.p.values)
    return df.dropna(subset=['z_stouffer']), km
def loco(fm, dsl, kfloor=3):
    folds = sorted({fold[d] for d in dsl}); res = {}
    genes, T, Z, Y, V, N = stack(fm, dsl)
    for f in folds:
        keep = [i for i, d in enumerate(dsl) if fold[d] != f]; km = kmin_rule(len(keep), kfloor)
        z, k = stouffer(Z[:, keep], N[keep], km); res[FOLDNAME[f]] = pd.Series(z, index=genes)
    return pd.DataFrame(res)
hdr = lambda t: f'# {t}; run {datetime.date.today()}; script scripts/s3a_axes.py; inputs: ph_axis_v2/out/<mode>/*.tsv.gz (limma fits, scripts/s2_fit.R), ref/ (v1 tables from srep-mauney-integrated confounds64 out/tables_H)\n'
def write(df, name, title, index=True):
    with open(f'{OUT}/{name}', 'w') as fh:
        fh.write(hdr(title)); df.to_csv(fh, sep='\t', index=index, index_label=('gene' if index else None), float_format='%.5g')
info = {}
# ---------- v1 reproduction check
fm = fits['v1repro']['ctrl']; a, km = build_axis(fm, V1DS, 4)
j = a.join(v1c[['Z']], how='inner'); print('v1 reproduction: genes', len(j), 'corr', round(j.z_stouffer.corr(j.Z), 5), 'max|diff|', round((j.z_stouffer - j.Z).abs().max(), 3)); info['v1_repro_corr'] = float(j.z_stouffer.corr(j.Z))
# ---------- axes
AX = {}
for var in VARIANTS:
    mode = VARIANTS[var][0]
    for model, lab in (('ctrl', 'controls'), ('case', 'cases'), ('both', 'both')):
        dsl = variant_ds(var, model)
        if len(dsl) < 3: continue
        a, km = build_axis(fits[mode][model], dsl, 4); AX[(var, model)] = (a, dsl, km)
        cols = [(d, int(fits[mode][model][d].n.iloc[0])) for d in dsl]
        print(var, model, 'datasets', len(dsl), 'kmin', km, 'genes with z', len(a), 'persons', sum(c[1] for c in cols))
        if var == 'v2_bulk':
            write(a, f'axis_v2_{lab}.tsv', f'pH axis v2, {lab}; Stouffer z (weights sqrt(n)) over {len(dsl)} datasets [{",".join(dsl)}], genes in at least {km}; beta_re_sdunits = DerSimonian-Laird random-effects mean slope in SD of expression per pH unit'); 
        write(a, f'axis_{var}_{lab}.tsv.gz' if False else f'axis_{var}_{lab}.tsv', f'pH axis {var}, {lab}; datasets [{",".join(dsl)}]') if var != 'v2_bulk' else None
for var, model, lab in [('v2_bulk', 'ctrl', 'controls'), ('v2_bulk', 'both', 'both'), ('v2_bulk', 'case', 'cases')]:
    dsl = AX[(var, model)][1]; L = loco(fits[VARIANTS[var][0]][model], dsl)
    write(L, f'loco_v2_{lab}.tsv', f'pH axis v2 {lab}, leave-collection-out: each column is the Stouffer z without the datasets of that collection (fold); genes in at least 80 percent of the remaining datasets, minimum 3')
# ---------- datasets table
rows = []
for var in ['v2_bulk', 'v2_bulk_ofc_lcm']:
    for model in ['ctrl', 'case', 'both']:
        for d in variant_ds(var, model):
            x = fits[VARIANTS[var][0]][model][d]; rows.append(dict(variant=var, model=model, dataset=d, fold=FOLDNAME[fold[d]], tier=tier[d], n_persons=int(x.n.iloc[0]), n_genes=len(x)))
write(pd.DataFrame(rows), 'datasets_by_model.tsv', 'datasets entering each model; n_persons = persons in the model', index=False)
# ---------- v1 vs v2 gene-level
rows = []
v1z = v1c.Z
for (var, model), (a, dsl, km) in AX.items():
    z2 = a.z_stouffer; sh = z2.index.intersection(v1z.index)
    x, y = z2.loc[sh].values, v1z.loc[sh].values
    rng = np.random.default_rng(20261014); bs = []
    for _ in range(200):
        ii = rng.integers(0, len(sh), len(sh)); bs.append(sp(x[ii], y[ii]))
    rows.append(dict(variant=var, v2_model=model, n_datasets_v2=len(dsl), n_genes_v2=len(z2), n_genes_v1=len(v1z), n_shared=len(sh), spearman=sp(x, y), spearman_lo_gene_boot=np.percentile(bs, 2.5), spearman_hi_gene_boot=np.percentile(bs, 97.5), pearson=pr(x, y),
                     n_abs_z_gt_3_v1=int((v1z.abs() > 3).sum()), n_abs_z_gt_3_v2=int((z2.abs() > 3).sum()), n_q05_v1=int((v1c.q < 0.05).sum()), n_q05_v2=int((a.q_bh < 0.05).sum())))
C1 = pd.DataFrame(rows); write(C1, 'comparison_v1_v2_gene_level.tsv', 'v1 (stored combined_controls_pH.tsv) against v2 axes over shared genes; intervals bootstrap genes (genes treated as independent)', index=False)
print(C1[['variant','v2_model','n_datasets_v2','n_genes_v2','n_shared','spearman','spearman_lo_gene_boot','spearman_hi_gene_boot','n_q05_v1','n_q05_v2']].round(3).to_string())
# ---------- collection-level agreement (between-collection reliability, point estimates)
def coll_z(fm, dsl, kfloor=1):
    genes, T, Z, Y, V, N = stack(fm, dsl); out = {}
    for f in sorted({fold[d] for d in dsl}):
        idx = [i for i, d in enumerate(dsl) if fold[d] == f]; z, k = stouffer(Z[:, idx], N[idx], 1); out[FOLDNAME[f]] = pd.Series(z, index=genes)
    return pd.DataFrame(out)
def between(CZ):
    cols = list(CZ.columns); k = len(cols); rr = []; prs = {}
    for i in range(k):
        for j in range(i + 1, k):
            r = sp(CZ[cols[i]].values, CZ[cols[j]].values); rr.append(r); prs[(cols[i], cols[j])] = r
    return np.nanmean(rr), prs
rows = []; CZs = {}
CZ1 = pd.DataFrame({FOLDNAME[fold[d]]: pd.Series(zfromtp(fits['v1repro']['ctrl'][d].t.values, fits['v1repro']['ctrl'][d].p.values), index=fits['v1repro']['ctrl'][d].index) for d in V1DS}); CZs['v1 (five datasets, v1 persons)'] = CZ1
for var in VARIANTS:
    CZs[var] = coll_z(fits[VARIANTS[var][0]]['ctrl'], variant_ds(var, 'ctrl'))
for nm, CZ in CZs.items():
    rbar, prs = between(CZ); k = CZ.shape[1]
    rows.append(dict(axis=nm, n_collections=k, mean_pairwise_spearman=rbar, SB_all=sb(rbar, k), SB_4=sb(rbar, 4), SB_3=sb(rbar, 3), single=rbar, **{f'rho_{a}_{b}': v for (a, b), v in prs.items()}))
    # each collection against the Stouffer of the others
    for c in CZ.columns:
        pass
R1 = pd.DataFrame(rows); write(R1, 'reliability_between_collections_point.tsv', 'mean pairwise Spearman among collection-level control pH z (one Stouffer per collection), Spearman-Brown to k collections; point estimates (intervals: reliability_bootstrap.tsv)', index=False)
print(R1[['axis','n_collections','mean_pairwise_spearman','SB_all','SB_4','SB_3']].round(3).to_string())
json.dump(info, open(f'{OUT}/s3a_info.json', 'w'))
