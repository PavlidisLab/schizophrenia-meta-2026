"""S3b: donor bootstrap and split-half reliability. Usage: python s3b_boot.py <task> <variant> [B]
 task boot  : B donor bootstraps (controls model; also cases / both for v2_bulk): headline statistics, between-collection reliability, leave-collection-out r
 task halves: B random donor halves (controls), per-dataset halves combined by Stouffer; Spearman between halves, Spearman-Brown
 task cc    : B donor bootstraps of the control-versus-case agreement (per dataset and pooled)
 variants: v1repro | v2_bulk | v2_bulk_ofc | v2_bulk_ofc_lcm | v2_bulk_v1first"""
import sys, os, datetime, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from hl import *
task, var = sys.argv[1], sys.argv[2]; BREP = int(sys.argv[3]) if len(sys.argv) > 3 else 200
OUT = f'{B}/out/tables'; os.makedirs(OUT, exist_ok=True)
VAR = dict(VARIANTS); VAR['v1repro'] = ('v1repro', {0})
mode, tiers = VAR[var]
reg = registry('primary'); fold = reg.fold.to_dict(); tier = reg.tier.to_dict()
fits = load_fits(mode)
DSL = V1DS if var == 'v1repro' else [d for d in fits['ctrl'] if tier[d] in tiers]
MD = f'{B}/out/{mode}/matrices'
def load(d):
    Y = pd.read_csv(f'{MD}/{d}.expr.tsv.gz', sep='\t', index_col=0); C = pd.read_csv(f'{MD}/{d}.cov.tsv', sep='\t')
    Y = Y[C.person.values] if set(C.person) <= set(Y.columns) else Y
    C['sexM'] = (C.sex == 'M').astype(float); return Y, C
DATA = {d: load(d) for d in set(DSL) | set(fits.get('case', {})) & set(DSL)}
master = pd.Index(sorted(set().union(*[set(DATA[d][0].index) for d in DSL]))); G = len(master)
POS = {d: master.get_indexer(DATA[d][0].index) for d in DSL}
VALS = {d: DATA[d][0].values.T.astype(float) for d in DSL}   # persons x genes
def design(C, terms):
    cols = []
    for t in terms:
        if t == 'sexM' and min((C.sexM == 0).sum(), (C.sexM == 1).sum()) < 2: continue
        cols.append(C[t].values.astype(float))
    return np.column_stack([np.ones(len(C))] + cols)
def fit_group(d, idx, terms, j=1):
    """OLS on the persons idx (indices into the dataset's person rows); returns (t, df, std_beta, v) over the dataset genes"""
    Y = VALS[d][idx]; C = DATA[d][1].iloc[idx]; X = design(C, terms)
    if np.linalg.matrix_rank(X) < X.shape[1]: return None
    n, p = X.shape; XtXi = np.linalg.inv(X.T @ X); Bt = XtXi @ X.T @ Y; R = Y - X @ Bt; s2 = (R ** 2).sum(0) / (n - p)
    se = np.sqrt(s2 * XtXi[j, j]); t = Bt[j] / np.where(se > 0, se, np.nan); sd = Y.std(0, ddof=1)
    return t, n - p, Bt[j] / np.where(sd > 0, sd, np.nan), (se / np.where(sd > 0, sd, np.nan)) ** 2, n
def to_master(d, v): 
    out = np.full(G, np.nan); out[POS[d]] = v; return out
MODEL_TERMS = {'ctrl': ['pH', 'age', 'sexM'], 'case': ['pH', 'age', 'sexM'], 'both': ['pH', 'age', 'sexM', 'sz']}
def groups(d, model):
    C = DATA[d][1]
    return np.where(C.sz.values == 0)[0] if model == 'ctrl' else np.where(C.sz.values == 1)[0] if model == 'case' else np.arange(len(C))
def axis_from_Z(Zm, N, dsl, excl=(), kfloor=4):
    keep = [i for i, d in enumerate(dsl) if fold[d] not in excl]
    z, k = stouffer(Zm[:, keep], N[keep], kmin_rule(len(keep), kfloor)); return z
def coll_between(Zm, N, dsl):
    cols = {}
    for f in sorted({fold[d] for d in dsl}):
        idx = [i for i, d in enumerate(dsl) if fold[d] == f]; z, _ = stouffer(Zm[:, idx], N[idx], 1); cols[f] = z
    fs = list(cols); rr = [sp(cols[a], cols[b]) for i, a in enumerate(fs) for b in fs[i + 1:]]
    return float(np.nanmean(rr)), len(fs)
SZ = {}
for f in os.listdir(f'{B}/ref'):
    if f.startswith('sz_adjust_') and not f.endswith('noshared.tsv'):
        x = pd.read_csv(f'{B}/ref/{f}', sep='\t'); SZ[f[10:-4]] = pd.Series(zfromtp(x.t_base.values, x.p_base.values), index=x.gene.values)
def one_replicate(model, dsl, rng, resample=True):
    Zm = np.full((G, len(dsl)), np.nan); N = np.zeros(len(dsl))
    for j, d in enumerate(dsl):
        g = groups(d, model)
        if model == 'both':
            c0 = np.where(DATA[d][1].sz.values == 0)[0]; c1 = np.where(DATA[d][1].sz.values == 1)[0]
            idx = np.concatenate([rng.choice(c0, len(c0)), rng.choice(c1, len(c1))]) if resample else g
        else: idx = rng.choice(g, len(g)) if resample else g
        r = fit_group(d, idx, MODEL_TERMS[model])
        if r is None: continue
        t, df, _, _, n = r; Zm[POS[d], j] = t2z(t, df); N[j] = n
    return Zm, N
t0 = time.time(); rng = np.random.default_rng(20261014 + abs(hash(var)) % 1000)
hdr = lambda t: f'# {t}; run {datetime.date.today()}; script scripts/s3b_boot.py {task} {var} {BREP}; seed 20261014 + variant offset; OLS t (unmoderated), df n-p; persons resampled with replacement within dataset (and within group for the both model)\n'
if task == 'boot':
    rows = []; models = ['ctrl'] + (['case', 'both'] if var != 'v1repro' else [])
    for model in models:
        dsl = [d for d in DSL if d in fits[model]]
        if len(dsl) < 3: continue
        for b in range(BREP + 1):
            Zm, N = one_replicate(model, dsl, rng, resample=(b > 0))      # replicate 0 = the OLS point estimate on the original persons
            axz = pd.Series(axis_from_Z(Zm, N, dsl), index=master).dropna(); st = stats_for_axis(axz); st.update(model=model, rep=b, n_genes_z=len(axz))
            rb, kc = coll_between(Zm, N, dsl); st.update(between_rbar=rb, between_SB=sb(rb, kc), between_SB4=sb(rb, 4), between_SB3=sb(rb, 3), n_collections=kc)
            # leave-collection-out r with each v1 collection's own SZ z
            for nm, S in SETS.items():
                rs_ = []
                for d0 in V1DS:
                    if d0 not in SZ or fold[d0] not in {fold[x] for x in dsl}: continue
                    xl = pd.Series(axis_from_Z(Zm, N, dsl, excl=(fold[d0],), kfloor=3), index=master); g = S.index.intersection(xl.dropna().index).intersection(SZ[d0].index)
                    rs_.append(pr(xl.loc[g].values, SZ[d0].loc[g].values))
                st[f'{nm}_loco_median_r'] = float(np.median(rs_)) if rs_ else np.nan
            rows.append(st)
    R = pd.DataFrame(rows); R.to_csv(f'{OUT}/boot_{var}.tsv', sep='\t', index=False, float_format='%.6g')
elif task == 'perm':
    # pH permuted among the controls of each dataset (the same permutation scheme as the stored variance-inflation check); null distribution of the number of set genes at q < 0.05 in the low-pH direction
    rows = []; dsl = [d for d in DSL if d in fits['ctrl']]
    for b in range(BREP):
        Zm = np.full((G, len(dsl)), np.nan); N = np.zeros(len(dsl))
        for j, d in enumerate(dsl):
            g = groups(d, 'ctrl'); C = DATA[d][1].iloc[g].copy(); C['pH'] = rng.permutation(C.pH.values); X = design(C, MODEL_TERMS['ctrl'])
            Y = VALS[d][g]; n, p = X.shape; XtXi = np.linalg.inv(X.T @ X); Bt = XtXi @ X.T @ Y; R_ = Y - X @ Bt; s2 = (R_ ** 2).sum(0) / (n - p); se = np.sqrt(s2 * XtXi[1, 1]); t = Bt[1] / np.where(se > 0, se, np.nan)
            Zm[POS[d], j] = t2z(t, n - p); N[j] = n
        axz = pd.Series(axis_from_Z(Zm, N, dsl), index=master).dropna(); st = stats_for_axis(axz); st['rep'] = b; rows.append(st)
    pd.DataFrame(rows).to_csv(f'{OUT}/perm_{var}.tsv', sep='\t', index=False, float_format='%.6g')
elif task == 'halves':
    rows = []; dsl = [d for d in DSL if d in fits['ctrl'] and len(groups(d, 'ctrl')) >= 12]
    for b in range(BREP):
        Z1 = np.full((G, len(dsl)), np.nan); Z2 = Z1.copy(); N1 = np.zeros(len(dsl)); N2 = N1.copy()
        for j, d in enumerate(dsl):
            g = rng.permutation(groups(d, 'ctrl')); h = len(g) // 2; a, c = g[:h], g[h:]
            for arr, Nn, idx in ((Z1, N1, a), (Z2, N2, c)):
                r = fit_group(d, idx, MODEL_TERMS['ctrl'])
                if r is None: continue
                arr[POS[d], j] = t2z(r[0], r[1]); Nn[j] = r[4]
        km = kmin_rule(len(dsl), 4); z1, _ = stouffer(Z1, N1, km); z2, _ = stouffer(Z2, N2, km); r = sp(z1, z2)
        rows.append(dict(rep=b, rho_halves=r, SB=2 * r / (1 + r), n_datasets=len(dsl), n_genes=int((~np.isnan(z1) & ~np.isnan(z2)).sum())))
    R = pd.DataFrame(rows); R.to_csv(f'{OUT}/halves_{var}.tsv', sep='\t', index=False, float_format='%.6g')
    # single-dataset halves (control pH t within one dataset with >= 30 controls), as the stored panel
    rows = []
    for d in [x for x in dsl if len(groups(x, 'ctrl')) >= 30]:
        for b in range(BREP):
            g = rng.permutation(groups(d, 'ctrl')); h = len(g) // 2; ra = fit_group(d, g[:h], MODEL_TERMS['ctrl']); rc = fit_group(d, g[h:], MODEL_TERMS['ctrl'])
            r = sp(ra[0], rc[0]); rows.append(dict(dataset=d, rep=b, rho_halves=r, SB=2 * r / (1 + r), n_ctrl=len(g)))
    pd.DataFrame(rows).to_csv(f'{OUT}/halves_single_{var}.tsv', sep='\t', index=False, float_format='%.6g')
elif task == 'cc':
    dsl = [d for d in DSL if d in fits['ctrl'] and d in fits.get('case', {})]
    sets = {'all genes': None, '5plus': set(five.index), '209': set(nine.index)}
    def cc_once(resample):
        per = {}; pooled = {'c': [], 'k': []}
        for d in dsl:
            gc, gk = groups(d, 'ctrl'), groups(d, 'case')
            ic = rng.choice(gc, len(gc)) if resample else gc; ik = rng.choice(gk, len(gk)) if resample else gk
            both_idx = np.concatenate([ic, ik]); sd = VALS[d][both_idx].std(0, ddof=1)
            rc = fit_group(d, ic, MODEL_TERMS['ctrl']); rk = fit_group(d, ik, MODEL_TERMS['case'])
            if rc is None or rk is None: continue
            bc = rc[2] * (rc[2] * 0 + 1); 
            # standardize both groups by the SD over both groups of the dataset
            Yc = VALS[d][ic]; Yk = VALS[d][ik]
            def beta(Y, idx):
                C = DATA[d][1].iloc[idx]; X = design(C, MODEL_TERMS['ctrl']); XtXi = np.linalg.inv(X.T @ X); Bt = XtXi @ X.T @ Y; R_ = Y - X @ Bt; s2 = (R_ ** 2).sum(0) / (len(C) - X.shape[1]); return Bt[1], s2 * XtXi[1, 1]
            b_c, v_c = beta(Yc, ic); b_k, v_k = beta(Yk, ik); sdd = np.where(sd > 0, sd, np.nan)
            sc, vc = b_c / sdd, v_c / sdd ** 2; sk, vk = b_k / sdd, v_k / sdd ** 2
            tc, tk = rc[0], rk[0]; genes = master[POS[d]]
            per[d] = dict(sc=pd.Series(sc, index=genes), sk=pd.Series(sk, index=genes), vc=pd.Series(vc, index=genes), vk=pd.Series(vk, index=genes), rho_t=sp(tc, tk))
        return per
    def summarize(per):
        out = {}
        for d, x in per.items():
            ok = x['sc'].notna() & x['sk'].notna(); a, b = x['sc'][ok].values, x['sk'][ok].values
            lam = max(1 - x['vc'][ok].mean() / a.var(), 1e-3)
            out[(d, 'all genes')] = dict(rho_t=x['rho_t'], rho_slope=sp(a, b), ratio_raw=float((a * b).sum() / (a * a).sum()), ratio_corrected=float((a * b).sum() / (a * a).sum() / lam), reliability_ctrl=lam, n_genes=int(ok.sum()))
        # pooled: inverse-variance (fixed-effect) mean of the standardized slope per gene, separately in controls and cases
        for sname, ss in sets.items():
            num_c = num_k = den_c = den_k = 0
            gl = master
            C_ = pd.DataFrame({d: x['sc'] for d, x in per.items()}).reindex(gl); K_ = pd.DataFrame({d: x['sk'] for d, x in per.items()}).reindex(gl)
            VC = pd.DataFrame({d: x['vc'] for d, x in per.items()}).reindex(gl); VK = pd.DataFrame({d: x['vk'] for d, x in per.items()}).reindex(gl)
            okm = C_.notna() & K_.notna() & (VC > 0) & (VK > 0); wc = (1 / VC).where(okm, 0); wk = (1 / VK).where(okm, 0)
            pc = (C_.fillna(0) * wc).sum(1) / wc.sum(1); pk = (K_.fillna(0) * wk).sum(1) / wk.sum(1); pvc = 1 / wc.sum(1)
            ok = (okm.sum(1) >= max(2, int(0.5 * len(per)))) & pc.notna() & pk.notna()
            if ss is not None: ok &= pd.Series(master.isin(list(ss)), index=master)
            a, b = pc[ok].values, pk[ok].values
            if len(a) > 20:
                lam = max(1 - pvc[ok].mean() / a.var(), 1e-3); out[('pooled', sname)] = dict(rho_t=np.nan, rho_slope=sp(a, b), ratio_raw=float((a * b).sum() / (a * a).sum()), ratio_corrected=float((a * b).sum() / (a * a).sum() / lam), reliability_ctrl=lam, n_genes=int(ok.sum()))
        return out
    rows = []
    for b in range(BREP + 1):
        o = summarize(cc_once(resample=(b > 0)))
        for (d, s), v in o.items(): rows.append(dict(rep=b, unit=d, geneset=s, **v))
    pd.DataFrame(rows).to_csv(f'{OUT}/cc_boot_{var}.tsv', sep='\t', index=False, float_format='%.6g')
print(task, var, 'done', round(time.time() - t0), 's')
