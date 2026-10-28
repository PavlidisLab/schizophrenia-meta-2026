"""fig07 axis-dependent files for pH axis v2 (primary) with v1 (fork). Validation mode (--v1check) reproduces the stored 375 files from v1 inputs."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fdlib import *
W = f'{RUN}/work'; os.makedirs(f'{W}/axis_v2_tables', exist_ok=True)
c1 = load_v1(); c2 = load_v2(); S = load_sets(); RG = load_rg()
# ---- v2 table in the layout of combined_controls_pH.tsv (pH rows)
t2 = pd.DataFrame({'gene': c2.index, 'factor': 'pH', 'n_datasets': c2.n_datasets.values, 'Z': c2.z_stouffer.values,
                   'datasets': 'axis_v2:10_dataset_regions', 'expr_percentile': c1.expr_percentile.reindex(c2.index).values,
                   'p': c2.p.values, 'q': c2.q_bh.values})
def check_q():
    return float(np.nanmax(np.abs(bh(t2.p.values) - t2.q.values)))
def rfmt(x, sig=None):
    """R fwrite-like formatting of a float"""
    if x is None or (isinstance(x, float) and math.isnan(x)): return ''
    x = float(x)
    if sig: x = float(f'{x:.{sig}g}')
    fx = repr(x); 
    if fx.endswith('.0'): fx = fx[:-2]
    m, e = f'{x:e}'.split('e') if x != 0 else ('0', '0')
    # shortest sci with minimal digits
    s = f'{x:.15g}'
    for d in range(0, 16):
        t = f'{x:.{d}e}'
        if float(t) == x: break
    mant, ex = t.split('e'); sci = f'{mant}e{ex[0]}{abs(int(ex)):02d}'
    fixed = np.format_float_positional(x, trim='-')
    return sci if len(sci) < len(fixed) else fixed
def null_test(tab, setdf, nrep=10000, seed=SEED):
    """port of combine_datasets.R expression-matched test. tab: index gene, Z, q, expr_percentile (rows with percentile only)."""
    tab = tab[tab.expr_percentile.notna()]
    bins = bins_from(tab.expr_percentile); by_bin = {b: np.where(bins.values == b)[0] for b in range(10)}
    Zs = np.sign(tab.Z.values); hit = (tab.q.values < 0.05)
    out = {}
    rng = np.random.default_rng(seed)
    for st, r in setdf.groupby('set', sort=False):
        r = r.drop_duplicates('gene'); keep = r.gene.isin(tab.index).values; r = r[keep]
        idx = tab.index.get_indexer(r.gene); sg = r.sz_sign.values
        obs_idx = int((hit[idx] & (Zs[idx] == -sg)).sum()); any_idx = int(hit[idx].sum())
        need = np.array([(bins.values[idx] == b).sum() for b in range(10)])
        # exact expectation with R's positional pairing (draws in bin order, sg in set-file order)
        pos_bin = np.repeat(np.arange(10), need)
        fneg = np.array([(hit[by_bin[b]] & (Zs[by_bin[b]] == -1)).mean() if len(by_bin[b]) else 0 for b in range(10)])
        fpos = np.array([(hit[by_bin[b]] & (Zs[by_bin[b]] == 1)).mean() if len(by_bin[b]) else 0 for b in range(10)])
        exp_exact = float(np.where(sg > 0, fneg[pos_bin], np.where(sg < 0, fpos[pos_bin], 0)).sum())
        exp_any = float(sum(need[b] * hit[by_bin[b]].mean() for b in range(10) if len(by_bin[b])))
        # Monte Carlo
        n1 = np.zeros(nrep); n2 = np.zeros(nrep); off = 0
        for b in range(10):
            k = need[b]
            if k == 0: continue
            pool = by_bin[b]
            perm = np.argsort(rng.random((nrep, len(pool))), axis=1)[:, :k]
            drawn = pool[perm]                       # nrep x k
            sgb = sg[off:off + k]; off += k
            n1 += (hit[drawn] & (Zs[drawn] == -sgb[None, :])).sum(1); n2 += hit[drawn].sum(1)
        out[st] = dict(n_set=len(setdf[setdf.set == st].drop_duplicates('gene')), n_measured=len(idx), observed_idx=obs_idx, observed_any_sign=any_idx,
                       expected=round(float(n1.mean()), 2), expected_exact=exp_exact, p=signif((1 + int((n1 >= obs_idx).sum())) / (nrep + 1), 3),
                       any_mean=float(n2.mean()), any_p=(1 + int((n2 >= any_idx).sum())) / (nrep + 1), nrep=nrep)
    return out
if __name__ == '__main__':
    print('q check v2 table BH:', check_q())
    r1 = null_test(c1.reset_index().set_index('gene'), S)
    st = rd(f'{RES}/combined_summary.tsv'); st = st[st.test.str.startswith('combined control pH q < 0.05, expression-matched')]
    for _, r in st.iterrows():
        m = r1[r['set']]; print(r['set'], 'stored obs', r.observed, 'exp', r.expected, 'p', r.combined_p_one_sided, '| mine obs(idx)', m['observed_idx'], 'exp MC', m['expected'], 'exact', round(m['expected_exact'], 2), 'p', m['p'], 'any', m['observed_any_sign'])
