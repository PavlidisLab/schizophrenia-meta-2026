#!/usr/bin/env python3
"""Agonal-axis analogue of the pH rows of Fig. 3b / Fig. 3d (RA, 2026-10-22).
(a) Pearson r, GTEx slow-minus-fast axis (mean t over BA9, BA24, Cortex) vs the schizophrenia effect z_sz_raw (age+sex; the statistic used for the pH r).
(b) set genes with axis q < 0.05 whose sign matches the SZ direction (SZ-up genes higher after slow death, SZ-down lower) vs 10,000 expression-decile-matched random sets (null_test of f07_core.py, seed 20261015)."""
import os, sys, math, gzip, hashlib
import numpy as np, pandas as pd
from scipy import stats
H = os.path.expanduser("~"); MS = "."  # run from the repository root
REPO = os.environ.get("ANALYSIS_REPO", f"{H}/Dev/schizophrenia-transcriptome-replication")  # analysis repository checkout
SC = f"{MS}/figures/data/fig07/fig07_a2_a4_gene_scatter.tsv"
SC5 = f"{MS}/figures/data/fig07/fig07_a2_a4_gene_scatter_ge5.tsv"
HARDY = f"{REPO}/analyses/gtex_ischemia_2026-10-08/hardy_t_by_region.tsv.gz"
CNT = f"{REPO}/figure_data/gtex_death_speed_all_tissues/donor_counts_by_tissue.tsv"
PCT = f"{REPO}/analyses/confounds/results/combined_controls_pH.tsv"
REG = ["Brain - Frontal Cortex (BA9)", "Brain - Anterior cingulate cortex (BA24)", "Brain - Cortex"]
NREP, SEED = 10000, 20261015
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]; return np.minimum(q, 1)
def rd(p): return pd.read_csv(p, sep="\t", comment="#", dtype=str, keep_default_na=False)
def bins_from(pct):
    e = np.nanquantile(pct.values, np.arange(11) / 10)
    return np.clip(np.searchsorted(e, pct.values, side="left") - 1, 0, 9)
# axis
cnt = rd(CNT); cnt = cnt[cnt.tissue.isin(REG)].set_index("tissue"); df = {r: int(cnt.loc[r, "n_hardy_4"]) + int(cnt.loc[r, "n_hardy_2"]) - 5 for r in REG}
hd = pd.read_csv(HARDY, sep="\t", index_col=0)[REG].dropna()
P = pd.DataFrame({r: 2 * stats.t.sf(hd[r].abs(), df[r]) for r in REG}, index=hd.index)
same = (hd > 0).all(axis=1) | (hd < 0).all(axis=1)
pa = np.where(same, P.max(axis=1), 1.0)
ax = pd.DataFrame({"axis_t": hd.mean(axis=1), "p_iu": pa, "q": bh(pa)}, index=hd.index)
pc = rd(PCT); pc = pc[pc.factor == "pH"].set_index("gene"); pc["expr_percentile"] = pd.to_numeric(pc.expr_percentile, errors="coerce")
ax["pct"] = pc.expr_percentile.reindex(ax.index)
pool = ax[ax.pct.notna()].copy(); pool["bin"] = bins_from(pool.pct)
hit = (pool.q.values < 0.05); Zs = np.sign(pool.axis_t.values)
by_bin = {b: np.where(pool.bin.values == b)[0] for b in range(10)}
# sets
sc = rd(SC).drop_duplicates("gene").set_index("gene")
for c in ("z_ph_controls", "z_sz_raw"): sc[c] = pd.to_numeric(sc[c], errors="coerce")
sc["sgn"] = sc.sz_direction.map({"up": 1, "down": -1}).fillna(0).astype(int)
sc["axis_t"] = ax.axis_t.reindex(sc.index)
sets = {"107 genes in 10 or more lists": sc[sc.in_ge10 == "yes"], "383 recurrent genes": sc[sc.in_176 == "yes"]}
sets["105 genes in 10 or more lists with a pH measure"] = sets["107 genes in 10 or more lists"].dropna(subset=["z_ph_controls"])
sets["371 of the 383 with a pH measure"] = sets["383 recurrent genes"].dropna(subset=["z_ph_controls"])
def fisher(r, n):
    z = np.arctanh(r); se = 1 / math.sqrt(n - 3); return math.tanh(z - 1.96 * se), math.tanh(z + 1.96 * se)
rng = np.random.default_rng(SEED); rows_a = []; rows_b = []
for name, s in sets.items():
    # (a)
    d = s.dropna(subset=["axis_t", "z_sz_raw"]); r = stats.pearsonr(d.axis_t, d.z_sz_raw)[0]; lo, hi = fisher(r, len(d))
    bs = []; ar = d.axis_t.values; zr = d.z_sz_raw.values; bg = np.random.default_rng(1)
    for _ in range(2000):
        i = bg.integers(0, len(d), len(d)); bs.append(np.corrcoef(ar[i], zr[i])[0, 1])
    sp = stats.spearmanr(d.axis_t, d.z_sz_raw)[0]
    rows_a.append(dict(set=name, n_set=len(s), n_with_axis_and_sz_effect=len(d), pearson_r=round(r, 4), ci95_fisher=f"{lo:.3f} to {hi:.3f}", ci95_gene_bootstrap=f"{np.percentile(bs,2.5):.3f} to {np.percentile(bs,97.5):.3f}", spearman_rho=round(sp, 4)))
    # (b)
    m = s[s.index.isin(pool.index)]; idx = pool.index.get_indexer(m.index); sg = m.sgn.values
    obs = int((hit[idx] & (Zs[idx] == sg) & (sg != 0)).sum()); anysign = int(hit[idx].sum())
    need = np.array([(pool.bin.values[idx] == b).sum() for b in range(10)])
    n1 = np.zeros(NREP); off = 0
    for b in range(10):
        k = need[b]
        if k == 0: continue
        pl = by_bin[b]; drawn = pl[np.argsort(rng.random((NREP, len(pl))), axis=1)[:, :k]]
        sgb = sg[off:off + k]; off += k   # positional pairing as in null_test (set order within bins is not sorted; see note)
        n1 += (hit[drawn] & (Zs[drawn] == sgb[None, :]) & (sgb[None, :] != 0)).sum(1)
    rows_b.append(dict(set=name, n_set=len(s), n_measured_axis=len(m), observed_q05_in_slow_death_direction=obs, expected_matched_mean=round(float(n1.mean()), 2), ratio=round(obs / float(n1.mean()), 2), p=round((1 + int((n1 >= obs).sum())) / (NREP + 1), 5), observed_q05_any_sign=anysign, n_up=int((sg > 0).sum()), n_down=int((sg < 0).sum()), n_no_direction=int((sg == 0).sum()), n_draws=NREP))
hdr = ["# figure: manuscript_ph.md Results, 'The recurrent genes are enriched for the agonal state-associated pattern' (analogue of the pH rows, Fig. 3b / 3d)",
       "# script: agonal_enrichment_383_2026-10-22/scripts/agonal_enrichment.py, run 2026-10-22",
       f"# input: {SC} sha256 {sha(SC)} (gene sets: in_ge10, in_176 = the 383; z_sz_raw = schizophrenia effect, age+sex, combined over the ten dataset-regions; sz_direction)",
       f"# input: analyses/gtex_ischemia_2026-10-08/hardy_t_by_region.tsv.gz sha256 {sha(HARDY)}; donor_counts_by_tissue.tsv sha256 {sha(CNT)}; combined_controls_pH.tsv sha256 {sha(PCT)} (expr_percentile for deciles)",
       "# axis: mean over BA9, BA24, Cortex of the t of slow death (Hardy 4) minus fast natural death (Hardy 2); positive = higher after slow death. Axis q: intersection-union p (largest regional two-sided p, 1 if signs differ), BH over all genes measured in the three regions (conservative; same as agonal_axis_per_gene_383.tsv).",
       f"# axis genes: {len(ax)} measured in all three regions; {len(pool)} of them with an expression percentile (null pool)"]
for fn, rows, note in (("agonal_t1_correlation.tsv", rows_a, "# (a) Pearson r of the axis value with z_sz_raw, the statistic used for the pH r (not direction-signed). Interval: Fisher, genes independent (as the pH r); gene bootstrap 2,000 resamples alongside."),
                       ("agonal_t2_observed_vs_expected.tsv", rows_b, f"# (b) observed = set genes at axis q < 0.05 whose axis sign equals the schizophrenia direction (SZ-up genes higher after slow death, SZ-down lower); expected = mean over {NREP} random sets matched gene by gene on expression decile, same SZ directions, seed {SEED}; p = (1 + n_ge_observed)/(draws + 1). Genes without a direction cannot count as observed.")):
    with open(fn, "w") as f:
        f.write("\n".join(hdr + [note]) + "\n"); pd.DataFrame(rows).to_csv(f, sep="\t", index=False)
    print(pd.DataFrame(rows).to_string())
