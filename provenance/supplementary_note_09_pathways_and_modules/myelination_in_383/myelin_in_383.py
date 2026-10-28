#!/usr/bin/env python3
"""Myelination and oxidative phosphorylation genes among the 383 recurrent genes and the 107 genes in >= 10 lists, against chance; list concentration of the myelination hits.
Inputs: recurrence_genes.tsv (manuscript provenance, scenario 'published, brain': 15,844 tested genes; recurrent = q_bh < 0.1 (383); ge10 = n_lists_observed >= 10 (107); columns de_prior_hitlists, expr_decile), term_group_genes.tsv (M09 results_A, groups myelin and oxphos, versions terms and hits, per-gene list directions in 'lists'),
recurrence_units.tsv (cohort cluster of each list). Chance: 20,000 random gene sets (exact hypergeometric draw per stratum) with the same number of genes in each stratum of (a) none, (b) DE-prior quintile (hitlist-count cuts at the 0.2/0.4/0.6/0.8 quantiles of the 15,844), (c) DE-prior quintile x expression quintile; drawn without replacement from the tested genes
(for the 'hits' sets also from tested genes with >= 1 list, since list genes are in >= 1 list by construction). p = share of draws at least (upper) / at most (lower) as many as observed, +1 smoothed.
Usage: myelin_in_383.py RECURRENCE_GENES TERM_GROUP_GENES RECURRENCE_UNITS OUTDIR
"""
import sys, hashlib, re, collections
import numpy as np, pandas as pd
RG, TG, UN, OUT = sys.argv[1:5]
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
rng = np.random.default_rng(20261023)
d = pd.read_csv(RG, sep="\t", comment="#"); d = d[d.scenario == "published, brain"].drop_duplicates("gene").set_index("gene")
d["rec383"] = (d.q_bh < 0.1).astype(int); d["ge10"] = (d.n_lists_observed >= 10).astype(int)
assert d.rec383.sum() == 383 and d.ge10.sum() == 107
cuts = d.de_prior_hitlists.quantile([.2, .4, .6, .8]).values
d["dq"] = np.where(d.de_prior_hitlists.isna(), -1, np.digitize(d.de_prior_hitlists.fillna(-1), cuts))
d["eq"] = np.digitize(d.expr_decile, [3, 5, 7, 9])      # expression deciles 1-10 -> quintiles 0-4
tg = pd.read_csv(TG, sep="\t", comment="#", dtype=str); tg = tg[tg.group.isin(["myelin", "oxphos"])]
NS = 20000; rows = []; gene_rows = []
def draw(strata_key, gset, pool_df):
    """counts of rec383 / ge10 in NS stratum-matched random sets: per stratum the count is hypergeometric (m genes drawn without replacement from the stratum's pool)"""
    need = gset.groupby(strata_key).size()
    rec = np.zeros(NS, int); g10 = np.zeros(NS, int)
    for k, m in need.items():
        pool = pool_df[pool_df[strata_key] == k] if strata_key != "none" else pool_df
        N = len(pool); m = min(m, N)
        rec += rng.hypergeometric(int(pool.rec383.sum()), N - int(pool.rec383.sum()), m, size=NS)
        g10 += rng.hypergeometric(int(pool.ge10.sum()), N - int(pool.ge10.sum()), m, size=NS)
    return rec, g10
d["none"] = 0; d["dq_eq"] = d["dq"] * 10 + d["eq"]
for g in ("myelin", "oxphos"):
    for v in ("terms", "hits"):
        sub = tg[(tg.group == g) & (tg.version == v)].drop_duplicates("gene")
        inu = sub[sub.gene.isin(d.index)]; gs = d.loc[inu.gene]
        obs = (int(gs.rec383.sum()), int(gs.ge10.sum()))
        pools = {"all tested genes": d, "tested genes with >= 1 list": d[d.n_lists_observed >= 1]} if v == "hits" else {"all tested genes": d}
        for pn, pool in pools.items():
            for mname, key in (("uniform", "none"), ("DE-prior quintile", "dq"), ("DE-prior x expression quintile", "dq_eq")):
                gset = gs.copy()
                if key != "none":
                    # genes whose stratum is empty in the pool are dropped from the draw count (reported in n_drawn)
                    gset = gset[gset[key].isin(pool[key].unique())]
                rec, g10 = draw(key, gset, pool)
                for lab, o, nd in (("383", obs[0], rec), ("107", obs[1], g10)):
                    rows.append(dict(group=g, version=v, n_term_file=len(sub), n_tested=len(inu), n_drawn=len(gset), background=pn, matched_on=mname, set=lab, observed=o, expected=float(nd.mean()), null_lo=float(np.percentile(nd, 2.5)), null_hi=float(np.percentile(nd, 97.5)),
                                     ratio=o / nd.mean() if nd.mean() > 0 else np.nan, p_upper=float((1 + (nd >= o).sum()) / (NS + 1)), p_lower=float((1 + (nd <= o).sum()) / (NS + 1))))
        for gene in inu.gene:
            r = d.loc[gene]
            gene_rows.append(dict(group=g, version=v, gene=gene, in_383=int(r.rec383), in_107=int(r.ge10), n_lists_observed=int(r.n_lists_observed), q_bh=r.q_bh, de_prior_hitlists=r.de_prior_hitlists, expr_decile=int(r.expr_decile), sz_direction=r.sz_direction if "sz_direction" in r else "", up=r.up, down=r.down))
R = pd.DataFrame(rows); G = pd.DataFrame(gene_rows)
# ---- list concentration of the myelination hits
un = pd.read_csv(UN, sep="\t", comment="#", header=None)
cl = {}
for r in un.itertuples(index=False):
    if str(r[3]).startswith("yes") and "published, brain" in str(r[0]) or str(r[0]) == "published only, all tissues" and str(r[3]).startswith("yes"): cl[str(r[2])] = str(r[4])
mh = tg[(tg.group == "myelin") & (tg.version == "hits")].drop_duplicates("gene")
per = collections.defaultdict(lambda: [0, 0, 0, 0]); nlist_per_gene = []
for r in mh.itertuples():
    toks = [t for t in str(r.lists).split(";") if ":" in t]; names = [t.rsplit(":", 1)[0].strip() for t in toks]; nlist_per_gene.append(len(names))
    for t in toks:
        nm, dr = t.rsplit(":", 1); per[nm.strip()][0] += 1; per[nm.strip()][1] += dr.strip() == "up"; per[nm.strip()][2] += dr.strip() == "down"
    rr = d.loc[r.gene] if r.gene in d.index else None
L = pd.DataFrame([dict(list=k, cohort_cluster=cl.get(k, ""), myelin_hit_genes=v[0], up=v[1], down=v[2]) for k, v in per.items()]).sort_values("myelin_hit_genes", ascending=False)
L["share_of_gene_list_reports"] = L.myelin_hit_genes / L.myelin_hit_genes.sum()
# genes recurrent only because of the six lists: lists from outside the six
six = set(per)
studs = d.studies.fillna("")
mg = G[(G.group == "myelin") & (G.version == "hits")].copy()
def n_outside(g):
    names = [s.strip() for s in studs.loc[g].split(";") if s.strip()]
    return sum(1 for s in names if re.sub(r" \(own list\)", "", s) not in six)
mg["n_lists_outside_the_six"] = [n_outside(g) for g in mg.gene]
mg["n_genes_lists_per_gene_in_six"] = nlist_per_gene[:len(mg)] if len(nlist_per_gene) == len(mg) else np.nan
hdr = f"# myelin_in_383.py (myelination_in_383); recurrence_genes.tsv sha256 {sha(RG)}; term_group_genes.tsv {sha(TG)}; recurrence_units.tsv {sha(UN)}; {NS} draws, seed 20261023\n"
for nm, df in (("enrichment_vs_chance.tsv", R), ("genes.tsv", G), ("myelin_hits_by_list.tsv", L), ("myelin_hits_recurrent_detail.tsv", mg[mg.in_383 + mg.in_107 > 0] if len(mg) else mg)):
    with open(f"{OUT}/{nm}", "w") as fh: fh.write(hdr); df.round(5).to_csv(fh, sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(R[(R.matched_on != "uniform") | (R.background == "all tested genes")].round(3).to_string())
print(L.to_string()); print("genes in six-list hits:", len(mh), "lists per gene:", collections.Counter(nlist_per_gene))
print(mg[(mg.in_383 == 1)][["gene", "in_107", "n_lists_observed", "n_lists_outside_the_six"]].to_string())
