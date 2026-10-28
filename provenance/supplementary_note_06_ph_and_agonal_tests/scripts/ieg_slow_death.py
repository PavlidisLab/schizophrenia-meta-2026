#!/usr/bin/env python3
"""Are immediate early genes lower with slow death (GTEx Hardy slow vs fast) or a high agonal factor (Pritzker GSE92538 GPL17027)? 2026-10-09.
Axes (t > 0 = higher in slow death / higher with more agonal): analyses/gtex_ischemic_time/hardy_t_by_region.tsv.gz (slow death Hardy 4 vs fast natural 2,
adjusted for ischemic time, age, sex; cortical columns Cortex, Frontal Cortex (BA9), ACC (BA24); axis = mean of non-missing), analyses/axes_pritzker_agonal_2026-10-09/axis_pritzker_agonal_GPL17027.tsv.gz
(t_agonal, and t_agonal_with_pH). Sets: Tyssowski 2018 rapid PRG, delayed PRG, secondary response genes as mapped to human in ph_sensitivity_gene_features_2026-10-09/results_laptop/activity_pH_mean.tsv
(gene lists in its 'genes' column); expression decile from combined_controls_pH.tsv expr_percentile. Set score = mean t; null = 10,000 random sets of the same size matched on expression decile; two-sided permutation p.
Outputs: gene_t.tsv (named genes), set_scores.tsv."""
import numpy as np, pandas as pd, os
R = "schizophrenia-transcriptome-replication/analyses"; H = os.path.dirname(os.path.abspath(__file__)); rng = np.random.default_rng(20261009); NP = 10000
g = pd.read_csv(f"{R}/gtex_ischemia_2026-10-08/hardy_t_by_region.tsv.gz", sep="\t", index_col=0, compression="gzip")
cols = [c for c in g.columns if c in ("Brain - Cortex", "Brain - Frontal Cortex (BA9)", "Brain - Anterior cingulate cortex (BA24)")]
hardy = g[cols].mean(axis=1, skipna=True).where(g[cols].notna().sum(axis=1) >= 2)
p = pd.read_csv(f"{R}/axes_pritzker_agonal_2026-10-09/axis_pritzker_agonal_GPL17027.tsv.gz", sep="\t", comment="#", index_col="gene", compression="gzip")
ph = pd.read_csv(f"{R}/confounds/results/combined_controls_pH.tsv", sep="\t", comment="#"); ph = ph[ph.factor == "pH"].set_index("gene")
D = pd.DataFrame({"hardy": hardy, "agonal": p.t_agonal.groupby(level=0).mean(), "agonal_pH": p.t_agonal_with_pH.groupby(level=0).mean(), "pct": ph.expr_percentile}); D["dec"] = np.minimum((D.pct * 10).astype(int, errors="ignore"), 9) if D.pct.notna().all() else np.minimum((D.pct.fillna(-1) * 10).astype(int), 9)
named = ["EGR1", "NR4A1", "NR4A2", "DUSP1", "FOS", "ARC", "JUNB", "NPAS4"]
rows = []
for gname in named:
    r = {"gene": gname}
    for c in cols: r[c.replace("Brain - ", "GTEx ")] = round(g.loc[gname, c], 2) if gname in g.index and pd.notna(g.loc[gname, c]) else None
    r["GTEx_cortical_mean_t"] = round(D.hardy.get(gname, np.nan), 2); r["Pritzker_t_agonal"] = round(D.agonal.get(gname, np.nan), 2); r["Pritzker_t_agonal_with_pH"] = round(D.agonal_pH.get(gname, np.nan), 2)
    rows.append(r)
pd.DataFrame(rows).to_csv(f"{H}/gene_t.tsv", sep="\t", index=False); print(pd.DataFrame(rows).to_string(index=False))
A = pd.read_csv(f"{R}/ph_sensitivity_gene_features_2026-10-09/results_laptop/activity_pH_mean.tsv", sep="\t", comment="#")
sets = {r.activity_class: r.genes.split(",") for r in A.itertuples()}
Dn = D.dropna(subset=["pct"]); by = {k: Dn.index[Dn.dec == k].values for k in range(10)}
out = []
for sname, G in sets.items():
    for axis in ("hardy", "agonal", "agonal_pH"):
        Gg = [x for x in G if x in D.index and pd.notna(D.loc[x, axis]) and x in Dn.index]
        if len(Gg) < 5: continue
        need = Dn.loc[Gg].dec.value_counts().to_dict(); obs = D.loc[Gg, axis].mean()
        null = np.array([D.loc[np.concatenate([rng.choice(by[k], c, replace=False) for k, c in need.items()]), axis].mean() for _ in range(NP)])
        out.append({"set": sname, "axis": axis, "n_genes": len(Gg), "mean_t": round(obs, 3), "matched_null_mean": round(null.mean(), 3), "effect": round(obs - null.mean(), 3), "p_two_sided": (1 + np.sum(np.abs(null - null.mean()) >= abs(obs - null.mean()))) / (NP + 1)})
pd.DataFrame(out).to_csv(f"{H}/set_scores.tsv", sep="\t", index=False); print(pd.DataFrame(out).to_string(index=False))
