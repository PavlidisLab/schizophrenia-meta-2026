#!/usr/bin/env python3
"""Tables from out_<listset>/ (observed.json, observed_pairs.npz, reps_*.npz). usage: python3 summarize.py LISTSET"""
import sys, os, json, hashlib, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ls = sys.argv[1]; D = os.path.join(HERE, "out_" + ls)
obs = json.load(open(os.path.join(D, "observed.json"))); op = np.load(os.path.join(D, "observed_pairs.npz"))
run = json.load(open(os.path.join(D, "run.json"))); inc = obs["inc"]; unit = obs["unit"]; L = len(inc)
import itertools; P = np.array(list(itertools.combinations(range(L), 2)))
M = ["mean_pathway_jaccard", "shared_pathways_total", "mean_gene_jaccard", "shared_genes_total"]; KS = [2, 3, 5, 10]
def summ(jt, st, jg, sg, m):
    return np.array([np.nanmean(jt[..., m], -1), st[..., m].sum(-1), np.nanmean(jg[..., m], -1), sg[..., m].sum(-1)])
NAMES = {("indep", "uniform"): "uniform_independent", ("indep", "deprior"): "de_prior_independent", ("rot", "R1", "S0"): "rot_R1_S0_lists_independently",
         ("rot", "R2", "S0"): "rot_R2_S0_clusters_together", ("rot", "R1", "S3"): "rot_R1_DEprior40_lists_independently", ("rot", "R2", "S3"): "rot_R2_DEprior40_clusters_together", ("strat", "S3"): "control_DEprior40_strata_independent_no_coexpression"}
def ci(obs_v, v):
    lo, hi = np.nanpercentile(v, [2.5, 97.5]); mu = np.nanmean(v)
    return dict(observed=obs_v, null_mean=mu, null_lo95=lo, null_hi95=hi, obs_over_null=obs_v / mu, log2_obs_over_null=np.log2(obs_v / mu),
                log2_lo95=np.log2(obs_v / hi) if hi > 0 else np.nan, log2_hi95=np.log2(obs_v / lo) if lo > 0 else np.nan, p_upper=(1 + (v >= obs_v - 1e-12).sum()) / (1 + len(v)), n_replicates=len(v))
# unit masks for jackknife
lists_units = [("list:" + inc[i], np.array([i not in (a, b) for a, b in P]) if False else None) for i in range(L)]
def mask_without(idx): s = set(idx); return np.array([a not in s and b not in s for a, b in P])
units = {}
for i in range(L): units[("list", inc[i])] = [i]
for u in sorted(set(unit)): units[("cluster", u)] = [i for i in range(L) if unit[i] == u]
rows = []; jrows = []; sizes = {}
run["cfgs"] = [list(k) for k in NAMES]
for cfg in run["cfgs"]:
    key = tuple(cfg); name = NAMES[key]; f = os.path.join(D, "reps_" + "_".join(cfg[1:]) + ".npz")
    if not os.path.exists(f): continue
    z = np.load(f); jt, st, jg, sg, ge = z["jt"].astype(float), z["st"], z["jg"].astype(float), z["sg"], z["ge"]; sizes[name] = float(z["size"].mean())
    allm = np.ones(len(P), bool); sobs = summ(op["jt"].astype(float), op["st"], op["jg"].astype(float), op["sg"], allm); srep = summ(jt, st, jg, sg, allm)
    for mi, m in enumerate(M):
        r = dict(level="pathway" if "pathway" in m else "gene", metric=m, null=name); r.update(ci(sobs[mi], srep[mi])); r.update(n_lists=L, n_pairs=len(P)); rows.append(r)
        th = {}
        for (ut, un), idx in units.items():
            mk = mask_without(idx); so = summ(op["jt"].astype(float), op["st"], op["jg"].astype(float), op["sg"], mk)[mi]; sr = summ(jt, st, jg, sg, mk)[mi].mean()
            th.setdefault(ut, []).append(np.log(so / sr))
        row = dict(null=name, metric=m, ratio=r["obs_over_null"])
        for ut, v in th.items():
            v = np.array(v); n = len(v); se = np.sqrt((n - 1) / n * ((v - v.mean()) ** 2).sum()); row[f"jack_{ut}_lo"] = r["obs_over_null"] * np.exp(-1.96 * se); row[f"jack_{ut}_hi"] = r["obs_over_null"] * np.exp(1.96 * se)
            row[f"jack_{ut}_n_units"] = n
        jrows.append(row)
    for qi, q in enumerate(KS):
        r = dict(level="gene", metric=f"genes_in_{q}_or_more_lists", null=name); r.update(ci(obs["genes_ge_k"][str(q)], ge[:, qi])); r.update(n_lists=L, n_pairs=len(P)); rows.append(r)
T = pd.DataFrame(rows); J = pd.DataFrame(jrows)
Tm = T.merge(J[["null", "metric"] + [c for c in J.columns if c.startswith("jack_")]], on=["null", "metric"], how="left")
json.dump(sizes, open(os.path.join(D, "mean_size_ratio.json"), "w"))
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
head = [f"# figure: Supplementary Note 5 Figure 1a, coexpression-preserving nulls (this folder replaces only the null)",
        f"# script: analyses/coexpression_null_models/coexpr_null.py (remote) + summarize.py; list set '{ls}': {L} lists, {len(P)} pairs, scenario 'published, brain', rule A primary lists (analyses/primary_lists.json sha256 {sha(os.path.join(HERE,'..','primary_lists.json'))})",
        "# statistics as fairness_lists/a5b_go.py (GO BP 5-200 genes; top pathways = up to 20 terms with p < 0.01 and >= 2 hits; mean pathway/gene Jaccard and shared totals over pairs); genes_in_k_or_more_lists = annotated tested genes that are hits in >= k lists",
        f"# nulls: uniform_independent and de_prior_independent = a5b_go.py's draws (re-run here); rot_* = each list's hits moved by one circular shift along the GTEx v8 cortex co-expression ordering (rotation_lib.py of addendum 4a), within strata S0 (none) or 40 equal-count DE-prior bins; R1 lists shifted independently, R2 lists of one cohort cluster shifted together; {run['nrep']} draws each, seed {run['seed0']}",
        "# log2_lo95 / log2_hi95 = log2 of observed over the null's upper / lower 95% limit (as fig12 panel a); p_upper = (1 + replicates >= observed) / (1 + n); jack_list / jack_cluster = ratio x exp(+-1.96 SE) of a delete-one jackknife over lists (units) / cohort clusters (multi-bank studies singly), pair metrics only"]
out = os.path.join(HERE, f"suppnote5_a_coexpr_null_{ls}.tsv")
with open(out, "w") as fh:
    for h in head: fh.write(h + "\n")
    Tm.to_csv(fh, sep="\t", index=False, float_format="%.5g")
print(Tm[["level", "metric", "null", "observed", "null_mean", "null_lo95", "null_hi95", "obs_over_null", "log2_lo95", "log2_hi95"]].to_string())
print(J.to_string())
