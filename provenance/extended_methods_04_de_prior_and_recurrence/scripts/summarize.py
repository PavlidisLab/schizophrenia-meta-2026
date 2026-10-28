#!/usr/bin/env python3
"""Summaries and figure from the TSVs of stability.py (laptop, light). Writes summary_subsets.tsv, summary_bootstrap.tsv, thresholds.tsv, fig_prior_stability.svg."""
import os, sys, math, json, hashlib, datetime
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(H, "out")
rd = lambda f: pd.read_csv(os.path.join(O, f), sep="\t", comment="#")
sh = lambda f: hashlib.sha256(open(os.path.join(O, f), "rb").read()).hexdigest()[:12]
def head(src): return [f"# script: analyses/de_prior_stability_2026-10-12/summarize.py, run {datetime.date.today()}", f"# source: {src}"]
def put(fn, df, src):
    with open(os.path.join(O, fn), "w") as fh:
        fh.write("\n".join(head(src)) + "\n"); df.to_csv(fh, sep="\t", index=False, float_format="%.5g")
q = lambda x: (np.nanmedian(x), np.nanpercentile(x, 2.5), np.nanpercentile(x, 97.5))
S = rd("subset_runs.tsv"); rows = []
for (pool, n), g in S.groupby(["pool", "n"], sort=False):
    for m in ("rho_full_count", "rho_disjoint_count", "rho_full_freq", "rho_disjoint_freq", "auroc_count", "auroc_freq"):
        v = g[m].values
        if np.isnan(v).all(): rows.append(dict(pool=pool, n=n, metric=m, median=np.nan, lo=np.nan, hi=np.nan, n_reps=len(g), note="not available (2n > pool)" if "disjoint" in m else "")); continue
        md, lo, hi = q(v); rows.append(dict(pool=pool, n=n, metric=m, median=md, lo=lo, hi=hi, n_reps=len(g), note="n = pool: identical to the full prior" if n == g["n"].max() and pool and (n in (8712, 242)) and "full" in m else ""))
SS = pd.DataFrame(rows).sort_values(["pool", "metric", "n"]); put("summary_subsets.tsv", SS, f"out/subset_runs.tsv sha256 {sh('subset_runs.tsv')}; median and 2.5-97.5% over subsets")
B = rd("bootstrap_draws.tsv"); G = rd("bootstrap_genes.tsv"); br = []
for m in ("rho_pair_count", "rho_full_count", "auroc_count", "rho_pair_freq", "rho_full_freq", "auroc_freq", "n_contrasts_distinct"):
    md, lo, hi = q(B[m].values); br.append(dict(quantity=m, median=md, lo=lo, hi=hi, mean=B[m].mean()))
f = G.freq.values; sd = G.boot_sd_freq.values; ok = np.isfinite(sd) & np.isfinite(f)
lo45, hi45 = np.nanpercentile(f[ok], [45, 55]); band = ok & (f >= lo45) & (f <= hi45)
ntest = G.n_tested.values; binom = f * (1 - f) / ntest
br.append(dict(quantity="boot_sd_freq, all universe genes (median; 2.5-97.5% over genes)", median=np.median(sd[ok]), lo=np.percentile(sd[ok], 2.5), hi=np.percentile(sd[ok], 97.5), mean=sd[ok].mean()))
br.append(dict(quantity=f"boot_var_freq, genes with full frequency in the 45-55th percentile ({lo45:.4f} to {hi45:.4f}; n={int(band.sum())}); median, IQR-free 2.5-97.5% over genes", median=np.median(sd[band] ** 2), lo=np.percentile(sd[band] ** 2, 2.5), hi=np.percentile(sd[band] ** 2, 97.5), mean=(sd[band] ** 2).mean()))
br.append(dict(quantity="boot_sd_freq, same band", median=np.median(sd[band]), lo=np.percentile(sd[band], 2.5), hi=np.percentile(sd[band], 97.5), mean=sd[band].mean()))
br.append(dict(quantity="binomial p(1-p)/n_tested, same band (reference, contrasts independent)", median=np.median(binom[band]), lo=np.percentile(binom[band], 2.5), hi=np.percentile(binom[band], 97.5), mean=binom[band].mean()))
br.append(dict(quantity="median full frequency, universe genes", median=np.median(f[ok]), lo=np.percentile(f[ok], 2.5), hi=np.percentile(f[ok], 97.5), mean=f[ok].mean()))
br.append(dict(quantity="boot_sd_freq / freq (CV), band", median=np.median(sd[band] / f[band]), lo=np.percentile(sd[band] / f[band], 2.5), hi=np.percentile(sd[band] / f[band], 97.5), mean=(sd[band] / f[band]).mean()))
br.append(dict(quantity="genes with a draw in which they were untested (frequency undefined)", median=int((G.draws_untested > 0).sum()), lo=np.nan, hi=np.nan, mean=np.nan))
put("summary_bootstrap.tsv", pd.DataFrame(br), f"out/bootstrap_draws.tsv sha256 {sh('bootstrap_draws.tsv')}; out/bootstrap_genes.tsv sha256 {sh('bootstrap_genes.tsv')}; median and 2.5-97.5% over draws (or genes where stated)")
# thresholds: n at which median rho (vs full prior and vs disjoint subset) first reaches 0.9, log-linear interpolation; level-off by increments per doubling
tr = []
for pool in S.pool.unique():
    for m in ("rho_full_count", "rho_disjoint_count", "rho_full_freq", "rho_disjoint_freq"):
        d = SS[(SS.pool == pool) & (SS.metric == m)].dropna(subset=["median"]).sort_values("n"); x, y = d.n.values.astype(float), d["median"].values
        first = next((int(a) for a, b in zip(x, y) if b >= 0.9), None); interp = None
        for i in range(1, len(x)):
            if y[i] >= 0.9 > y[i - 1]: interp = float(np.exp(np.log(x[i - 1]) + (0.9 - y[i - 1]) / (y[i] - y[i - 1]) * (np.log(x[i]) - np.log(x[i - 1])))); break
        incs = [(x[i - 1], x[i], (y[i] - y[i - 1]) / math.log2(x[i] / x[i - 1])) for i in range(1, len(x))]
        lev = next((f"{int(a)}->{int(b)}" for a, b, v in incs if v < 0.01), None)
        tr.append(dict(pool=pool, metric=m, first_grid_n_median_ge_0p9=first if first else "not reached", interpolated_n_for_0p9=round(interp) if interp else "", first_step_with_increment_lt_0p01_per_doubling=lev or "none",
                       increments_per_doubling="; ".join(f"{int(a)}->{int(b)}:{v:.3f}" for a, b, v in incs), max_n=int(x.max())))
put("thresholds.tsv", pd.DataFrame(tr), "summary_subsets.tsv")
# figure
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "sans-serif"], "svg.fonttype": "none", "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "axes.grid.axis": "y", "grid.color": "#e5e7eb", "grid.linewidth": 0.6, "axes.edgecolor": "#6b7280", "xtick.color": "#6b7280", "ytick.color": "#6b7280",
                     "axes.labelcolor": "#1f2937", "axes.titleweight": "normal", "axes.titlelocation": "left", "font.size": 9, "axes.titlesize": 10})
fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.3))
col = {"all experiments": "#2563eb", "postmortem brain": "#f59e0b"}
for ax, (m_full, m_dis, ttl) in zip(axs[:2], (("rho_full_count", "rho_disjoint_count", "Rank agreement of the prior"), (None, None, None))):
    pass
ax = axs[0]
for pool, c in col.items():
    for m, ls, lab in (("rho_full_count", "-", "vs full prior"), ("rho_disjoint_count", "--", "vs disjoint subset")):
        d = SS[(SS.pool == pool) & (SS.metric == m)].dropna(subset=["median"]).sort_values("n")
        ax.plot(d.n, d["median"], ls, color=c, marker="o", ms=3, lw=1.3, label=f"{pool}, {lab}"); ax.fill_between(d.n, d.lo, d.hi, color=c, alpha=0.12, lw=0)
ax.axhline(0.9, color="#6b7280", lw=0.8, ls=":"); ax.set_xscale("log"); ax.set_xlabel("contrasts in the subset (n)"); ax.set_ylabel("Spearman rho (n_de_total)"); ax.set_title("Rank agreement"); ax.set_ylim(0, 1.02)
ax.legend(frameon=False, fontsize=6.5, loc="lower right")
ax = axs[1]
for pool, c in col.items():
    d = SS[(SS.pool == pool) & (SS.metric == "auroc_count")].sort_values("n"); ax.plot(d.n, d["median"], "-", color=c, marker="o", ms=3, lw=1.3, label=pool); ax.fill_between(d.n, d.lo, d.hi, color=c, alpha=0.12, lw=0)
ax.axhline(0.5, color="#6b7280", lw=0.8, ls=":"); ax.set_xscale("log"); ax.set_xlabel("contrasts in the subset (n)"); ax.set_ylabel("AUROC of the prior for the 209"); ax.set_title("AUROC for the 209"); ax.legend(frameon=False, fontsize=7)
ax = axs[2]
for v, c in ((B.auroc_count.values, "#f59e0b"),):
    ax.hist(v, bins=30, color=c, alpha=0.85)
ax.axvline(0.9245, color="#1f2937", lw=1); ax.set_xlabel("AUROC of the 209 (bootstrap priors)"); ax.set_ylabel("draws"); ax.set_title("Postmortem-brain bootstrap")
fig.text(0.01, 0.005, "Source: de_prior_stability_2026-10-12/out/summary_subsets.tsv, bootstrap_draws.tsv. Bands: 2.5-97.5% over 200 subsets; line at rho = 0.9 and observed AUROC 0.924.", fontsize=7, color="#6b7280")
fig.tight_layout(rect=(0, 0.04, 1, 1)); fig.savefig(os.path.join(O, "fig_prior_stability.svg")); print("figure ok")
