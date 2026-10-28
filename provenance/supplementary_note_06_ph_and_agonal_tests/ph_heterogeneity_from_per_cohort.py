#!/usr/bin/env python3
"""Cochran's Q, I-squared, DerSimonian-Laird mean and 95% prediction interval over the cohort rows of
figures/data/ph_by_cohort/per_cohort_80.tsv (set A: A_est with its 95% t interval; the standard error is the interval
half-width divided by the t quantile on n_SZ + n_CTRL - 2 degrees of freedom).
Usage: python3 ph_heterogeneity_from_per_cohort.py per_cohort_80.tsv [previous_per_cohort_80.tsv] > out.tsv"""
import sys
import numpy as np
import pandas as pd
from scipy import stats


def het(path):
    df = pd.read_csv(path, sep="\t", comment="#")
    d = df[~df["cohort"].str.startswith("all")]
    y = d["A_est"].values
    tq = stats.t.ppf(0.975, (d["A_n_SZ"] + d["A_n_CTRL"] - 2).clip(lower=1))
    se = ((d["A_hi"] - d["A_lo"]) / (2 * tq)).values
    k = len(y)
    w = 1 / se**2
    m = (w * y).sum() / w.sum()
    q = (w * (y - m) ** 2).sum()
    c = w.sum() - (w**2).sum() / w.sum()
    tau2 = max(0.0, (q - (k - 1)) / c)
    i2 = max(0.0, (q - (k - 1)) / q)
    ws = 1 / (se**2 + tau2)
    mu = (ws * y).sum() / ws.sum()
    sem = ws.sum() ** -0.5
    half = stats.t.ppf(0.975, k - 2) * np.sqrt(tau2 + sem**2)
    return dict(n_cohorts=k, Q=q, df=k - 1, I2=i2, tau2=tau2, re_mean=mu, re_lo=mu - 1.96 * sem, re_hi=mu + 1.96 * sem,
                pred_lo=mu - half, pred_hi=mu + half)


rows = [("current", het(sys.argv[1]))]
if len(sys.argv) > 2:
    rows.append(("previous grouping", het(sys.argv[2])))
print("grouping\t" + "\t".join(rows[0][1].keys()))
for name, r in rows:
    print(name + "\t" + "\t".join(f"{v:.4f}" if isinstance(v, float) else str(v) for v in r.values()))
