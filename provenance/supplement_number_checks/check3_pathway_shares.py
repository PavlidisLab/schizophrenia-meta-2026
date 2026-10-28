#!/usr/bin/env python3
"""Check 3: shares of the excess of shared top pathways over the DE-prior null. Inputs: provenance/supplementary_note_09_pathways_and_modules/shared_pathways_by_group.tsv and a5b_agreement.tsv.
excess = observed shared pathway pairs (710) - DE-prior null mean (265.6, 95% 211.4 to 344.3) = 444.4; share of the excess = pairs of the group / excess, assuming the null contributes about 0 pairs to these terms (deprior_null_mean 0.01 to 0.03 lists per term).
Usage: check3_pathway_shares.py SHARED_BY_GROUP A5B_AGREEMENT OUT"""
import sys, hashlib, pandas as pd
G, A, OUT = sys.argv[1:4]; sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
g = pd.read_csv(G, sep="\t", comment="#"); a = pd.read_csv(A, sep="\t", comment="#"); r = a[a.metric == "shared_pathways_total"].iloc[0]
obs, nm, lo, hi = r.observed, r.deprior_null_mean, r.deprior_null_lo95, r.deprior_null_hi95; ex = obs - nm
rows = []
for _, x in g.iterrows():
    p = float(x.shared_pathway_pairs); rows.append(dict(group=x.group, shared_pathway_pairs=int(p), share_of_all_710_pairs=p / obs, share_of_excess=p / ex, share_of_excess_if_null_is_mean_hi95=p / (obs - hi), share_of_excess_if_null_is_mean_lo95=p / (obs - lo)))
R = pd.DataFrame(rows).round(4); open(OUT, "w").write(f"# check3_pathway_shares.py; shared_pathways_by_group.tsv {sha(G)}; a5b_agreement.tsv {sha(A)}; observed {obs:.0f}, DE-prior null mean {nm} (95% {lo} to {hi}), excess {ex:.1f}\n"); R.to_csv(OUT, sep="\t", index=False, mode="a"); print(R.to_string())
