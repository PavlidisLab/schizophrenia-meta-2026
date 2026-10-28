#!/usr/bin/env python3
"""Check 6: which people-counts belong to which pH sentence. Reads per_cohort_80.tsv (expr14 and other columns, pooled row), fig07_a1_ph_difference.tsv and fig07_a1_ph_per_person.tsv.
Usage: check6_ph_people.py PER_COHORT FIG07_DIFF FIG07_PERSON OUT"""
import sys, hashlib, pandas as pd
PC, DF, PP, OUT = sys.argv[1:5]; sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
c = pd.read_csv(PC, sep="\t", comment="#"); p = c[c.cohort.str.startswith("all cohorts")].iloc[0]; df = pd.read_csv(DF, sep="\t", comment="#"); df = df[df.in_14_dataset_regions == "yes"]; pp = pd.read_csv(PP, sep="\t", comment="#")
rows = [("people with pH linked to the 80 studies (set A; Note 6 pooled)", int(p.A_n_SZ), int(p.A_n_CTRL), f"{p.A_est:.3f} ({p.A_lo:.3f} to {p.A_hi:.3f})", "per_cohort_80.tsv pooled row, A_*"),
        ("  of them in the 14 datasets with expression data (expr14)", int(p.expr14_n_SZ), int(p.expr14_n_CTRL), f"{p.expr14_est:.3f} ({p.expr14_lo:.3f} to {p.expr14_hi:.3f})", "per_cohort_80.tsv pooled row, expr14_*"),
        ("  of them not in those datasets (other)", int(p.other_n_SZ), int(p.other_n_CTRL), f"{p.other_est:.3f} ({p.other_lo:.3f} to {p.other_hi:.3f})", "per_cohort_80.tsv pooled row, other_*"),
        ("samples with pH in the 14 dataset-regions (a person can be in 2 or 3 regions)", int(df.n_sz.sum()), int(df.n_ctl.sum()), "", "fig07_a1_ph_difference.tsv, in_14_dataset_regions = yes"),
        ("person rows in the per-person table of the 14 dataset-regions", int((pp.group != "control").sum()), int((pp.group == "control").sum()), "", "fig07_a1_ph_per_person.tsv"),
        ("ten datasets of the pH-in-controls analysis, each person once (Fig. 3b, Extended Methods 7)", 237, 279, "", "Figures_ph.md; manuscript_ph.md Methods (279 controls, 237 cases; 516 people)")]
R = pd.DataFrame(rows, columns=["count", "n_cases", "n_controls", "case minus control pH (95% interval)", "source"]); R["n_total"] = R.n_cases + R.n_controls
open(OUT, "w").write(f"# check6_ph_people.py; per_cohort_80.tsv {sha(PC)}; fig07_a1_ph_difference.tsv {sha(DF)}; fig07_a1_ph_per_person.tsv {sha(PP)}\n"); R.to_csv(OUT, sep="\t", index=False, mode="a"); print(R.to_string())
