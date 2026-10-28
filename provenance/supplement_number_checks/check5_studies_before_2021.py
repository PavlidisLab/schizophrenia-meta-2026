#!/usr/bin/env python3
"""Check 5: of the 80 base studies, how many were published before 2021, never pooled in a meta-analysis, and never reanalysed.
Year of publication: the four-digit year in the study label ('Hakak 2001'), which equals the lineage-node year wherever both exist (0 differences in 71 studies). fig01_c_lineage_nodes_80.tsv has no row (by study_id = node id) for 9 of the 80 studies, so a merge on it drops them.
never pooled = ever_pooled == no; never reanalysed (of those) = ever_reanalysed == no (equals ever_reused_any == no in this table).
Usage: check5_studies_before_2021.py FIG01_D FIG01_C OUT"""
import sys, hashlib, pandas as pd
D, C, OUT = sys.argv[1:4]; sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
d = pd.read_csv(D, sep="\t", comment="#", dtype=str); d["year_from_label"] = d.label.str.extract(r"(\d{4})")[0].astype(int); n = pd.read_csv(C, sep="\t", comment="#", dtype=str)
m = d.merge(n[["id", "year"]], left_on="study_id", right_on="id", how="left"); m["year_from_node"] = pd.to_numeric(m.year, errors="coerce")
def counts(y):
    b = m[y < 2021]; nev = b[b.ever_pooled == "no"]; return len(b), len(nev), int((nev.ever_reanalysed == "no").sum())
rows = [dict(year_source="label year (all 80)", before_2021=counts(m.year_from_label)[0], never_pooled=counts(m.year_from_label)[1], never_reanalysed=counts(m.year_from_label)[2], studies_without_year=0),
        dict(year_source="lineage-node year (merge on study_id = node id)", before_2021=counts(m.year_from_node)[0], never_pooled=counts(m.year_from_node)[1], never_reanalysed=counts(m.year_from_node)[2], studies_without_year=int(m.year_from_node.isna().sum()))]
R = pd.DataFrame(rows); open(OUT, "w").write(f"# check5_studies_before_2021.py; fig01_d_neglected_datasets_80.tsv {sha(D)}; fig01_c_lineage_nodes_80.tsv {sha(C)}\n# studies with no node-year row: " + "; ".join(f"{r.label} ({r.year_from_label}; pooled {r.ever_pooled}, reanalysed {r.ever_reanalysed})" for r in m[m.year_from_node.isna()].itertuples()) + "\n"); R.to_csv(OUT, sep="\t", index=False, mode="a"); print(R.to_string()); print(open(OUT).read().split("\n")[1])
