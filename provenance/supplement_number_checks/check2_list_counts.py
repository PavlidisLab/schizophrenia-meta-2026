#!/usr/bin/env python3
"""Check 2: stated and matching list counts. Reads the delivered a1_listcount_candidates.tsv rows and recomputes stated/match from two analyses/summary.tsv versions (expected_n filled = count stated in the paper; matches_paper = yes).
summary 9ef302baa9a0 is the file behind main_text_values M4 (668 lists, 328 studies); d189eb75c687 is the file behind a1_listcount_candidates.tsv; the 2026-10-08 file behind workflow_counts.tsv is not on disk.
Usage: check2_list_counts.py A1_TSV SUMMARY_9EF SUMMARY_D189 OVERLAY_JSON OUT"""
import sys, json, hashlib, pandas as pd
A1, S1, S2, OV, OUT = sys.argv[1:6]; sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
N = json.load(open(OV))["nodes"]; cat = lambda n: str(n.get("category"))
c665 = {str(n["pmid"]) for n in N if cat(n) not in ("2", "review", "1m") and not str(n.get("source")).startswith("Cited source") and n.get("pmid")}
rows = []
a1 = pd.read_csv(A1, sep="\t", comment="#"); r = a1[(a1.source == "lineage_srep") & (a1.list_status == "extracted+gemma+reproduced") & (a1.role == "any")]
for _, x in r.iterrows(): rows.append(dict(source=f"a1_listcount_candidates.tsv, {x.population}", population=x.population, lists_stated=int(x.stated), lists_match=int(x.match)))
for name, p in (("summary 9ef302baa9a0 (M4)", S1), ("summary d189eb75c687 (a1 input)", S2)):
    S = pd.read_csv(p, sep="\t", comment="#", dtype=str); E = S[S.status == "extracted"].copy(); E["pmid"] = E.list.str.split("__").str[0]
    for lab, sub in (("all extracted lists", E), ("studies matched by PMID to the 665 corpus", E[E.pmid.isin(c665)])):
        st = sub[sub.expected_n.notna() & (sub.expected_n.str.strip() != "")]; rows.append(dict(source=name, population=f"{lab}: {len(sub)} lists, {sub.pmid.nunique()} studies", lists_stated=len(st), lists_match=int((st.matches_paper == "yes").sum())))
rows.append(dict(source="provenance/sources/workflow_counts.tsv (built 2026-10-08)", population="all extracted lists of analyses/summary.tsv at that date (708 rows with status extracted)", lists_stated=460, lists_match=401))
R = pd.DataFrame(rows); open(OUT, "w").write(f"# check2_list_counts.py; a1 {sha(A1)}; summary_9ef302baa9a0 {sha(S1)}; summary_d189eb75c687 {sha(S2)}; lineage overlay {sha(OV)}\n"); R.to_csv(OUT, sep="\t", index=False, mode="a"); print(R.to_string())
