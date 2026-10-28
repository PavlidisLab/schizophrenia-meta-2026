#!/usr/bin/env python3
"""fig07 panel c, v2: signed pH score by the pH handling call from the full-text re-read (ph_agonal_handling_64.tsv).
The 2026-10-09 files were built inline; this script reproduces them (checked against the 2026-10-09 inputs: 64 of 64 rows identical).
Reads figure_data/fig07/fig07_c_ph_score_by_adjustment.tsv (export_fig02_fig07.py) and
analyses/ph_agonal_handling_by_study/ph_agonal_handling_64.tsv. No model calls.
    python3 analyses/figure_data_exports/export_fig07_handling_v2.py"""
import csv, hashlib, os
import numpy as np
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ADJ = os.path.join(REPO, "figure_data/fig07/fig07_c_ph_score_by_adjustment.tsv")
PH = os.path.join(REPO, "analyses/ph_agonal_handling_by_study/ph_agonal_handling_64.tsv")
OUT = os.path.join(REPO, "figure_data/fig07")
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
rd = lambda p: list(csv.DictReader((l for l in open(p, encoding="utf8") if not l.startswith("#")), delimiter="\t"))
ph = {r["pmid"]: r for r in rd(PH)}
rows = []
for a in rd(ADJ):
    p = ph.get(a["study_id"])
    if p:
        rows.append({"study_id": a["study_id"], "label": a["label"], "ph_handling_h2_old": a["ph_handling"], "signed_ph_score": a["signed_ph_score"],
                     "n_genes": a["n_genes"], "n_scored": a["n_scored"], "scored_ge10": a["scored_ge10"], "ph_status_in_h2": a["ph_status_in_h2"],
                     "pmid": p["pmid"], "ph_call": p["ph_call"], "ph_basis": p["ph_basis"], "h2_cov_ph": p["h2_cov_ph"]})
cols = list(rows[0])
with open(os.path.join(OUT, "fig07_c_ph_score_by_handling_v2.tsv"), "w") as fh:
    fh.write("# figure: fig07 panel c (v2, pH handling re-read from full text, 2026-10-09; re-exported 2026-10-11)\n")
    fh.write(f"# script: analyses/figure_data_exports/export_fig07_handling_v2.py; signed_ph_score as in fig07_c_ph_score_by_adjustment.tsv sha256 {sha(ADJ)}; ph_call from analyses/ph_agonal_handling_by_study/ph_agonal_handling_64.tsv sha256 {sha(PH)}\n")
    fh.write(f"# {len(rows)} primary brain studies; ph_handling_h2_old is the earlier h2 classification\n")
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
by = {}
for r in rows:
    by.setdefault(r["ph_call"], []).append(r)
with open(os.path.join(OUT, "fig07_c_ph_score_by_handling_v2_summary.tsv"), "w") as fh:
    fh.write("# figure: fig07 panel c (v2) summary; median and quartiles of signed_ph_score over studies with scored_ge10 = yes; n_all = all studies with the call\n")
    fh.write("ph_call\tn\tmedian\tq1\tq3\tn_all\n")
    for c in sorted(by):
        s = [float(r["signed_ph_score"]) for r in by[c] if r["scored_ge10"] == "yes" and r["signed_ph_score"] not in ("", "NA")]
        f = lambda x: "" if not s else repr(float(np.round(x, 3)))
        fh.write(f"{c}\t{len(s)}\t{f(np.median(s) if s else 0)}\t{f(np.quantile(s, .25) if s else 0)}\t{f(np.quantile(s, .75) if s else 0)}\t{len(by[c])}\n")
print(len(rows), "rows;", len(by), "calls")
