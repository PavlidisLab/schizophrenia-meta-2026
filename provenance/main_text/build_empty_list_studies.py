#!/usr/bin/env python3
"""n_tested for the 8 studies whose derived list is empty at FDR 5% (RA, 2026-10-19). Same definition as n_tested in fig02_b_study_annotation.tsv: n_k of the study's universe in universes.json."""
import csv, json, hashlib, os, datetime
H = os.path.expanduser("~"); B = f"{H}/"
U = B + "ref_B80_inputs/universes.json"; S = B + "base80_studies.tsv"
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
u = json.load(open(U))["lists"]
base = [r for r in csv.DictReader(open(S), delimiter="\t") if r["empty_list_dropped"] == "True"]
assert len(base) == 8
rows = []
for r in base:
    ks = [k for k, v in u.items() if v.get("study") == r["pmid"] and "derived" in k]
    nks = {u[k]["n_k"] for k in ks}
    assert len(nks) == 1 and int(r["universe_n_k"]) in nks, (r["label"], nks)
    rows.append((r["pmid"], r["label"], nks.pop(), "derived_limma_fit_after_filter", 0, r["universe_type"], len(ks)))
hd = ["# figure: Fig. 2c interactive version, grey 'genes tested' bar for the 8 studies with an empty list at FDR 5%",
 "# script: fig2c_interactive_data_2026-10-19/build_empty_list_studies.py, run %s" % datetime.date.today(),
 f"# input: recurrence_on_80_studies/ref_B80_inputs/universes.json sha256 {sha(U)}; recurrence_on_80_studies/base80_studies.tsv sha256 {sha(S)}",
 "# n_tested = n_k of the study's universe as recurrence.py uses it (universes.json, 'derived' entries of the study; the same value as universe_n_k in base80_studies.tsv; asserted equal). It is the union over the study's derived units of every gene in our limma fit after the project expression filter, not the authors' tested set; n_tested_class = derived_limma_fit_after_filter (the 72 listed studies use the classes of fig02_b_study_annotation.tsv)",
 "# n_list_genes = 0: our derived list of the study is empty at FDR 5%, so the study is dropped from the recurrence analysis (72 of the 80 studies contribute a list); n_derived_units = derived units pooled into the universe (Maschietto 2015 and de Baumont 2015 share one GEO series, GSE62191, and the same universe)"]
with open("empty_list_studies_n_tested.tsv", "w") as f:
    f.write("\n".join(hd) + "\n" + "\t".join(["study_id", "label", "n_tested", "n_tested_class", "n_list_genes", "universe_type", "n_derived_units"]) + "\n")
    for r in rows: f.write("\t".join(map(str, r)) + "\n")
print(rows)
