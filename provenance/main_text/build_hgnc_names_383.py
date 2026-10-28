#!/usr/bin/env python3
"""Approved HGNC name and id for the 383 recurrent genes (RA, 2026-10-19), from the local HGNC complete set."""
import csv, hashlib, os, datetime
H = os.path.expanduser("~")
HG = f"{H}/Dev/schizophrenia-transcriptome-replication/data/reference/hgnc_complete_set.txt"
REC = f"{H}/out/recurrence_genes.tsv"
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
rows = [r for r in csv.DictReader((l for l in open(REC) if not l.startswith("#")), delimiter="\t") if r["scenario"] == "published, brain"]
S = [r["gene"] for r in rows if float(r["q_bh"]) < 0.1]; assert len(S) == 383
appr, prev, alias = {}, {}, {}
for r in csv.DictReader(open(HG), delimiter="\t"):
    if r["status"] != "Approved": continue
    appr[r["symbol"]] = r
    for s in r["prev_symbol"].split("|"):
        if s: prev.setdefault(s, []).append(r["symbol"])
    for s in r["alias_symbol"].split("|"):
        if s: alias.setdefault(s, []).append(r["symbol"])
mt = datetime.date.fromtimestamp(os.path.getmtime(HG)).isoformat()
hd = ["# figure: Fig. 2c interactive version, hover card (gene full names)",
 "# script: fig2c_interactive_data_2026-10-19/build_hgnc_names_383.py, run %s" % datetime.date.today(),
 f"# input: data/reference/hgnc_complete_set.txt sha256 {sha(HG)}; file modification date {mt} (the download date was not recorded in the file; the modification date is the date the file was written locally)",
 f"# input: recurrence_on_80_studies/out/recurrence_genes.tsv sha256 {sha(REC)} (scenario 'published, brain', q_bh < 0.1 = 383 genes)",
 "# match: approved_symbol = exact HGNC approved symbol; previous_symbol / alias_symbol = the gene symbol is not currently approved and maps to the one gene named in hgnc_symbol_now; ambiguous = several candidates, none chosen; not_found = no HGNC record"]
n = {"approved_symbol": 0}
with open("hgnc_names_383.tsv", "w") as f:
    f.write("\n".join(hd) + "\n" + "\t".join(["gene", "match", "hgnc_symbol_now", "approved_name", "hgnc_id", "locus_type"]) + "\n")
    for g in S:
        if g in appr: m, s = "approved_symbol", g
        elif len(prev.get(g, [])) == 1: m, s = "previous_symbol", prev[g][0]
        elif len(alias.get(g, [])) == 1 and g not in prev: m, s = "alias_symbol", alias[g][0]
        elif g in prev or g in alias: m, s = "ambiguous", ""
        else: m, s = "not_found", ""
        n[m] = n.get(m, 0) + 1
        r = appr.get(s, {})
        f.write("\t".join([g, m, s, r.get("name", ""), r.get("hgnc_id", ""), r.get("locus_type", "")]) + "\n")
print(n)
