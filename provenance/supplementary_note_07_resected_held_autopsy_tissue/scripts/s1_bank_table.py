#!/usr/bin/env python3
"""Task 1: bank-level overlap of the 80 base studies with the Liharska 2025 PM banks (HBTRC, NYBB, UMBEB).
Reads lineage.json, base80_studies.tsv, recurrence_units.tsv, recurrence_genes.tsv (B80 out), donor_overlap tables/INDEX.tsv.
Writes bank_evidence_80.tsv. Evidence strings are curated from the study text in data/raw/<pmid>/ (see README); `python3 -I s1_bank_table.py`."""
import csv, json, os, re, glob, hashlib, sys
sys.dont_write_bytecode = True
R = os.path.expanduser("schizophrenia-transcriptome-replication"); A = os.path.expanduser("agent-mail/sz-transcriptome")
B80 = f"{A}/recurrence_on_80_studies"; OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
def tsv(p): return list(csv.DictReader((l for l in open(p, encoding="utf8") if not l.startswith("#")), delimiter="\t"))
L = json.load(open(f"{R}/corpus/lineage.json")); node = {x["pmid"]: x for x in L["nodes"]}
base = tsv(f"{B80}/base80_studies.tsv")
units = {r["representative"]: r for r in tsv(f"{B80}/out/recurrence_units.tsv") if r["scenario"] == "published, brain"}
genes = [r for r in tsv(f"{B80}/out/recurrence_genes.tsv") if r["scenario"] == "published, brain" and float(r["q_bh"]) < 0.1]
assert len(genes) == 383
cnt = {}
for g in genes:
    for s in g["studies"].split(";"): cnt[s] = cnt.get(s, 0) + 1
idx = {r["study"]: r for r in tsv(f"{R}/analyses/donor_overlap/tables/INDEX.tsv")}
# curated bank evidence (study text quoted in lineage.json / data/raw/<pmid>/); PM bank: HBTRC | UMBEB | NYBB
EV = {
 "24702465": ("HBTRC", "Harvard Brain Tissue Resource Center (HBTRC), McLean Hospital (lineage.json parent; data/raw/24702465 abstract/GEO text); same 16 subjects as the PV series (donor_overlap report)", "McLean 18 (GSE37981 table)"),
 "24628518": ("HBTRC", "HBTRC at McLean Hospital; same subjects as the pyramidal-neuron study (lineage.json parent; data/raw/24628518)", "McLean 16 (GSE46509 table)"),
 "34440415": ("HBTRC", "'international brain banks (Harvard Brain Tissue Resource Center, MRC Edinburgh Brain Bank, the University of Mississippi Medical Center and the London Brain Bank)' (data/raw/34440415/pmc_fulltext.txt)", "none"),
 "38448582": ("HBTRC", "191 donors from the Harvard Brain Tissue Resource Center (McLean) / NIH NeuroBioBank (lineage.json lane_by; data/raw/38448582)", "McLean 191"),
 "38781388": ("HBTRC", "McLean/Harvard Brain Tissue Resource Center cohort, 24 SZ + 24 control (+ Mount Sinai 41 SZ + 51 control) (lineage.json bank_names; data/raw/38781388)", "McLean 48, MSSM 92"),
 "31288836": ("HBTRC+UMBEB", "'NIH NeuroBioBank (The Harvard Brain Tissue Resource Center ...; ... University of Miami Brain Endowment Bank ...)' (data/raw/31288836); donor table McLean 9, Miami 2, CT Lab 37, HSB 13", "McLean 9, Miami 2 (controls)"),
 "19255580": ("HBTRC (reused data)", "second cohort = DLPFC BA9, 16 SZ vs 26 controls, 'Harvard Brain Bank' = Glatt 2005 raw data re-analysed (lineage.json tissue_n, lane_by; list hbb_ba9); first cohort Charing Cross Hospital (UK)", "none for Maycox; Glatt 2005 table McLean 46 (Glatt is not a base study)"),
 "17822540": ("HBTRC", "three banks: Stanley (30), Maudsley (48), Harvard Brain Tissue Resource Centre (32) (lineage.json parent; data/raw/17822540)", "none"),
 "31053723": ("HBTRC+UMBEB", "NIH NeuroBioBank Pittsburgh, Harvard BTRC, Sepulveda, Miami (lineage.json tissue_n; data/raw/31053723). Derived list holds 0 genes, so the study is dropped from the 72 lists", "McLean 24, Maryland 34, Pittsburgh 1, Sepulveda 24 (no Miami rows)"),
}
rows = []
for b in base:
    p = b["pmid"]; n = node.get(p, {}); u = units.get(p)
    ev = EV.get(p)
    d = ""
    key = {"24702465": "pietersen2014_GSE37981", "24628518": "pietersen2014_GSE46509", "38448582": "ling2024", "38781388": "ruzicka2024", "31288836": "mendizabal2019", "31053723": "pai2019"}.get(p)
    rows.append(dict(pmid=p, study=b["label"], lineage_lane=n.get("lane", ""), lineage_banks="+".join(n.get("banks") or []),
        shares_PM_bank="yes" if ev else "no bank named in the text checked",
        pm_bank=ev[0] if ev else "", bank_evidence=ev[1] if ev else "", donor_table_in_repo=ev[2] if ev else "",
        cohort_cluster=(u["cluster"] if u and u["cluster"] else ("left out of cohort tier (multi-bank)" if u else "")),
        contributes_list_among_72="no (empty list, dropped)" if b["empty_list_dropped"] == "True" else "yes",
        list_genes=b["list_genes"], genes_of_383_named=cnt.get(b["label"], 0) if u else 0))
w = ["# script: scripts/s1_bank_table.py, run 2026-10-18",
     f"# inputs (sha256 first 12): lineage.json {sha(R+'/corpus/lineage.json')}; base80_studies.tsv {sha(B80+'/base80_studies.tsv')}; recurrence_units.tsv {sha(B80+'/out/recurrence_units.tsv')}; recurrence_genes.tsv {sha(B80+'/out/recurrence_genes.tsv')}",
     "# genes_of_383_named = genes of the 383 (published, brain, q_bh < 0.1) whose 'studies' column names the study; 0 also for studies that are not the representative of their same-data group",
     "# shares_PM_bank is 'no bank named in the text checked' where the available text (several studies: abstract only) names none of Harvard/HBTRC/McLean, Columbia/NYBB, Miami/UMBEB; it is not proof of absence"]
with open(f"{OUT}/bank_evidence_80.tsv", "w") as fh:
    fh.write("\n".join(w) + "\n"); wr = csv.DictWriter(fh, list(rows[0]), delimiter="\t"); wr.writeheader(); wr.writerows(rows)
sh = [r for r in rows if r["shares_PM_bank"] == "yes"]
print(len(sh), "sharing;", sum(r["contributes_list_among_72"] == "yes" for r in sh), "with a list;", [(r["study"], r["genes_of_383_named"]) for r in sh])
# text coverage of the studies that name no bank
cov = []
for b in base:
    fs = [f for f in glob.glob(f"{R}/data/raw/{b['pmid']}/**/*", recursive=True) if os.path.isfile(f) and re.search(r"\.(txt|xml|html)$", f)]
    cov.append((b["label"], max([os.path.getsize(f) for f in fs] or [0])))
print("studies with no text file or <20 kB text:", sum(c[1] < 20000 for c in cov), "of", len(cov))
open(f"{OUT}/logs/text_coverage.tsv", "w").write("# max text-file size (bytes) under data/raw/<pmid>/ per base study (proxy for whether Methods text was readable)\nstudy\tmax_text_bytes\n" + "\n".join(f"{a}\t{c}" for a, c in cov) + "\n")
