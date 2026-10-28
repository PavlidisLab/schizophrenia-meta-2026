#!/usr/bin/env python3
"""Long table of literature coexpression modules that have a gene list and a stated direction in schizophrenia, with each module's discovery bank lineage.

Sources (analysis repository, REPO):
  - analyses/coexpression_modules_2026-10-08/modules_batch{1-4}.tsv: module gene lists read from supplements of 79 postmortem-brain coexpression papers; kept: direction in {higher, lower} (mixed and none_stated dropped) and present in modules_axes_all.tsv (the set Supplementary Note 5 scores on the axes)
  - Gandal 2018 (PMID 30545856) Table S5 (data/raw/30545856/zotero/aat8127_table_s5.xlsx): geneModules and Module-DiseaseAssociation (Network == gene); direction = sign of the SCZ t; stated_significant = SCZ fdr < 0.05
  - Kang and Oldham 2026 CoPA Table S14 (side-studies/coexpression_claims_survey/copa_kang_oldham_2026/biorxiv_2026.08.25.747130_media-14.xlsx): modules called Lower/Higher in SCZ, per dataset (CMC, SZBD); stated_significant = yes (listed as dysregulated)
Discovery lineage: bank labels of side-studies/coexpression_claims_survey/bank_overlap_claims.tsv (rules in bank_overlap.py), union over the paper's claims; Gandal = its claim row; CoPA: CMC -> CMC lineage, SZBD -> McLean/Harvard.
Output: modules_long.tsv (module x gene), modules_index.tsv (one row per module).
    python3 build_modules.py REPO OUTDIR
"""
import sys, re, hashlib, csv, zipfile
import xml.etree.ElementTree as ET
from pathlib import Path
import pandas as pd
REPO = Path(sys.argv[1]); OUT = Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True)
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
def sheets(path):
    z = zipfile.ZipFile(path)
    ss = [''.join(t.text or '' for t in si.iter(NS + 't')) for si in ET.parse(z.open('xl/sharedStrings.xml')).getroot()]
    rels = {r.get('Id'): r.get('Target') for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
    out = {}
    for s in ET.fromstring(z.read('xl/workbook.xml')).find(NS + 'sheets'):
        tgt = rels[s.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')].lstrip('/'); tgt = tgt if tgt.startswith('xl/') else 'xl/' + tgt
        rows = []
        for r in ET.parse(z.open(tgt)).getroot().iter(NS + 'row'):
            row = {}
            for c in r.findall(NS + 'c'):
                v = c.find(NS + 'v')
                if v is not None: row[re.match(r'[A-Z]+', c.get('r')).group()] = ss[int(v.text)] if c.get('t') == 's' else v.text
            rows.append(row)
        out[s.get('name')] = rows
    return out
inputs = {}
def inp(p): inputs[str(p.relative_to(REPO))] = sha(p)[:12]; return p
# paper -> lineage
bo = pd.read_csv(inp(REPO / "side-studies/coexpression_claims_survey/bank_overlap_claims.tsv"), sep="\t", dtype=str).fillna("")
inv = pd.read_csv(inp(REPO / "side-studies/coexpression_claims_survey/coexpression_inventory.tsv"), sep="\t", comment="#", dtype=str, quoting=csv.QUOTE_NONE).fillna("")
inv_p = inv.drop_duplicates("paper").set_index("paper")["pmid"].to_dict()
paper_lineage = {}
for r in bo.itertuples():
    paper_lineage.setdefault(r.paper, set()).update(r.discovery.split(";"))
pmid_lineage = {}
for p, s in paper_lineage.items():
    pm = inv_p.get(p)
    if pm: pmid_lineage.setdefault(pm, set()).update(s)
# the manuscript's module set: modules_axes_all.tsv (>= 20 symbol genes with a pH axis value; Ensembl-id genes dropped), direction stated
ax = pd.read_csv(inp(REPO / "analyses/coexpression_modules_2026-10-08/modules_axes_all.tsv"), sep="\t", dtype=str)
keep = set(map(tuple, ax[ax.direction.isin(["higher", "lower"])][["pmid", "module"]].values))
mods = []
for i in range(1, 5):
    d = pd.read_csv(inp(REPO / f"analyses/coexpression_modules_2026-10-08/modules_batch{i}.tsv"), sep="\t", dtype=str, quoting=csv.QUOTE_NONE).fillna("")
    d = d[d.direction.str.startswith(("higher", "lower")) & ~d.gene.str.startswith("ENSG")].drop_duplicates(["pmid", "module_id", "gene"])
    for (pmid, paper, mid), g in d.groupby(["pmid", "paper", "module_id"], sort=False):
        if (pmid, mid) not in keep or paper.startswith("Gandal"): continue   # Gandal 2018 enters once, from its Table S5 (below)
        mods.append(dict(source="paper_supplement", pmid=pmid, paper=paper, module=mid, direction=1 if g.direction.iloc[0].startswith("higher") else -1, sz_linked=g.sz_linked.iloc[0], stated_significant="",
                         lineage=";".join(sorted(pmid_lineage.get(pmid, {"unnamed"}))), region=g.region.iloc[0], genes=sorted(set(g.gene) - {""})))
# Gandal
S5 = sheets(inp(REPO / "data/raw/30545856/zotero/aat8127_table_s5.xlsx"))
gm = S5["geneModules"]; h = {v: k for k, v in gm[0].items()}
mg = {}
for r in gm[1:]: mg.setdefault(r.get(h["Module"]), []).append(r.get(h["gene_name"]))
da = S5["Module-DiseaseAssociation"]; hd = {v: k for k, v in da[0].items()}
assoc = {r[hd["ModulName"]]: (float(r[hd["t"]]), float(r[hd["fdr"]])) for r in da[1:] if r.get(hd["Network"]) == "gene" and r.get(hd["Group"]) == "SCZ"}
g_lin = ";".join(sorted(paper_lineage.get("Gandal 2018 (Science, cross-disorder)", {"unnamed"})))
for m, gl in mg.items():
    if m == "geneM0" or m not in assoc: continue
    t, f = assoc[m]
    mods.append(dict(source="Gandal2018_S5", pmid="30545856", paper="Gandal 2018", module=m, direction=1 if t > 0 else -1, sz_linked="yes" if f < 0.05 else "no", stated_significant="yes" if f < 0.05 else "no",
                     lineage=g_lin, region="cortex", genes=sorted(set(x for x in gl if x)), stat=f"SCZ t {t:.2f} fdr {f:.3g}"))
# CoPA
S14 = sheets(inp(REPO / "side-studies/coexpression_claims_survey/copa_kang_oldham_2026/biorxiv_2026.08.25.747130_media-14.xlsx"))
recs = []
for nm, rr in S14.items():
    if len(rr) > 1000 and rr[0].get("A") == "Dataset":
        hh = list(rr[0].values()); kk = list(rr[0].keys())
        for r in rr[1:]:
            d = dict(zip(hh, [r.get(k) for k in kk]))
            if d.get("Disease") == "SCZ": recs.append(d)
C = pd.DataFrame(recs).drop_duplicates()
LIN = {"CMC": "CMC lineage", "SZBD": "McLean/Harvard"}
for (ds, mt, m, dr), g in C.groupby(["Dataset", "Module type", "Module", "Direction"]):
    if dr.startswith("Mixed") or not (dr.startswith("Lower") or dr.startswith("Higher")) or g.Gene.nunique() < 30: continue   # >= 30 genes, as in module_axes.py
    mods.append(dict(source="CoPA_Kang_Oldham_2026_S14", pmid="42800954", paper="Kang and Oldham 2026", module=f"{ds}|{mt}|{m}", direction=1 if dr.startswith("Higher") else -1, sz_linked="yes", stated_significant="yes",
                     lineage=LIN.get(ds, "unnamed"), region="snRNA", genes=sorted(set(g.Gene.dropna()))))
idx = pd.DataFrame([{k: v for k, v in m.items() if k != "genes"} | {"n_genes": len(m["genes"])} for m in mods])
idx.insert(0, "module_uid", [f"M{i:04d}" for i in range(len(idx))])
long = pd.DataFrame([(u, g) for u, m in zip(idx.module_uid, mods) for g in m["genes"]], columns=["module_uid", "gene"])
hdr = ["# build_modules.py (module_retest); inputs sha256 prefixes: " + "; ".join(f"{k} {v}" for k, v in inputs.items()),
       "# one row per module with a gene list and a stated direction in schizophrenia (direction +1 = higher, -1 = lower in SZ, as the source states); lineage = discovery bank lineage labels (bank_overlap.py rules)"]
for name, df in (("modules_index.tsv", idx), ("modules_long.tsv", long)):
    with open(OUT / name, "w") as fh:
        fh.write("\n".join(hdr) + "\n"); df.to_csv(fh, sep="\t", index=False)
print(idx.groupby("source").size(), idx.lineage.value_counts().head(20), idx.n_genes.describe())
