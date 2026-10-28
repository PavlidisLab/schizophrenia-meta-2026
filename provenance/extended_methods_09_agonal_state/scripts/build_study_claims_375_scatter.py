#!/usr/bin/env python3
"""One row per published brain study of the analysis population (62 units on the 515-study corpus; 64 before 2026-10-11): what the study claims, how it handled pH, and how its own hit list lines up with the dying-process axes.

Columns: label, pmid, year, cluster (cohort cluster of recurrence.py), tissue lane, claimed themes with the verbatim quote each was assigned from (corpus/themes.json, Sonnet read of title and abstract),
pH handling (h2/h2_studies.tsv cov_ph, with the hand corrections), n genes with a direction,
  ph_score       mean over the list's genes of -(SZ sign) x control pH Z   (positive = the list moves with low pH the way it moves in SZ; ph_by_adjustment.py convention)
  agonal_score   mean over the list's genes of (SZ sign) x GTEx agonal t (slow vs fast natural death, cortical mean; positive = the list moves with the agonal state the way it moves in SZ)
  dachet_share   share of the list's genes whose SZ direction equals their direction in dying fresh cortex at 24 h (GSE153669)
  agonal_share   share of the list's genes whose SZ direction equals the sign of the agonal t.
Chance level of the two shares is the share expected from the gene-level base rates; each is also given against a shuffled-sign null (200 shuffles of the SZ signs within the study's list; z and one-sided p).

    python3 analyses/study_claims_2026-10-08/build_study_claims.py /path/to/GSE153669_RAW.tar
"""
import csv, io, json, os, re, sys, tarfile, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET

import os as _os
HERE = Path(_os.environ["OUTDIR"]); HERE.mkdir(parents=True, exist_ok=True)   # copy: output folder and analyses dir from the environment
GL = Path(_os.environ["ANALYSES"])
sys.path.insert(0, str(GL))
import recurrence as R  # noqa: E402

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
rng = np.random.default_rng(20261008)


def rd(p):
    return list(csv.DictReader((l for l in open(GL / p, encoding="utf-8") if not l.startswith("#")), delimiter="\t"))


def xl_rows(data):
    z = zipfile.ZipFile(io.BytesIO(data))
    ss = [''.join(t.text or '' for t in si.iter(NS + 't')) for si in ET.parse(z.open('xl/sharedStrings.xml')).getroot()] if 'xl/sharedStrings.xml' in z.namelist() else []
    rels = {r.get('Id'): r.get('Target') for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
    s = ET.fromstring(z.read('xl/workbook.xml')).find(NS + 'sheets')[0]
    tgt = rels[s.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')].lstrip('/'); tgt = tgt if tgt.startswith('xl/') else 'xl/' + tgt
    o = {}
    for r in ET.parse(z.open(tgt)).getroot().iter(NS + 'row'):
        row = {}
        for c in r.findall(NS + 'c'):
            v = c.find(NS + 'v')
            if v is not None: row[re.match(r'[A-Z]+', c.get('r')).group()] = ss[int(v.text)] if c.get('t') == 's' else v.text
        if 'A' in row and 'B' in row:
            try: o[row['A']] = float(row['B'])
            except ValueError: pass
    return o


E_ = {}
with (tarfile.open(sys.argv[1]) if len(sys.argv) > 1 else __import__("contextlib").nullcontext(None)) as tf:   # copy: GSE153669_RAW.tar is not available; without it the Dachet columns stay blank
  if tf is not None:
    for m in tf.getmembers():
        t = re.search(r'_(\d+)H\.xlsx$', m.name)
        if t and int(t.group(1)) in (0, 24): E_[int(t.group(1))] = xl_rows(tf.extractfile(m).read())
import csv as _c
_sc=os.environ["DACHET_SCATTER"]   # copy: 24 h change read from the stored gene scatter (fig08_a_dachet_gene_scatter.tsv, change_24h = log2(x24+1)-log2(x0+1)) instead of the tar
dachet={r["gene"]:float(r["change_24h"]) for r in _c.DictReader((l for l in open(_sc) if not l.startswith("#")),delimiter="\t") if r["change_24h"]!=""}
z = {r["gene"]: float(r["Z"]) for r in rd("confounds/results/combined_controls_pH.tsv") if r["factor"] == "pH" and r["Z"]}
H = pd.read_csv(GL / "gtex_ischemia_2026-10-08/hardy_t_by_region.tsv.gz", sep="\t", index_col=0)
agonal = H[["Brain - Frontal Cortex (BA9)", "Brain - Anterior cingulate cortex (BA24)", "Brain - Cortex"]].mean(axis=1, skipna=False).dropna().to_dict()
themes = json.load(open(GL.parent / "corpus/themes.json"))["studies"]
nodes = {str(n["pmid"]): n for n in json.load(open(GL.parent / "corpus/lineage.json"))["nodes"] if n.get("pmid")}
h2 = {x["pmid"]: x for x in rd("h2/h2_studies.tsv")}
units = [x for x in rd("recurrence_units.tsv") if x["scenario"] == "published, brain" and x["modelled"].lower().startswith("y")]
pl = json.load(open(GL / "primary_lists.json"))["studies"]


def share_stat(signs, vals):
    s = np.array(signs); v = np.array(vals)
    obs = float(np.mean(np.sign(v) == s))
    null = np.array([np.mean(np.sign(v) == rng.permutation(s)) for _ in range(200)])
    return obs, float((obs - null.mean()) / null.std()) if null.std() > 0 else np.nan, float((1 + (null >= obs).sum()) / 201)


rows = []
for u in units:
    p = u["representative"]; genes = {}
    for ln in pl[p]["primary"]:
        for g, d in R.read_list(ln).items(): genes[g] = d
    sg = {g: {"up": 1, "down": -1}.get(d, 0) for g, d in genes.items()}
    node = nodes.get(p, {}); th = themes.get(node.get("id", p), {}).get("themes", [])
    gs_ph = [(sg[g], z[g]) for g in sg if sg[g] and g in z]
    gs_ag = [(sg[g], agonal[g]) for g in sg if sg[g] and g in agonal]
    gs_da = [(sg[g], dachet[g]) for g in sg if sg[g] and g in dachet]
    row = {"label": u["label"], "pmid": p, "year": node.get("year", ""), "cluster": u["cluster"], "lane": u["lane"],
           "themes": "; ".join(x["theme"] for x in th), "theme_quotes": " | ".join(f"[{x['theme']}] {x['quote']}" for x in th),
           "ph_status": h2.get(p, {}).get("cov_ph", "") or "not_mentioned", "n_signed_genes": sum(1 for g in sg if sg[g])}
    row["ph_score"] = round(float(np.mean([-s * zz for s, zz in gs_ph])), 3) if len(gs_ph) >= 10 else ""
    row["agonal_score"] = round(float(np.mean([s * a for s, a in gs_ag])), 3) if len(gs_ag) >= 10 else ""
    for nm, gl in (("agonal", gs_ag), ("dachet", gs_da)):
        if len(gl) >= 10:
            o, zz, pp = share_stat([s for s, _ in gl], [v for _, v in gl]); row[f"{nm}_share"] = round(o, 3); row[f"{nm}_share_z"] = round(zz, 2); row[f"{nm}_share_p"] = round(pp, 3)
        else:
            row[f"{nm}_share"] = row[f"{nm}_share_z"] = row[f"{nm}_share_p"] = ""
    rows.append(row)
with open(HERE / "study_claims.tsv", "w", newline="", encoding="utf-8") as fh:
    import hashlib as _h
    _pl = _h.sha256(open(GL / "primary_lists.json", "rb").read()).hexdigest()[:12]
    fh.write(f"# build_study_claims.py, run 2026-10-11; {len(rows)} published brain study units (recurrence_units.tsv, one representative list per unit; primary_lists.json sha256 {_pl}); axes: control pH Z, GTEx agonal t (cortical mean), Dachet 24 h change; shares tested against 200 sign shuffles within the list\n")
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
df = pd.DataFrame(rows)
for c in ("ph_score", "agonal_score", "agonal_share", "dachet_share"): df[c] = pd.to_numeric(df[c], errors="coerce")
t = df.dropna(subset=["agonal_score"])
print(f"{len(df)} studies; {len(t)} with >= 10 scored genes")
print(f"ph_score median {t.ph_score.median():.2f}, positive in {(t.ph_score > 0).sum()}; agonal_score median {t.agonal_score.median():.3f}, positive in {(t.agonal_score > 0).sum()}; agonal_share median {t.agonal_share.median():.3f} (>0.5 in {(t.agonal_share > 0.5).sum()}); dachet_share median {t.dachet_share.median():.3f} (>0.5 in {(t.dachet_share > 0.5).sum()})")
print("studies with agonal_share_p < 0.05:", int((pd.to_numeric(t.agonal_share_p) < 0.05).sum()), "; dachet_share_p < 0.05:", int((pd.to_numeric(t.dachet_share_p) < 0.05).sum()))
