#!/usr/bin/env python3
"""GTEx v8 brain: association of every gene with total ischemic time and with the Hardy death scale, region by region, against the control pH axis and the Dachet 24 h change.

Inputs (open access, GTEx portal): GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz, GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt (SMTSISCH, SMRIN, SMTSD, SMAFRZE),
GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt (AGE band, SEX, DTHHRDY). analyses/ph_and_confound_analyses/results/combined_controls_pH.tsv; Dachet 24 h change is recomputed from GSE153669_RAW.tar.

Per region (13 brain regions, RNA-seq samples with an ischemic time), per gene with median TPM >= 1 in that region: OLS of log2(TPM + 1) on ischemic time (hours), age (band midpoint), sex and Hardy
scale as a factor (dummies for 0, 1, 3, 4 against 2); the t statistic of ischemic time is the gene's ischemic-time statistic. For the Hardy scale: 'slow death' (4) against 'fast natural death' (2), same covariates with
ischemic time. Outputs: Spearman rho of each statistic with the pH Z, the PMI Z and the Dachet change; correlation of the ischemic-time statistics between regions.

    python3 analyses/gtex_ischemic_time/gtex_ischemic_time.py GTEX_TPM.gct.gz ATTR.txt SUBJ.txt GSE153669_RAW.tar
"""
import csv, io, re, sys, tarfile, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"

attr = pd.read_csv(sys.argv[2], sep="\t", usecols=["SAMPID", "SMTSD", "SMTS", "SMAFRZE", "SMTSISCH", "SMRIN"])
subj = pd.read_csv(sys.argv[3], sep="\t")
brain = attr[(attr.SMTS == "Brain") & (attr.SMAFRZE == "RNASEQ") & attr.SMTSISCH.notna()].copy()
brain["SUBJID"] = brain.SAMPID.str.split("-").str[:2].str.join("-")
brain = brain.merge(subj, on="SUBJID")
brain["age"] = brain.AGE.map(lambda a: np.mean([int(x) for x in a.split("-")]))
brain["isch_h"] = brain.SMTSISCH / 60.0
brain["male"] = (brain.SEX == 1).astype(float)
hdr = pd.read_csv(sys.argv[1], sep="\t", skiprows=2, nrows=0).columns
use = [s for s in brain.SAMPID if s in hdr]
brain = brain[brain.SAMPID.isin(use)]
tpm = pd.read_csv(sys.argv[1], sep="\t", skiprows=2, usecols=["Name", "Description"] + use, dtype={c: np.float32 for c in use})
print("brain samples", len(use), "genes", len(tpm))
sym = tpm.Description.values
X = np.log2(tpm[use].values + 1).astype(np.float32)
idx = {s: i for i, s in enumerate(use)}


def axis(factor):
    d = {}
    for l in open(ROOT / "analyses/ph_and_confound_analyses/results/combined_controls_pH.tsv"):
        if l.startswith('#') or l.startswith('gene\t'): continue
        f = l.rstrip('\n').split('\t')
        if f[1] == factor:
            try: d[f[0]] = float(f[3])
            except ValueError: pass
    return d


pH, PMI = axis("pH"), axis("PMI")
# Dachet 24 h change
def rows_of(data):
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
D = {}
with tarfile.open(sys.argv[4]) as tf:
    for m in tf.getmembers():
        t = re.search(r'_(\d+)H\.xlsx$', m.name)
        if t: D[int(t.group(1))] = rows_of(tf.extractfile(m).read())
dachet = {g: np.log2(D[24][g] + 1) - np.log2(D[0][g] + 1) for g in set(D[0]) & set(D[24])}


def ols_t(Y, M, col):
    """t statistic of column `col` of design M for every row of Y (genes x samples)."""
    Q, R = np.linalg.qr(M)
    beta = np.linalg.solve(R, Q.T @ Y.T)                  # p x genes
    res = Y.T - M @ beta
    dof = M.shape[0] - M.shape[1]
    s2 = (res ** 2).sum(0) / dof
    cov_diag = np.diag(np.linalg.inv(R.T @ R))[col]
    return beta[col] / np.sqrt(s2 * cov_diag), beta[col]


results, tmaps, hmaps = [], {}, {}
for region, sub in brain.groupby("SMTSD"):
    cols = [idx[s] for s in sub.SAMPID]
    Y = X[:, cols].astype(np.float64)
    keep = np.median(tpm[sub.SAMPID].values, axis=1) >= 1
    h = sub.DTHHRDY.values
    dummies = np.column_stack([(h == k).astype(float) for k in (0, 1, 3, 4)])
    M = np.column_stack([np.ones(len(sub)), sub.isch_h.values, sub.age.values, sub.male.values, dummies])
    t_isch, b_isch = ols_t(Y[keep], M, 1)
    # Hardy: slow (4) vs fast natural (2) among those two groups plus ischemic time, age, sex
    sel = np.isin(h, (2, 4))
    M2 = np.column_stack([np.ones(sel.sum()), (h[sel] == 4).astype(float), sub.isch_h.values[sel], sub.age.values[sel], sub.male.values[sel]])
    t_hardy, _ = ols_t(Y[keep][:, sel], M2, 1)
    g = sym[keep]
    tmaps[region] = pd.Series(t_isch, index=g).groupby(level=0).first()
    hmaps[region] = pd.Series(t_hardy, index=g).groupby(level=0).first()
    for nm, tv in (("ischemic time", t_isch), ("slow vs fast natural death", t_hardy)):
        ser = pd.Series(tv, index=g).groupby(level=0).first()
        row = {"region": region, "samples": len(sub), "genes": int(keep.sum()), "statistic": nm}
        for an, ad in (("pH Z", pH), ("PMI Z", PMI), ("Dachet 24 h change", dachet)):
            k = [x for x in ser.index if x in ad]
            row[f"rho_{an}"] = round(spearmanr(ser[k].values, [ad[x] for x in k])[0], 3)
            row[f"n_{an}"] = len(k)
        results.append(row)
    print(region, len(sub), round(results[-2]["rho_pH Z"], 3), round(results[-2]["rho_Dachet 24 h change"], 3))
pd.DataFrame(results).to_csv(HERE / "results_by_region.tsv", sep="\t", index=False)
reg = list(tmaps)
M_ = pd.DataFrame({r: tmaps[r] for r in reg})
C = M_.corr(method="spearman")
C.round(3).to_csv(HERE / "ischemic_t_between_regions.tsv", sep="\t")
off = C.values[np.triu_indices(len(reg), 1)]
pd.DataFrame({r: hmaps[r] for r in reg}).round(4).to_csv(HERE / "hardy_t_by_region.tsv.gz", sep="\t")
pd.DataFrame({r: tmaps[r] for r in reg}).round(4).to_csv(HERE / "ischemic_t_by_region.tsv.gz", sep="\t")
Ch = pd.DataFrame({r: hmaps[r] for r in reg}).corr(method="spearman")
Ch.round(3).to_csv(HERE / "hardy_t_between_regions.tsv", sep="\t")
offh = Ch.values[np.triu_indices(len(reg), 1)]
print("between-region rho of slow-vs-fast-death t: median", round(float(np.median(offh)), 3), "range", round(float(offh.min()), 3), "to", round(float(offh.max()), 3))
print("between-region rho of ischemic-time t: median", round(float(np.median(off)), 3), "range", round(float(off.min()), 3), "to", round(float(off.max()), 3))
