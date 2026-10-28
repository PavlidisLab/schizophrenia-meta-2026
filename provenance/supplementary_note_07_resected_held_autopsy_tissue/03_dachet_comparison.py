#!/usr/bin/env python3
"""Step 3: the Dachet et al. 2021 ex vivo time course (GSE153669; one patient, room temperature, 0 to 24 h) on the same
footing as the Yaqubi contrasts: Spearman rho of the change from 0 h (log2(x + 1), loader copied from
analyses/dachet_live_tissue_2026-10-08/build_report.py) with the pH axis v2 and the agonal axis, over all genes and within the 383 and 107 recurrent genes
(bootstrap over genes, matched null of 2,000 sets as in 02_axes_and_sets.py), and the direct similarity of the Dachet change with the Yaqubi effects.
Usage: python3 03_dachet_comparison.py <out dir> <GSE153669_RAW.tar>      Seed 20261018.
"""
import collections, hashlib, io, os, random, re, sys, tarfile, zipfile, datetime
sys.dont_write_bytecode = True
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from scipy import stats

OUT, RAW = sys.argv[1], sys.argv[2]
AN = os.path.expanduser("<analysis repo>/analyses")
sys.path.insert(0, AN)
A = os.path.expanduser("~/analysis-delivery")
B80 = f"{A}/base80_recurrence_2026-10-16/out"
PHAX, GTEX = f"{A}/ph_axis_v2_2026-10-14/axis_v2_controls.tsv.gz", f"{AN}/gtex_ischemia_2026-10-08/hardy_t_by_region.tsv.gz"
SEED, NPERM, NBOOT = 20261018, 2000, 1000
import extract
extract.HGNC = extract.load_hgnc()
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
HOURS = [0, 1, 2, 4, 8, 12, 24]


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()[:12]


def xlsx_sheet(data):
    z = zipfile.ZipFile(io.BytesIO(data))
    ss = [''.join(t.text or '' for t in si.iter(NS + 't')) for si in ET.parse(z.open('xl/sharedStrings.xml')).getroot()] if 'xl/sharedStrings.xml' in z.namelist() else []
    rels = {r.get('Id'): r.get('Target') for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
    sheets = ET.fromstring(z.read('xl/workbook.xml')).find(NS + 'sheets')
    s = [x for x in sheets][0]
    tgt = rels[s.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')].lstrip('/')
    tgt = tgt if tgt.startswith('xl/') else 'xl/' + tgt
    rows = []
    for r in ET.parse(z.open(tgt)).getroot().iter(NS + 'row'):
        row = {}
        for c in r.findall(NS + 'c'):
            v = c.find(NS + 'v')
            if v is not None:
                row[re.match(r'[A-Z]+', c.get('r')).group()] = ss[int(v.text)] if c.get('t') == 's' else v.text
        rows.append(row)
    return rows


Ex = {}
with tarfile.open(RAW) as tf:
    for m in tf.getmembers():
        t = re.search(r'_(\d+)H\.xlsx$', m.name)
        if t:
            rows = xlsx_sheet(tf.extractfile(m).read())
            Ex[int(t.group(1))] = {r['A']: float(r['B']) for r in rows[1:] if 'A' in r and 'B' in r and re.match(r'^-?[\d.eE+-]+$', r['B'])}
genes = sorted(set.intersection(*[set(Ex[h]) for h in HOURS]))
L = np.array([[np.log2(Ex[h][g] + 1) for h in HOURS] for g in genes])
D = L - L[:, [0]]
HGc = {}
def hg(s):
    if s not in HGc:
        HGc[s] = extract.symbol_to_hgnc(str(s))[0] or None
    return HGc[s]
dch = pd.DataFrame(D, index=[hg(g) for g in genes], columns=HOURS)
dch = dch[dch.index.notna()]
dch = dch[~dch.index.duplicated(keep="first")]
print("Dachet genes mapped:", len(dch))

ph = pd.read_csv(PHAX, sep="\t", comment="#", index_col=0)["z_stouffer"]
Hd = pd.read_csv(GTEX, sep="\t", index_col=0)
ag = Hd[["Brain - Frontal Cortex (BA9)", "Brain - Anterior cingulate cortex (BA24)", "Brain - Cortex"]].mean(axis=1, skipna=False).dropna()
AX = {"pH axis v2": ph, "agonal axis (GTEx slow-minus-fast)": ag}

# recurrent sets and strata (as ad_383.py)
rec = pd.read_csv(f"{B80}/recurrence_genes.tsv", sep="\t", comment="#")
rec = rec[rec.scenario == "published, brain"]
sign, dec, prior = {}, {}, {}
for r in rec.itertuples():
    d = int(r.up) - int(r.down)
    sign[r.gene] = (1 if d > 0 else -1) if d else 0
    dec[r.gene] = r.expr_decile
    prior[r.gene] = float(r.de_prior_hitlists) if not pd.isna(r.de_prior_hitlists) else np.nan
med = float(np.nanmedian(list(prior.values())))
prior = {g: (med if np.isnan(v) else v) for g, v in prior.items()}
S383 = list(rec[rec.q_bh < 0.1].gene)
sets = pd.read_csv(f"{B80}/recurrent_sets.tsv", sep="\t", comment="#")
S107 = list(sets[sets.set == "study_brain_ge10"].gene)
pool = [g for g in sign if sign[g] != 0 and g not in set(S383)]
qc = np.quantile([prior[g] for g in pool], [.2, .4, .6, .8])
st = {g: (dec[g], int(np.searchsorted(qc, prior[g], side="right"))) for g in sign}
by = collections.defaultdict(list)
for g in pool:
    by[st[g]].append(g)
rng = np.random.default_rng(SEED)


def draws(genes_, seed):
    need = collections.Counter(st[g] for g in genes_)
    r = random.Random(seed)
    out = []
    for _ in range(NPERM):
        row = []
        for s, k in need.items():
            c = by[s]
            row += r.sample(c, k) if len(c) >= k else r.choices(c, k=k)
        out.append(row)
    return out


DR = {"383 recurrent (q_bh<0.1)": (S383, draws(S383, SEED)), "107 in >=10 lists": (S107, draws(S107, SEED + 1))}


def rk(x):
    return stats.rankdata(x)


def rho_boot(x, y):
    r = stats.spearmanr(x, y)[0]
    n = len(x)
    b = [np.corrcoef(rk(x[i]), rk(y[i]))[0, 1] for i in (rng.integers(0, n, n) for _ in range(NBOOT))]
    return r, np.quantile(b, [.025, .975])


rows = []
for h in HOURS[1:]:
    e = dch[h]
    for an, ax in AX.items():
        j = pd.concat([e, ax], axis=1, join="inner").dropna()
        r, ci = rho_boot(j.iloc[:, 0].values, j.iloc[:, 1].values)
        rows.append(dict(scope="all genes", hours=h, axis=an, n_genes=len(j), spearman=round(r, 3), boot95_low=round(ci[0], 3), boot95_high=round(ci[1], 3), null_mean="", null_q025="", null_q975=""))
        for sname, (sg, dr) in DR.items():
            sg_ok = [g for g in sg if g in e.index and g in ax.index and not np.isnan(e[g]) and not np.isnan(ax[g])]
            ev, av = e.reindex(sg_ok).values, ax.reindex(sg_ok).values
            r, ci = rho_boot(ev, av)
            nl = []
            for row in dr:
                gg = [g for g in row if g in e.index and g in ax.index]
                gg = [g for g in gg if not np.isnan(e[g]) and not np.isnan(ax[g])]
                nl.append(stats.spearmanr(e.reindex(gg).values, ax.reindex(gg).values)[0])
            nl = np.array(nl)
            rows.append(dict(scope=sname, hours=h, axis=an, n_genes=len(sg_ok), spearman=round(r, 3), boot95_low=round(ci[0], 3), boot95_high=round(ci[1], 3),
                             null_mean=round(nl.mean(), 3), null_q025=round(np.quantile(nl, .025), 3), null_q975=round(np.quantile(nl, .975), 3)))
R = pd.DataFrame(rows)

# direct similarity: Dachet change vs Yaqubi effects
E = pd.read_csv(f"{OUT}/tables/contrast_effects.tsv", sep="\t", comment="#", index_col=0)
rows2 = []
pairs = [(4, "S_6h20C_vs_0h"), (8, "S_6h20C_vs_0h"), (8, "S_24h20C_vs_0h"), (24, "S_24h20C_vs_0h"), (24, "S_36h20C_vs_0h"), (24, "S_6h20C_vs_0h"), (24, "S_24h4C_vs_0h")]
for h, c in pairs:
    for kind in ("lfc", "t"):
        j = pd.concat([dch[h].rename("d"), E[f"{kind}_{c}"].rename("y")], axis=1, join="inner").dropna()
        r, ci = rho_boot(j.d.values, j.y.values)
        rows2.append(dict(dachet_hours=h, yaqubi_contrast=c, yaqubi_effect=kind, n_genes=len(j), spearman=round(r, 3), boot95_low=round(ci[0], 3), boot95_high=round(ci[1], 3)))
R2 = pd.DataFrame(rows2)
H = [f"# script: 03_dachet_comparison.py, run {datetime.date.today().isoformat()}; NEW analysis; seed {SEED}; {NPERM} matched sets; {NBOOT} bootstrap resamples over genes",
     f"# inputs: GSE153669_RAW.tar sha256 {sha(RAW)}; axis_v2_controls.tsv.gz {sha(PHAX)}; hardy_t_by_region.tsv.gz {sha(GTEX)}; contrast_effects.tsv {sha(f'{OUT}/tables/contrast_effects.tsv')}; recurrence_genes.tsv {sha(f'{B80}/recurrence_genes.tsv')}",
     "# effect = log2(x + 1) at the hour minus at 0 h, single patient, room temperature; pH axis positive = higher at higher pH; agonal axis positive = higher in slow death"]
for df, nm in ((R, "dachet_axis_correlation.tsv"), (R2, "dachet_vs_yaqubi_effects.tsv")):
    with open(f"{OUT}/tables/{nm}", "w") as f:
        f.write("\n".join(H) + "\n")
        df.to_csv(f, sep="\t", index=False)
print(R.to_string())
print(R2.to_string())
