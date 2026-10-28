#!/usr/bin/env python3
"""GTEx v8, every tissue (SMTSD): Hardy 3 vs Hardy 2 and Hardy 4 vs Hardy 3 contrasts, scored against the same axes as gtex_death_speed_all_tissues.py (copied from it; same OLS, covariates, gene filter, axes, null; the original is unchanged).

Per tissue (RNA-seq samples with an ischemic time; one sample per donor per tissue), per gene with median TPM >= 1 over all the tissue's samples: OLS of log2(TPM + 1) on contrast group (dummy), ischemic time (hours),
age (band midpoint) and sex, restricted to the two Hardy groups compared. Contrasts: 3v2 = Hardy 3 (intermediate death) vs Hardy 2 (fast natural death); 4v3 = Hardy 4 (slow death) vs Hardy 3. t > 0 = higher in the first-named group.
A contrast is run in a tissue only if both groups have >= MIN_PER_GROUP (10) samples; constant covariates are dropped.
Axes as in the original: control pH Z, Dachet 24 h change, brain 4-vs-2 axis (mean t of Cortex, BA9, BA24 in hardy_t_by_region.tsv.gz; circular for those three tissues only in the 4v2 sense), SZ direction of the 176 recurrent genes, IEG and rPRG sets vs expression-decile-matched null.
Additionally writes brain axes for the new contrasts (mean t of Cortex, BA9, BA24 of this script's 3v2 and 4v3, all three regions required, as decomposition_axes.py builds GTEX) for the pair-agreement analysis (hardy3_axis_removal.py).

    python gtex_hardy3_contrasts.py TPM.gct.gz ATTR.txt SUBJ.txt GSE153669_RAW.tar ROOT OUTDIR
"""
import hashlib, io, re, sys, tarfile, zipfile, csv, datetime
from pathlib import Path
import numpy as np, pandas as pd
import xml.etree.ElementTree as ET
from scipy.stats import spearmanr

TPM, ATTR, SUBJ, RAW, ROOT, OUT = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], Path(sys.argv[5]), Path(sys.argv[6])
OUT.mkdir(parents=True, exist_ok=True)
MIN_PER_GROUP, NPERM, SEED = 10, 10000, 20261010
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
SCRIPT = "analyses/gtex_death_speed_all_tissues/gtex_hardy3_contrasts.py"
IEG = ["FOS", "FOSB", "JUN", "JUNB", "EGR1", "EGR2", "EGR3", "NR4A1", "NR4A2", "NPAS4", "ARC", "IER2", "DUSP1"]  # dachet_live_tissue_2026-10-08/build_report.py line 155


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""): h.update(b)
    return h.hexdigest()


P = {"pH": ROOT / "analyses/ph_and_confound_analyses/results/combined_controls_pH.tsv", "rec": ROOT / "analyses/recurrence_genes.tsv",
     "act": ROOT / "analyses/ph_sensitivity_gene_features_2026-10-09/results_laptop/activity_pH_mean.tsv", "hardy": ROOT / "analyses/gtex_ischemic_time/hardy_t_by_region.tsv.gz"}
SRC = {"gtex_tpm": TPM, "gtex_sample_attributes": ATTR, "gtex_subject_phenotypes": SUBJ, "GSE153669_RAW.tar": RAW, **{k: str(v) for k, v in P.items()}}
SHAS = {k: sha(v) for k, v in SRC.items()}
HEADER = [f"# script: {SCRIPT}", f"# date: {datetime.date.today().isoformat()}"] + [f"# source {k}: {Path(v).name if k.startswith('gtex') or k.endswith('.tar') else v.replace(str(ROOT) + '/', '')} sha256 {SHAS[k]}" for k, v in SRC.items()]
HEADER.append(f"# parameters: min per group {MIN_PER_GROUP}; permutations {NPERM}; seed {SEED}; covariates ischemic time (h), age band midpoint, sex; gene filter median TPM >= 1 within tissue")


def write_tsv(name, df, note, index=False, gz=False):
    p = OUT / name
    import gzip
    op = gzip.open if gz else open
    with op(p, "wt") as f:
        f.write("\n".join(HEADER) + f"\n# {note}\n")
        df.to_csv(f, sep="\t", index=index)
    print("wrote", p)


# ---- samples
attr = pd.read_csv(ATTR, sep="\t", usecols=["SAMPID", "SMTSD", "SMTS", "SMAFRZE", "SMTSISCH"])
subj = pd.read_csv(SUBJ, sep="\t")
s = attr[(attr.SMAFRZE == "RNASEQ") & attr.SMTSISCH.notna()].copy()
s["SUBJID"] = s.SAMPID.str.split("-").str[:2].str.join("-")
s = s.merge(subj, on="SUBJID")
s["age"] = s.AGE.map(lambda a: np.mean([int(x) for x in a.split("-")]))
s["isch_h"] = s.SMTSISCH / 60.0
s["male"] = (s.SEX == 1).astype(float)
hdr = pd.read_csv(TPM, sep="\t", skiprows=2, nrows=0).columns
s = s[s.SAMPID.isin(set(hdr))].reset_index(drop=True)
print("samples", len(s), "tissues", s.SMTSD.nunique())

# ---- donor counts (all categories) and which contrasts run
rows = []
for t, sub in s.groupby("SMTSD"):
    c = sub.DTHHRDY.value_counts()
    r = {"tissue": t, "SMTS": sub.SMTS.iloc[0], "n_samples": len(sub)}
    for k in (0, 1, 2, 3, 4): r[f"n_hardy_{k}"] = int(c.get(k, 0))
    r["n_hardy_missing"] = int(sub.DTHHRDY.isna().sum())
    for nm, a, b in (("3v2", 3, 2), ("4v3", 4, 3)):
        ok = r[f"n_hardy_{a}"] >= MIN_PER_GROUP and r[f"n_hardy_{b}"] >= MIN_PER_GROUP
        r[f"run_{nm}"] = ok
        r[f"skip_reason_{nm}"] = "" if ok else f"fewer than {MIN_PER_GROUP} samples in Hardy {a if r[f'n_hardy_{a}'] < MIN_PER_GROUP else b}"
    rows.append(r)
counts = pd.DataFrame(rows)
write_tsv("donor_counts_hardy3_contrasts.tsv", counts, "samples (one per donor per tissue; RNA-seq samples with an ischemic time) per Hardy category; DTHHRDY from SubjectPhenotypesDS", False)

# ---- expression
t0 = pd.read_csv(TPM, sep="\t", skiprows=2, usecols=["Name", "Description"] + list(s.SAMPID), dtype={c: np.float32 for c in s.SAMPID})
sym = t0.Description.values
Xraw = t0[list(s.SAMPID)].values
del t0
X = np.log2(Xraw + 1).astype(np.float32)
print("matrix", X.shape)

# ---- axes
ct = pd.read_csv(P["pH"], sep="\t", comment="#")
pH = ct[ct.factor == "pH"].set_index("gene")
pHZ = pH.iloc[:, 2].astype(float).to_dict() if False else None
# column 4 of the raw file (index 3), as the brain script reads it
pHd = {}
for l in open(P["pH"]):
    if l.startswith("#") or l.startswith("gene\t"): continue
    f = l.rstrip("\n").split("\t")
    if f[1] == "pH":
        try: pHd[f[0]] = float(f[3])
        except ValueError: pass
pct = pH.expr_percentile
dec = np.minimum((pct.dropna() * 10).astype(int), 9)

def rows_of(data):
    z = zipfile.ZipFile(io.BytesIO(data))
    ss = [''.join(t.text or '' for t in si.iter(NS + 't')) for si in ET.parse(z.open('xl/sharedStrings.xml')).getroot()] if 'xl/sharedStrings.xml' in z.namelist() else []
    rels = {r.get('Id'): r.get('Target') for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
    sh = ET.fromstring(z.read('xl/workbook.xml')).find(NS + 'sheets')[0]
    tgt = rels[sh.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')].lstrip('/'); tgt = tgt if tgt.startswith('xl/') else 'xl/' + tgt
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
with tarfile.open(RAW) as tf:
    for m in tf.getmembers():
        t = re.search(r'_(\d+)H\.xlsx$', m.name)
        if t: D[int(t.group(1))] = rows_of(tf.extractfile(m).read())
dachet = {g: np.log2(D[24][g] + 1) - np.log2(D[0][g] + 1) for g in set(D[0]) & set(D[24])}

hg = pd.read_csv(P["hardy"], sep="\t", index_col=0, compression="gzip")
bcols = [c for c in hg.columns if c in ("Brain - Cortex", "Brain - Frontal Cortex (BA9)", "Brain - Anterior cingulate cortex (BA24)")]
assert len(bcols) == 3, bcols
brain_axis = hg[bcols].mean(axis=1, skipna=True).where(hg[bcols].notna().sum(axis=1) >= 2).dropna().to_dict()

R = [r for r in csv.DictReader((l for l in open(P["rec"]) if not l.startswith("#")), delimiter="\t") if r["scenario"] == "published, brain"]
rec176 = [r for r in R if float(r["q_bh"]) < 0.1]
szdir = {r["gene"]: (1.0 if float(r["up"]) > float(r["down"]) else -1.0) for r in rec176 if float(r["up"]) != float(r["down"])}
print("recurrent q<0.1:", len(rec176), "with direction:", len(szdir))
act = pd.read_csv(P["act"], sep="\t", comment="#")
rprg = act[act.activity_class == "rPRG"].genes.iloc[0].split(",")
SETS = {"IEG": IEG, "rPRG_Tyssowski": rprg}
CONT = {"pH_Z": pHd, "Dachet_24h_change": dachet, "brain_slow_death_axis": brain_axis, "SZ_direction_176": szdir}
by = {k: dec.index[dec == k].values for k in range(10)}
rng = np.random.default_rng(SEED)


def ols_t(Y, M, col):
    Q, Rm = np.linalg.qr(M)
    beta = np.linalg.solve(Rm, Q.T @ Y.T)
    res = Y.T - M @ beta
    s2 = (res ** 2).sum(0) / (M.shape[0] - M.shape[1])
    cd = np.diag(np.linalg.inv(Rm.T @ Rm))[col]
    return beta[col] / np.sqrt(s2 * cd), beta[col]

def set_stat(ser, genes):
    G = [g for g in genes if g in ser.index and g in dec.index]
    if len(G) < 5: return None
    need = dec.loc[G].value_counts().to_dict()
    obs = ser[G].mean()
    pool = {k: ser.reindex(by[k]).dropna().values for k in need}
    null = np.array([np.mean(np.concatenate([rng.choice(pool[k], c, replace=False) for k, c in need.items()])) for _ in range(NPERM)])
    return len(G), obs, obs - null.mean(), (1 + np.sum(np.abs(null - null.mean()) >= abs(obs - null.mean()))) / (NPERM + 1)

res, tmat = [], {"3v2": {}, "4v3": {}}
for t, sub in s.groupby("SMTSD"):
    ix = sub.index.values
    cr = counts[counts.tissue == t].iloc[0]
    keep = np.median(Xraw[:, ix], axis=1) >= 1
    h = sub.DTHHRDY.values
    for nm, a, b in (("3v2", 3, 2), ("4v3", 4, 3)):
        if not cr[f"run_{nm}"]: continue
        sel = np.isin(h, (a, b))
        covs = {"ischemic_time": sub.isch_h.values[sel], "age": sub.age.values[sel], "sex": sub.male.values[sel]}
        dropped = [k for k, v in covs.items() if np.ptp(v) == 0]  # constant covariate (e.g. sex in single-sex tissues) makes the design singular
        M = np.column_stack([np.ones(sel.sum()), (h[sel] == a).astype(float)] + [v for k, v in covs.items() if k not in dropped])
        Y = X[:, ix[sel]][keep].astype(np.float64)
        tv, _ = ols_t(Y, M, 1)
        ser = pd.Series(tv, index=sym[keep]).groupby(level=0).first()
        tmat[nm][t] = ser
        r = {"tissue": t, "SMTS": cr.SMTS, "contrast": "hardy3_vs_hardy2" if nm == "3v2" else "hardy4_vs_hardy3", "n_group_a": int((h[sel] == a).sum()), "n_group_b": int((h[sel] == b).sum()), "n_genes": int(keep.sum()), "covariates_dropped_constant": ",".join(dropped)}
        for an, ad in CONT.items():
            k = [x for x in ser.index if x in ad]
            rho = spearmanr(ser[k].values, [ad[x] for x in k])[0] if len(k) > 10 else np.nan
            r[f"rho_{an}"] = round(float(rho), 4); r[f"n_{an}"] = len(k)
        for sn, gl in SETS.items():
            o = set_stat(ser, gl)
            r.update({f"{sn}_n": o[0], f"{sn}_mean_t": round(o[1], 3), f"{sn}_effect_vs_matched_null": round(o[2], 3), f"{sn}_p": o[3]} if o else {f"{sn}_n": 0})
        res.append(r)
        print(t, nm, r["n_group_a"], r["n_group_b"], r["rho_pH_Z"], r["rho_brain_slow_death_axis"], flush=True)
rho = pd.DataFrame(res)
write_tsv("rho_hardy3_contrasts_by_tissue.tsv", rho, "per tissue and contrast: Spearman rho of gene t (first-named group higher = positive) with each axis; set mean t and effect over expression-decile-matched null. brain_slow_death_axis is the 4-vs-2 axis from Cortex, BA9, BA24 (circular for those three tissues: their own 4v2 samples include Hardy 2 donors shared with 3v2 and Hardy 4 donors shared with 4v3)", False)
for nm, fn in (("3v2", "hardy3_vs_hardy2_t_by_tissue.tsv.gz"), ("4v3", "hardy4_vs_hardy3_t_by_tissue.tsv.gz")):
    if tmat[nm]:
        write_tsv(fn, pd.DataFrame(tmat[nm]).round(4), f"gene x tissue t statistic, contrast {nm}; rows = gene symbol (first of duplicates); blank = gene not tested in that tissue", True, True)

# ---- brain axes for the new contrasts (all three regions required, as decomposition_axes.py builds GTEX)
B3 = ["Brain - Cortex", "Brain - Frontal Cortex (BA9)", "Brain - Anterior cingulate cortex (BA24)"]
for nm, fn in (("3v2", "axis_brain_hardy3_vs_hardy2.tsv.gz"), ("4v3", "axis_brain_hardy4_vs_hardy3.tsv.gz")):
    assert all(c in tmat[nm] for c in B3), (nm, [c for c in B3 if c not in tmat[nm]])
    ax = pd.DataFrame({c: tmat[nm][c] for c in B3})
    out = pd.DataFrame({"mean_t_3_regions": ax.mean(axis=1, skipna=False)}).dropna()
    write_tsv(fn, out.round(5), f"gene -> mean over Cortex, BA9, BA24 of the {nm} t (all three required); {len(out)} genes; same construction as the 4v2 GTEX axis in decomposition_axes.py", True, True)

# ---- agreement between the two contrasts' gene t vectors, per tissue, and summaries
ag = []
for t_ in sorted(set(tmat["3v2"]) & set(tmat["4v3"])):
    a_, b_ = tmat["3v2"][t_], tmat["4v3"][t_]
    k = a_.index.intersection(b_.index)
    ag.append({"tissue": t_, "n_genes": len(k), "rho_3v2_vs_4v3": round(float(spearmanr(a_[k], b_[k])[0]), 4)})
write_tsv("rho_3v2_vs_4v3_by_tissue.tsv", pd.DataFrame(ag), "Spearman rho between the gene t vectors of the two contrasts within a tissue (tissues where both ran)", False)
rho["group"] = np.where(rho.SMTS == "Brain", "brain", "non_brain")
cols = [c for c in rho.columns if c.startswith("rho_") or c.endswith("_mean_t") or c.endswith("_effect_vs_matched_null")]
sm = []
for (c_, g_), d in rho.groupby(["contrast", "group"]):
    r_ = {"contrast": c_, "group": g_, "n_tissues": len(d), "median_n_group_a": float(d.n_group_a.median()), "median_n_group_b": float(d.n_group_b.median())}
    for c in cols:
        r_[f"median_{c}"] = round(float(d[c].median()), 4)
        r_[f"n_pos_{c}"] = int((d[c] > 0).sum())
    sm.append(r_)
write_tsv("summary_medians_hardy3_contrasts.tsv", pd.DataFrame(sm), "median over tissues (brain = SMTS Brain, 13 regions; non_brain = all others that ran) of each rho / set mean t / set effect; n_pos_* = tissues with value > 0", False)
