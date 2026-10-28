#!/usr/bin/env python3
"""Immediate early genes (IEGs) across the death, processing and living-versus-postmortem contrasts.

Gene sets
  S7     the 7 Tyssowski rapid primary response genes (rPRG) inside the 383 recurrent genes of the core-80 base
         (EGR1, NR4A1, DUSP1, FOS, MAFF, PPP1R15A, GADD45G); membership is DERIVED here (rPRG class x 383) and compared with the
         N15 class file and qc_recompute C item 10e (they agree; see log).
  rPRG / dPRG / SRG   Tyssowski et al. 2018 Table S6 classes (mouse KCl-stimulated neurons), mapped to human by symbol exactly as
         analyses/ph_sensitivity_gene_features_2026-10-09/activity_overlap.py (upper-case; HGNC approved, or an alias / previous symbol
         that points to one approved gene). Per contrast the set is restricted to genes with a value in that contrast.
SZ direction of a gene = sign(up - down) of base80_recurrence_2026-10-16/out/recurrence_genes.tsv ('published, brain'); 0 = no direction.
Effect per contrast: t statistic (positive = higher in the later / dead / slow / postmortem / warmer / case condition); Dachet = log2FC
(single patient, no t, no P). Null: 10,000 random sets matched on strata, observed set mean minus null mean, 95% range of the null means.
  Variant A  expression decile only (decile from analyses/confounds/results/combined_controls_pH.tsv expr_percentile, pool = all genes with a
             value in the contrast and a decile; the set's own genes are in the pool) = convention of ieg_slow_death_2026-10-09 / Note 13.
  Variant B  expression decile x DE-prior quintile, strata and pool exactly as ad_383_2026-10-17/ad_383.py (pool = recurrence_genes 'published,
             brain' genes with an SZ direction, excluding the 383; prior quintile cuts from the pool). Set genes absent from the recurrence
             table have no stratum and are dropped from the matched mean (counted in n_matched).
Seeds: numpy default_rng([20261018, contrast index, set index, variant]) per combination. No model calls. Reads only; writes only into the parent folder.
Run:  nice -n 19 python3 ieg_across_contrasts.py
"""
import collections, csv, datetime, functools, gzip, hashlib, io, json, multiprocessing as mp, os, re, sys, tarfile, zipfile
sys.dont_write_bytecode = True
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from scipy import stats

try:
    os.nice(19)
except Exception:
    pass
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.dirname(HERE)
REPO = os.path.expanduser("<analysis repo>"); AN = f"{REPO}/analyses"
A = os.path.expanduser("~/analysis-delivery")
sys.path.insert(0, AN)
SEED, NPERM, NPROC = 20261018, 10000, 6
TODAY = datetime.date.today().isoformat()
IEG7 = ["EGR1", "NR4A1", "DUSP1", "FOS", "MAFF", "PPP1R15A", "GADD45G"]

B80 = f"{A}/base80_recurrence_2026-10-16/out"
P = dict(
    rec=f"{B80}/recurrence_genes.tsv", sets=f"{B80}/recurrent_sets.tsv",
    prior=f"{AN}/de_prior_lists/our_prior/1_global_our_prior.tsv", ph=f"{AN}/confounds/results/combined_controls_pH.tsv",
    tys=f"{REPO}/data/raw/tyssowski2018/NIHMS958209-supplement-6.xlsx", hgnc=f"{REPO}/data/reference/hgnc_complete_set.txt",
    hl=f"{AN}/ieg_halflife_2026-10-07/raw/all_HLs_human_PC1.txt",
    n15=f"{A}/rerun_core80_2026-10-16/L14/N15_ieg_halflife/classes/ieg_classes_q01_recurrence_genes.tsv",
    yaq=f"{A}/yaqubi_2026_axes_2026-10-18/tables/contrast_effects.tsv", yaqdef=f"{A}/yaqubi_2026_axes_2026-10-18/tables/contrast_definitions.tsv",
    dachet=f"{A}/GSE153669_download_2026-10-15/GSE153669_RAW.tar",
    st1=f"{A}/liharska_living_vs_pm_2026-10-18/data/liv_pm_de_ST1.tsv", ct4=f"{A}/liharska_living_vs_pm_2026-10-18/data/collado_torres/TableS4.csv",
    hardy=f"{AN}/gtex_ischemia_2026-10-08/hardy_t_by_region.tsv.gz", gtexall=f"{AN}/gtex_death_speed_all_tissues_2026-10-10/slow_vs_fast_t_by_tissue.tsv.gz",
    gtexn=f"{AN}/gtex_death_speed_all_tissues_2026-10-10/donor_counts_by_tissue.tsv",
    prit=f"{AN}/axes_pritzker_agonal_2026-10-09/axis_pritzker_agonal_GPL17027.tsv.gz",
    bde=f"{AN}/fairness_sample_level/results/b_de", bind=f"{AN}/fairness_sample_level/results/b_independent_set.tsv", blev=f"{AN}/fairness_sample_level/results/b_levels.tsv",
    sc=f"{AN}/single_cell/results/sc_cross_long.tsv.gz",
)


def sha(p, n=12):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()[:n]


def tsv(p):
    return list(csv.DictReader((l for l in open(p, encoding="utf8") if not l.startswith("#")), delimiter="\t"))


LOG = open(f"{OUT}/logs/run.log", "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); LOG.write(s + "\n"); LOG.flush()


import extract
extract.HGNC = extract.load_hgnc()


@functools.lru_cache(maxsize=None)
def hg(s):
    r = extract.symbol_to_hgnc(str(s))
    return r[0] or None


ENS2SYM = {e: s for e, (s, _) in extract.HGNC[3].items()}

# ------------------------------------------------------------------ recurrence population (as ad_383.py)
rec = [r for r in tsv(P["rec"]) if r["scenario"] == "published, brain"]
sz_sign, dec_rec, prior, nup, ndown = {}, {}, {}, {}, {}
for r in rec:
    d = int(r["up"]) - int(r["down"])
    sz_sign[r["gene"]] = (1 if d > 0 else -1) if d else 0
    dec_rec[r["gene"]] = r["expr_decile"]; nup[r["gene"]] = int(r["up"]); ndown[r["gene"]] = int(r["down"])
    prior[r["gene"]] = float(r["de_prior_hitlists"]) if r["de_prior_hitlists"] != "" else np.nan
med = float(np.nanmedian(list(prior.values())))
prior = {g: (med if np.isnan(v) else v) for g, v in prior.items()}
S383 = [r["gene"] for r in rec if float(r["q_bh"]) < 0.1]
assert len(S383) == 383
assert {r["gene"] for r in tsv(P["sets"]) if r["set"] == "study_brain"} == set(S383)
nlists = {r["gene"]: int(r["n_lists_observed"]) for r in rec}
poolB = [g for g in sz_sign if sz_sign[g] != 0 and g not in set(S383)]
qcuts = np.quantile([prior[g] for g in poolB], [0.2, 0.4, 0.6, 0.8])
stratumB = {g: (dec_rec[g], int(np.searchsorted(qcuts, prior[g], side="right"))) for g in sz_sign}
poolB_set = set(poolB)
say("recurrence universe", len(sz_sign), "pool B", len(poolB), "383:", len(S383))

# decile for variant A
ph = pd.read_csv(P["ph"], sep="\t", comment="#"); ph = ph[ph.factor == "pH"].set_index("gene")
pct = ph.expr_percentile.dropna()
decA = np.minimum((pct * 10).astype(int), 9)

# ------------------------------------------------------------------ Tyssowski classes (activity_overlap.py mapping)
import openpyxl
H = pd.read_csv(P["hgnc"], sep="\t", dtype=str, usecols=["symbol", "alias_symbol", "prev_symbol", "status"]); H = H[H.status == "Approved"]
appr = set(H.symbol); alias = {}
for _, r in H.iterrows():
    for col in ("alias_symbol", "prev_symbol"):
        for a in (str(r[col]).split("|") if isinstance(r[col], str) else []):
            alias.setdefault(a.upper(), set()).add(r.symbol)


def to_human(s):
    u = s.upper()
    if u in appr:
        return u
    c = alias.get(u, set())
    return next(iter(c)) if len(c) == 1 else None


wb = openpyxl.load_workbook(P["tys"], read_only=True)
CLASSES, MAPINFO = {}, {}
for sheet, name in (("rapid PRGs", "rPRG"), ("delayed PRGs", "dPRG"), ("SRGs", "SRG")):
    mouse = [str(r[0]) for r in wb[sheet].iter_rows(min_row=2, values_only=True) if r[0]]
    hum = sorted({h for h in (to_human(m) for m in mouse) if h})
    CLASSES[name] = hum; MAPINFO[name] = (len(mouse), sum(1 for m in mouse if to_human(m)), len(hum))
say("Tyssowski mapping (mouse symbols, mapped, distinct human):", MAPINFO)
G7 = sorted(set(CLASSES["rPRG"]) & set(S383))
n15 = pd.read_csv(P["n15"], sep="\t", comment="#")
n15_r = sorted(n15[n15["class"] == "rPRG"].gene)
say("S7 derived:", G7, "| equals N15 rPRG list:", G7 == n15_r, "| equals the 7 named in the brief:", sorted(IEG7) == G7)
assert sorted(IEG7) == G7
SETS = {"S7 (rPRG in the 383)": G7, "rPRG (Tyssowski, all mapped)": CLASSES["rPRG"], "dPRG (Tyssowski)": CLASSES["dPRG"], "SRG (Tyssowski)": CLASSES["SRG"]}
SETNAMES = list(SETS)

# ------------------------------------------------------------------ half-life (N15 convention)
hl = pd.read_csv(P["hl"], sep="\t").dropna(subset=["external_gene_name"]).drop_duplicates("external_gene_name").set_index("external_gene_name")["human_PC1"].dropna()
allv = hl.values
hl_pct = {g: 100 * float((allv < v).mean()) for g, v in hl.items()}
say("half-life genes (distinct symbols with PC1):", len(hl))
HLSUM = {}
for sn_, G in SETS.items():
    v = [hl_pct[g] for g in G if g in hl_pct]
    HLSUM[sn_] = dict(n_set=len(G), n_halflife=len(v), pct_median=float(np.median(v)) if v else np.nan, pct_min=min(v) if v else np.nan, pct_max=max(v) if v else np.nan,
                      n_in_shortest_16pct=sum(x < 16.0 for x in v), pct_mean=float(np.mean(v)) if v else np.nan)
say("half-life summary:", json.dumps(HLSUM, default=float))

# ------------------------------------------------------------------ contrasts
CON = collections.OrderedDict()   # name -> dict(group, unit, stat Series, p Series or None, note)
SRC = {}


def addc(name, group, unit, stat, p=None, note=""):
    stat = stat[~stat.index.duplicated(keep="first")].dropna()
    if p is not None:
        p = p[~p.index.duplicated(keep="first")].reindex(stat.index)
    CON[name] = dict(group=group, unit=unit, stat=stat, p=p, note=note)


def tp(t, df):
    return pd.Series(2 * stats.t.sf(np.abs(t.values), df), index=t.index)


# (1) and (2) Yaqubi
yq = pd.read_csv(P["yaq"], sep="\t", comment="#", index_col=0)
ydef = pd.read_csv(P["yaqdef"], sep="\t", comment="#").set_index("contrast")
SRC["yaqubi_contrast_effects"] = P["yaq"]; SRC["yaqubi_contrast_definitions"] = P["yaqdef"]
YL = {"S_6h4C_vs_0h": "6 h at 4 C vs frozen", "S_24h4C_vs_0h": "24 h at 4 C vs frozen", "S_36h4C_vs_0h": "36 h at 4 C vs frozen",
      "S_6h20C_vs_0h": "6 h at 20 C vs frozen", "S_24h20C_vs_0h": "24 h at 20 C vs frozen", "S_36h20C_vs_0h": "36 h at 20 C vs frozen",
      "S_20C_minus_4C_6h": "20 C minus 4 C, 6 h", "S_20C_minus_4C_24h": "20 C minus 4 C, 24 h", "S_20C_minus_4C_36h": "20 C minus 4 C, 36 h",
      "A_SPMI_vs_CTRL": "autopsy short PMI vs frozen", "A_LPMI_vs_CTRL": "autopsy long PMI vs frozen"}
for k, lab in YL.items():
    t = yq[f"t_{k}"]; df = float(ydef.loc[k, "resid_df"])
    addc("Yaqubi: " + lab, "Yaqubi surgical holds" if k.startswith("S_") else "Yaqubi autopsy", "t (OLS, own refit, not DESeq2)", t, tp(t.dropna(), df),
         f"P from t with {int(df)} residual df; positive = higher in held / autopsy than immediately frozen (20 C minus 4 C: warmer higher)")

# (4) Dachet time course (loader copied from yaqubi_2026_axes_2026-10-18/scripts/03_dachet_comparison.py)
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
HOURS = [0, 1, 2, 4, 8, 12, 24]


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
with tarfile.open(P["dachet"]) as tf:
    for m in tf.getmembers():
        t = re.search(r'_(\d+)H\.xlsx$', m.name)
        if t:
            rows = xlsx_sheet(tf.extractfile(m).read())
            Ex[int(t.group(1))] = {r['A']: float(r['B']) for r in rows[1:] if 'A' in r and 'B' in r and re.match(r'^-?[\d.eE+-]+$', r['B'])}
dgenes = sorted(set.intersection(*[set(Ex[h]) for h in HOURS]))
Ld = np.array([[np.log2(Ex[h][g] + 1) for h in HOURS] for g in dgenes])
Dd = pd.DataFrame(Ld - Ld[:, [0]], index=[hg(g) for g in dgenes], columns=HOURS)
Dd = Dd[Dd.index.notna()]; Dd = Dd[~Dd.index.duplicated(keep="first")]
SRC["dachet_GSE153669_RAW.tar"] = P["dachet"]
say("Dachet genes in all 7 time points:", len(dgenes), "mapped:", len(Dd))
for h in HOURS[1:]:
    addc(f"Dachet: {h} h vs 0 h", "Dachet time course (1 patient)", "log2FC (log2(x+1); no t, no P)", Dd[h], None, "one resected-tissue sample per time point; positive = higher at the later time")

# (3) Liharska / Collado-Torres
SRC["liharska_ST1"] = P["st1"]; SRC["collado_torres_TableS4"] = P["ct4"]


def best(df, sym, pc):
    df = df.assign(_s=sym).dropna(subset=["_s", pc]).sort_values(pc).drop_duplicates("_s")
    return df.set_index("_s")


l = pd.read_csv(P["st1"], sep="\t", comment="#"); l["sym"] = l.gene.map(ENS2SYM); l1 = best(l, l.sym, "P.Value")
addc("Liharska: living vs postmortem (original ST1)", "Living vs postmortem", "t (moderated)", l1.t, l1["P.Value"], "positive = higher in postmortem")
c = pd.read_csv(P["ct4"], low_memory=False); c = c[c.isExprs]
c["sym"] = c.ensemblID.map(ENS2SYM); c["sym"] = c.sym.where(c.sym.notna(), c.Symbol.map(lambda s: hg(s) if isinstance(s, str) else None))
for tag, lab in (("original", "Collado-Torres reproduction (original)"), ("without_qSVA", "Collado-Torres without qSVA"), ("with_qSVA", "Collado-Torres with qSVA")):
    ci = best(c, c.sym, f"{tag}_P.Value").dropna(subset=[f"{tag}_t"])
    addc(f"Liharska: living vs postmortem ({lab})", "Living vs postmortem", "t (limma-voom)", ci[f"{tag}_t"], ci[f"{tag}_P.Value"], "positive = higher in postmortem")

# (5) GTEx slow (Hardy 4) vs fast natural (2)
SRC["gtex_hardy_t_by_region"] = P["hardy"]; SRC["gtex_slow_vs_fast_t_by_tissue"] = P["gtexall"]; SRC["gtex_donor_counts"] = P["gtexn"]
hd = pd.read_csv(P["hardy"], sep="\t", index_col=0)
ga = pd.read_csv(P["gtexall"], sep="\t", comment="#", index_col=0)
gn = pd.read_csv(P["gtexn"], sep="\t", comment="#").set_index("tissue")
for tcol in hd.columns:
    n = int(gn.loc[tcol, "n_hardy_4"] + gn.loc[tcol, "n_hardy_2"])
    addc(f"GTEx slow vs fast death: {tcol.replace('Brain - ', '')}", "GTEx brain regions", "t (OLS; ischemic time, age, sex)", hd[tcol], tp(hd[tcol].dropna(), n - 5),
         f"P approximate: t distribution, df = {n}-5 (n Hardy 4 + n Hardy 2, 5 model terms)")
cort = ["Brain - Cortex", "Brain - Frontal Cortex (BA9)", "Brain - Anterior cingulate cortex (BA24)"]
cm = hd[cort].mean(axis=1, skipna=True).where(hd[cort].notna().sum(axis=1) >= 2)
addc("GTEx slow vs fast death: cortical mean (Cortex, BA9, BA24)", "GTEx brain regions", "t (mean of 3 regions)", cm, None, "Note 6 convention; no P (mean of three t)")
chk = max(float((ga[t].sub(hd[t]).abs().max())) for t in hd.columns if t in ga.columns)
say("max |difference| between hardy_t_by_region and the all-tissue table for brain regions:", chk)
for tcol in ga.columns:
    if tcol.startswith("Brain - "):
        continue
    n = int(gn.loc[tcol, "n_hardy_4"] + gn.loc[tcol, "n_hardy_2"])
    addc(f"GTEx slow vs fast death: {tcol}", "GTEx other tissues", "t (OLS; ischemic time, age, sex)", ga[tcol], tp(ga[tcol].dropna(), n - 5), f"P approximate: t distribution, df = {n}-5")
GTEX_TISSUES_ALL = [t for t in ga.columns]
say("GTEx tissues with a slow-vs-fast table:", len(GTEX_TISSUES_ALL), "(brain", sum(t.startswith('Brain') for t in GTEX_TISSUES_ALL), ")")

# (6) Pritzker agonal factor
SRC["pritzker_agonal_axis"] = P["prit"]
with gzip.open(P["prit"], "rt") as f:
    hdr = f.readline()
dfs = [int(x) for x in re.findall(r"df (\d+)", hdr)]
pa = pd.read_csv(P["prit"], sep="\t", comment="#", index_col="gene")
pa = pa.groupby(level=0).mean(numeric_only=True)
addc("Pritzker agonal factor (GSE92538, n 104)", "Pritzker agonal score", "t (OLS: dx + age + sex + PMI + agonal)", pa.t_agonal, tp(pa.t_agonal.dropna(), dfs[0]), f"P from t, df {dfs[0]}; positive = higher with a higher agonal score")
addc("Pritzker agonal factor + pH (n 92)", "Pritzker agonal score", "t (OLS, with tissue pH)", pa.t_agonal_with_pH, tp(pa.t_agonal_with_pH.dropna(), dfs[1]), f"P from t, df {dfs[1]}")

# (7) SZ case-control: 22 bulk datasets of the sample-level reanalysis (b_de), L0 unadjusted and L4 all measured covariates
SRC["b_independent_set"] = P["bind"]
bind = pd.read_csv(P["bind"], sep="\t")["dataset"].tolist()
assert len(bind) == 22
TD = {"L0": {}, "L4": {}}
for ds in bind:
    p_ = f"{P['bde']}/{ds}.tsv.gz"; SRC[f"b_de/{ds}"] = p_
    d = pd.read_csv(p_, sep="\t", index_col=0)
    d = d[~d.index.duplicated(keep="first")]
    for lv in TD:
        TD[lv][ds] = d[f"t_{lv}"]
SZDS = {}
for lv, lab in (("L0", "unadjusted"), ("L4", "all measured covariates")):
    M = pd.DataFrame(TD[lv]); SZDS[lv] = M
    mt = M.mean(axis=1)[M.notna().sum(axis=1) >= 11]
    addc(f"SZ case-control: mean t over the 22 bulk datasets ({lv}, {lab})", "SZ case-control comparator", "mean t over datasets", mt, None,
         "positive = higher in cases; genes in >= 11 of the 22 datasets; the SZ comparator is not independent of the recurrence directions")

# (7b) single-nucleus SZ case-control, class level: Endresz 2026 meta and Lee 2025 (PsychAD)
SRC["sc_cross_long"] = P["sc"]
scd = pd.read_csv(P["sc"], sep="\t", comment="#", usecols=["source", "gemma_class", "cell_type", "gene", "t"])
for src_, lab in (("Endresz 2026 meta", "Endresz 2026 meta-analysis"), ("Lee 2025 (PsychAD)", "PsychAD (Lee 2025)")):
    x = scd[(scd.source == src_) & scd.gemma_class.notna()]
    for cls, xx in x.groupby("gemma_class"):
        s_ = xx.groupby("gene").t.mean()
        addc(f"SZ single nucleus: {lab}, {cls}", "SZ single-nucleus comparator", "mean t over the class's cell types", s_, None,
             "positive = higher in SZ; genes not measured in the class are absent (expression filter)")
del scd
CNAMES = list(CON)
say("contrasts:", len(CNAMES))
SRC.update(rec=P["rec"], sets=P["sets"], prior=P["prior"], ph=P["ph"], tyssowski=P["tys"], hgnc=P["hgnc"], halflife=P["hl"], n15=P["n15"])


# ------------------------------------------------------------------ null engine
def draw_sums(rng, arr, k, D=NPERM):
    N = len(arr)
    if k == 1:
        return arr[rng.integers(0, N, D)]
    if N < k:
        return arr[rng.integers(0, N, (D, k))].sum(1)
    if k * k / (2 * N) < 2.0:
        idx = rng.integers(0, N, (D, k))
        while True:
            s = np.sort(idx, 1); bad = (s[:, 1:] == s[:, :-1]).any(1)
            if not bad.any():
                break
            idx[bad] = rng.integers(0, N, (int(bad.sum()), k))
        return arr[idx].sum(1)
    out = np.zeros(D)
    for a in range(0, D, 1000):
        b = min(D, a + 1000)
        idx = np.argpartition(rng.random((b - a, N)), k - 1, axis=1)[:, :k]
        out[a:b] = arr[idx].sum(1)
    return out


def run_combo(args):
    ci, si, variant = args
    name = CNAMES[ci]; c = CON[name]; v = c["stat"]; G = SETS[SETNAMES[si]]
    rng = np.random.default_rng([SEED, ci, si, 0 if variant == "A" else 1])
    gv = [g for g in G if g in v.index]
    if variant == "A":
        strat = {g: int(decA[g]) for g in gv if g in decA.index}
        pool = collections.defaultdict(list)
        for g in v.index.intersection(decA.index):
            pool[int(decA[g])].append(g)
    else:
        strat = {g: stratumB[g] for g in gv if g in stratumB}
        pool = collections.defaultdict(list)
        for g in v.index:
            if g in poolB_set:
                pool[stratumB[g]].append(g)
    need = collections.Counter(strat.values())
    matched = [g for g in strat if len(pool[strat[g]]) > 0]
    need = collections.Counter(strat[g] for g in matched)
    if len(matched) == 0:
        return (ci, si, variant, 0, np.nan, np.nan, np.nan, np.nan, 0)
    tot = np.zeros(NPERM); short = 0
    for s, k in need.items():
        arr = v.reindex(pool[s]).values
        short += int(len(arr) < k)
        tot += draw_sums(rng, arr, k)
    nm = tot / len(matched)
    obs = float(v.reindex(matched).mean())
    return (ci, si, variant, len(matched), obs, float(nm.mean()), float(np.percentile(nm, 2.5)), float(np.percentile(nm, 97.5)), short)


if __name__ == "__main__":
    combos = [(ci, si, vr) for ci in range(len(CNAMES)) for si in range(len(SETNAMES)) for vr in ("A", "B")]
    say("combinations:", len(combos), "draws each:", NPERM)
    with mp.get_context("fork").Pool(NPROC) as pool_:
        res = pool_.map(run_combo, combos, chunksize=8)
    R = {(r[0], r[1], r[2]): r for r in res}

    rows = []
    for ci, name in enumerate(CNAMES):
        c = CON[name]; v = c["stat"]; pv = c["p"]
        for si, sn_ in enumerate(SETNAMES):
            G = SETS[sn_]; gv = [g for g in G if g in v.index]
            vals = v.reindex(gv)
            up = int((vals > 0).sum()); dn = int((vals < 0).sum())
            if pv is not None:
                pp = pv.reindex(gv)
                up5 = int(((vals > 0) & (pp < 0.05)).sum()); dn5 = int(((vals < 0) & (pp < 0.05)).sum())
            else:
                up5 = dn5 = None
            szs = [g for g in gv if sz_sign.get(g, 0) != 0 and vals[g] != 0]
            conc = sum(1 for g in szs if np.sign(vals[g]) == sz_sign[g])
            h = HLSUM[sn_]
            for vr, vlab in (("A", "A: expression decile (Note 13 convention)"), ("B", "B: expression decile x DE-prior quintile (ad_383.py)")):
                r = R[(ci, si, vr)]
                rows.append(dict(contrast_group=c["group"], contrast=name, set=sn_, stat_unit=c["unit"], null_variant=vlab, n_set_genes=len(G), n_with_value=len(gv), n_matched=r[3],
                                 mean_stat=r[4], matched_null_mean=r[5], null_lo95=r[6], null_hi95=r[7], difference=(r[4] - r[5]) if r[3] else np.nan,
                                 n_up=up, n_down=dn, n_up_p05=up5, n_down_p05=dn5, n_with_sz_direction=len(szs), n_concordant_with_sz=conc,
                                 frac_concordant=(conc / len(szs)) if szs else np.nan, strata_short_of_pool=r[8],
                                 halflife_n=h["n_halflife"], halflife_pct_median=h["pct_median"], halflife_pct_min=h["pct_min"], halflife_pct_max=h["pct_max"],
                                 halflife_n_in_shortest_16pct=h["n_in_shortest_16pct"]))
    T = pd.DataFrame(rows)

    # per-gene table for the 7
    hl_ = {}
    pg = []
    for name in CNAMES:
        c = CON[name]; v = c["stat"]; pv = c["p"]
        for g in G7:
            e = float(v[g]) if g in v.index else np.nan
            p_ = float(pv[g]) if (pv is not None and g in v.index and pd.notna(pv[g])) else np.nan
            pg.append(dict(contrast_group=c["group"], contrast=name, gene=g, effect=e, stat_unit=c["unit"], P=p_,
                           direction=("up" if e > 0 else "down" if e < 0 else ("not measured" if np.isnan(e) else "zero")),
                           sz_direction={1: "up", -1: "down", 0: "none"}[sz_sign[g]], sz_lists_up=nup[g], sz_lists_down=ndown[g],
                           concordant_with_sz=(None if np.isnan(e) or e == 0 or sz_sign[g] == 0 else bool(np.sign(e) == sz_sign[g])),
                           P_below_0p05=(None if np.isnan(p_) else bool(p_ < 0.05)), halflife_pct=hl_pct.get(g, np.nan)))
    for lv in ("L0", "L4"):
        M = SZDS[lv]
        for ds in bind:
            for g in G7:
                e = float(M.loc[g, ds]) if g in M.index and pd.notna(M.loc[g, ds]) else np.nan
                pg.append(dict(contrast_group="SZ per dataset (22 bulk)", contrast=f"SZ case-control {ds} ({lv})", gene=g, effect=e, stat_unit="t", P=np.nan,
                               direction=("up" if e > 0 else "down" if e < 0 else "not measured"), sz_direction={1: "up", -1: "down", 0: "none"}[sz_sign[g]], sz_lists_up=nup[g], sz_lists_down=ndown[g],
                               concordant_with_sz=(None if np.isnan(e) or e == 0 or sz_sign[g] == 0 else bool(np.sign(e) == sz_sign[g])), P_below_0p05=None, halflife_pct=hl_pct.get(g, np.nan)))
    PG = pd.DataFrame(pg)

    SHAS = {k: sha(v) for k, v in SRC.items() if os.path.isfile(v)}
    hdr = [f"# script: scripts/ieg_across_contrasts.py (NEW analysis), run {TODAY}; seeds: numpy default_rng([{SEED}, contrast index, set index, variant]); {NPERM} matched draws per combination; no model calls",
           f"# gene sets: S7 = {','.join(G7)} (rPRG x 383 of core-80 base); Tyssowski 2018 Table S6 mapped by symbol; mapping (mouse symbols, mapped, distinct human): {json.dumps(MAPINFO)}",
           "# sign convention: positive effect = higher in the later / dead / slow-death / postmortem / warmer / case condition; Dachet effect = log2FC; null 95% range = 2.5 to 97.5 percentile of the null set means",
           "# inputs sha256 (first 12): " + "; ".join(f"{k} {v}" for k, v in SHAS.items())]
    json.dump(dict(date=TODAY, seed=SEED, nperm=NPERM, inputs_sha256_12=SHAS, mapinfo=MAPINFO, halflife_summary=HLSUM, n_contrasts=len(CNAMES), python=sys.version.split()[0],
                   numpy=np.__version__, pandas=pd.__version__), open(f"{OUT}/run_meta.json", "w"), indent=1, default=float)
    for fn, df in (("ieg_table.tsv", T), ("per_gene_7.tsv", PG)):
        with open(f"{OUT}/{fn}", "w") as f:
            f.write("\n".join(hdr) + "\n")
            df.to_csv(f, sep="\t", index=False, float_format="%.5g")
    # half-life per gene for the sets (all genes of the four sets)
    hlrows = [dict(set=sn_, gene=g, class_membership=("rPRG" if g in CLASSES["rPRG"] else "dPRG" if g in CLASSES["dPRG"] else "SRG"), in_383=(g in set(S383)), halflife_PC1=hl.get(g, np.nan),
                   halflife_pct_genomewide=hl_pct.get(g, np.nan)) for sn_, G in SETS.items() if sn_ != "S7 (rPRG in the 383)" for g in G]
    hlrows += [dict(set="S7 (rPRG in the 383)", gene=g, class_membership="rPRG", in_383=True, halflife_PC1=hl.get(g, np.nan), halflife_pct_genomewide=hl_pct.get(g, np.nan)) for g in G7]
    with open(f"{OUT}/halflife_by_gene.tsv", "w") as f:
        f.write("\n".join(hdr[:3] + [hdr[3], "# half-life: Agarwal and Kelley 2022 human PC1 (lower = shorter-lived); percentile = percent of the " + str(len(hl)) + " genome-wide genes with a lower PC1 (N15 convention)"]) + "\n")
        pd.DataFrame(hlrows).to_csv(f, sep="\t", index=False, float_format="%.5g")
    say("wrote ieg_table.tsv", len(T), "rows; per_gene_7.tsv", len(PG), "rows")
