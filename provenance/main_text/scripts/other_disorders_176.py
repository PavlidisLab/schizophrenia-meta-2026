#!/usr/bin/env python3
"""Analysis B: do the 176 study-tier recurrent genes (recurrence_genes.tsv, scenario 'published, brain', q_bh < 0.1) recur in
postmortem-brain DE results for bipolar disorder (BD), major depression (MDD) and autism (ASD)?

Sources (all on disk; nothing fetched):
  G18-array  Gandal et al. 2018 Science 359:693 (PMID 29439242), microarray meta-analysis, Table S1 (github CSV): ASD, SCZ, BD, MDD beta, P, FDR.
  G18-RNAseq Gandal et al. 2018 Science 362:eaat8127 (PMID 30545856), PsychENCODE RNA-seq gene-level DGE (Table S1, sheet DGE): ASD, SCZ, BD log2FC, t, p, fdr.
  Gemma      Gemma postmortem-brain BD and MDD result sets already downloaded by de_prior_lists/checks_2026-09-29 (53 lists, GET only),
             and 6 autism result sets fetched by fetch_gemma_asd.py (gemma_asd_cache/, GET only).

(1) Over-representation of the 176 among DE genes (FDR < 0.05, and < 0.10), by four nulls:
    hypergeometric; DE-prior-adjusted expectation (logistic DE ~ log1p(prior) + log1p(prior)^2 + expr_decile, fit on the pool without the 176, Poisson-binomial tail);
    matched permutation (decile x DE-prior quintile strata).
(2) Direction: sign(disorder effect) vs the SZ direction (sign(up - down) of recurrence_genes.tsv): binomial vs 50% and vs the matched null.
(3) Mean SZ-aligned signed z, and Spearman of the disorder's signed z with pH and agonal axes among the 176 (null: matched sets).
(4) Hagihara 2018 Table S1 per-subject pH by diagnosis, recomputed.

    PYTHONPATH=<dir containing pyxlsb> python3 analyses/other_disorders_recurrent_genes/other_disorders_176.py
Fixed seed 20261009. Pool for the null: published-brain genes with a SZ sign, minus the 176.
"""
import collections, csv, glob, gzip, hashlib, io, math, os, random, sys, zipfile
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
AN = os.path.dirname(HERE)
REPO = os.path.dirname(AN)
sys.path.insert(0, AN)
SEED, NPERM_MAIN, NPERM_LIST, NPERM_AXIS = 20261009, 10000, 2000, 5000

REC = os.path.join(AN, "recurrence_genes.tsv")
PRIOR = os.path.join(AN, "de_prior_lists/our_prior/1_global_our_prior.tsv")
G18A_S1 = os.path.join(REPO, "data/raw/29439242/github/TableS1_Microarray_Meta-Analyses.csv")
G18B_S1 = os.path.join(REPO, "data/raw/30545856/browser/aat8127_table_s1.xlsx")
GEM_DIR = os.path.join(AN, "de_prior_lists/checks_2026-09-29")
GEM_LISTS = os.path.join(GEM_DIR, "check1_lists.tsv")
PH_Z = os.path.join(AN, "confounds/results/combined_controls_pH.tsv")
GTEX = os.path.join(AN, "gtex_ischemia_2026-10-08/hardy_t_by_region.tsv.gz")
HAT = os.path.join(REPO, "data/raw/40921414/supplement/NPR2-45-e70053-s001.zip")
DAI = os.path.join(REPO, "data/raw/33841074/supp.zip")
HAG = os.path.join(REPO, "data/raw/28776581/supplement/npp2017167x2.xlsx")
SCORE176 = os.path.join(AN, "agonal_axes_score_176_2026-10-09/score_176.tsv")


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()[:12]


def tsv(p):
    return list(csv.DictReader((l for l in open(p, encoding="utf8") if not l.startswith("#")), delimiter="\t"))


try:
    import extract, functools
    extract.HGNC = extract.load_hgnc()
    @functools.lru_cache(maxsize=None)
    def hg(s):
        r = extract.symbol_to_hgnc(s)
        return r[0] or None
except Exception as e:  # pragma: no cover
    print("extract.symbol_to_hgnc unavailable, symbols used as printed:", e)
    def hg(s):
        return s


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


# ------------------------------------------------------------------ population
rec = [r for r in tsv(REC) if r["scenario"] == "published, brain"]
sign, dec, prior = {}, {}, {}
for r in rec:
    d = int(r["up"]) - int(r["down"])
    sign[r["gene"]] = (1 if d > 0 else -1) if d else 0
    dec[r["gene"]] = r["expr_decile"]
    prior[r["gene"]] = float(r["de_prior_hitlists"]) if r["de_prior_hitlists"] != "" else np.nan
med = float(np.nanmedian(list(prior.values())))
prior = {g: (med if np.isnan(v) else v) for g, v in prior.items()}
S176 = [r["gene"] for r in rec if float(r["q_bh"]) < 0.1]
assert len(S176) > 0
  # size of the recurrent set is data-dependent (176 before the 2026-10-11 corpus change)
S176set = set(S176)
pool = [g for g in sign if sign[g] != 0 and g not in S176set]
# strata: expr_decile x prior quintile (quintiles over the pool)
pv = np.array([prior[g] for g in pool])
qcuts = np.quantile(pv, [0.2, 0.4, 0.6, 0.8])
def pq(g):
    return int(np.searchsorted(qcuts, prior[g], side="right"))
stratum = {g: (dec[g], pq(g)) for g in list(sign)}
by_st = collections.defaultdict(list)
for g in pool:
    by_st[stratum[g]].append(g)
need = collections.Counter(stratum[g] for g in S176)
short = {s: k for s, k in need.items() if len(by_st[s]) < k}
print("strata needed", len(need), "| strata with fewer pool genes than needed:", short)

ALL = sorted(set(sign))                      # master gene index
GI = {g: i for i, g in enumerate(ALL)}
s176_idx = np.array([GI[g] for g in S176])
sgn = np.array([sign[g] for g in ALL], float)
prior_arr = np.array([prior[g] for g in ALL])
dec_arr = np.array([dec[g] for g in ALL])

rng = random.Random(SEED)
def draws(n):
    out = np.zeros((n, len(S176)), dtype=np.int32)
    for i in range(n):
        row = []
        for s, k in need.items():
            cand = by_st[s]
            row += rng.sample(cand, k) if len(cand) >= k else rng.choices(cand, k=k)
        out[i] = [GI[g] for g in row]
    return out

print("drawing matched null sets ...")
DR = draws(NPERM_MAIN)


# ------------------------------------------------------------------ disorder tables
class Table:
    """per gene (master index): signed z, p, fdr; NaN where unmeasured"""
    def __init__(self, name, source, disorder, z, p, q, extra=""):
        self.name, self.source, self.disorder, self.extra = name, source, disorder, extra
        self.z, self.p, self.q = z, p, q
        self.meas = ~np.isnan(z)


def vec(d, key):
    a = np.full(len(ALL), np.nan)
    for g, v in d.items():
        if g in GI:
            a[GI[g]] = v[key]
    return a


def zfrom(p, beta):
    p = max(min(float(p), 1.0), 1e-300)
    return float(np.sign(beta)) * float(stats.norm.isf(p / 2))


tables = []
# --- Gandal 2018 array
df = pd.read_csv(G18A_S1)
for dis in ("SCZ", "BD", "MDD", "ASD"):
    best = {}
    for sym, b, p, f in zip(df["hgnc_symbol"].fillna(df["external_gene_id"]), df[f"{dis}.beta_log2FC"], df[f"{dis}.P.value"], df[f"{dis}.FDR"]):
        if not isinstance(sym, str) or not sym or np.isnan(p):
            continue
        g = hg(sym)
        if g and g in GI and (g not in best or p < best[g]["p"]):
            best[g] = dict(z=zfrom(p, b), p=p, q=f)
    tables.append(Table(f"G18-array {dis}", "Gandal 2018 Science (PMID 29439242) microarray meta-analysis Table S1", dis, vec(best, "z"), vec(best, "p"), vec(best, "q")))

# --- Gandal 2018 PsychENCODE RNA-seq
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
def read_dge(path):
    z = zipfile.ZipFile(path)
    ss = []
    for ev, el in ET.iterparse(z.open("xl/sharedStrings.xml")):
        if el.tag == NS + "si":
            ss.append("".join(t.text or "" for t in el.iter(NS + "t")))
            el.clear()
    rels = {r.get("Id"): r.get("Target") for r in ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))}
    for s in ET.fromstring(z.read("xl/workbook.xml")).find(NS + "sheets"):
        if s.get("name") == "DGE":
            tgt = rels[s.get("{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id")].lstrip("/")
            tgt = tgt if tgt.startswith("xl/") else "xl/" + tgt
            rows = []
            for ev, el in ET.iterparse(z.open(tgt)):
                if el.tag == NS + "row":
                    row = {}
                    for c in el.findall(NS + "c"):
                        v = c.find(NS + "v")
                        if v is not None:
                            col = "".join(ch for ch in c.get("r") if ch.isalpha())
                            row[col] = ss[int(v.text)] if c.get("t") == "s" else v.text
                    rows.append(row)
                    el.clear()
            return rows
rows = read_dge(G18B_S1)
hdr = {v: k for k, v in rows[0].items()}
for dis in ("SCZ", "BD", "ASD"):
    best = {}
    for r in rows[1:]:
        try:
            sym = r[hdr["gene_name"]]
            b, p, f = float(r[hdr[f"{dis}.log2FC"]]), float(r[hdr[f"{dis}.p.value"]]), float(r[hdr[f"{dis}.fdr"]])
        except (KeyError, ValueError):
            continue
        g = hg(sym)
        if g and g in GI and (g not in best or p < best[g]["p"]):
            best[g] = dict(z=zfrom(p, b), p=p, q=f)
    tables.append(Table(f"G18-RNAseq {dis}", "Gandal 2018 Science (PMID 30545856) PsychENCODE RNA-seq DGE (Table S1)", dis, vec(best, "z"), vec(best, "p"), vec(best, "q")))

# --- Gemma BD / MDD per list
lists = tsv(GEM_LISTS)
for L in lists:
    L["cache"] = os.path.join(GEM_DIR, "cache")
ASD_SEL = os.path.join(HERE, "gemma_asd_selection.tsv")
for r in tsv(ASD_SEL):
    lists.append(dict(list_id=f"ASD__{r['short_name']}__" + ("all" if r["subset"] == "all" else r["subset"].replace(" ", "_")), disorder="ASD", short_name=r["short_name"],
                      result_set=r["result_set"], level_id=r["level_id"], n_case="", n_control="", genes_q05="",
                      cache=os.path.join(HERE, "gemma_asd_cache")))
gem_rows, gem_tabs = [], collections.defaultdict(list)
for L in lists:
    f = glob.glob(os.path.join(L["cache"], L["short_name"], L["result_set"] + ".tsv.gz"))
    if not f:
        continue
    d = pd.read_csv(f[0], sep="\t", comment="#", low_memory=False)
    lv = L["level_id"]
    pc, lc = f"contrast_{lv}_pvalue", f"contrast_{lv}_log2fc"
    if pc not in d.columns:
        continue
    d = d[d["gene_official_symbol"].notna() & ~d["gene_official_symbol"].astype(str).str.contains(r"\|")]
    d = d[d[pc].notna()]
    best = {}
    for sym, p, b in zip(d["gene_official_symbol"], d[pc], d[lc]):
        g = hg(str(sym))
        if g and g in GI and (g not in best or p < best[g]["p"]):
            best[g] = dict(p=float(p), b=float(b))
    if not best:
        continue
    genes = list(best)
    q = bh([best[g]["p"] for g in genes])
    dd = {g: dict(z=zfrom(best[g]["p"], best[g]["b"]), p=best[g]["p"], q=qq) for g, qq in zip(genes, q)}
    t = Table(f"Gemma {L['list_id']}", f"Gemma result set {L['result_set']} ({L['short_name']})", L["disorder"], vec(dd, "z"), vec(dd, "p"), vec(dd, "q"),
              extra=f"n_case {L['n_case']}, n_control {L['n_control']}")
    gem_tabs[L["disorder"]].append(t)
    gem_rows.append((L, t))

# pooled Gemma per disorder: Stouffer z over lists, union of DE; with and without lists holding > 18.3% DE of tested genes
def pooled(dis, drop_big):
    ts = [t for t in gem_tabs[dis]]
    keep = []
    for t in ts:
        frac = np.nansum(t.q < 0.05) / max(t.meas.sum(), 1)
        if drop_big and frac > 0.183:
            continue
        keep.append(t)
    Z = np.vstack([t.z for t in keep])
    k = (~np.isnan(Z)).sum(0)
    zs = np.where(k > 0, np.nansum(Z, 0) / np.sqrt(np.maximum(k, 1)), np.nan)
    Q = np.vstack([t.q for t in keep])
    de_any = np.nanmin(np.where(np.isnan(Q), 9, Q), 0)
    de_any = np.where(k > 0, de_any, np.nan)
    nde = np.nansum(Q < 0.05, 0).astype(float)
    nde = np.where(k > 0, nde, np.nan)
    # q-like column for union: min q over lists (conservative "any list"), p proxy unused
    return Table(f"Gemma {dis} pooled" + (" (lists with <= 18.3% DE)" if drop_big else " (all lists)"), "Gemma lists in de_prior_lists/checks_2026-09-29",
                 dis, zs, np.full(len(ALL), np.nan), de_any, extra=f"{len(keep)} lists, signed z = Stouffer over lists; DE = q < threshold in any list"), keep

pooled_tabs = []
for dis in ("BD", "MDD", "ASD"):
    for drop in (False, True):
        t, keep = pooled(dis, drop)
        pooled_tabs.append(t)
tables_main = tables + pooled_tabs


# ------------------------------------------------------------------ statistics
def pb_tail(ps, k):
    dist = np.zeros(len(ps) + 1)
    dist[0] = 1
    for p in ps:
        dist[1:] = dist[1:] * (1 - p) + dist[:-1] * p
        dist[0] *= (1 - p)
    return float(dist[k:].sum())


def fit_expected(t, thr):
    """DE ~ log1p(prior) + log1p(prior)^2 + C(decile), fit on the pool; predicted p for the 176 measured in t"""
    pool_idx = np.array([GI[g] for g in pool])
    m = t.meas[pool_idx]
    idx = pool_idx[m]
    y = (t.q[idx] < thr).astype(int)
    if y.sum() < 10 or y.sum() > len(y) - 10:
        return None
    lp = np.log1p(prior_arr[idx])
    D = pd.DataFrame(dict(y=y, lp=lp, lp2=lp ** 2, dec=dec_arr[idx]))
    try:
        mod = smf.glm("y ~ lp + lp2 + C(dec)", D, family=sm.families.Binomial()).fit()
    except Exception:
        return None
    s = s176_idx[t.meas[s176_idx]]
    T = pd.DataFrame(dict(lp=np.log1p(prior_arr[s]), dec=dec_arr[s]))
    T["lp2"] = T.lp ** 2
    try:
        pr = mod.predict(T)
    except Exception:
        return None
    return np.asarray(pr)


def enrich_row(t, thr, DRm, use_logit=True):
    meas176 = t.meas[s176_idx]
    n176 = int(meas176.sum())
    obs = int(((t.q[s176_idx] < thr) & meas176).sum())
    poolm = np.array([GI[g] for g in pool])
    pm = t.meas[poolm]
    n_bg = int(pm.sum())
    n_de_bg = int((t.q[poolm][pm] < thr).sum())
    hyper = float(stats.hypergeom.sf(obs - 1, n_bg + n176, n_de_bg + obs, n176))
    de = np.where(t.meas[DRm], t.q[DRm] < thr, False)
    cnt = de.sum(1)
    nm = t.meas[DRm].sum(1)
    # null counts are over drawn genes measured in t; rescale to the observed number measured
    null = cnt * (n176 / np.maximum(nm, 1))
    p_perm = (1 + (null >= obs - 1e-9).sum()) / (1 + len(null))
    row = dict(table=t.name, disorder=t.disorder, fdr_thr=thr, n_176_measured=n176, n_pool_measured=n_bg, n_pool_DE=n_de_bg,
               pool_DE_frac=round(n_de_bg / max(n_bg, 1), 4), obs_176_DE=obs, frac_176_DE=round(obs / max(n176, 1), 4),
               hypergeom_p=hyper, null_mean=round(float(null.mean()), 2), null_sd=round(float(null.std()), 2),
               fold_vs_null=round(obs / null.mean(), 2) if null.mean() > 0 else "", perm_p=round(p_perm, 4))
    ex = fit_expected(t, thr) if use_logit else None
    if ex is not None:
        row.update(logit_expected=round(float(ex.sum()), 2), logit_p=pb_tail(ex, obs))
    else:
        row.update(logit_expected="", logit_p="")
    return row


def dir_rows(t, thr, DRm):
    """sign concordance with the SZ direction, among the 176 that are DE in t (thr) and among all measured"""
    out = []
    s = np.array([i for i in s176_idx if sgn[i] != 0])
    for label, sel_fn in (("DE q<%.2f" % thr, lambda idx: (t.q[idx] < thr) & t.meas[idx]),
                          ("nominal p<0.05", lambda idx: (t.p[idx] < 0.05) & t.meas[idx]) if not np.isnan(t.p).all() else (None, None),
                          ("all measured", lambda idx: t.meas[idx])):
        if label is None:
            continue
        if label in ("nominal p<0.05", "all measured") and thr != 0.05:
            continue
        m = sel_fn(s)
        n = int(m.sum())
        conc = int(((np.sign(t.z[s][m]) == sgn[s][m])).sum())
        # null
        M = sel_fn(DRm)
        C = (np.sign(t.z[DRm]) == sgn[DRm]) & M
        n_null = M.sum(1)
        ok = n_null > 0
        frac = np.where(ok, C.sum(1) / np.maximum(n_null, 1), np.nan)
        fobs = conc / n if n else float("nan")
        pn = (1 + (frac[ok] >= fobs - 1e-9).sum()) / (1 + ok.sum()) if n else float("nan")
        bt_two = stats.binomtest(conc, n, 0.5).pvalue if n else float("nan")
        bt_one = stats.binomtest(conc, n, 0.5, alternative="greater").pvalue if n else float("nan")
        out.append(dict(table=t.name, disorder=t.disorder, set=label, n_176_with_SZ_sign_in_set=n, concordant=conc, discordant=n - conc,
                        frac_concordant=round(fobs, 3) if n else "", binom_p_two_sided_vs_50=bt_two, binom_p_one_sided_greater=bt_one,
                        null_frac_mean=round(float(np.nanmean(frac)), 3), null_frac_sd=round(float(np.nanstd(frac)), 3),
                        perm_p_frac_ge_obs=round(pn, 4) if n else ""))
    return out


# --- (1) enrichment and (2) direction for the main tables
enr, drr = [], []
for t in tables_main:
    for thr in (0.05, 0.10):
        enr.append(enrich_row(t, thr, DR))
        drr += dir_rows(t, thr, DR)
# per-list (Gemma), fewer permutations
DRl = DR[:NPERM_LIST]
lst_rows = []
for L, t in gem_rows:
    r = enrich_row(t, 0.05, DRl, use_logit=False)
    r.update(list_id=L["list_id"], n_case=L["n_case"], n_control=L["n_control"], genes_q05_in_list_file=L["genes_q05"])
    lst_rows.append(r)

# --- (3) SZ-aligned mean z and axes
def zalign(t, DRm):
    s = np.array([i for i in s176_idx if sgn[i] != 0 and t.meas[i]])
    obs = float(np.mean(sgn[s] * t.z[s]))
    Zd = t.z[DRm] * sgn[DRm]
    nn = np.nanmean(Zd, 1)
    return len(s), obs, nn


axes = {}
axes["control pH Z (low pH = negative)"] = ({}, -1)
for l in open(PH_Z):
    if l.startswith("#") or l.startswith("gene\t"):
        continue
    f = l.rstrip("\n").split("\t")
    if f[1] == "pH" and f[0] in GI:
        try:
            axes["control pH Z (low pH = negative)"][0][f[0]] = float(f[3])
        except ValueError:
            pass
Hd = pd.read_csv(GTEX, sep="\t", index_col=0)
ag = Hd[["Brain - Frontal Cortex (BA9)", "Brain - Anterior cingulate cortex (BA24)", "Brain - Cortex"]].mean(axis=1, skipna=False).dropna().to_dict()
axes["GTEx agonal t (slow vs fast death)"] = ({g: v for g, v in ag.items() if g in GI}, +1)
# Hatano and Dai, if pyxlsb is importable
try:
    from pyxlsb import open_workbook
    tmp = os.path.join(HERE, "_tmp_hatano")
    os.makedirs(tmp, exist_ok=True)
    with zipfile.ZipFile(HAT) as zz:
        zz.extractall(tmp)
    for f, name in (("0001", "BD"), ("0004", "SZ"), ("0003", "control")):
        p = glob.glob(os.path.join(tmp, f"*{f}*"))[0]
        with open_workbook(p) as wb, wb.get_sheet(1) as sh:
            rr = [[c.v for c in r] for r in sh.rows()]
        hh = rr[1]
        gi, pi = hh.index("GeneSymbol"), hh.index("pH")
        acc = collections.defaultdict(list)
        for r in rr[2:]:
            if len(r) > pi and r[gi] and isinstance(r[pi], (int, float)):
                acc[str(r[gi]).strip()].append(float(r[pi]))
        axes[f"Hatano {name} brains pH r"] = ({g: float(np.mean(v)) for g, v in acc.items() if g in GI}, -1)
    import shutil
    shutil.rmtree(tmp)
    HAT_OK = True
except Exception as e:
    print("Hatano axes skipped:", e)
    HAT_OK = False

rng2 = np.random.default_rng(SEED)
def axis_rho(t, ax, orient, DRm, nperm):
    a = np.full(len(ALL), np.nan)
    for g, v in ax.items():
        a[GI[g]] = v
    a = a * orient          # positive = disorder-like direction if gene is up in the disorder
    s = s176_idx[t.meas[s176_idx] & ~np.isnan(a[s176_idx])]
    if len(s) < 10:
        return None
    obs = stats.spearmanr(t.z[s], a[s])[0]
    nulls = []
    for i in range(min(nperm, len(DRm))):
        d = DRm[i]
        d = d[t.meas[d] & ~np.isnan(a[d])]
        nulls.append(stats.spearmanr(t.z[d], a[d])[0] if len(d) > 10 else np.nan)
    nulls = np.array(nulls)
    nulls = nulls[~np.isnan(nulls)]
    return len(s), obs, float(nulls.mean()), float(nulls.std()), (obs - nulls.mean()) / nulls.std(), (1 + (nulls >= obs).sum()) / (1 + len(nulls))

axrows, zrows = [], []
for t in tables_main:
    n, obs, nn = zalign(t, DR)
    zrows.append(dict(table=t.name, disorder=t.disorder, n_176_with_sign_and_value=n, mean_SZ_aligned_z=round(obs, 3),
                      null_mean=round(float(np.nanmean(nn)), 3), null_sd=round(float(np.nanstd(nn)), 3),
                      z_vs_null=round((obs - np.nanmean(nn)) / np.nanstd(nn), 2), p_one_sided=round((1 + (nn >= obs).sum()) / (1 + len(nn)), 4)))
    for an, (ax, orient) in axes.items():
        r = axis_rho(t, ax, orient, DR, NPERM_AXIS)
        if r:
            axrows.append(dict(table=t.name, disorder=t.disorder, axis=an, n_genes=r[0], spearman_rho=round(r[1], 3), null_mean=round(r[2], 3),
                               null_sd=round(r[3], 3), z=round(r[4], 2), p_one_sided_rho_ge_obs=round(r[5], 4)))


# ------------------------------------------------------------------ Hagihara Table S1
raw = pd.read_excel(HAG, sheet_name="TableS1", header=None, skiprows=3)
raw = raw.iloc[:, :9]
raw.columns = ["dataset", "id", "dx_raw", "sex", "pH", "pmi", "age", "ap_life", "ap_final"]
raw = raw[raw.dataset.notna() & raw.dx_raw.notna()].copy()  # Dean 2016 rows have no ID
for c in ("pH", "pmi", "age"):
    raw[c] = pd.to_numeric(raw[c], errors="coerce")
mp = {"Control": "CTL", "CON": "CTL", "Schizophrenia": "SZ", "SCZ": "SZ", "SCZ-A": "SZ", "Schizophrenia (MRDS)": "SZ", "Schizophrenia (non-MRDS)": "SZ", "Bipolar": "BD"}
raw["dx"] = raw.dx_raw.astype(str).str.strip().map(mp)
H = raw[raw.dx.notna() & raw.pH.notna()].copy()
hrows = []
def hag_model(d, label):
    d = d.copy()
    d["dx"] = pd.Categorical(d.dx, categories=["CTL", "SZ", "BD"] if "BD" in set(d.dx) else ["CTL", "SZ"])
    m = smf.ols("pH ~ C(dx) + C(dataset)", d).fit()
    sd = math.sqrt(m.scale)
    for dx in ("SZ", "BD"):
        k = f"C(dx)[T.{dx}]"
        if k in m.params:
            ci = m.conf_int().loc[k]
            hrows.append(dict(model=label, dx=dx, n_dx=int((d.dx == dx).sum()), n_control=int((d.dx == "CTL").sum()), n_datasets=d.dataset.nunique(),
                              coef_pH_units=round(m.params[k], 3), ci_low=round(ci[0], 3), ci_high=round(ci[1], 3), t=round(m.tvalues[k], 2), p=m.pvalues[k],
                              residual_sd=round(sd, 3), coef_in_residual_sd=round(m.params[k] / sd, 2)))
    return m
hag_model(H, "all rows, dataset fixed effect")
hag_model(H[H.dataset != "Dean et al, 2016"], "without Dean 2016")
both_ds = [ds for ds in H.dataset.unique() if {"SZ", "BD", "CTL"} <= set(H[H.dataset == ds].dx)]
Hb = H[H.dataset.isin(both_ds)]
mb = hag_model(Hb, "datasets with SZ, BD and controls (" + ", ".join(sorted(both_ds)) + ")")
# SZ vs BD in those shared-control datasets
wt = mb.t_test("C(dx)[T.SZ] - C(dx)[T.BD] = 0")
hrows.append(dict(model="SZ minus BD in the shared-control datasets", dx="SZ-BD", n_dx="", n_control="", n_datasets=len(both_ds),
                  coef_pH_units=round(float(np.squeeze(wt.effect)), 3), ci_low=round(float(wt.conf_int()[0][0]), 3), ci_high=round(float(wt.conf_int()[0][1]), 3),
                  t=round(float(np.squeeze(wt.tvalue)), 2), p=float(np.squeeze(wt.pvalue)), residual_sd="", coef_in_residual_sd=""))
# per-dataset raw differences
for ds, g in H.groupby("dataset"):
    c = g[g.dx == "CTL"].pH
    for dx in ("SZ", "BD"):
        p = g[g.dx == dx].pH
        if len(c) > 2 and len(p) > 2:
            hrows.append(dict(model=f"per dataset: {ds}", dx=dx, n_dx=len(p), n_control=len(c), n_datasets=1, coef_pH_units=round(p.mean() - c.mean(), 3),
                              ci_low="", ci_high="", t=round(stats.ttest_ind(p, c, equal_var=False).statistic, 2), p=stats.ttest_ind(p, c, equal_var=False).pvalue,
                              residual_sd="", coef_in_residual_sd=""))

# Hatano BD row from score_176.tsv (read, not recomputed)
sc = tsv(SCORE176)
hat_bd = [r for r in sc if r["axis"].startswith("Hatano BD")]


# ------------------------------------------------------------------ write
def write(name, rows, header_lines):
    cols = list(rows[0])
    for r in rows:
        for c in r:
            if c not in cols:
                cols.append(c)
    with open(os.path.join(HERE, name), "w") as f:
        for h in header_lines:
            f.write("# " + h + "\n")
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(("%.3g" % r.get(c) if isinstance(r.get(c), float) and c.endswith("_p") else str(r.get(c, ""))) for c in cols) + "\n")

common = [
    "script: analyses/other_disorders_recurrent_genes/other_disorders_176.py, run 2026-10-09; seed %d" % SEED,
    f"population: recurrence_genes.tsv sha256 {sha(REC)}, scenario 'published, brain', q_bh < 0.1 = 176 genes; SZ sign = sign(up - down) (175 with a sign); DE prior = de_prior_hitlists column (n_de_total of de_prior_lists/our_prior/1_global_our_prior.tsv sha256 {sha(PRIOR)})",
    f"null pool: {len(pool)} published-brain genes with a SZ sign, excluding the 176; strata = expr_decile x DE-prior quintile (quintile cuts over the pool); {NPERM_MAIN} matched sets (per-list rows: {NPERM_LIST})",
    f"Gemma autism result sets: gemma_asd_cache/ (fetched 2026-10-10 GET only; gemma_asd_selection.tsv sha256 {sha(ASD_SEL)}); "
    f"inputs: Gandal2018 array TableS1 csv sha256 {sha(G18A_S1)}; Gandal2018 PsychENCODE aat8127_table_s1.xlsx sha256 {sha(G18B_S1)}; Gemma lists check1_lists.tsv sha256 {sha(GEM_LISTS)} (result-set tables in checks_2026-09-29/cache/)",
]
write("enrichment_176_by_disorder.tsv", enr, common + ["hypergeom_p: background = pool genes measured in the table; logit_*: DE ~ log1p(prior) + log1p(prior)^2 + C(expr_decile) fit on the pool, Poisson-binomial upper tail; perm_p: matched null counts rescaled to the number of the 176 measured"])
write("enrichment_176_gemma_per_list.tsv", lst_rows, common + ["one row per Gemma BD or MDD postmortem-brain list at q < 0.05 (BH over genes, minimum-p probe per gene); perm_p from the first %d matched sets; tests are not corrected across lists" % NPERM_LIST])
write("direction_176_by_disorder.tsv", drr, common + ["concordant = sign of the disorder's effect equals the SZ sign; binomial against 50%; perm_p = proportion of matched sets with concordant fraction >= observed (null genes carry their own SZ sign)"])
write("szaligned_z_176_by_disorder.tsv", zrows, common + ["mean_SZ_aligned_z = mean over the 176 (with a SZ sign) of SZ sign x the disorder's signed z (z from the table's p and effect sign)"])
write("axes_within_disorder_176.tsv", axrows, common + [
    f"inputs: combined_controls_pH.tsv sha256 {sha(PH_Z)}; hardy_t_by_region.tsv.gz sha256 {sha(GTEX)}; Hatano zip sha256 {sha(HAT)}",
    "spearman_rho between the disorder's signed z and the axis oriented so positive = the gene is higher in the disorder when pH is lower (pH axes) or when death is faster (agonal); null = same rho in matched sets (%d)" % NPERM_AXIS])
write("hagihara_table_s1_pH_by_diagnosis.tsv", hrows, [
    "script: other_disorders_176.py, run 2026-10-09", f"input: data/raw/28776581/supplement/npp2017167x2.xlsx sha256 {sha(HAG)}, sheet TableS1; SCZ-A and Dean 2016 SZ arms counted as SZ; NOS (1 row) dropped",
    "pH ~ diagnosis + dataset (OLS); coef = difference from controls in pH units"])
if hat_bd:
    write("hatano_bd_ph_score_from_score_176.tsv", hat_bd, ["copied from analyses/agonal_axes_score_176_2026-10-09/score_176.tsv (not recomputed); file sha256 %s" % sha(SCORE176)])

for name, rows_ in (("ENRICH", enr), ("DIR", drr), ("Z", zrows), ("AXES", axrows), ("HAG", hrows)):
    print("=====", name)
    print(pd.DataFrame(rows_).to_string())
print("Gemma lists:", len(gem_rows), {d: len(v) for d, v in gem_tabs.items()})
