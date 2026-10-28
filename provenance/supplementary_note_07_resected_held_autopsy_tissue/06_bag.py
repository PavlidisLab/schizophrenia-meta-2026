#!/usr/bin/env python3
"""Brain Artifact Genes (BAGs; Supplementary Data 4 = genes common to SPMI and LPMI vs CTRL; 857 up, 425 down as released) against the
recurrent genes and the axes, plus the GO themes the paper names (oxidative phosphorylation / ATP synthesis; response to LPS).
Set-up (population, strata, matched sets, axes) copied from 02_axes_and_sets.py. Usage: python3 06_bag.py <out dir> <dir with moesm3.xlsx>   Seed 20261018.
"""
import collections, gzip, hashlib, os, random, sys, datetime, math
sys.dont_write_bytecode = True
import numpy as np
import pandas as pd
from scipy import stats

OUT, SUPP = sys.argv[1], sys.argv[2]
REPO = os.path.expanduser("<analysis repo>")
AN = os.path.join(REPO, "analyses")
sys.path.insert(0, AN)
A = os.path.expanduser("~/analysis-delivery")
B80 = f"{A}/base80_recurrence_2026-10-16/out"
REC, SETS = f"{B80}/recurrence_genes.tsv", f"{B80}/recurrent_sets.tsv"
PHAX = f"{A}/ph_axis_v2_2026-10-14/axis_v2_controls.tsv.gz"
GTEX = f"{AN}/gtex_ischemia_2026-10-08/hardy_t_by_region.tsv.gz"
SPEC = f"{AN}/interneuron_ph_2026-10-01/results/specificity_HBCC.tsv.gz"
EFF = f"{OUT}/tables/contrast_effects.tsv"
SEED, NPERM, NPERM_RHO, NBOOT = 20261018, 10000, 5000, 1000
GOF = os.path.join(AN, "fairness_lists", "go_bp_sets.tsv.gz")
TGG = os.path.expanduser("M09/results_seedA_noextra_new/term_group_genes.tsv")
EFF = ""
import extract
extract.HGNC = extract.load_hgnc()
HGc = {}


def hg(s):
    if s not in HGc:
        r = extract.symbol_to_hgnc(str(s))
        HGc[s] = r[0] or None
    return HGc[s]


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()[:12]


rec = pd.read_csv(REC, sep="\t", comment="#")
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
assert len(S383) == 383
sets = pd.read_csv(SETS, sep="\t", comment="#")
assert set(sets[sets.set == "study_brain"].gene) == set(S383)
S107 = list(sets[sets.set == "study_brain_ge10"].gene)
assert len(S107) == 107
ALL = sorted(sign)
GI = {g: i for i, g in enumerate(ALL)}
NG = len(ALL)
sgn = np.array([sign[g] for g in ALL], float)
pool = [g for g in sign if sign[g] != 0 and g not in set(S383)]
qcuts = np.quantile([prior[g] for g in pool], [0.2, 0.4, 0.6, 0.8])
stratum = {g: (dec[g], int(np.searchsorted(qcuts, prior[g], side="right"))) for g in sign}
by_st = collections.defaultdict(list)
for g in pool:
    by_st[stratum[g]].append(g)


class GS:
    def __init__(self, name, genes, seed):
        self.name, self.genes = name, genes
        self.idx = np.array([GI[g] for g in genes])
        need = collections.Counter(stratum[g] for g in genes)
        rng = random.Random(seed)
        out = np.zeros((NPERM, len(genes)), dtype=np.int32)
        for i in range(NPERM):
            row = []
            for s, k in need.items():
                c = by_st[s]
                row += rng.sample(c, k) if len(c) >= k else rng.choices(c, k=k)
            out[i] = [GI[g] for g in row]
        self.DR = out
        self.short = {s: k for s, k in need.items() if len(by_st[s]) < k}


GSETS = [GS("383 recurrent (q_bh<0.1)", S383, SEED), GS("107 in >=10 lists", S107, SEED + 1)]
print("sets built; strata short of pool genes:", [g.short for g in GSETS])

# ---------------------------------------------------------------- axes (HGNC symbols already)
ph = pd.read_csv(PHAX, sep="\t", comment="#", index_col=0)["z_stouffer"].to_dict()
Hd = pd.read_csv(GTEX, sep="\t", index_col=0)
ag = Hd[["Brain - Frontal Cortex (BA9)", "Brain - Anterior cingulate cortex (BA24)", "Brain - Cortex"]].mean(axis=1, skipna=False).dropna().to_dict()
AXES = {"pH axis v2": ph, "agonal axis (GTEx slow-minus-fast)": ag}


def full(d):
    a = np.full(NG, np.nan)
    for g, v in d.items():
        if g in GI:
            a[GI[g]] = v
    return a


AX_ALL = {k: pd.Series(v) for k, v in AXES.items()}
AX_U = {k: full(v) for k, v in AXES.items()}


X3 = os.path.join(SUPP, "moesm3.xlsx")
s4 = pd.read_excel(X3, sheet_name="Supplementary Data 4", header=1, engine="openpyxl")
up = s4.iloc[:, :3].dropna(subset=[s4.columns[0]]); up.columns = ["gene", "lfc_SPMI", "lfc_LPMI"]
dn = s4.iloc[:, 4:7].dropna(subset=[s4.columns[4]]); dn.columns = ["gene", "lfc_SPMI", "lfc_LPMI"]
up["dir"], dn["dir"] = 1, -1
bag = pd.concat([up, dn])
bag["g"] = [hg(x) for x in bag.gene]
bag = bag.dropna(subset=["g"]).drop_duplicates("g").set_index("g")
bag["lfc_mean"] = bag[["lfc_SPMI", "lfc_LPMI"]].astype(float).mean(axis=1)
print("BAG released rows", len(up), len(dn), "mapped unique", len(bag), "up", int((bag.dir > 0).sum()), "down", int((bag.dir < 0).sum()))
BU = np.zeros(NG, bool); BD = np.zeros(NG, bool)
for g, r in bag.iterrows():
    if g in GI:
        (BU if r.dir > 0 else BD)[GI[g]] = True
rows = []
def wil(k, n, z=1.96):
    p = k / n; c = (p + z*z/(2*n)) / (1 + z*z/n); h = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / (1 + z*z/n); return c-h, c+h
for gs in GSETS:
    for nm, M, d_ in (("BAG up", BU, 1), ("BAG down", BD, -1), ("BAG all", BU | BD, 0)):
        obs = int(M[gs.idx].sum()); nl = M[gs.DR].sum(1)
        row = dict(gene_set=gs.name, bag_set=nm, bag_genes_in_universe=int(M.sum()), overlap=obs, null_mean=round(float(nl.mean()), 2), null_q025=float(np.quantile(nl, .025)), null_q975=float(np.quantile(nl, .975)),
                   fold=round(obs / nl.mean(), 2), fold_null_q025=round(float(np.quantile(nl, .025)) / nl.mean(), 2), fold_null_q975=round(float(np.quantile(nl, .975)) / nl.mean(), 2))
        # direction: SZ sign vs BAG direction among overlap genes (for up/down sets the BAG direction is fixed; for all, per-gene direction)
        dirv = np.where(BU, 1, np.where(BD, -1, 0))
        s = np.array([i for i in gs.idx if M[i] and sgn[i] != 0])
        conc = int((sgn[s] == dirv[s]).sum()); n = len(s)
        D = gs.DR; Mk = M[D] & (sgn[D] != 0); Ck = Mk & (sgn[D] == dirv[D])
        fr = Ck.sum(1) / np.maximum(Mk.sum(1), 1)
        lo, hi = wil(conc, n) if n else (np.nan, np.nan)
        row.update(n_with_SZ_sign=n, concordant=conc, frac_concordant=round(conc / n, 3) if n else "", wilson95_low=round(lo, 3), wilson95_high=round(hi, 3), null_frac_mean=round(float(fr.mean()), 3),
                   null_frac_q025=round(float(np.quantile(fr, .025)), 3), null_frac_q975=round(float(np.quantile(fr, .975)), 3))
        rows.append(row)
BAGOV = pd.DataFrame(rows)

# BAG log2FC pattern vs axes
rng = np.random.default_rng(SEED)
rk = stats.rankdata
rows = []
for an, ax in AXES.items():
    a = full(ax)
    for scope, genes_ in (("all BAG genes", list(bag.index)), ("BAG genes in the 383", [g for g in bag.index if g in set(S383)]), ("BAG genes in the 107", [g for g in bag.index if g in set(S107)])):
        for lab in ("lfc_mean", "lfc_SPMI", "lfc_LPMI"):
            gg = [g for g in genes_ if g in GI and not np.isnan(a[GI[g]])]
            x = bag.loc[gg, lab].astype(float).values; y = a[[GI[g] for g in gg]]
            if len(gg) < 8: continue
            rs, rp = stats.spearmanr(x, y)[0], np.corrcoef(x, y)[0, 1]
            bs = [np.corrcoef(rk(x[i]), rk(y[i]))[0, 1] for i in (rng.integers(0, len(x), len(x)) for _ in range(NBOOT))]
            bp = [np.corrcoef(x[i], y[i])[0, 1] for i in (rng.integers(0, len(x), len(x)) for _ in range(NBOOT))]
            rows.append(dict(axis=an, scope=scope, bag_lfc=lab, n_genes=len(gg), spearman=round(rs, 3), spearman_boot95_low=round(np.quantile(bs, .025), 3), spearman_boot95_high=round(np.quantile(bs, .975), 3),
                             pearson=round(rp, 3), pearson_boot95_low=round(np.quantile(bp, .025), 3), pearson_boot95_high=round(np.quantile(bp, .975), 3)))
BAGAX = pd.DataFrame(rows)

# GO themes
goterms = {"oxidative phosphorylation / ATP synthesis (14 GO BP terms of the M09 oxphos group)": ["GO:0009142", "GO:0009145", "GO:0022900", "GO:0015986", "GO:0042776", "GO:0042775", "GO:0042773",
           "GO:0006754", "GO:0009206", "GO:0006119", "GO:0022904", "GO:0009201", "GO:0019646", "GO:0006120"],
           "response to lipopolysaccharide (GO:0032496)": ["GO:0032496"]}
gm = {}
for l in gzip.open(GOF, "rt"):
    if l.startswith("#") or l.startswith("term\t"): continue
    f = l.rstrip("\n").split("\t"); gm[f[0]] = {hg(x) for x in f[2].split(",") if hg(x)}
phs = pd.Series(ph).dropna()
top_hi = set(phs.sort_values(ascending=False).index[:200]); top_lo = set(phs.sort_values().index[:200])
tg = pd.read_csv(TGG, sep="\t", comment="#")
ox_terms = set(tg[(tg.group == "oxphos") & (tg.version == "terms")].gene); ox_hits = set(tg[(tg.group == "oxphos") & (tg.version == "hits")].gene)
bag_up, bag_dn = {g for g in bag.index if bag.loc[g, "dir"] > 0}, {g for g in bag.index if bag.loc[g, "dir"] < 0}
bg = set(phs.index)
def ov(A, B):
    A, B = A & bg, B & bg
    k = len(A & B); e = len(A) * len(B) / len(bg)
    return len(A), len(B), k, round(e, 2), round(k / e, 2) if e > 0 else ""
rows = []
themes = {}
for tn, ids in goterms.items():
    mem = set().union(*[gm[i] for i in ids])
    themes[tn] = mem
cases = {"BAG up (all)": bag_up, "BAG down (all)": bag_dn}
for tn, mem in themes.items():
    cases[f"BAG up and annotated: {tn}"] = bag_up & mem
    cases[f"BAG down and annotated: {tn}"] = bag_dn & mem
for cn, A in cases.items():
    for bn, B in (("pH-axis top 200 (highest z: higher at higher pH)", top_hi), ("pH-axis bottom 200 (lowest z: higher at lower pH)", top_lo),
                  ("oxphos term group, GO members (M09 'terms')", ox_terms), ("oxphos term group, genes in at least one list (M09 'hits')", ox_hits),
                  ("383 recurrent", set(S383)), ("107 recurrent", set(S107))):
        a, b, k, e, f_ = ov(A, B)
        rows.append(dict(bag_subset=cn, compared_with=bn, n_subset_in_axis_universe=a, n_other_in_axis_universe=b, overlap=k, expected_if_independent=e, fold=f_, background_genes=len(bg)))
GOT = pd.DataFrame(rows)
H = [f"# script: 06_bag.py, run {datetime.date.today().isoformat()}; NEW analysis (addition); seed {SEED}; {NPERM} matched sets; {NBOOT} bootstrap resamples over genes",
     f"# inputs: MOESM3_ESM.xlsx Supplementary Data 4 sha256 {sha(X3)}; axis_v2_controls.tsv.gz {sha(PHAX)}; hardy_t_by_region.tsv.gz {sha(GTEX)}; recurrence_genes.tsv {sha(REC)}; recurrent_sets.tsv {sha(SETS)}; go_bp_sets.tsv.gz {sha(GOF)}; term_group_genes.tsv {sha(TGG)}",
     "# BAG direction as released (up/down in SPMI and LPMI vs CTRL); expected overlap for GO/axis rows = hypergeometric expectation with the genes of the pH axis as background (not DE-prior matched); fold for the recurrent sets in bag_vs_recurrent.tsv is over the matched null"]
def w(df, name):
    with open(f"{OUT}/tables/{name}", "w") as f:
        f.write("\n".join(H) + "\n"); df.to_csv(f, sep="\t", index=False)
w(BAGOV, "bag_vs_recurrent.tsv"); w(BAGAX, "bag_pattern_axis_correlation.tsv"); w(GOT, "bag_go_themes_overlap.tsv")
for n_, d_ in (("mapping", None),):
    pass
pd.set_option("display.width", 250, "display.max_columns", 40)
print(BAGOV.to_string()); print(BAGAX.to_string()); print(GOT.to_string())
