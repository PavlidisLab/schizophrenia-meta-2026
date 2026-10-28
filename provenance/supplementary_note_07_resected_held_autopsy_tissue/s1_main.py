#!/usr/bin/env python3
"""Is the living-versus-postmortem (LIV-PM) prefrontal-cortex expression difference of Liharska et al. 2025
aligned with the recurrent schizophrenia (SZ) pattern (the 383 genes of the core-80 base; the 107 in >= 10 lists), the bulk pH axis v2 and
the GTEx agonal axis? Same matched-null design as ad_383.py (expression decile x DE-prior quintile strata).

Signatures (stat = per-gene statistic; sign convention: positive = higher in postmortem for LIV-PM variants, higher in cases for disease signatures):
  LIV-PM  : Liharska ST1 (MOESM8, moderated t), Collado-Torres 2023 Table S4 (original / without qSVA / with qSVA; limma-voom t)
  disease : Liharska ST5 (MOESM12) ALZ1, ALZ2, SCZ1, SCZ2, PD1, PD2 (logFC and P only: stat = sign(logFC) * qnorm(1 - P/2))
Axes (raw, not re-oriented): pH axis v2 controls z_stouffer (positive = higher at higher pH); GTEx agonal = mean over BA9, BA24, Cortex of the
slow(Hardy 4)-versus-fast-natural(2) death t (positive = higher after slow death).
Fixed seed 20261018. Writes only into the parent folder.   python3 -I s1_main.py
"""
import sys, os, csv, gzip, hashlib, math, random, collections, datetime
sys.dont_write_bytecode = True
import numpy as np, pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.dirname(HERE); DATA = os.path.join(OUT, "data")
REPO = os.path.expanduser("<analysis repo>"); AN = os.path.join(REPO, "analyses")
sys.path.insert(0, AN)
SEED, NSETS, NBOOT_SET, NBOOT_GW = 20261018, 10000, 2000, 1000
rng = np.random.default_rng(SEED)
B80 = os.path.expanduser("out")
REC, SETS = f"{B80}/recurrence_genes.tsv", f"{B80}/recurrent_sets.tsv"
PRIOR = f"{AN}/de_prior_lists/our_prior/1_global_our_prior.tsv"
PHV2 = os.path.expanduser("axis_v2_controls.tsv.gz")
GTEX = f"{AN}/gtex_ischemia_2026-10-08/hardy_t_by_region.tsv.gz"
SPEC = f"{AN}/interneuron_ph_2026-10-01/results/specificity_HBCC.tsv.gz"
ST1, ST5 = f"{DATA}/liv_pm_de_ST1.tsv", f"{DATA}/disease_de_ST5.tsv"
CT4 = f"{DATA}/collado_torres/TableS4.csv"
TODAY = "2026-10-18"


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()[:12]


def tsv(p):
    return list(csv.DictReader((l for l in open(p, encoding="utf8") if not l.startswith("#")), delimiter="\t"))


import extract
extract.HGNC = extract.load_hgnc()
ENS2SYM = {e: s for e, (s, _) in extract.HGNC[3].items()}
_hgc = {}


def hg(s):
    if s not in _hgc:
        _hgc[s] = extract.symbol_to_hgnc(s)[0] or None
    return _hgc[s]


# ---------------------------------------------------------------- our recurrence population (identical to ad_383.py)
rec = [r for r in tsv(REC) if r["scenario"] == "published, brain"]
sign, dec, prior = {}, {}, {}
for r in rec:
    d = int(r["up"]) - int(r["down"])
    sign[r["gene"]] = (1 if d > 0 else -1) if d else 0
    dec[r["gene"]] = r["expr_decile"]
    prior[r["gene"]] = float(r["de_prior_hitlists"]) if r["de_prior_hitlists"] != "" else np.nan
med = float(np.nanmedian(list(prior.values())))
prior = {g: (med if np.isnan(v) else v) for g, v in prior.items()}
S383 = [r["gene"] for r in rec if float(r["q_bh"]) < 0.1]
assert len(S383) == 383
sets_rows = tsv(SETS)
assert {r["gene"] for r in sets_rows if r["set"] == "study_brain"} == set(S383)
S107 = [r["gene"] for r in sets_rows if r["set"] == "study_brain_ge10"]
assert len(S107) == 107
ALL = sorted(sign); GI = {g: i for i, g in enumerate(ALL)}
sgn = np.array([sign[g] for g in ALL], float)
prior_arr = np.array([prior[g] for g in ALL]); dec_arr = np.array([dec[g] for g in ALL])
pool = [g for g in sign if sign[g] != 0 and g not in set(S383)]
qcuts = np.quantile([prior[g] for g in pool], [0.2, 0.4, 0.6, 0.8])
stratum = {g: (dec[g], int(np.searchsorted(qcuts, prior[g], side="right"))) for g in sign}
by_st = collections.defaultdict(list)
for g in pool:
    by_st[stratum[g]].append(g)
pool_idx = np.array([GI[g] for g in pool])
log = open(os.path.join(OUT, "logs_s1_main.txt"), "w")


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s); log.write(s + "\n"); log.flush()


P("383:", len(S383), "107:", len(S107), "signs 0 in 383:", sum(sign[g] == 0 for g in S383), "universe:", len(ALL), "pool:", len(pool))


class GeneSet:
    def __init__(self, name, genes, seed):
        self.name, self.genes = name, genes
        self.idx = np.array([GI[g] for g in genes])
        need = collections.Counter(stratum[g] for g in genes)
        short = {s: k for s, k in need.items() if len(by_st[s]) < k}
        r = random.Random(seed)
        out = np.zeros((NSETS, len(genes)), dtype=np.int32)
        for i in range(NSETS):
            row = []
            for s, k in need.items():
                c = by_st[s]
                row += r.sample(c, k) if len(c) >= k else r.choices(c, k=k)
            out[i] = [GI[g] for g in row]
        self.DR = out
        P(name, len(genes), "strata", len(need), "strata short of pool genes:", len(short))


def sub(genes, s):
    return [g for g in genes if sign[g] == s]


GSETS = {}
for k, (nm, genes) in enumerate([("383", S383), ("107", S107)]):
    GSETS[nm] = GeneSet(nm, [g for g in genes], SEED + 10 * k)
    GSETS[nm + "_up"] = GeneSet(nm + "_up", sub(genes, 1), SEED + 10 * k + 1)
    GSETS[nm + "_down"] = GeneSet(nm + "_down", sub(genes, -1), SEED + 10 * k + 2)
P({k: len(v.genes) for k, v in GSETS.items()})

# ---------------------------------------------------------------- signatures
class Sig:
    def __init__(self, name, kind, stat, q, effect=None, src=""):
        self.name, self.kind, self.src = name, kind, src
        self.stat = stat          # Series by approved symbol
        self.q = q                # Series of BH q (same index)
        self.effect = effect      # Series (logFC) or None
        self.ALL_stat = stat.reindex(ALL).values
        self.ALL_q = q.reindex(ALL).values
        self.meas = ~np.isnan(self.ALL_stat)


def zfrom(p, beta):
    p = np.clip(np.asarray(p, float), 1e-300, 1.0)
    return np.sign(beta) * stats.norm.isf(p / 2)


def best_per_symbol(df, sym, p):
    df = df.assign(_s=sym).dropna(subset=["_s", p])
    df = df.sort_values(p).drop_duplicates("_s")
    return df.set_index("_s")


sigs = []
meta = {}
# Liharska ST1
l = pd.read_csv(ST1, sep="\t", comment="#")
l["sym"] = l.gene.map(ENS2SYM)
meta["ST1_rows"], meta["ST1_unmapped_ENSG"] = len(l), int(l.sym.isna().sum())
l1 = best_per_symbol(l, l.sym, "P.Value")
meta["ST1_symbols"] = len(l1)
sigs.append(Sig("LIV-PM Liharska (ST1)", "livpm", l1.t, l1["adj.P.Val"], l1.logFC, "MOESM8 sheet Data"))
# Collado-Torres S4
c = pd.read_csv(CT4, low_memory=False)
c = c[c.isExprs]
c["sym"] = c.ensemblID.map(ENS2SYM)
c["sym"] = c.sym.where(c.sym.notna(), c.Symbol.map(lambda s: hg(str(s)) if isinstance(s, str) else None))
meta["CT_rows_expressed"], meta["CT_unmapped"] = len(c), int(c.sym.isna().sum())
for tag, label in (("original", "LIV-PM Collado-Torres reproduction of Liharska"), ("without_qSVA", "LIV-PM Collado-Torres, no qSVA"), ("with_qSVA", "LIV-PM Collado-Torres, with qSVA")):
    pc = f"{tag}_P.Value"
    ci = best_per_symbol(c, c.sym, pc)
    ci = ci.dropna(subset=[f"{tag}_t"])
    sigs.append(Sig(label, "livpm", ci[f"{tag}_t"], ci[f"{tag}_adj.P.Val"], ci[f"{tag}_logFC"], f"Collado-Torres Table S4 columns {tag}_*"))
# disease ST5
d = pd.read_csv(ST5, sep="\t", comment="#")
d["sym"] = d.gene.map(ENS2SYM)
meta["ST5_rows"], meta["ST5_unmapped_ENSG"] = len(d), int(d.sym.isna().sum())
for nm in ("SCZ1", "SCZ2", "ALZ1", "ALZ2", "PD1", "PD2"):
    x = best_per_symbol(d, d.sym, f"{nm} P.Value").dropna(subset=[f"{nm} logFC"])
    z = pd.Series(zfrom(x[f"{nm} P.Value"].values, x[f"{nm} logFC"].values), index=x.index)
    sigs.append(Sig(nm, "disease", z, x[f"{nm} adj.P.Val"], x[f"{nm} logFC"], f"MOESM12 {nm} columns; stat = signed z from P"))
for s in sigs:
    P(s.name, "genes", int(s.meas.sum()), "of universe;", "genes with q<0.05:", int((s.q < 0.05).sum()), "(up", int(((s.q < 0.05) & (s.stat > 0)).sum()), "down", int(((s.q < 0.05) & (s.stat < 0)).sum()), ")")
P(meta)

# ---------------------------------------------------------------- axes
ph = pd.read_csv(PHV2, sep="\t", comment="#", index_col=0).z_stouffer
ph.index = [hg(s) or s for s in ph.index]
ph = ph.groupby(level=0).first()
Hd = pd.read_csv(GTEX, sep="\t", index_col=0)
ag = Hd[["Brain - Frontal Cortex (BA9)", "Brain - Anterior cingulate cortex (BA24)", "Brain - Cortex"]].mean(axis=1, skipna=False).dropna()
ag.index = [hg(s) or s for s in ag.index]
ag = ag.groupby(level=0).first()
AXES = {"pH axis v2 (positive = higher at higher pH)": ph, "agonal axis (positive = higher after slow death)": ag}
AXSHORT = {k: ("pH" if k.startswith("pH") else "agonal") for k in AXES}
P("pH genes", len(ph), "agonal genes", len(ag))

# axes-removed variants: normal-score of ranks, OLS on both axes over all genes with all three, residual = stat with axes removed
def nscore(v):
    r = stats.rankdata(v); return stats.norm.ppf((r - 0.5) / len(v))


for s in list(sigs):
    if s.name in ("LIV-PM Liharska (ST1)", "LIV-PM Collado-Torres, no qSVA", "LIV-PM Collado-Torres, with qSVA", "SCZ1", "SCZ2"):
        df = pd.concat([s.stat.rename("y"), ph.rename("ph"), ag.rename("ag")], axis=1, join="inner").dropna()
        Y, X = nscore(df.y.values), np.column_stack([np.ones(len(df)), nscore(df.ph.values), nscore(df.ag.values)])
        beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
        res = pd.Series(Y - X @ beta, index=df.index)
        r2 = 1 - res.var() / pd.Series(Y).var()
        meta[f"R2_axes_on_{s.name}"] = round(float(r2), 4)
        P(f"{s.name}: R2 of normal-scored stat on pH + agonal normal scores (n={len(df)}): {r2:.4f}")
        sigs.append(Sig(s.name + " | pH and agonal axes removed", "residual", res, s.q.reindex(res.index), None, f"residual of normal scores on pH and agonal normal scores, n={len(df)}"))
sigd = {s.name: s for s in sigs}

# ---------------------------------------------------------------- stats helpers
def spearman_boot(x, y, nb, rg):
    ok = ~(np.isnan(x) | np.isnan(y)); x, y = x[ok], y[ok]
    n = len(x); rho = stats.spearmanr(x, y)[0]; pr = stats.pearsonr(x, y)[0]
    bs, bp = [], []
    for _ in range(nb):
        i = rg.integers(0, n, n)
        bs.append(stats.spearmanr(x[i], y[i])[0]); bp.append(stats.pearsonr(x[i], y[i])[0])
    return n, rho, np.percentile(bs, [2.5, 97.5]), pr, np.percentile(bp, [2.5, 97.5])


def null_spearman(vals_by_set, axis_by_set):
    """vals, axis: (n_sets, k) arrays; Spearman per row, NaN rows via loop."""
    n = vals_by_set.shape[0]
    out = np.full(n, np.nan)
    bad = np.isnan(vals_by_set).any(1) | np.isnan(axis_by_set).any(1)
    g = ~bad
    if g.any():
        rx = stats.rankdata(vals_by_set[g], axis=1); ry = stats.rankdata(axis_by_set[g], axis=1)
        rx -= rx.mean(1, keepdims=True); ry -= ry.mean(1, keepdims=True)
        den = np.sqrt((rx ** 2).sum(1) * (ry ** 2).sum(1))
        out[g] = np.where(den > 0, (rx * ry).sum(1) / np.where(den > 0, den, 1), np.nan)
    for i in np.where(bad)[0]:
        m = ~(np.isnan(vals_by_set[i]) | np.isnan(axis_by_set[i]))
        if m.sum() > 10:
            out[i] = stats.spearmanr(vals_by_set[i][m], axis_by_set[i][m])[0]
    return out


def q(v):
    v = v[~np.isnan(v)]
    return float(np.mean(v)), float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


def woolf(a, b, c_, d_):
    a, b, c_, d_ = [x + 0.5 for x in (a, b, c_, d_)]
    orr = a * d_ / (b * c_); se = math.sqrt(1 / a + 1 / b + 1 / c_ + 1 / d_)
    return orr, orr * math.exp(-1.96 * se), orr * math.exp(1.96 * se)


def axis_vec(ax):
    return ax.reindex(ALL).values


rows_gw, rows_dir, rows_ax, rows_ov, rows_valid, rows_part = [], [], [], [], [], []

# ---------------------------------------------------------------- 1. genome-wide correlations with the axes and with each other
for s in sigs:
    for an, ax in AXES.items():
        x = s.stat.reindex(ph.index.union(ag.index)).values
        n, rho, rci, pr, pci = spearman_boot(x, ax.reindex(ph.index.union(ag.index)).values, NBOOT_GW, rng)
        rows_gw.append(dict(signature=s.name, kind=s.kind, axis=an, n_genes=n, spearman=round(rho, 4), spearman_boot95_low=round(rci[0], 4), spearman_boot95_high=round(rci[1], 4),
                            pearson=round(pr, 4), pearson_boot95_low=round(pci[0], 4), pearson_boot95_high=round(pci[1], 4)))
        P(s.name, "|", an[:5], "n", n, "rho", round(rho, 3), rci.round(3), "r", round(pr, 3))

# ---------------------------------------------------------------- 2. validation against the paper's own correlations
paper = {"ALZ1": 0.40, "ALZ2": 0.27, "SCZ1": 0.20, "SCZ2": 0.17, "PD1": 0.64, "PD2": 0.39}
paper_pd = "paper PD1 0.64 / PD2 0.39 use LIV-PM signatures from PM2 (PD1) or PM1 (PD2) plus half of the LIV samples, not in the supplement; primary signature used here"
lv = l1
for nm, pv in paper.items():
    s = sigd[nm]
    for stat_label, xl, yl in (("logFC", lv.logFC, s.effect), ("signed z / t", lv.t, s.stat)):
        j = pd.concat([xl.rename("a"), yl.rename("b")], axis=1, join="inner").dropna()
        n, rho, rci, _, _ = spearman_boot(j.a.values, j.b.values, 500, rng)
        rows_valid.append(dict(disease_signature=nm, liv_pm_statistic=stat_label, n_genes=n, spearman_here=round(rho, 4), boot95_low=round(rci[0], 4), boot95_high=round(rci[1], 4),
                               paper_spearman=pv, difference=round(rho - pv, 4), note=paper_pd if nm.startswith("PD") else ""))
# Collado-Torres 'original' vs Liharska ST1 (identity check of mapping and sign)
j = pd.concat([sigd["LIV-PM Liharska (ST1)"].stat.rename("a"), sigd["LIV-PM Collado-Torres reproduction of Liharska"].stat.rename("b")], axis=1, join="inner").dropna()
rows_valid.append(dict(disease_signature="Collado-Torres 'original' t versus Liharska ST1 t", liv_pm_statistic="t", n_genes=len(j), spearman_here=round(stats.spearmanr(j.a, j.b)[0], 6), note="same results re-exported by Collado-Torres (sign and gene mapping check)"))
for nm in ("LIV-PM Collado-Torres, no qSVA", "LIV-PM Collado-Torres, with qSVA"):
    j = pd.concat([sigd["LIV-PM Liharska (ST1)"].stat.rename("a"), sigd[nm].stat.rename("b")], axis=1, join="inner").dropna()
    rows_valid.append(dict(disease_signature=f"{nm} t versus Liharska ST1 t", liv_pm_statistic="t", n_genes=len(j), spearman_here=round(stats.spearmanr(j.a, j.b)[0], 4)))
# SCZ1 vs SCZ2 etc reproduce
for a_, b_, pv in (("ALZ1", "ALZ2", 0.49), ("SCZ1", "SCZ2", 0.26), ("PD1", "PD2", 0.63)):
    j = pd.concat([sigd[a_].effect.rename("a"), sigd[b_].effect.rename("b")], axis=1, join="inner").dropna()
    rows_valid.append(dict(disease_signature=f"{a_} vs {b_} (logFC)", liv_pm_statistic="-", n_genes=len(j), spearman_here=round(stats.spearmanr(j.a, j.b)[0], 4), paper_spearman=pv, difference=round(stats.spearmanr(j.a, j.b)[0] - pv, 4)))
# DEG counts vs paper (LIV DEG / PM DEG)
for lab, ncode in (("PM DEG", 1), ("LIV DEG", -1)):
    nn = int(((l1["adj.P.Val"] < 0.05) & (np.sign(l1.t) == ncode)).sum())
    nlab = int((l1["DE status"] == lab).sum())
    rows_valid.append(dict(disease_signature=f"count of {lab} (adj.P<0.05 and sign) vs 'DE status' column", liv_pm_statistic="t", n_genes=nn, spearman_here=np.nan, note=f"genes in DE status column (after min-p per symbol) = {nlab}"))

# ---------------------------------------------------------------- 3. set analyses
def set_analyses(s):
    st, qq, meas = s.ALL_stat, s.ALL_q, s.meas
    for nm in ("383", "107"):
        gs = GSETS[nm]
        idx = gs.idx[(sgn[gs.idx] != 0) & meas[gs.idx]]
        if len(idx) < 10:
            continue
        DR = gs.DR
        # direction: concordance and Spearman with SZ sign and aligned stat
        v, sg = st[idx], sgn[idx]
        conc = int((np.sign(v) == sg).sum()); n = len(idx)
        Vn, Sn = st[DR], sgn[DR]
        Mn = ~np.isnan(Vn)
        null_frac = ((np.sign(Vn) == Sn) & Mn).sum(1) / np.maximum(Mn.sum(1), 1)
        rho = stats.spearmanr(v, sg)[0]
        null_rho = null_spearman(Vn, Sn)
        aligned = float(np.mean(v * sg)); null_al = np.nanmean(np.where(Mn, Vn * Sn, np.nan), 1)
        # bootstrap over set genes
        bf, br, ba = [], [], []
        for _ in range(NBOOT_SET):
            i = rng.integers(0, n, n)
            bf.append((np.sign(v[i]) == sg[i]).mean()); br.append(stats.spearmanr(v[i], sg[i])[0]); ba.append(np.mean(v[i] * sg[i]))
        nf, nf_lo, nf_hi = q(null_frac); nr, nr_lo, nr_hi = q(null_rho); na_, na_lo, na_hi = q(null_al)
        # restricted to DEGs of the signature
        deg = qq[idx] < 0.05
        concd = int((np.sign(v[deg]) == sg[deg]).sum()) if deg.any() else 0
        rows_dir.append(dict(signature=s.name, set=nm, n_genes_with_SZ_sign_measured=n, n_concordant=conc, frac_concordant=round(conc / n, 4),
                             frac_boot95_low=round(np.percentile(bf, 2.5), 4), frac_boot95_high=round(np.percentile(bf, 97.5), 4),
                             null_frac_mean=round(nf, 4), null_frac_q025=round(nf_lo, 4), null_frac_q975=round(nf_hi, 4),
                             frac_minus_null=round(conc / n - nf, 4), perm_p_frac_ge_obs=round((1 + (null_frac >= conc / n - 1e-12).sum()) / (1 + NSETS), 4),
                             spearman_stat_vs_SZ_sign=round(rho, 4), spearman_boot95_low=round(np.percentile(br, 2.5), 4), spearman_boot95_high=round(np.percentile(br, 97.5), 4),
                             null_spearman_mean=round(nr, 4), null_spearman_q025=round(nr_lo, 4), null_spearman_q975=round(nr_hi, 4),
                             mean_SZ_aligned_stat=round(aligned, 4), aligned_boot95_low=round(np.percentile(ba, 2.5), 4), aligned_boot95_high=round(np.percentile(ba, 97.5), 4),
                             null_aligned_mean=round(na_, 4), null_aligned_q025=round(na_lo, 4), null_aligned_q975=round(na_hi, 4),
                             n_signature_DEG_in_set=int(deg.sum()), n_DEG_concordant=concd))
        # axes within set (across set genes)
        for an, ax in AXES.items():
            a = axis_vec(ax)
            ii = idx[~np.isnan(a[idx])]
            if len(ii) < 10:
                continue
            r_ = stats.spearmanr(st[ii], a[ii])[0]
            Vg, Ag = st[DR], a[DR]
            nrr = null_spearman(Vg, Ag)
            bb = []
            for _ in range(NBOOT_SET):
                k = rng.integers(0, len(ii), len(ii)); bb.append(stats.spearmanr(st[ii][k], a[ii][k])[0])
            nm_, lo, hi = q(nrr)
            rows_ax.append(dict(signature=s.name, set=nm, axis=an, n_genes=len(ii), spearman=round(r_, 4), boot95_low=round(np.percentile(bb, 2.5), 4), boot95_high=round(np.percentile(bb, 97.5), 4),
                                null_mean=round(nm_, 4), null_q025=round(lo, 4), null_q975=round(hi, 4), perm_p_two_sided=round((1 + (np.abs(nrr[~np.isnan(nrr)] - nm_) >= abs(r_ - nm_) - 1e-12).sum()) / (1 + (~np.isnan(nrr)).sum()), 4)))
        # SZ sign vs axes within the set (the SZ sign's own alignment, for reference; once per signature is redundant, kept small)
    # overlap of PM-up / LIV-up DEGs with 383 / 107 up / down / all
    pm = ~np.isnan(qq)
    for nm in ("383", "107"):
        for part in ("", "_up", "_down"):
            gs = GSETS[nm + part]
            idx = gs.idx[meas[gs.idx]]
            n = len(idx)
            if n < 5:
                continue
            Dn = gs.DR
            Mn = meas[Dn]
            for dlab, dsign in (("higher in postmortem / cases (stat>0)", 1), ("higher in living / controls (stat<0)", -1)):
                f = lambda ix: (qq[ix] < 0.05) & (np.sign(st[ix]) == dsign) & meas[ix]
                obs = int(f(idx).sum())
                null_f = (f(Dn) & Mn).sum(1) / np.maximum(Mn.sum(1), 1)
                pi = pool_idx[meas[pool_idx]]
                bg_de = int(f(pi).sum()); bg_n = len(pi)
                orr, olo, ohi = woolf(obs, n - obs, bg_de, bg_n - bg_de)
                bb = []
                for _ in range(NBOOT_SET):
                    i = rng.integers(0, n, n); bb.append(f(idx[i]).mean())
                nfm = null_f.mean()
                rows_ov.append(dict(signature=s.name, set=nm + (part.replace("_", " ") if part else " all"), deg_direction=dlab, n_set_measured=n, n_set_DEG=obs,
                                    frac_set_DEG=round(obs / n, 4), frac_boot95_low=round(np.percentile(bb, 2.5), 4), frac_boot95_high=round(np.percentile(bb, 97.5), 4),
                                    matched_null_frac_mean=round(nfm, 4), matched_null_q025=round(np.percentile(null_f, 2.5), 4), matched_null_q975=round(np.percentile(null_f, 97.5), 4),
                                    fold_over_matched=round(obs / n / nfm, 3) if nfm > 0 else np.nan,
                                    fold_boot95_low=round(np.percentile(bb, 2.5) / nfm, 3) if nfm > 0 else np.nan, fold_boot95_high=round(np.percentile(bb, 97.5) / nfm, 3) if nfm > 0 else np.nan,
                                    table_set_DEG=obs, table_set_notDEG=n - obs, table_pool_DEG=bg_de, table_pool_notDEG=bg_n - bg_de,
                                    fisher_OR_vs_pool=round(orr, 3), OR_woolf95_low=round(olo, 3), OR_woolf95_high=round(ohi, 3)))


for s in sigs:
    set_analyses(s)
    P("done", s.name)

# ---------------------------------------------------------------- 4. partial rank correlations
def partial_rank(y, x, zs):
    """partial Spearman of y,x given zs (list), via OLS residuals of normal scores"""
    Y, Xv = nscore(y), nscore(x)
    Z = np.column_stack([np.ones(len(y))] + [nscore(z) for z in zs])
    ry = Y - Z @ np.linalg.lstsq(Z, Y, rcond=None)[0]; rx = Xv - Z @ np.linalg.lstsq(Z, Xv, rcond=None)[0]
    return float(np.corrcoef(ry, rx)[0, 1])


for liv in ("LIV-PM Liharska (ST1)", "LIV-PM Collado-Torres, no qSVA", "LIV-PM Collado-Torres, with qSVA"):
    for dis in ("SCZ1", "SCZ2", "ALZ1", "ALZ2", "PD1", "PD2"):
        j = pd.concat([sigd[liv].stat.rename("l"), sigd[dis].stat.rename("d"), ph.rename("ph"), ag.rename("ag")], axis=1, join="inner").dropna()
        raw = stats.spearmanr(j.l, j.d)[0]
        full = partial_rank(j.l.values, j.d.values, [j.ph.values, j.ag.values])
        bs = []
        for _ in range(NBOOT_GW // 2):
            i = rng.integers(0, len(j), len(j)); jj = j.iloc[i]
            bs.append((stats.spearmanr(jj.l, jj.d)[0], partial_rank(jj.l.values, jj.d.values, [jj.ph.values, jj.ag.values])))
        bs = np.array(bs)
        rows_part.append(dict(liv_pm_signature=liv, disease_signature=dis, n_genes=len(j), spearman_raw=round(raw, 4), raw_boot95_low=round(np.percentile(bs[:, 0], 2.5), 4), raw_boot95_high=round(np.percentile(bs[:, 0], 97.5), 4),
                              partial_given_pH_and_agonal=round(full, 4), partial_boot95_low=round(np.percentile(bs[:, 1], 2.5), 4), partial_boot95_high=round(np.percentile(bs[:, 1], 97.5), 4)))
    # within the 383: LIV-PM vs SZ sign controlling for axes
    idx = GSETS["383"].idx
    idx = idx[(sgn[idx] != 0) & sigd[liv].meas[idx] & ~np.isnan(axis_vec(ph)[idx]) & ~np.isnan(axis_vec(ag)[idx])]
    y, x = sigd[liv].ALL_stat[idx], sgn[idx]
    zs = [axis_vec(ph)[idx], axis_vec(ag)[idx]]
    rawr = stats.spearmanr(y, x)[0]; pr = partial_rank(y, x, zs)
    bs = []
    for _ in range(NBOOT_SET):
        i = rng.integers(0, len(idx), len(idx)); bs.append((stats.spearmanr(y[i], x[i])[0], partial_rank(y[i], x[i], [z[i] for z in zs])))
    bs = np.array(bs)
    rows_part.append(dict(liv_pm_signature=liv, disease_signature="383-gene SZ sign (within the 383)", n_genes=len(idx), spearman_raw=round(rawr, 4), raw_boot95_low=round(np.percentile(bs[:, 0], 2.5), 4), raw_boot95_high=round(np.percentile(bs[:, 0], 97.5), 4),
                          partial_given_pH_and_agonal=round(pr, 4), partial_boot95_low=round(np.percentile(bs[:, 1], 2.5), 4), partial_boot95_high=round(np.percentile(bs[:, 1], 97.5), 4)))
    # SZ sign vs axes inside the 383 (context)
P("partials done")

# ---------------------------------------------------------------- 5. cell-type view
sp = pd.read_csv(SPEC, sep="\t", comment="#")
sp["sym"] = sp.gene.map(lambda s: hg(str(s)))
sp = sp.dropna(subset=["sym"]).drop_duplicates("sym").set_index("sym")
CLASSES = ["Astro", "EN", "IN_other", "LAMP5", "Micro", "OPC", "Oligo", "PVALB", "SST", "VIP", "Vascular"]
dec_sp = pd.qcut(sp.log2_mean_cpm_tissue.rank(method="first"), 10, labels=False)
rows_ct = []
NNULL_CT, NTOP = 10000, 100
for s in [sigd[n] for n in ("LIV-PM Liharska (ST1)", "LIV-PM Collado-Torres, no qSVA", "LIV-PM Collado-Torres, with qSVA", "LIV-PM Liharska (ST1) | pH and agonal axes removed", "SCZ1", "SCZ2")]:
    st = s.stat.reindex(sp.index)
    ok = st.notna()
    # restrict to the bulk-measured ones
    for cl in CLASSES:
        sc = sp[f"spec_{cl}"].where(ok & (sp[f"cpm_{cl}"] >= 1))
        both = sc.notna()
        rho, p_ = stats.spearmanr(st[both], sc[both])
        top = sc.nlargest(NTOP).index
        tv = st[top].values
        # expression-decile-matched random sets (decile of mean tissue log2 CPM), drawn from the measured genes not in the top set
        pool_ct = [g for g in sp.index[ok] if g not in set(top)]
        pdec = dec_sp[pool_ct].values
        byd = {k: np.array(pool_ct)[pdec == k] for k in range(10)}
        need = collections.Counter(dec_sp[top].values)
        rg = np.random.default_rng(SEED + CLASSES.index(cl))
        stv = st.to_dict()
        nm_mean, nm_pm = [], []
        stv_arr = {k: np.array([stv[g] for g in byd[k]]) for k in byd}
        for _ in range(NNULL_CT):
            vals = np.concatenate([rg.choice(stv_arr[k], size=n, replace=len(stv_arr[k]) < n) for k, n in need.items()])
            nm_mean.append(vals.mean()); nm_pm.append((vals > 0).mean())
        bm = [tv[rg.integers(0, len(tv), len(tv))].mean() for _ in range(2000)]
        bp = [(tv[rg.integers(0, len(tv), len(tv))] > 0).mean() for _ in range(2000)]
        rows_ct.append(dict(signature=s.name, cell_class=cl, n_markers=len(top), spearman_stat_vs_specificity=round(rho, 4), n_genes_spearman=int(both.sum()),
                            markers_mean_stat=round(float(tv.mean()), 4), mean_boot95_low=round(float(np.percentile(bm, 2.5)), 4), mean_boot95_high=round(float(np.percentile(bm, 97.5)), 4),
                            matched_null_mean=round(float(np.mean(nm_mean)), 4), matched_null_q025=round(float(np.percentile(nm_mean, 2.5)), 4), matched_null_q975=round(float(np.percentile(nm_mean, 97.5)), 4),
                            effect_over_null=round(float(tv.mean() - np.mean(nm_mean)), 4),
                            frac_markers_stat_positive=round(float((tv > 0).mean()), 4), frac_boot95_low=round(float(np.percentile(bp, 2.5)), 4), frac_boot95_high=round(float(np.percentile(bp, 97.5)), 4),
                            matched_null_frac_mean=round(float(np.mean(nm_pm)), 4)))
    P("celltype done", s.name)

# ---------------------------------------------------------------- write
inputs = [f"inputs (sha256 prefix): Liharska ST1 {sha(ST1)} (extracted from MOESM8 sha256 prefix in the file header), ST5 {sha(ST5)} (MOESM12), Collado-Torres Table S4 {sha(CT4)} (github LieberInstitute/living_brain_reanalysis processed-data/SupplementaryTables/TableS4.csv, downloaded 2026-10-18)",
          f"recurrence_genes.tsv {sha(REC)}; recurrent_sets.tsv {sha(SETS)}; DE prior {sha(PRIOR)}; pH axis v2 controls {sha(PHV2)}; GTEx hardy_t_by_region {sha(GTEX)}; HBCC specificity {sha(SPEC)}",
          f"script: scripts/s1_main.py, run {TODAY}; seed {SEED}; matched sets {NSETS}; bootstrap over genes: {NBOOT_SET} (set analyses), {NBOOT_GW} (genome-wide); symbols: ENSG -> HGNC approved symbol (HGNC ensembl_gene_id column via extract.load_hgnc), minimum-P row per symbol; NEW analysis",
          f"gene mapping: {meta}",
          "intervals: percentile bootstrap over genes (genes are not independent: co-expression); matched null = 10,000 sets drawn from published-brain genes with an SZ sign outside the 383, same expression decile x DE-prior quintile strata as ad_383.py; null columns give mean and 2.5%-97.5% range; perm P descriptive only"]


def write(name, rows, extra=()):
    df = pd.DataFrame(rows)
    with open(os.path.join(OUT, name), "w") as f:
        for h in inputs + list(extra):
            f.write("# " + h + "\n")
        df.to_csv(f, sep="\t", index=False)
    P("wrote", name, len(df))


write("genomewide_axis_correlations.tsv", rows_gw, ["positive stat = higher in postmortem (LIV-PM) / higher in cases (disease); pH axis positive = higher at higher pH; agonal positive = higher after slow death; n = genes with both values"])
write("validation_vs_paper.tsv", rows_valid, ["paper values from Liharska et al. 2025 Results (preprint text: ALZ1/ALZ2/SCZ1/SCZ2 vs primary LIV-PM; PD1/PD2 vs LIV-PM signatures from PM2/PM1 plus half of LIV, not in the supplement); statistic used by the paper (logFC or t) is not stated in the text, both given"])
write("sets_direction.tsv", rows_dir, ["SZ sign = sign(up - down) over lists (recurrence_genes.tsv); concordant = sign of the signature statistic equals the SZ sign (LIV-PM: positive = higher in postmortem; so concordant means the gene is higher in postmortem where it is higher in SZ cases)"])
write("sets_axes.tsv", rows_ax, ["Spearman across set genes between the signature statistic and the axis"])
write("sets_overlap.tsv", rows_ov, ["DEG = signature adj.P < 0.05 (BH within the signature, as published for ST1/ST5; Collado-Torres adj.P.Val); direction by sign of the statistic; set 'up'/'down' = 383 or 107 genes whose SZ sign is up/down in cases; fold over matched = set fraction / mean fraction of matched sets (strata of that subset); table_* = 2x2 against all pool genes measured (OR with Woolf interval, descriptive)"])
write("partial_correlations.tsv", rows_part, ["partial = Spearman-type partial correlation through normal scores, given pH and agonal axes (OLS residuals); genome-wide rows use genes with all four values"])
write("celltype_view.tsv", rows_ct, ["markers = top 100 genes by PsychAD HBCC control specificity per class (class mean CPM >= 1) among genes with a statistic; matched null = random genes of the same decile of mean tissue log2 CPM (10,000 draws); stat > 0 = higher in postmortem (LIV-PM) / cases (disease)"])
json_meta = os.path.join(OUT, "run_meta.json")
import json
json.dump(dict(date=TODAY, seed=SEED, nsets=NSETS, meta=meta, python=sys.version.split()[0], numpy=np.__version__, pandas=pd.__version__), open(json_meta, "w"), indent=1)
