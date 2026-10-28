#!/usr/bin/env python3
"""Step 2: do the gene-level effects of the Yaqubi et al. 2026 processing conditions look like the bulk pH axis v2, the
agonal axis, or the recurrent schizophrenia genes (383 = set study_brain, 107 = study_brain_ge10 of base80_recurrence_2026-10-16)?

Adapted from ad_383.py (same population, same strata = expression decile x DE-prior quintile,
same pool = published-brain genes with a SZ sign outside the 383, same 10,000 matched sets drawn per gene set).
Signs: effect = processed/autopsy minus 0 h (CTRL); pH axis v2 (z_stouffer of axis_v2_controls) positive = higher at higher pH;
agonal axis = GTEx slow-minus-fast-death t (mean of BA9, BA24, Cortex) positive = higher in slow death. SZ sign = sign(up - down) of recurrence_genes.tsv.
Usage: python3 02_axes_and_sets.py <out dir> <dir with MOESM3_ESM.xlsx>
Seed 20261018.
"""
import collections, hashlib, os, random, sys, datetime, math
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

# ---------------------------------------------------------------- effects
E = pd.read_csv(EFF, sep="\t", comment="#", index_col=0)
tcols = [c for c in E.columns if c.startswith("t_")]
rcols = [c for c in E.columns if c.startswith("relLFC_")]
contrasts = [c[2:] for c in tcols]
EFFV = {}      # label -> pd.Series (all genes)
for c in tcols:
    EFFV["t:" + c[2:]] = E[c].dropna()
for c in rcols:
    EFFV["relLFC:" + c[7:]] = E[c].dropna()


def rank(x):
    return stats.rankdata(x)


def spearman_pearson_boot(x, y, rng):
    r_s, r_p = stats.spearmanr(x, y)[0], np.corrcoef(x, y)[0, 1]
    n = len(x)
    bs, bp = [], []
    for _ in range(NBOOT):
        i = rng.integers(0, n, n)
        bs.append(np.corrcoef(rank(x[i]), rank(y[i]))[0, 1])
        bp.append(np.corrcoef(x[i], y[i])[0, 1])
    return r_s, np.quantile(bs, [.025, .975]), r_p, np.quantile(bp, [.025, .975])


rng = np.random.default_rng(SEED)
# ---- (1) all genes
rows = []
for lab, s in EFFV.items():
    for an, ax in AX_ALL.items():
        j = pd.concat([s, ax], axis=1, join="inner").dropna()
        x, y = j.iloc[:, 0].values, j.iloc[:, 1].values
        rs, cs, rp, cp = spearman_pearson_boot(x, y, rng)
        rows.append(dict(effect=lab, axis=an, n_genes=len(j), spearman=round(rs, 3), spearman_boot95_low=round(cs[0], 3), spearman_boot95_high=round(cs[1], 3),
                         pearson=round(rp, 3), pearson_boot95_low=round(cp[0], 3), pearson_boot95_high=round(cp[1], 3)))
AXALL = pd.DataFrame(rows)

# ---- (2) within recurrent sets with matched null
def rho_null(gs, e, a, which):
    """Spearman or Pearson between e and a within the set and within each matched set."""
    ok = ~np.isnan(e) & ~np.isnan(a)
    s = gs.idx[ok[gs.idx]]
    f = (lambda x, y: stats.spearmanr(x, y)[0]) if which == "spearman" else (lambda x, y: np.corrcoef(x, y)[0, 1])
    obs = f(e[s], a[s])
    nl = []
    for i in range(NPERM_RHO):
        d = gs.DR[i]
        d = d[ok[d]]
        nl.append(f(e[d], a[d]) if len(d) > 10 else np.nan)
    nl = np.array(nl)
    nl = nl[~np.isnan(nl)]
    # bootstrap over the set's genes
    bs = []
    for _ in range(NBOOT):
        i = rng.integers(0, len(s), len(s))
        bs.append(f(e[s][i], a[s][i]))
    return len(s), obs, np.quantile(bs, [.025, .975]), nl


rows = []
for gs in GSETS:
    for lab, ser in EFFV.items():
        e = full(ser.to_dict())
        for an, a in AX_U.items():
            for which in ("spearman", "pearson"):
                n, obs, ci, nl = rho_null(gs, e, a, which)
                rows.append(dict(gene_set=gs.name, effect=lab, axis=an, statistic=which, n_genes=n, value=round(obs, 3), boot95_low=round(ci[0], 3), boot95_high=round(ci[1], 3),
                                 null_mean=round(nl.mean(), 3), null_sd=round(nl.std(), 3), null_q025=round(np.quantile(nl, .025), 3), null_q975=round(np.quantile(nl, .975), 3),
                                 diff_from_null_mean=round(obs - nl.mean(), 3), z_vs_null=round((obs - nl.mean()) / nl.std(), 2),
                                 perm_p_two_sided=round((1 + (np.abs(nl - nl.mean()) >= abs(obs - nl.mean())).sum()) / (1 + len(nl)), 4)))
AXSET = pd.DataFrame(rows)
print("axis correlations done")

# ---- (3) direction concordance of the SZ sign with each effect sign
def wilson(k, n, z=1.96):
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return c - h, c + h


rows = []
refs = dict(EFFV)
# reference "effects": the axes themselves in the orientation 'low pH' (-pH z) and 'slow death' (+t), signed so positive = SZ-like direction hypothesis
refs["ref: low-pH direction (minus pH axis z)"] = pd.Series({g: -v for g, v in ph.items()})
refs["ref: slow-death direction (GTEx t)"] = pd.Series(ag)
for gs in GSETS:
    for lab, ser in refs.items():
        e = full(ser.to_dict())
        s = np.array([i for i in gs.idx if sgn[i] != 0 and not np.isnan(e[i])])
        n = len(s)
        conc = int((np.sign(e[s]) == sgn[s]).sum())
        DR = gs.DR
        M = ~np.isnan(e[DR]) & (sgn[DR] != 0)
        C = (np.sign(e[DR]) == sgn[DR]) & M
        frac = C.sum(1) / np.maximum(M.sum(1), 1)
        lo, hi = wilson(conc, n)
        rows.append(dict(gene_set=gs.name, effect=lab, n_genes=n, concordant=conc, frac_concordant=round(conc / n, 3), wilson95_low=round(lo, 3), wilson95_high=round(hi, 3),
                         null_mean=round(frac.mean(), 3), null_sd=round(frac.std(), 3), null_q025=round(np.quantile(frac, .025), 3), null_q975=round(np.quantile(frac, .975), 3),
                         diff_from_null_mean=round(conc / n - frac.mean(), 3), perm_p_two_sided=round((1 + (np.abs(frac - frac.mean()) >= abs(conc / n - frac.mean())).sum()) / (1 + len(frac)), 4)))
CONC = pd.DataFrame(rows)
print("concordance done")

# ---------------------------------------------------------------- gene sets from the supplementary tables
X3 = os.path.join(SUPP, "moesm3.xlsx")
def hgl(series):
    return {hg(s) for s in series.dropna().astype(str) if hg(s)}

s5 = pd.read_excel(X3, sheet_name="Supplementary Data 5", header=1, engine="openpyxl")
blocks = {"S_36h20C_vs_0h": 0, "S_24h20C_vs_0h": 7, "S_24h4C_vs_0h": 14, "S_36h4C_vs_0h": 21, "S_6h20C_vs_0h": 28}
deg = {}
for nm, c0 in blocks.items():
    b = s5.iloc[:, c0:c0 + 6].dropna(how="all")
    b.columns = ["id", "sym", "baseMean", "lfc", "p", "padj"]
    b = b.dropna(subset=["sym"])
    deg[nm] = b
    print(nm, len(b), "up", int((b.lfc > 0).sum()), "down", int((b.lfc < 0).sum()), "min|lfc|", round(b.lfc.abs().min(), 2), "max padj", b.padj.max())
up = {nm: hgl(b[b.lfc > 0].sym) for nm, b in deg.items()}
dn = {nm: hgl(b[b.lfc < 0].sym) for nm, b in deg.items()}
art_up_all5 = set.intersection(*up.values())
print("up in all 5 DEG-yielding conditions:", len(art_up_all5), "(paper: 237 'artifact-responsive genes')")
s2 = pd.read_excel(X3, sheet_name="Supplementary Data 2", header=1, engine="openpyxl").iloc[:, :7].dropna(subset=["symbol"])
s3 = pd.read_excel(X3, sheet_name="Supplementary Data 3", header=1, engine="openpyxl").iloc[:, :6].dropna(subset=["symbol"])
s4 = pd.read_excel(X3, sheet_name="Supplementary Data 4", header=1, engine="openpyxl")
bag_up, bag_dn = hgl(s4.iloc[:, 0]), hgl(s4.iloc[:, 4])
SETS_X = {
    "artifact-responsive (up in all 5 DEG-yielding conditions, recomputed from Supp. Data 5)": art_up_all5,
    "processing DEG up in any condition": set.union(*up.values()), "processing DEG down in any condition": set.union(*dn.values()),
    "LPMI vs CTRL DEG up": hgl(s2[s2.log2FoldChange > 0].symbol), "LPMI vs CTRL DEG down": hgl(s2[s2.log2FoldChange < 0].symbol),
    "SPMI vs CTRL DEG up": hgl(s3[s3.log2FoldChange > 0].symbol), "SPMI vs CTRL DEG down": hgl(s3[s3.log2FoldChange < 0].symbol),
    "BAG up (shared SPMI and LPMI, Supp. Data 4)": bag_up, "BAG down (shared SPMI and LPMI, Supp. Data 4)": bag_dn,
}
for nm in blocks:
    SETS_X[f"{nm} DEG up"] = up[nm]
    SETS_X[f"{nm} DEG down"] = dn[nm]
# the paper's own overlap fractions as a check on the set definitions
allproc = set.union(*up.values()) | set.union(*dn.values())
chk = dict(spmi_deg=len(SETS_X["SPMI vs CTRL DEG up"] | SETS_X["SPMI vs CTRL DEG down"]), lpmi_deg=len(SETS_X["LPMI vs CTRL DEG up"] | SETS_X["LPMI vs CTRL DEG down"]))
chk["frac_spmi_deg_in_processing_union"] = round(len((SETS_X["SPMI vs CTRL DEG up"] | SETS_X["SPMI vs CTRL DEG down"]) & allproc) / chk["spmi_deg"], 3)
chk["frac_lpmi_deg_in_processing_union"] = round(len((SETS_X["LPMI vs CTRL DEG up"] | SETS_X["LPMI vs CTRL DEG down"]) & allproc) / chk["lpmi_deg"], 3)
print("check (paper: 15.8% SPMI, 11.3% LPMI):", chk)

U = set(ALL)
rows = []
for gs in GSETS:
    gmask_set = set(gs.genes)
    for nm, st in SETS_X.items():
        inU = np.zeros(NG, bool)
        for g in st:
            if g in GI:
                inU[GI[g]] = True
        obs = int(inU[gs.idx].sum())
        nl = inU[gs.DR].sum(1)
        ov = [g for g in gs.genes if inU[GI[g]]]
        sz_conc = None
        # direction of the overlap: for 'up' sets the processing effect is +, for 'down' it is -
        dirn = -1 if "down" in nm else (1 if (" up" in nm or "(up in" in nm) else 0)
        if dirn and ov:
            sg = [sign[g] for g in ov if sign[g] != 0]
            sz_conc = int(sum(1 for s_ in sg if s_ == dirn)), len(sg)
        pool_frac = float(inU[[GI[g] for g in pool]].mean())
        rows.append(dict(gene_set=gs.name, processing_set=nm, n_processing_set_in_universe=int(inU.sum()), n_set_genes=len(gs.genes), overlap=obs,
                         frac_of_set=round(obs / len(gs.genes), 3), pool_frac=round(pool_frac, 3), null_mean=round(float(nl.mean()), 2), null_q025=float(np.quantile(nl, .025)), null_q975=float(np.quantile(nl, .975)),
                         fold_vs_matched_null=round(obs / nl.mean(), 2) if nl.mean() > 0 else "", fold_null_q025=round(float(np.quantile(nl, .025)) / nl.mean(), 2) if nl.mean() > 0 else "",
                         fold_null_q975=round(float(np.quantile(nl, .975)) / nl.mean(), 2) if nl.mean() > 0 else "",
                         perm_p_upper=round((1 + (nl >= obs).sum()) / (1 + len(nl)), 4),
                         overlap_with_SZ_sign_and_direction_known=sz_conc[1] if sz_conc else "", overlap_SZ_direction_same_as_processing_direction=sz_conc[0] if sz_conc else ""))
OVL = pd.DataFrame(rows)
print("overlap done")

# ---------------------------------------------------------------- cell-type specificity view
sp = pd.read_csv(SPEC, sep="\t", comment="#")
sp["g"] = [hg(s) for s in sp.gene]
sp = sp.dropna(subset=["g"]).sort_values("log2_mean_cpm_tissue", ascending=False).drop_duplicates("g").set_index("g")
speccols = [c for c in sp.columns if c.startswith("spec_")]
rows = []
effs_ct = {k: v for k, v in EFFV.items()}
effs_ct["axis: pH axis v2"] = AX_ALL["pH axis v2"]
effs_ct["axis: agonal axis (GTEx slow-minus-fast)"] = AX_ALL["agonal axis (GTEx slow-minus-fast)"]
for lab, ser in effs_ct.items():
    for c in speccols:
        j = pd.concat([ser.rename("e"), sp[c].rename("s")], axis=1, join="inner").dropna()
        x, y = j.e.values, j.s.values
        r = stats.spearmanr(x, y)[0]
        bs = []
        for _ in range(300):
            i = rng.integers(0, len(x), len(x))
            bs.append(np.corrcoef(rank(x[i]), rank(y[i]))[0, 1])
        rows.append(dict(effect=lab, specificity=c[5:], n_genes=len(j), spearman=round(r, 3), boot95_low=round(np.quantile(bs, .025), 3), boot95_high=round(np.quantile(bs, .975), 3)))
CT = pd.DataFrame(rows)
print("cell-type done")

# ---------------------------------------------------------------- write
today = datetime.date.today().isoformat()
H = [f"# script: 02_axes_and_sets.py, run {today}; NEW analysis; seed {SEED}; {NPERM} matched sets (axis correlations use the first {NPERM_RHO}); {NBOOT} bootstrap resamples over genes (300 for cell-type rows)",
     f"# inputs: contrast_effects.tsv sha256 {sha(EFF)}; axis_v2_controls.tsv.gz {sha(PHAX)} (z_stouffer; positive = higher at higher pH); hardy_t_by_region.tsv.gz {sha(GTEX)} (mean of BA9, BA24, Cortex; positive = higher in slow death);",
     f"# recurrence_genes.tsv {sha(REC)}; recurrent_sets.tsv {sha(SETS)}; specificity_HBCC.tsv.gz {sha(SPEC)}; MOESM3_ESM.xlsx (Supplementary Data 1-7) {sha(X3)}",
     "# effect sign: processed or autopsy minus 0 h CTRL (positive = higher after processing); 't:' = own OLS t from released counts, 'relLFC:' = released whole-transcriptome log2FC (baseMean >= 10)",
     "# matched null: published-brain genes with a SZ sign outside the 383, strata = expression decile x DE-prior quintile (ad_383.py); genes are not independent, intervals over genes are descriptive"]


def w(df, name, extra=""):
    with open(f"{OUT}/tables/{name}", "w") as f:
        f.write("\n".join(H) + "\n" + (extra and "# " + extra + "\n"))
        df.to_csv(f, sep="\t", index=False)


w(AXALL, "axis_correlation_all_genes.tsv", "Spearman and Pearson over all genes measured in both the effect table and the axis; bootstrap over genes")
w(AXSET, "axis_correlation_recurrent_sets.tsv", "within the recurrent set; null = same statistic in matched gene sets")
w(CONC, "direction_concordance.tsv", "share of set genes whose SZ sign (up - down of recurrence_genes.tsv) equals the sign of the effect; null = same share in matched sets; ref rows give the same test for the axes in the direction a lower pH or slower death would predict")
w(OVL, "overlap_processing_sets.tsv", f"overlap of processing/autopsy DEG sets with the recurrent sets; null = overlap of the processing set with matched sets; set definition check: {chk}")
w(CT, "celltype_specificity_correlation.tsv", "Spearman over all genes between the effect and the HBCC control-class specificity (PsychAD)")
pd.Series(chk).to_csv(f"{OUT}/tables/set_definition_check.tsv", sep="\t", header=False)
print(AXALL.to_string())
print(CONC.to_string())
print(OVL.to_string())
