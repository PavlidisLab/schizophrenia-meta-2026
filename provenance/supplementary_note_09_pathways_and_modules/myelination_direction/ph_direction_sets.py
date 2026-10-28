#!/usr/bin/env python3
"""Direction of the control-pH association of the myelination and oxidative phosphorylation gene sets (ten dataset-regions, controls only; frink).

Per dataset-region (matrices of fig3b_ten/out/primary/matrices, the inputs of pathways_ten_2026-10-20/ph_pathways_ten.R): controls (sz == 0) with pH, age, sex; per gene OLS  expression ~ pH + age + male  (as ph_pathways_ten.R).
slope = change in the gene's expression, in SD of that gene, per SD of pH (positive = HIGHER at higher pH, i.e. lower at low pH; negative = higher at low pH). Set value per region = mean over the set's genes of the slope; interval = percentile
interval (genes constant in a resample skipped) of 500 bootstrap resamples of the controls (genes kept together, so gene-gene correlation is retained). Pooled = mean of the region values; interval = t interval over the 10 regions (9 df).
Sets: term_group_genes.tsv (M09 results_A, sha256 c8c03204dffb), groups myelin and oxphos, versions terms (all GO genes of the group's terms) and hits (genes on the lists); sz_direction and the per-list directions are read from the same file.
Usage: ph_direction_sets.py MATRICES TERM_GROUP_GENES AXIS_V2 OUTDIR
"""
import sys, os, gzip, io, hashlib, collections
import numpy as np, pandas as pd
from scipy import stats
MAT, TG, AX, OUT = sys.argv[1:5]
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
TEN = ["GSE208338", "GSE35977", "GSE202537_caudate", "GSE21138", "GSE17612", "GSE80655_AnCg", "GSE53987_striatum", "GSE224683", "GSE92538_GPL96", "GSE92538_GPL570"]
tg = pd.read_csv(TG, sep="\t", comment="#", dtype=str)
print(tg.columns.tolist())
tg = tg[tg.group.isin(["myelin", "oxphos"])]
sets = {f"{g}_{v}": sorted(set(tg[(tg.group == g) & (tg.version == v)].gene)) for g in ("myelin", "oxphos") for v in ("terms", "hits")}
rng = np.random.default_rng(20261023)
def slopes(Y, p, age, male):
    X = np.column_stack([np.ones(len(p)), (p - p.mean()) / p.std(ddof=1), (age - age.mean()) / age.std(ddof=1), male])
    Ys = (Y - Y.mean(0)) / Y.std(0, ddof=1)
    return np.linalg.lstsq(X, Ys, rcond=None)[0][1]
rows = []; gene_rows = []
for D in TEN:
    E = pd.read_csv(f"{MAT}/{D}.expr.tsv.gz", sep="\t", index_col=0); cv = pd.read_csv(f"{MAT}/{D}.cov.tsv", sep="\t").set_index("person").loc[E.columns]
    cv["male"] = (cv.sex == "M").astype(float)
    c = cv[(cv.sz == 0) & cv[["pH", "age", "male"]].notna().all(1)]
    Y = E[c.index].T; p, a, m = c.pH.values.astype(float), c.age.values.astype(float), c.male.values
    genes = list(Y.columns); Yv = Y.values.astype(float); gi = {g: i for i, g in enumerate(genes)}
    for S, gl in sets.items():
        ix = [gi[g] for g in gl if g in gi]
        if len(ix) < 3: continue
        obs = slopes(Yv[:, ix], p, a, m)
        bs = []
        for _ in range(500):
            j = rng.integers(0, len(p), len(p)); bs.append(np.nanmean(slopes(Yv[j][:, ix], p[j], a[j], m[j])))
        rows.append(dict(dataset=D, set=S, n_ctl=len(p), n_genes_measured=len(ix), mean_slope=np.nanmean(obs), lo95=np.nanpercentile(bs, 2.5), hi95=np.nanpercentile(bs, 97.5),
                         median_gene_slope=float(np.nanmedian(obs)), share_genes_negative=float((obs[np.isfinite(obs)] < 0).mean())))
        if S.endswith("terms"):
            for g, s in zip([genes[i] for i in ix], obs): gene_rows.append(dict(dataset=D, set=S, gene=g, slope=s))
R = pd.DataFrame(rows)
pool = []
for S, h in R.groupby("set"):
    k = len(h); mu = h.mean_slope.mean(); se = h.mean_slope.std(ddof=1) / np.sqrt(k); tc = stats.t.ppf(0.975, k - 1)
    pool.append(dict(set=S, n_regions=k, n_regions_negative=int((h.mean_slope < 0).sum()), n_ctl_total=int(h.n_ctl.sum()), pooled_mean_slope=mu, lo95=mu - tc * se, hi95=mu + tc * se,
                     median_of_region_medians=float(h.median_gene_slope.median())))
P = pd.DataFrame(pool)
# gene-level pooled slope (mean over regions where measured) and the axis file's own sign
G = pd.DataFrame(gene_rows).groupby(["set", "gene"]).slope.mean().reset_index()
ax = pd.read_csv(AX, sep="\t", comment="#", index_col=0)
# per-gene list directions
def tally(sub):
    up = dn = 0; lists = set()
    for gcell in sub["lists"].dropna():
        for tok in str(gcell).split(";"):
            if ":" not in tok: continue
            nm, d = tok.rsplit(":", 1); lists.add(nm)
            up += d.strip() == "up"; dn += d.strip() == "down"
    return up, dn, lists
lc = [c for c in tg.columns if c not in ("group", "version", "gene", "sz_direction")]
print("columns after gene:", lc)
dirs = []
for S in sets:
    g, v = S.split("_"); sub = tg[(tg.group == g) & (tg.version == v)].drop_duplicates("gene")
    up_g = int((sub.sz_direction == "up").sum()); dn_g = int((sub.sz_direction == "down").sum()); mx = int((~sub.sz_direction.isin(["up", "down"])).sum())
    listcol = [c for c in lc if sub[c].astype(str).str.contains(":(?:up|down)").any()]
    up = dn = 0; names = set()
    if listcol:
        u, d, n = tally(sub.rename(columns={listcol[0]: "lists"})); up, dn, names = u, d, n
    gs = G[G.set == (S if S.endswith("terms") else S)].set_index("gene").slope if S.endswith("terms") else None
    gsl = G[G.set == f"{g}_terms"].set_index("gene").slope
    inset = sub.gene[sub.gene.isin(gsl.index)]
    axz = ax.reindex(sub.gene)
    dirs.append(dict(set=S, n_genes=len(sub), genes_reported_up=up_g, genes_reported_down=dn_g, genes_mixed=mx, gene_list_reports_up=up, gene_list_reports_down=dn, n_lists=len(names),
                     genes_in_axis_file=int(axz.z_stouffer.notna().sum()), mean_axis_z_stouffer=float(axz.z_stouffer.mean()), mean_axis_beta_re_sdunits=float(axz.beta_re_sdunits.mean()),
                     share_axis_z_negative=float((axz.z_stouffer.dropna() < 0).mean()), mean_pooled_slope=float(gsl.reindex(inset).mean())))
Dd = pd.DataFrame(dirs)
inputs = f"matrices expr sha256 prefixes: " + " ".join(f"{d}={sha(f'{MAT}/{d}.expr.tsv.gz')}" for d in TEN)
hdr = f"# ph_direction_sets.py (myelination_direction); term_group_genes.tsv sha256 {sha(TG)}; axis_v2_controls.tsv.gz {sha(AX)}; {inputs}\n# slope = SD of gene expression per SD of pH, controls, adjusted for age and sex; POSITIVE = higher at higher pH (lower at low pH), NEGATIVE = higher at low pH\n"
os.makedirs(OUT, exist_ok=True)
for nm, df in (("by_region.tsv", R), ("pooled.tsv", P), ("list_direction_and_axis.tsv", Dd)):
    with open(f"{OUT}/{nm}", "w") as fh: fh.write(hdr); df.round(4).to_csv(fh, sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(R.round(3).to_string()); print(P.round(3).to_string()); print(Dd.round(3).to_string())
