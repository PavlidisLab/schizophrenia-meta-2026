"""Summaries of 02_sample_level_blood.R output, GTEx and pH-axis lookups for hemoglobin / blood-score genes."""
import os, sys, gzip, math, csv
import numpy as np, pandas as pd
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from common import *
A = os.path.join(REPO, "analyses")
def rd(fn): return pd.read_csv(os.path.join(HERE, fn), sep="\t", comment="#")
def put(df, fn, extra, inputs):
    with open(os.path.join(HERE, fn), "w") as f:
        f.write(hdr("03_summarize.py", inputs)); f.write(extra.rstrip("\n") + "\n"); df.to_csv(f, sep="\t", index=False)
def dl(est, se):
    est, se = np.asarray(est, float), np.asarray(se, float); w = 1/se**2; fe = (w*est).sum()/w.sum()
    Q = (w*(est-fe)**2).sum(); k = len(est); tau2 = max(0, (Q-(k-1))/(w.sum()-(w**2).sum()/w.sum())) if k > 1 else 0
    ws = 1/(se**2+tau2); m = (ws*est).sum()/ws.sum(); s = math.sqrt(1/ws.sum())
    I2 = max(0, (Q-(k-1))/Q) if Q > 0 else 0
    return m, s, m-1.96*s, m+1.96*s, 2*(1-stats.norm.cdf(abs(m/s))), tau2, I2
SL = ["samplelevel_blood_score_case_control.tsv","samplelevel_assoc_hb_and_scores.tsv","samplelevel_hb_sz_effect_by_adjustment.tsv","samplelevel_209_signed_effect_by_adjustment.tsv"]
slp = [os.path.join("analyses/hemoglobin_gene_analyses", s) for s in SL]
# (a) case-control blood scores
cc = rd(SL[0]); cc = cc[cc.adj == "age+sex"]
rows = []
for sc in cc.score.unique():
    for label, sub in [("bulk, all datasets", cc[(cc.score == sc) & (cc.compartment == "bulk")]),
                       ("bulk, one dataset per cohort cluster (largest n)", None),
                       ("bulk, cluster-level (fixed effect within cluster, random effects across clusters)", "CL"),
                       ("non-bulk (LCM / sorted nuclei)", cc[(cc.score == sc) & (cc.compartment != "bulk")])]:
        if isinstance(sub, str):
            b = cc[(cc.score == sc) & (cc.compartment == "bulk")]; rr = []
            for c, g in b.groupby("cluster"):
                w = 1/g.se**2; rr.append((c, (w*g.coef).sum()/w.sum(), (1/w.sum())**0.5, g.p.min()))
            sub = pd.DataFrame(rr, columns=["cluster","coef","se","p"])
        elif sub is None:
            b = cc[(cc.score == sc) & (cc.compartment == "bulk")].copy(); b["n"] = b.n_sz + b.n_ctl
            b["cl"] = b.cluster.where(~b.cluster.isin(["pritzker"]), "pritzker"); sub = b.sort_values("n", ascending=False).drop_duplicates("cl")
        if len(sub) < 2: continue
        m, s, lo, hi, p, tau2, I2 = dl(sub.coef, sub.se)
        rows.append([sc, label, len(sub), m, lo, hi, p, I2, int((sub.coef > 0).sum()), int((sub.p < 0.05).sum())])
meta = pd.DataFrame(rows, columns=["score","set","n_datasets","random_effects_SD_units","ci_lo","ci_hi","p","I2","n_positive","n_p_lt_0.05"])
put(meta, "blood_score_case_control_meta.tsv", "# SZ minus control in SD of the control score, DerSimonian-Laird random effects over datasets; datasets of one cohort cluster share donors, so the 'one per cluster' row is the independent-donor version", slp)
print(meta.round(3).to_string())
# (b) associations
asc = rd(SL[1]); asc = asc[asc.compartment == "bulk"]
rows = []
for it in ["HBB","HBA1","HBA2","HBD","HBG2","ALAS2","SLC4A1","CA1","ery_mean_z","nonHb_ery_PC1","leuk_mean_z"]:
    for tr in ["pH","PMI","RIN","agonal","age","male"]:
        for scope in ["controls", "all"]:
            d = asc[(asc["item"] == it) & (asc.trait == tr) & (asc.scope == scope)]
            if len(d) == 0: continue
            rows.append([it, tr, scope, len(d), d.t.mean(), (d.t/np.sqrt(len(d))).sum() if False else d.t.sum()/math.sqrt(len(d)), int((d.t > 0).sum()), int((d.t < 0).sum()), int((d.p < 0.05).sum()), int(((d.p < 0.05) & (d.t > 0)).sum()), int(((d.p < 0.05) & (d.t < 0)).sum())])
asum = pd.DataFrame(rows, columns=["item","trait","scope","n_datasets","mean_t","stouffer_Z","n_t_pos","n_t_neg","n_p_lt_0.05","n_sig_pos","n_sig_neg"])
put(asum, "samplelevel_assoc_summary.tsv", "# bulk datasets only; t of lm(y ~ trait + age + sex [+ dx]); Stouffer Z is unweighted and treats datasets as independent (several share donors). positive t = higher expression with higher pH/PMI/RIN/age, in males", slp)
print(asum[asum.scope == "controls"].round(2).to_string())
# (c) hb SZ effect by adjustment
ad = rd(SL[2]); ad = ad[ad.compartment == "bulk"]; ad["se"] = (ad.logfc/ad.t).abs()
rows = []
for g in ["HBB","HBA1","HBA2","HBD","HBG2"]:
    for mdl in ad.model.unique():
        d = ad[(ad.gene == g) & (ad.model == mdl)]
        if len(d) < 2: continue
        m, s, lo, hi, p, tau2, I2 = dl(d.logfc, d.se)
        rows.append([g, mdl, len(d), d.logfc.median(), m, lo, hi, p, I2, int((d.logfc > 0).sum()), int((d.p < 0.05).sum()), int(((d.p < 0.05) & (d.logfc > 0)).sum())])
adm = pd.DataFrame(rows, columns=["gene","model","n_datasets","median_log2FC","re_log2FC","ci_lo","ci_hi","p","I2","n_log2FC_pos","n_p_lt_0.05","n_sig_pos"])
put(adm, "hb_sz_effect_meta_by_adjustment.tsv", "# SZ vs control log2FC in bulk datasets, random-effects; models: base = age+sex; nonHb_PC1 / ery_PC1 = + blood-gene PC1; quality_base / quality_pHPMI / quality_pHPMI_nonHbPC1 on the complete-case samples of the 18 datasets with pH and PMI", slp)
print(adm.round(3).to_string())
# (d) 209
e = rd(SL[3]); e = e[e.compartment == "bulk"]
base = e[e.model == "base"].set_index("dataset"); rows = []
for ds, b in base.iterrows():
    r = {"dataset": ds, "n_genes": b.n_genes, "base_signed_sum": b.sum_signed_logfc}
    sm = lambda m: (lambda x: x.sum_signed_logfc.iloc[0] if len(x) else np.nan)(e[(e.dataset == ds) & (e.model == m)])
    r["ratio_nonHb_PC1"] = sm("nonHb_PC1")/b.sum_signed_logfc; r["ratio_ery_PC1"] = sm("ery_PC1")/b.sum_signed_logfc
    r["ratio_quality_base"] = sm("quality_base")/b.sum_signed_logfc
    r["ratio_quality_pHPMI"] = sm("quality_pHPMI")/sm("quality_base"); r["ratio_quality_pHPMI_nonHbPC1"] = sm("quality_pHPMI_nonHbPC1")/sm("quality_base")
    n = e[(e.dataset == ds) & (e.model == "null_randomPC1_ratio_quantiles")]
    r["null_median"], r["null_lo2.5"], r["null_hi97.5"] = n.mean_signed_t.iloc[0], n.n_same_dir_p05.iloc[0], n.n_opp_dir_p05.iloc[0]
    r["obs_below_null_2.5"] = bool(r.get("ratio_nonHb_PC1", np.nan) < r["null_lo2.5"])
    rows.append(r)
r209 = pd.DataFrame(rows)
put(r209, "recurrent209_signed_effect_ratios.tsv", "# ratio = signed sum of log2FC (209 study_brain genes, sign = SZ direction) after adjustment / before; ratio_quality_* columns are relative to quality_base (same complete-case samples), ratio_quality_base to base; null = 200 random 9-gene PC1 adjustments (2.5 / 50 / 97.5 percentiles)", slp)
print(r209.round(2).to_string())
print("median ratio nonHb_PC1:", r209.ratio_nonHb_PC1.median(), " ery_PC1:", r209.ratio_ery_PC1.median(), " n obs<null2.5:", r209["obs_below_null_2.5"].sum(), "of", len(r209))
# (e) GTEx brain: death-speed and ischemic-time t for hemoglobin and blood genes
def rdgz(p):
    rows = [l.rstrip("\n").split("\t") for l in gzip.open(os.path.join(A, p), "rt") if not l.startswith("#")]
    return pd.DataFrame(rows[1:], columns=["gene"]+rows[0][1:]).set_index("gene").apply(pd.to_numeric, errors="coerce")
H = rdgz("gtex_ischemic_time/hardy_t_by_region.tsv.gz"); I = rdgz("gtex_ischemic_time/ischemic_t_by_region.tsv.gz")
sl = rdgz("gtex_death_speed_all_tissues/slow_vs_fast_t_by_tissue.tsv.gz"); vt = rdgz("gtex_death_speed_all_tissues/ventilator_vs_fast_t_by_tissue.tsv.gz")
genes = [g for g in dict.fromkeys(HB + ERY + LEUK)]
rows = []
for reg in H.columns:
    for g in genes:
        t_h = H.at[g, reg] if g in H.index else np.nan; t_i = I.at[g, reg] if g in I.index else np.nan
        col = H[reg].dropna(); rk = (col.rank(ascending=False)[g] if g in col.index else np.nan)
        rows.append([reg, g, t_h, (1 - (rk-1)/len(col)) if rk == rk else np.nan, len(col), t_i, sl.at[g, reg] if (g in sl.index and reg in sl.columns) else np.nan, vt.at[g, reg] if (g in vt.index and reg in vt.columns) else np.nan])
gt = pd.DataFrame(rows, columns=["region","gene","t_slow_vs_fast_death","percentile_of_genes_tested","n_genes_tested","t_ischemic_time","t_slow_vs_fast_all_tissue_file","t_ventilator_vs_fast"])
put(gt, "gtex_brain_hardy_ischemic_hemoglobin.tsv", "# GTEx v8, per brain region, t of OLS (log2 TPM+1 ~ group + ischemic time + age + sex); t > 0: higher in slow death (Hardy 4) than fast natural death (Hardy 2); genes with region median TPM < 1 absent; percentile 1 = most positive t among tested genes", ["analyses/gtex_ischemic_time/hardy_t_by_region.tsv.gz","analyses/gtex_ischemic_time/ischemic_t_by_region.tsv.gz","analyses/gtex_death_speed_all_tissues/slow_vs_fast_t_by_tissue.tsv.gz","analyses/gtex_death_speed_all_tissues/ventilator_vs_fast_t_by_tissue.tsv.gz"])
x = gt[gt.gene.isin(["HBB","HBA1","HBA2","ALAS2","SLC4A1","SLC25A37"])]
print(x.groupby("gene").agg(n_reg=("t_slow_vs_fast_death", "count"), mean_t=("t_slow_vs_fast_death", "mean"), min_t=("t_slow_vs_fast_death", "min"), max_t=("t_slow_vs_fast_death", "max"), pct=("percentile_of_genes_tested", "mean"), isch=("t_ischemic_time", "mean"), vent=("t_ventilator_vs_fast", "mean")).round(2))
# (f) pH Z
ph = pd.read_csv(os.path.join(A, "confounds/results/combined_controls_pH.tsv"), sep="\t", comment="#")
print(ph[ph.gene.isin(genes)].pivot(index="gene", columns="factor", values="Z").round(2).dropna(how="all"))
tot = ph[ph.factor == "pH"]; print("n genes pH", len(tot))
sel = ph[ph.gene.isin(genes)].copy()
sel["rank_pH_Z_low_first"] = sel.apply(lambda r: (ph[(ph.factor == r.factor)].Z < r.Z).sum()+1, axis=1); sel["n_genes"] = sel.factor.map(ph.factor.value_counts())
put(sel[["gene","factor","n_datasets","Z","datasets","expr_percentile","p","q","rank_pH_Z_low_first","n_genes"]], "pH_pmi_axis_hemoglobin.tsv", "# combined control-only Z (Stouffer over datasets); Z < 0: higher expression at lower pH / shorter PMI", ["analyses/ph_and_confound_analyses/results/combined_controls_pH.tsv"])
print(sel[sel.factor == "pH"][["gene","Z","rank_pH_Z_low_first","n_genes","q"]].round(3).to_string())
