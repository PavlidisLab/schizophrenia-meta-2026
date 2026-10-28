#!/usr/bin/env python3
"""NEW 2026-10-18 (composition_sim). Step 3: summarise sim_reps_long.tsv.gz into the result TSVs (provenance header lines on each).
Run: python summarise.py <run_dir> <final_dir>. Intervals: 2.5th and 97.5th percentiles over reps; mean = back-transformed mean of Fisher z."""
import sys, os, json, numpy as np, pandas as pd
RUN, FIN = sys.argv[1], sys.argv[2]; OUT = os.path.join(RUN, "out"); os.makedirs(FIN, exist_ok=True)
SHAS = open(os.path.join(OUT, "input_sha256.txt")).read().split("\n")
def sh(name): return [l.split("\t")[0] for l in SHAS if l.endswith(name)][0]
HDR = ["# NEW analysis, composition-only simulation of the bulk pH and agonal axes; run 2026-10-18 on the remote environment; seed 20261018 (SeedSequence [seed, config_id, rep]); 200 reps per config",
       "# scripts: scripts/prep.py, scripts/sim.py, scripts/summarise.py; inputs sha256(16): axis_v2_controls.tsv.gz %s; hardy_t_by_region.tsv.gz %s; recurrent_sets.tsv %s; nsw_donor_covariates_used.tsv %s; PsychAD HBCC counts_class.npz %s; MSSM counts_class.npz %s; NSW pb_cols.tsv %s; donors_used.tsv %s" % (
        sh("axis_v2_controls.tsv.gz"), sh("hardy_t_by_region.tsv.gz"), sh("recurrent_sets.tsv"), sh("nsw_donor_covariates_used.tsv"), sh("HBCC_Cohort/counts_class.npz"), sh("MSSM_Cohort/counts_class.npz"), sh("pb_cols.tsv"), sh("donors_used.tsv")),
       "# sign: pH axis z_stouffer positive = higher expression at higher pH; simulated t = OLS t of log2CPM on pH (positive = higher at higher pH); agonal = GTEx Hardy slow-minus-fast-death t (positive = higher in slow death), mean of BA9, BA24, Cortex"]
def write(df, name, extra=()):
    with open(os.path.join(FIN, name), "w") as f:
        for l in HDR + list(extra): f.write(l + "\n")
        df.to_csv(f, sep="\t", index=False, float_format="%.5g")
L = pd.read_csv(os.path.join(OUT, "sim_reps_long.tsv.gz"), sep="\t"); C = pd.read_csv(os.path.join(OUT, "sim_configs.tsv"), sep="\t")
def fz_mean(v): v = np.clip(v, -0.999999, 0.999999); return float(np.tanh(np.nanmean(np.arctanh(v))))
def summ(g):
    v = g.value.values; isr = g.stat.iloc[0].startswith("rho_")
    return pd.Series(dict(mean=fz_mean(v) if isr else float(np.nanmean(v)), median=float(np.nanmedian(v)), lo=float(np.nanpercentile(v, 2.5)), hi=float(np.nanpercentile(v, 97.5)), sd=float(np.nanstd(v)), n_reps=int(np.isfinite(v).sum())))
S = L.groupby(["config_id", "stat"]).apply(summ).reset_index().merge(C, on="config_id")
cols = ["config_id", "tag", "profile", "variant", "loco", "N", "rna", "scale", "noise", "boot", "pc_beta"]
main = S[S.stat.str.match(r"^rho_(pH|agonal)_(all|rec383|rec107)$")].copy()
main[["axis", "gene_set"]] = main.stat.str.extract(r"^rho_(pH|agonal)_(\w+)$")
main = main.rename(columns={"mean": "rho_mean", "lo": "rho_lo", "hi": "rho_hi"})
n_g = L[L.stat.str.startswith("n_")].groupby(["config_id", "stat"]).value.first().reset_index()
n_g[["axis", "gene_set"]] = n_g.stat.str.extract(r"^n_(pH|agonal)_(\w+)$"); n_g = n_g.rename(columns={"value": "n_genes"})[["config_id", "axis", "gene_set", "n_genes"]]
main = main.merge(n_g, on=["config_id", "axis", "gene_set"])
write(main[cols + ["axis", "gene_set", "n_genes", "rho_mean", "median", "rho_lo", "rho_hi", "sd", "n_reps"]].sort_values(["tag", "profile", "variant", "loco", "N", "rna", "scale", "pc_beta", "axis", "gene_set"]), "sim_rho_summary.tsv",
      ["# rho = Spearman across genes between the simulated per-gene t and the real axis; rho_lo/rho_hi = 2.5%/97.5% over 200 reps (each rep: donor bootstrap of the 32 NSW controls, refit of slopes, new pH draws, new composition residuals, new counts); mean = Fisher-z mean"])
# deciles
dec = S[S.stat.str.contains("decile")].copy(); dec[["axis", "decile"]] = dec.stat.str.extract(r"^rho_(pH|agonal)_decile(\d+)$")
write(dec[cols + ["axis", "decile", "mean", "lo", "hi"]].sort_values(["tag", "profile", "rna", "axis", "decile"]), "sim_rho_by_expression_decile.tsv", ["# decile 1 = lowest mean simulated log2CPM; genes with both axis values; rho within decile, mean and 2.5/97.5% over reps"])
# expected noise-free
E = pd.read_csv(os.path.join(OUT, "sim_expected_noise_free.tsv"), sep="\t")
keepc = [c for c in E.columns if c.startswith("rho_") and "decile" not in c or c.startswith("n_")]
write(E[["profile", "rna", "scale", "variant"] + keepc], "sim_expected_noise_free.tsv", ["# noise-free expectation: slope of log2 mixture CPM on the 32 observed control pH values, fitted slopes, no residual noise, no sampling noise. scale 0 has zero slope; its rho is from a 1e-12 jitter and is not meaningful"])
Ed = E[["profile", "rna", "scale", "variant"] + [c for c in E.columns if "decile" in c]]
# ceiling
R = pd.read_csv("%s/input/reliability_summary.tsv" % RUN, sep="\t", comment="#"); R = R[R.axis == "v2 bulk (primary assignment)"].iloc[0]
ceil = pd.DataFrame([
    dict(basis="between-collection Spearman-Brown reliability (5 collections)", reliability=R.between_SB_point, rel_lo=R.between_SB_lo, rel_hi=R.between_SB_hi),
    dict(basis="split-half Spearman-Brown reliability (9 datasets)", reliability=R.halves_SB_mean, rel_lo=R.halves_SB_lo, rel_hi=R.halves_SB_hi)])
for c in ["reliability", "rel_lo", "rel_hi"]: ceil["sqrt_" + c] = np.sqrt(ceil[c])
write(ceil, "ceiling_reliability.tsv", ["# source: ph_axis_v2_2026-10-14/reliability_summary.tsv row 'v2 bulk (primary assignment)' (sha256(16) %s). Ceiling = sqrt(reliability) = the largest correlation a model that equalled the true axis could have with the observed axis (classical attenuation); applies to the pH axis over all genes; it is not available for the 383/107 subsets or the agonal axis" % open(os.path.join(OUT, "rel_sha.txt")).read().strip()])
# fraction of ceiling for headline configs
head = main[(main.tag == "main_grid") & (main.N == 200) & (main.scale == 1.0) & (main.axis == "pH") & (main.gene_set == "all")][["profile", "rna", "rho_mean", "rho_lo", "rho_hi"]].copy()
for _, r in ceil.iterrows():
    k = "between" if r.basis.startswith("between") else "halves"
    head["frac_of_ceiling_" + k] = head.rho_mean / r.sqrt_reliability; head["frac_lo_" + k] = head.rho_lo / r.sqrt_reliability; head["frac_hi_" + k] = head.rho_hi / r.sqrt_reliability
write(head, "sim_fraction_of_ceiling.tsv", ["# headline = N=200, scale 1, pH axis, all genes; fraction = rho / sqrt(reliability); interval endpoints divided by the same constant (reliability uncertainty not propagated)"])
# loco
M0 = main[(main.tag == "main_grid") & (main.N == 200) & (main.rna == 1) & (main.scale == 1.0)].set_index(["profile", "axis", "gene_set"])
lo = main[main.tag == "leave_one_class_out"].copy()
rows = []
for _, r in lo.iterrows():
    m = M0.loc[(r.profile, r.axis, r.gene_set)]
    d = m.rho_mean - r.rho_mean; se = np.sqrt(m.sd ** 2 + r.sd ** 2)
    rows.append(dict(profile=r.profile, left_out=r.loco, axis=r.axis, gene_set=r.gene_set, rho_all_classes=m.rho_mean, rho_without=r.rho_mean, rho_without_lo=r.rho_lo, rho_without_hi=r.rho_hi, drop=d, drop_lo=d - 1.96 * se, drop_hi=d + 1.96 * se))
write(pd.DataFrame(rows), "sim_leave_one_class_out.tsv", ["# leave-one-class-out: the left-out class keeps its model fraction at the mean pH (no pH dependence, no noise); the other classes keep slopes and noise. drop = rho(all classes) - rho(without); interval = 1.96 x sqrt(SD1^2 + SD2^2) of the two rep distributions (configs are not paired)"])
# positive control
pc = S[S.config_id.isin(C[C.tag == "positive_control"].config_id) & S.stat.str.startswith("pc_")]
pcw = pc.pivot_table(index=["profile", "pc_beta", "scale"], columns="stat", values=["mean", "lo", "hi"])
pcw.columns = ["%s_%s" % (b, a) for a, b in pcw.columns]; pcw = pcw.reset_index()
rr = main[(main.tag == "positive_control") & (main.axis == "pH") & (main.gene_set == "all")][["profile", "pc_beta", "scale", "rho_mean", "rho_lo", "rho_hi"]].rename(columns={"rho_mean": "rho_with_real_pH_axis", "rho_lo": "rho_axis_lo", "rho_hi": "rho_axis_hi"})
write(pcw.merge(rr, on=["profile", "pc_beta", "scale"]), "sim_positive_control.tsv", ["# POSITIVE CONTROL (labelled): cell-intrinsic pH effect injected in a random 5% of genes (log2 fold change per pH unit = +/- pc_beta, random sign, all classes), on top of the composition model (scale 1) or alone (scale 0). auc = probability a random injected gene has larger |t| than a random other gene; sign_agree = fraction of injected genes with t of the injected sign; rho_t_vs_beta = Spearman of t with the injected effect over all genes"])
# fraction slopes (generated data) vs observed
ff = S[S.stat.str.startswith("fracslope_") & (S.tag == "main_grid") & (S.N == 200) & (S.rna == 1) & (S.scale == 1.0)].copy(); ff["class"] = ff.stat.str.replace("fracslope_", "")
nsw = pd.read_csv(os.path.join(OUT, "nsw_control_composition.tsv"), sep="\t"); K7 = ["EN", "IN", "Astro", "OPC", "Oligo", "Immune", "Endo"]
X = np.column_stack([np.ones(len(nsw)), nsw.ph, nsw.age, nsw.male]); Fr = nsw[K7].values / nsw[K7].values.sum(1, keepdims=True)
Bf = np.linalg.lstsq(X, Fr, rcond=None)[0][1]
obs = pd.DataFrame(dict(**{"class": K7}, nsw_fraction_slope_per_pH_adj_age_sex=Bf))
sl = pd.read_csv(os.path.join(OUT, "nsw_class_slopes_point.tsv"), sep="\t").rename(columns={"class_nsw": "class"})
cm = {"EN": "Exc_Neurons", "IN": "In_Neurons", "Astro": "Astrocytes", "OPC": "OPC", "Oligo": "Oligodendrocyte", "Immune": "Microglia", "Endo": "Endothelial"}
sl["nsw_class_label"] = sl["class"].map(cm)
ff2 = ff.pivot_table(index=["class", "profile"], values=["mean", "lo", "hi"]).reset_index().rename(columns={"mean": "simulated_fraction_slope_mean", "lo": "simulated_lo", "hi": "simulated_hi"})
out = sl.drop(columns="class_psychad").merge(obs, on="class").merge(ff2, on="class")
write(out, "nsw_class_composition_slopes.tsv", ["# per-class composition slopes on pH in the 32 NSW control donors (CLR ~ pH + age + male; fraction slope = OLS of class fraction on pH + age + male, per pH unit); class mapping NSW label -> PsychAD class in nsw_class_label/class; PsychAD Mural has no NSW counterpart (fixed at the PsychAD control mean fraction, no pH dependence); simulated_* = slope of the generated class fraction on pH over reps (N=200, RNA-weighted, scale 1; Mural row not in NSW)"])
# gene counts & profile summary
json.dump(json.load(open(os.path.join(OUT, "sim_gene_counts.json"))), open(os.path.join(FIN, "sim_gene_counts.json"), "w"), indent=1)
pd.read_csv(os.path.join(OUT, "profile_summary.tsv"), sep="\t").to_csv(os.path.join(FIN, "profile_summary.tsv"), sep="\t", index=False)
print("summaries written")
