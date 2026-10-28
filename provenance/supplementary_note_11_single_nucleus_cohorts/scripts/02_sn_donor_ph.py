#!/usr/bin/env python3
"""Single-nucleus donors with measured pH: does the donors' own pH explain the SZ effect (A) and detection vs loss (B).
# script: analyses/single_nucleus_donor_ph/scripts/02_sn_donor_ph.py   (run on remote, nice -n 10, one thread, ~minutes)
# data: snap_sz/pb/{cmc,gse254569}/{counts_class.npz,cols_class.tsv,genes.tsv,cache/*.meta.npz} (class-level pseudobulk counts and
#   per-nucleus class codes, built 2026-09-27 by the snap_sz pipeline); donors tables from 01_clean_ph.py (ph_clean);
#   Z_pH = analyses/ph_and_confound_analyses/results/combined_controls_pH.tsv (factor pH); gene sets = analyses/ph_and_confound_analyses/recurrent_sets.tsv.
# usage: 02_sn_donor_ph.py <pb_root> <donor_dir> <axis_dir> <out_dir> [B]
"""
import sys, os, json, hashlib, itertools
import numpy as np, pandas as pd
from scipy import sparse
from scipy.stats import spearmanr

PB, DON, AX, OUT = sys.argv[1:5]
B = int(sys.argv[5]) if len(sys.argv) > 5 else 1000
SEED = 20261012
NPERM = 300           # within-diagnosis pH permutations for the null of the change in r
MIN_NUC = 30          # pseudobulk needs >= 30 nuclei
DETECT = 0.8          # gene kept in a class when counts > 0 in >= 80% of that class's pseudobulks
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(SEED)

# ---------- axis and sets
z = pd.read_csv(f"{AX}/results/combined_controls_pH.tsv", sep="\t", comment="#")
zph = z[z.factor == "pH"].drop_duplicates("gene").set_index("gene").Z.dropna()
rs = pd.read_csv(f"{AX}/recurrent_sets.tsv", sep="\t", comment="#")
s209 = rs[rs.set == "study_brain"].drop_duplicates("gene").set_index("gene")
s5p = rs[rs.set == "study_brain_5plus"].drop_duplicates("gene")
sign209 = s209.sz_direction.map({"up": 1, "down": -1})  # mixed -> NaN, excluded from the signed mean
SETS = {"study_brain_209": list(s209.index), "study_brain_5plus": list(s5p.gene), "all_Z_genes": list(zph.index)}

CLS_LABEL = {"glutamatergic": "excitatory", "gabaergic": "interneuron", "astrocyte": "astrocyte", "oligodendrocyte": "oligodendrocyte",
             "polydendrocyte": "OPC", "microglia": "microglia", "endothelia": "vascular"}


def ols_t(Y, X, j):
    """Y: n x g, X: n x p, return coef and t of column j for each gene."""
    XtXi = np.linalg.inv(X.T @ X)
    beta = XtXi @ X.T @ Y
    res = Y - X @ beta
    s2 = (res ** 2).sum(0) / (X.shape[0] - X.shape[1])
    se = np.sqrt(s2 * XtXi[j, j])
    return beta[j], beta[j] / np.where(se > 0, se, np.nan)


def design(d, cols):
    return np.column_stack([np.ones(len(d))] + [d[c].to_numpy(float) for c in cols])


def pct(a):
    a = np.asarray(a, float)
    return np.nanpercentile(a, 2.5), np.nanpercentile(a, 97.5)


def corr(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    return np.corrcoef(a[m], b[m])[0, 1]


rowsA, rowsP, rowsF, rowsE, rowsD, donor_tabs, meta = [], [], [], [], [], [], {}
for cohort, key in [("CMC", "cmc"), ("GSE254569", "gse254569")]:
    don = pd.read_csv(f"{DON}/donors_{key}_ph_clean.tsv", sep="\t", dtype={"donor_id": str, "bioassay_id": str})
    cc = pd.read_csv(f"{PB}/{key}/cols_class.tsv", sep="\t", dtype={"donor_id": str, "bioassay_id": str})
    genes = pd.read_csv(f"{PB}/{key}/genes.tsv", sep="\t", dtype=str).gene_name.to_numpy()
    zz = np.load(f"{PB}/{key}/counts_class.npz")   # stored as CSC, genes x pseudobulks
    M = sparse.csc_matrix((zz["data"], zz["indices"], zz["indptr"]), shape=tuple(zz["shape"])).T.tocsr()
    assert M.shape[0] == len(cc), (M.shape, len(cc))
    # nucleus counts per donor x class code from the per-nucleus cache
    cnt = {}
    for _, r in cc.drop_duplicates("donor_id").iterrows():
        f = f"{PB}/{key}/cache/{r.bioassay_id}.meta.npz"
        cls = np.load(f)["cls"]
        cnt[r.donor_id] = np.bincount(cls[cls >= 0], minlength=7)
    code2name = {}
    for k in range(7):
        cand = None
        for name, g in cc.groupby("cell_type"):
            ok = all(cnt[r.donor_id][k] == r.n_nuclei for r in g.itertuples())
            if ok:
                cand = name if cand is None else "AMBIG"
        code2name[k] = cand
    assert None not in code2name.values() and "AMBIG" not in code2name.values(), code2name
    meta[cohort] = {"class_code_to_name": code2name, "genes": int(M.shape[1]), "pseudobulks": int(M.shape[0])}
    # donor-level composition
    frac = []
    for d, c in cnt.items():
        row = {"donor_id": d, "n_nuclei_total": int(c.sum())}
        for k, name in code2name.items():
            row["n_" + CLS_LABEL[name]] = int(c[k])
        frac.append(row)
    fr = pd.DataFrame(frac)
    for lab in CLS_LABEL.values():
        fr["frac_" + lab] = fr["n_" + lab] / fr.n_nuclei_total
    fr["frac_interneuron_of_neurons"] = fr.n_interneuron / (fr.n_interneuron + fr.n_excitatory)
    don = don.merge(fr, on="donor_id", how="left")
    # subtype shares among interneurons from labelled nuclei (GSE254569 only)
    if cohort == "GSE254569":
        o = pd.read_csv("Gemma/hdf5s/GSE254569_obs.tsv", sep="\t", usecols=["Donor", "celltypes_final"])
        o = o[o.celltypes_final.str.startswith("In_")]
        t = pd.crosstab(o.Donor, o.celltypes_final)
        sh = pd.DataFrame({"donor_id": t.index,
                           "share_SST_of_interneurons": (t["In_SST"] / t.sum(1)).values,
                           "share_PVALB_of_interneurons": ((t["In_PVALB_Ba"] + t["In_PVALB_Ch"]) / t.sum(1)).values})
        don = don.merge(sh, on="donor_id", how="left")
    donor_tabs.append(don.assign(cohort=cohort))
    dn = don.set_index("donor_id")

    lib = np.asarray(M.sum(1)).ravel()
    expr_all = np.log2(M.multiply(1e6 / lib[:, None]).tocsr().toarray() + 1.0)
    detect_all = (M > 0).toarray()
    # one column per symbol: the most expressed
    gser = pd.Series(np.arange(len(genes)), index=genes)
    mean_e = expr_all.mean(0)
    keep_idx = pd.Series(mean_e, index=range(len(genes))).groupby(genes).idxmax()
    sym2col = keep_idx.to_dict()

    # ---------------- B: composition vs pH
    def fit_slopes(d, ycol, covs, jname="ph_clean"):
        X = design(d, [jname] + covs)
        y = d[ycol].to_numpy(float)[:, None]
        b, t = ols_t(y, X, 1)
        return b[0], t[0]

    def boot(d, fn, nb=B):
        out = []
        n = len(d)
        for _ in range(nb):
            ix = rng.integers(0, n, n)
            try:
                out.append(fn(d.iloc[ix]))
            except np.linalg.LinAlgError:
                pass
        return np.array(out)

    valid = don[don.ph_clean.notna() & don.dx.isin(["control", "SZ"])].copy()
    valid["sz"] = (valid.dx == "SZ").astype(int)
    for grp, sub in [("controls", valid[valid.dx == "control"]), ("cases", valid[valid.dx == "SZ"]), ("all_dx_adjusted", valid)]:
        covs = ["age", "male"] + (["sz"] if grp == "all_dx_adjusted" else [])
        for ycol in ["frac_interneuron", "frac_interneuron_of_neurons", "frac_excitatory", "frac_astrocyte", "frac_oligodendrocyte",
                     "frac_OPC", "frac_microglia", "frac_vascular"]:
            s = sub.dropna(subset=[ycol, "age"])
            if len(s) < 8:
                continue
            est = fit_slopes(s, ycol, covs)[0]
            bs = boot(s, lambda d: fit_slopes(d, ycol, covs)[0])
            lo, hi = pct(bs)
            sdph = s.ph_clean.std()
            sdy = s[ycol].std()
            # unadjusted correlation as well
            rowsF.append(dict(cohort=cohort, group=grp, outcome=ycol, n=len(s), ph_min=s.ph_clean.min(), ph_max=s.ph_clean.max(),
                              ph_sd=sdph, slope_per_0p1_pH=est * 0.1, lo=lo * 0.1, hi=hi * 0.1,
                              slope_per_SD_pH_in_outcome_SD=est * sdph / sdy, lo_sd=lo * sdph / sdy, hi_sd=hi * sdph / sdy,
                              mean_outcome=s[ycol].mean(), spearman_unadjusted=spearmanr(s.ph_clean, s[ycol])[0]))
    # group difference in interneuron fraction before / after pH
    for ycol in ["frac_interneuron", "frac_interneuron_of_neurons"]:
        s = valid.dropna(subset=[ycol, "age"])
        def f2(d):
            b0 = ols_t(d[ycol].to_numpy(float)[:, None], design(d, ["sz", "age", "male"]), 1)[0][0]
            b1 = ols_t(d[ycol].to_numpy(float)[:, None], design(d, ["sz", "ph_clean", "age", "male"]), 1)[0][0]
            return b0, b1, b1 - b0
        est = f2(s); bs = boot(s, f2)
        for i, nm in enumerate(["sz_minus_control_age_sex", "sz_minus_control_age_sex_pH", "change_after_pH"]):
            lo, hi = pct(bs[:, i])
            rowsD.append(dict(cohort=cohort, outcome=ycol, quantity=nm, n=len(s), n_sz=int(s.sz.sum()), n_ctrl=int((1 - s.sz).sum()),
                              estimate=est[i], lo=lo, hi=hi, mean_control=s[s.sz == 0][ycol].mean()))

    # ---------------- per class pseudobulk analyses
    for cname, clab in CLS_LABEL.items():
        sel = (cc.cell_type == cname) & (cc.n_nuclei >= MIN_NUC)
        idx = np.where(sel.to_numpy())[0]
        d = cc.iloc[idx][["donor_id", "n_nuclei"]].merge(dn.reset_index(), on="donor_id", how="left")
        d["row"] = idx
        d = d[d.ph_clean.notna() & d.dx.isin(["control", "SZ"]) & d.age.notna()].reset_index(drop=True)
        d["sz"] = (d.dx == "SZ").astype(int)
        if len(d) < 15 or d.sz.sum() < 5 or (1 - d.sz).sum() < 5:
            continue
        E = expr_all[d.row.to_numpy()]
        det = detect_all[d.row.to_numpy()].mean(0) >= DETECT
        # restrict to representative symbol columns that pass detection
        cols = {g: c for g, c in sym2col.items() if det[c]}
        for sname, sg in SETS.items():
            gl = [g for g in sg if g in cols and g in zph.index]
            if len(gl) < 30:
                continue
            ci = np.array([cols[g] for g in gl]); Zv = zph.loc[gl].to_numpy(float)
            sg_sign = sign209.reindex(gl).to_numpy(float) if sname == "study_brain_209" else None
            Ysub = E[:, ci]
            for model, covs in [("age_sex", ["age", "male"]), ("age_sex_pmi" + ("_rin" if cohort == "GSE254569" else ""),
                                                              ["age", "male", "pmi"] + (["rin"] if cohort == "GSE254569" else []))]:
                dm = d.dropna(subset=covs).reset_index(drop=True)
                if len(dm) < 15 or dm.sz.sum() < 5 or (1 - dm.sz).sum() < 5:
                    continue
                keep = d.donor_id.isin(dm.donor_id).to_numpy()
                Y = Ysub[keep]
                def stats(Yv, dv):
                    out = {}
                    t0 = ols_t(Yv, design(dv, ["sz"] + covs), 1)[1]
                    t1 = ols_t(Yv, design(dv, ["sz", "ph_clean"] + covs), 1)[1]
                    r0, r1 = corr(t0, Zv), corr(t1, Zv)
                    out = [r0, r1, r1 - r0]
                    if sg_sign is not None:
                        m = np.isfinite(sg_sign)
                        out += [np.nanmean(sg_sign[m] * t0[m]), np.nanmean(sg_sign[m] * t1[m])]
                    else:
                        out += [np.nan, np.nan]
                    return out, t0, t1
                est, t0, t1 = stats(Y, dm)
                bs = []
                nm = len(dm)
                for _ in range(B):
                    ix = rng.integers(0, nm, nm)
                    try:
                        bs.append(stats(Y[ix], dm.iloc[ix].reset_index(drop=True))[0])
                    except np.linalg.LinAlgError:
                        pass
                bs = np.array(bs)
                # null for the change: pH permuted among donors within diagnosis group (keeps its group-mean difference, breaks gene links)
                nl = []
                for _ in range(NPERM):
                    dp = dm.copy()
                    for g_ in (0, 1):
                        m_ = (dp.sz == g_).to_numpy(); dp.loc[m_, "ph_clean"] = rng.permutation(dp.loc[m_, "ph_clean"].to_numpy())
                    t1p = ols_t(Y, design(dp, ["sz", "ph_clean"] + covs), 1)[1]
                    nl.append(corr(t1p, Zv) - est[0])
                row = dict(cohort=cohort, cell_class=clab, gene_set=sname, model=model, n_donors=nm, n_sz=int(dm.sz.sum()),
                           n_ctrl=int((1 - dm.sz).sum()), n_genes=len(gl),
                           ph_sd=dm.ph_clean.std(), spearman_noPH=spearmanr(t0, Zv, nan_policy="omit")[0], spearman_PH=spearmanr(t1, Zv, nan_policy="omit")[0])
                row["r_change_null_mean"] = float(np.mean(nl)); row["r_change_null_lo"], row["r_change_null_hi"] = pct(nl)
                for i, nme in enumerate(["r_noPH", "r_PH", "r_change", "mean_signed_t_noPH", "mean_signed_t_PH"]):
                    lo, hi = pct(bs[:, i]) if np.isfinite(bs[:, i]).any() else (np.nan, np.nan)
                    row[nme] = est[i]; row[nme + "_lo"] = lo; row[nme + "_hi"] = hi
                rowsA.append(row)
        # ---- pH signature within the class (t of pH given age, sex[, dx])
        for grp in ["controls", "cases", "all_dx_adjusted"]:
            sub = d if grp == "all_dx_adjusted" else d[d.sz == (1 if grp == "cases" else 0)]
            sub = sub.reset_index(drop=True)
            covs = ["age", "male"] + (["sz"] if grp == "all_dx_adjusted" else [])
            if len(sub) < 15:
                continue
            Esub = E[d.donor_id.isin(sub.donor_id).to_numpy()]
            for sname, sg in SETS.items():
                gl = [g for g in sg if g in cols and g in zph.index]
                if len(gl) < 30:
                    continue
                ci = np.array([cols[g] for g in gl]); Zv = zph.loc[gl].to_numpy(float)
                Y = Esub[:, ci]
                f = lambda Yv, dv: corr(ols_t(Yv, design(dv, ["ph_clean"] + covs), 1)[1], Zv)
                est = f(Y, sub)
                bs = []
                for _ in range(B):
                    ix = rng.integers(0, len(sub), len(sub))
                    try:
                        bs.append(f(Y[ix], sub.iloc[ix].reset_index(drop=True)))
                    except np.linalg.LinAlgError:
                        pass
                lo, hi = pct(bs)
                rowsP.append(dict(cohort=cohort, cell_class=clab, group=grp, gene_set=sname, n_donors=len(sub), n_genes=len(gl),
                                  ph_sd=sub.ph_clean.std(), r_tpH_vs_ZpH=est, lo=lo, hi=hi))
        # ---- B: marker expression within interneurons
        if clab == "interneuron":
            for gene in ["SST", "PVALB", "GAD1"]:
                if gene not in sym2col:
                    continue
                y = E[:, sym2col[gene]][:, None]
                for grp in ["controls", "cases", "all_dx_adjusted"]:
                    mask = (d.sz == (1 if grp == "cases" else 0)).to_numpy() if grp != "all_dx_adjusted" else np.ones(len(d), bool)
                    sub = d[mask].reset_index(drop=True); ys = y[mask]
                    covs = ["age", "male"] + (["sz"] if grp == "all_dx_adjusted" else [])
                    variants = [("age_sex", covs)]
                    shcol = {"SST": "share_SST_of_interneurons", "PVALB": "share_PVALB_of_interneurons"}.get(gene)
                    if cohort == "GSE254569" and shcol:
                        sub = sub.copy(); sub["logit_share"] = np.log(sub[shcol] / (1 - sub[shcol]))
                        variants.append(("age_sex_plus_subtype_share", covs + ["logit_share"]))
                    for vname, vcov in variants:
                        f = lambda Yv, dv: ols_t(Yv, design(dv, ["ph_clean"] + vcov), 1)[0][0]
                        est = f(ys, sub)
                        bs = []
                        for _ in range(B):
                            ix = rng.integers(0, len(sub), len(sub))
                            try:
                                bs.append(f(ys[ix], sub.iloc[ix].reset_index(drop=True)))
                            except np.linalg.LinAlgError:
                                pass
                        lo, hi = pct(bs)
                        rowsE.append(dict(cohort=cohort, gene=gene, group=grp, model=vname, n=len(sub), ph_sd=sub.ph_clean.std(),
                                          slope_log2CPM_per_0p1_pH=est * 0.1, lo=lo * 0.1, hi=hi * 0.1,
                                          mean_log2CPM=float(ys.mean()), sd_log2CPM=float(ys.std()),
                                          spearman_unadjusted=spearmanr(sub.ph_clean, ys[:, 0])[0]))
                # SZ difference before/after pH, same donors
                for vname, extra in [("age_sex", []), ("age_sex_plus_subtype_share", ["logit_share"])]:
                    if vname != "age_sex" and not (cohort == "GSE254569" and shcol):
                        continue
                    dd = d.copy()
                    if extra:
                        dd["logit_share"] = np.log(dd[shcol] / (1 - dd[shcol]))
                    def f2(Yv, dv):
                        b0 = ols_t(Yv, design(dv, ["sz", "age", "male"] + extra), 1)[0][0]
                        b1 = ols_t(Yv, design(dv, ["sz", "ph_clean", "age", "male"] + extra), 1)[0][0]
                        return [b0, b1, b1 - b0]
                    est = f2(y, dd); bs = []
                    for _ in range(B):
                        ix = rng.integers(0, len(dd), len(dd))
                        try:
                            bs.append(f2(y[ix], dd.iloc[ix].reset_index(drop=True)))
                        except np.linalg.LinAlgError:
                            pass
                    bs = np.array(bs)
                    for i, qn in enumerate(["sz_minus_control", "sz_minus_control_after_pH", "change_after_pH"]):
                        lo, hi = pct(bs[:, i])
                        rowsE.append(dict(cohort=cohort, gene=gene, group="sz_difference", model=vname + ":" + qn, n=len(dd),
                                          ph_sd=dd.ph_clean.std(), slope_log2CPM_per_0p1_pH=est[i], lo=lo, hi=hi,
                                          mean_log2CPM=float(y.mean()), sd_log2CPM=float(y.std()), spearman_unadjusted=np.nan))
    print(cohort, "done", flush=True)

pd.concat(donor_tabs).to_csv(f"{OUT}/donor_level_table.tsv", sep="\t", index=False)
for nm, rows in [("A_sz_effect_vs_ZpH_with_without_pH", rowsA), ("A_pH_signature_within_class", rowsP), ("B_composition_vs_pH", rowsF),
                 ("B_interneuron_markers_vs_pH", rowsE), ("B_interneuron_fraction_sz_difference", rowsD)]:
    pd.DataFrame(rows).to_csv(f"{OUT}/{nm}.tsv", sep="\t", index=False)
meta.update(dict(B=B, seed=SEED, min_nuclei=MIN_NUC, detect=DETECT))
json.dump(meta, open(f"{OUT}/run_meta.json", "w"), indent=1)
