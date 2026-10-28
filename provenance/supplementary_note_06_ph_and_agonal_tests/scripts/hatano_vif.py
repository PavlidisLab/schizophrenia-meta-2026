#!/usr/bin/env python3
"""STAT-012, Hatano set test: inter-gene correlation of the 209 recurrent genes (sign-aligned), VIF, corrected z and P.

The Hatano statistic (score_176.py, Supplementary Table 3): for each recurrent gene g with a Hatano pH correlation r_g,
s_g = -(SZ sign_g) x r_g; set score = mean of s_g over the m genes with a value; z = (score - null mean) / null SD against 10,000
expression-decile-matched random gene sets. Random gene sets carry no correlation among their genes, so the null SD is that of a mean of
independent genes. Under the null of no pH effect, cov(s_g, s_h) = sign_g sign_h cor(x_g, x_h) (x = expression across samples), so the
correlation of the score contributions is the sign-aligned expression correlation after removing pH, age, sex and diagnosis.
rho_bar = mean over gene pairs of sign_g sign_h cor(x_g, x_h);  VIF = 1 + (m - 1) rho_bar  (Wu and Smyth 2012, CAMERA);
corrected z = z / sqrt(VIF); corrected P = one-sided normal tail of the corrected z.

Reference expression (primary): the five independent collections of the project's pH analysis (GSE21138, GSE35977, GSE53987 BA46, GSE17612,
GSE92538 GPL96; cases and controls with pH, age, sex), preprocessed as analyses/within_group_ph_correlations/within_group_ph.R (Gemma processed
expression, HGNC symbol, one row per gene, arrays of one subject averaged, bottom tenth of mean expression dropped). Residuals within each
dataset (after pH, age, sex, diagnosis) are standardized and pooled. Sensitivity: each collection alone; GTEx BA9 and Cortex (RNA-seq;
residuals after ischemic time, age, sex); 95% interval from 1,000 sample bootstraps of the pooled matrix.
Also: a DE-prior-matched null for the same statistic (20 bins of the global DE prior in place of expression deciles).

Run:  PYTHONPATH=<dir with pyxlsb> python3 hatano_vif.py        (from this folder; needs gtex_extract/ from gtex_extract_remote.py)"""
import glob, gzip, json, os, re, sys, zipfile, shutil, tempfile
import numpy as np, pandas as pd
from pyxlsb import open_workbook
from scipy.stats import norm
from common import sha, mean_offdiag, vif, write_tsv, today

HERE = os.path.dirname(os.path.abspath(__file__))
AN = os.path.dirname(HERE); REPO = os.path.dirname(AN)
CONF = os.path.join(AN, "confounds")
rel = lambda p: os.path.relpath(p, REPO)
SUPP3 = os.path.join(AN, "recurrent_set_recomputations/3_ph_agonal_axes/supp_table3_ph_agonal_axes_209.tsv")
RSET = os.path.join(CONF, "recurrent_sets.tsv")
RECG = os.path.join(AN, "recurrence_genes.tsv")
HGNC = os.path.join(REPO, "data/reference/hgnc_complete_set.txt")
HAT = os.path.join(REPO, "data/raw/40921414/supplement/NPR2-45-e70053-s001.zip")
PRIOR = os.path.join(AN, "de_prior_lists/our_prior/1_global_our_prior.tsv")
GT = os.path.join(HERE, "gtex_extract")
SEED, NBOOT, NPERM, NSIM = 20261012, 1000, 10000, 5000
DS = [("GSE21138", 1880, "GSE21138", None), ("GSE35977", 5949, "GSE35977", None), ("GSE53987_BA46", 8359, "GSE53987", "Pre-frontal cortex (BA46)"),
      ("GSE17612", 1381, "GSE17612", None), ("GSE92538_GPL96", 18741, "GSE92538_GPL17027", None)]
AXES = {"Hatano SZ pH r": "0004", "Hatano control pH r": "0003", "Hatano pooled pH r": "0002", "Hatano BD pH r": "0001"}

# ---- gene set and SZ sign (recurrent_sets.tsv, set study_brain; sign = sign(up - down), as score_176.py)
rs = pd.read_csv(RSET, sep="\t", comment="#")
rs = rs[rs.set == "study_brain"]
assert len(rs) == 209
SIGN = {r.gene: int(np.sign(r.up - r.down)) for r in rs.itertuples() if r.up != r.down}
assert len(SIGN) == 208, len(SIGN)

hg = pd.read_csv(HGNC, sep="\t", usecols=["symbol", "entrez_id"], dtype=str).dropna()
ent2sym = dict(zip(hg.entrez_id.str.replace(r"\.0$", "", regex=True), hg.symbol))


def load_ds(name, gid, covname, region):
    with gzip.open(os.path.join(CONF, f"data/gemma/{gid}_processed.tsv.gz"), "rt") as fh:
        x = pd.read_csv(fh, sep="\t", comment="#", low_memory=False)
    x = x[x.NCBIid.notna() & ~x.NCBIid.astype(str).str.contains(r"\|")]
    x["sym"] = x.NCBIid.astype(float).astype(int).astype(str).map(ent2sym)
    x = x[x.sym.notna()]
    sc = list(x.columns[6:-1])  # last column is the symbol added above
    if all("GSM" in c for c in sc):
        gsm = {c: re.search(r"GSM[0-9]+", c).group() for c in sc}
    else:
        sj = json.load(open(os.path.join(CONF, f"data/gemma/{gid}_samples.json")))["data"]
        acc = {str(s["id"]): s["accession"]["accession"] for s in sj}
        gsm = {c: acc.get(re.sub(r".*BioAssayId=([0-9]+).*", r"\1", c)) for c in sc}
    M = x[sc].astype(float); M.columns = [gsm[c] for c in sc]; M.index = x.sym.values
    cv = pd.read_csv(os.path.join(CONF, f"data/covariates/{covname}.tsv"), sep="\t", comment="#", dtype=str)
    if "gemma_outlier" in cv: cv = cv[~cv.gemma_outlier.isin(["yes", "TRUE", "True", "1"])]
    if region: cv = cv[cv.region == region]
    cv = cv[cv.gsm.isin(M.columns) & cv.dx.isin(["control", "schizophrenia"])].drop_duplicates("gsm")
    for c in ("age", "pH"): cv[c] = pd.to_numeric(cv[c], errors="coerce")
    cv = cv.dropna(subset=["age", "pH", "sex"])
    E = M.loc[:, ~M.columns.duplicated()][cv.gsm.values]
    ok = (E.isna().mean() < 0.5).values; E = E.loc[:, ok]; cv = cv[ok]
    if "subject" in cv and cv.subject.fillna("").ne("").any() and cv.subject[cv.subject.fillna("") != ""].duplicated().any():
        E = E.T.groupby(cv.subject.values).mean().T
        cv = cv.drop_duplicates("subject"); E.columns = cv.gsm.values
    mu = E.mean(axis=1)
    E = E.loc[mu.sort_values(ascending=False).index]; E = E[~E.index.duplicated()]
    E = E.dropna(axis=0); mu = E.mean(axis=1); E = E[mu > mu.quantile(0.10)]
    return E, cv.set_index("gsm").loc[E.columns]


def resid_std(E, cols):
    X = np.column_stack([np.ones(E.shape[1])] + cols)
    Y = E.values.T
    R = Y - X @ np.linalg.lstsq(X, Y, rcond=None)[0]
    R = (R - R.mean(0)) / R.std(0, ddof=1)
    return pd.DataFrame(R, index=E.columns, columns=E.index)


def resid_collection(E, cv):
    cols = [cv.pH.values, cv.age.values]
    if cv.sex.nunique() > 1: cols.append((cv.sex == "M").astype(float).values)
    if cv.dx.nunique() > 1: cols.append((cv.dx == "schizophrenia").astype(float).values)
    return resid_std(E, cols)


def resid_gtex(path, smp):
    E = pd.read_csv(path, sep="\t", index_col=0)
    s = smp.set_index("SAMPID").loc[E.columns]
    return resid_std(E, [s.isch_h.values, s.age.values, s.male.values])


def rho_bar(R, genes):
    g = [x for x in genes if x in R.columns]
    s = np.array([SIGN[x] for x in g])
    with np.errstate(all="ignore"):
        C = R[g].corr(min_periods=10).values * np.outer(s, s)
    rb, npairs = mean_offdiag(C)
    return rb, npairs, len(g), C, g


def boot_rho(R, genes, rng, nb):
    g = [x for x in genes if x in R.columns]
    keep = R[g].isna().mean() < 0.2
    g = [x for x in g if keep[x]]
    s = np.array([SIGN[x] for x in g]); A = R[g].fillna(0.0).values
    out = []
    for _ in range(nb):
        B = A[rng.integers(0, len(A), len(A))]
        with np.errstate(all="ignore"):
            C = np.corrcoef(B, rowvar=False) * np.outer(s, s)
        out.append(mean_offdiag(C)[0])
    return np.percentile(out, [2.5, 97.5])


def hatano_r():
    tmp = tempfile.mkdtemp()
    out = {}
    with zipfile.ZipFile(HAT) as z: z.extractall(tmp)
    for name, f in AXES.items():
        p = glob.glob(os.path.join(tmp, f"*{f}*"))[0]
        with open_workbook(p) as wb, wb.get_sheet(1) as sh:
            rows = [[c.v for c in r] for r in sh.rows()]
        hdr = rows[1]; gi, pi = hdr.index("GeneSymbol"), hdr.index("pH")
        d = {}
        for r in rows[2:]:
            if len(r) > pi and r[gi] and isinstance(r[pi], (int, float)) and r[pi] == r[pi]:
                d.setdefault(str(r[gi]).strip(), []).append(float(r[pi]))
        out[name] = {g: float(np.mean(v)) for g, v in d.items()}
    shutil.rmtree(tmp)
    return out


def main():
    rng = np.random.default_rng(SEED)
    s3 = pd.read_csv(SUPP3, sep="\t", comment="#").set_index("axis")
    r_ax = hatano_r()
    # ---- reproduce the observed score of the committed table from the gene list and signs
    repro = {}
    for ax, d in r_ax.items():
        gs = [g for g in SIGN if g in d]
        sc = np.mean([-SIGN[g] * d[g] for g in gs])
        repro[ax] = (gs, sc)
        assert len(gs) == s3.loc[ax, "n_with_value"], (ax, len(gs), s3.loc[ax, "n_with_value"])
        assert abs(sc - s3.loc[ax, "mean_score"]) < 6e-5, (ax, sc, s3.loc[ax, "mean_score"])
        print("reproduced", ax, len(gs), round(sc, 4))
    # ---- reference expression
    R, info = {}, []
    for name, gid, cov, reg in DS:
        E, cv = load_ds(name, gid, cov, reg)
        R[name] = resid_collection(E, cv)
        info.append(dict(reference=name, kind="microarray collection", n_samples=len(cv), n_cases=int((cv.dx == "schizophrenia").sum()), n_controls=int((cv.dx == "control").sum()), genes_in_matrix=E.shape[0]))
    P = pd.concat(R.values(), axis=0, sort=False)
    smp = pd.read_csv(os.path.join(GT, "samples.tsv"), sep="\t")
    G = {}
    for reg, f in (("GTEx BA9", "expr_Brain__Frontal_Cortex_BA9.tsv.gz"), ("GTEx Cortex", "expr_Brain__Cortex.tsv.gz")):
        G[reg] = resid_gtex(os.path.join(GT, f), smp)
        info.append(dict(reference=reg, kind="RNA-seq GTEx v8", n_samples=G[reg].shape[0], n_cases=None, n_controls=None, genes_in_matrix=G[reg].shape[1]))
    refs = {"five collections pooled (primary)": P, "GTEx BA9": G["GTEx BA9"], "GTEx Cortex": G["GTEx Cortex"], **{k: v for k, v in R.items()}}
    rows = []
    for ax, (gs, sc) in repro.items():
        m = len(gs)
        z0, nsd, nmean = s3.loc[ax, "z"], s3.loc[ax, "null_sd"], s3.loc[ax, "null_mean"]
        zexact = (sc - nmean) / nsd
        for rname, Rm in refs.items():
            rb, npairs, mref, C, gref = rho_bar(Rm, gs)
            lo, hi = boot_rho(Rm, gs, rng, NBOOT) if rname.endswith("(primary)") else (np.nan, np.nan)
            v = vif(m, rb)
            zc = zexact / np.sqrt(v)
            rows.append(dict(axis=ax, reference=rname, m_genes_scored=m, genes_in_reference=mref, n_pairs=npairs, rho_bar=round(rb, 4),
                             rho_bar_lo=round(lo, 4), rho_bar_hi=round(hi, 4), VIF=round(v, 2), z_uncorrected_table=z0, z_uncorrected_exact=round(zexact, 3),
                             z_corrected=round(zc, 2), p_uncorrected_table=s3.loc[ax, "p_one_sided"], p_corrected_one_sided=float(norm.sf(zc)),
                             z_corrected_at_rho_hi=round(zexact / np.sqrt(vif(m, hi)), 2) if hi == hi else np.nan))
            if rname.endswith("(primary)"):
                print(ax, rname, "rho_bar", round(rb, 4), (round(lo, 4), round(hi, 4)), "VIF", round(v, 2), "z", z0, "->", round(zc, 2))
    # ---- sign-block decomposition of rho_bar in the primary reference (SZ-axis genes)
    gs = repro["Hatano SZ pH r"][0]
    rb, npairs, mref, C, gref = rho_bar(P, gs)
    sg = np.array([SIGN[g] for g in gref]); blocks = []
    for a, b, lab in ((1, 1, "up-up"), (-1, -1, "down-down"), (1, -1, "up-down")):
        ia, ib = np.where(sg == a)[0], np.where(sg == b)[0]
        v = [C[i, j] for i in ia for j in ib if (i < j if a == b else True) and np.isfinite(C[i, j])]
        raw = np.mean(v) * (1 if a == b else -1)
        blocks.append(dict(pair_type=lab, n_pairs=len(v), mean_sign_aligned_rho=round(float(np.mean(v)), 4), mean_raw_expression_rho=round(float(raw), 4)))
    # ---- DE-prior-matched null for the same statistic
    pr = pd.read_csv(PRIOR, sep="\t", usecols=["gene_symbol", "de_prior_norm"]).groupby("gene_symbol").de_prior_norm.max()
    rec = pd.read_csv(RECG, sep="\t", comment="#")
    rec = rec[rec.scenario == "published, brain"]
    pool_sign = {r.gene: int(np.sign(r.up - r.down)) for r in rec.itertuples() if r.up != r.down}
    pool_sign.update(SIGN)  # the 209's own signs from recurrent_sets.tsv take precedence
    dec = dict(zip(rec.gene, rec.expr_decile))
    de_rows = []
    for ax, (gs, sc) in repro.items():
        d = r_ax[ax]
        have = [g for g in pool_sign if g in d and g in pr.index and np.isfinite(pr[g])]
        vals = pd.DataFrame({"s": [-pool_sign[g] * d[g] for g in have], "pr": [pr[g] for g in have], "dec": [str(dec.get(g, "")) for g in have]}, index=have)
        qs = np.unique(np.quantile(pr[have].values, np.linspace(0, 1, 21)))
        vals["bin"] = np.clip(np.searchsorted(qs, vals.pr.values, side="right") - 1, 0, len(qs) - 2)
        gsin = [g for g in gs if g in vals.index]
        for lab, key in (("DE prior (20 quantile bins)", "bin"), ("expression decile (as score_176.py)", "dec")):
            need = vals.loc[gsin, key].value_counts()
            by = {k: vals.index[vals[key] == k].values for k in need.index}
            rg = np.random.default_rng(SEED + (0 if key == "bin" else 1))
            nul = np.array([vals.loc[np.concatenate([rg.choice(by[k], c, replace=False) for k, c in need.items()]), "s"].mean() for _ in range(NPERM)])
            obs = vals.loc[gsin, "s"].mean()
            m = len(gsin)
            rb = rho_bar(P, gsin)[0]
            v = vif(m, rb)
            z = (obs - nul.mean()) / nul.std(ddof=1)
            de_rows.append(dict(axis=ax, null=lab, m_genes=m, mean_score=round(obs, 4), null_mean=round(nul.mean(), 4), null_sd=round(nul.std(ddof=1), 4), z=round(z, 2),
                                p_one_sided=(1 + (nul >= obs).sum()) / (1 + NPERM), VIF_primary=round(v, 2), z_corrected=round(z / np.sqrt(v), 2),
                                p_corrected_one_sided=float(norm.sf(z / np.sqrt(v))), median_pool_de_prior_of_set=round(float(vals.loc[gsin, "pr"].median()), 3), median_pool_de_prior_of_pool=round(float(vals.pr.median()), 3)))
            print(de_rows[-1])
    # ---- direct check of the VIF: draw 25 samples from a reference set, attach a random standard-normal "pH" (no association), compute the
    # Hatano-style score over the SZ-axis genes; the SD of that score over draws is the null SD of the set score with the gene correlation kept.
    gs = repro["Hatano SZ pH r"][0]
    nsd_indep = float(s3.loc["Hatano SZ pH r", "null_sd"])
    sim = []
    for rname, Rm in refs.items():
        if rname.endswith("(primary)"): continue
        g = [x for x in gs if x in Rm.columns]
        A = Rm[g].fillna(0.0).values; sg_ = np.array([SIGN[x] for x in g])
        n = min(25, A.shape[0]); sc_ = []
        rg = np.random.default_rng(SEED + 7)
        for _ in range(NSIM):
            B = A[rg.choice(A.shape[0], n, replace=False)]; v = rg.standard_normal(n)
            Bc = (B - B.mean(0)); vc = v - v.mean()
            r = (Bc * vc[:, None]).sum(0) / np.sqrt((Bc ** 2).sum(0) * (vc ** 2).sum())
            sc_.append(np.mean(-sg_ * r))
        sd = float(np.std(sc_, ddof=1))
        rb = rho_bar(Rm, gs)[0]
        sim.append(dict(reference=rname, n_per_draw=n, draws=NSIM, m_genes=len(g), sd_of_score_simulated=round(sd, 4), sd_independent_genes_table=nsd_indep,
                        ratio_sq=round((sd / nsd_indep) ** 2, 1), VIF_formula=round(vif(len(g), rb), 1), z_observed_over_simulated_sd=round(float((repro["Hatano SZ pH r"][1] - s3.loc["Hatano SZ pH r", "null_mean"]) / sd), 2)))
        print(sim[-1])
    return pd.DataFrame(rows), pd.DataFrame(info), pd.DataFrame(blocks), pd.DataFrame(de_rows), pd.DataFrame(sim)


if __name__ == "__main__":
    t, info, blocks, de, sim = main()
    hdr = ["claim: STAT-012 - Hatano pH set-test z values (Supplementary Table 3) corrected for inter-gene correlation with a variance inflation factor",
           "script: analyses/gene_correlation_adjusted_tests/hatano_vif.py (common.py), run " + today(),
           f"source: {rel(SUPP3)} sha256 {sha(SUPP3)[:12]}; {rel(RSET)} sha256 {sha(RSET)[:12]} (set study_brain, 209 genes); {rel(HAT)} sha256 {sha(HAT)[:12]}",
           f"source: {rel(RECG)} sha256 {sha(RECG)[:12]} (pool signs and expr_decile of the DE-prior null only); {rel(PRIOR)} sha256 {sha(PRIOR)[:12]}",
           "source: Gemma processed expression of GSE21138 (1880), GSE35977 (5949), GSE53987 (8359), GSE17612 (1381), GSE92538 GPL96 (18741) under analyses/ph_and_confound_analyses/data/gemma/ with analyses/ph_and_confound_analyses/data/covariates/*.tsv",
           f"source: GTEx v8 extract analyses/gene_correlation_adjusted_tests/gtex_extract/ (extract_meta.json: TPM sha256 {json.load(open(os.path.join(GT, 'extract_meta.json')))['tpm_sha256'][:12]})",
           "method: rho_bar = mean over gene pairs of sign_g sign_h cor(x_g, x_h), x = residual expression after pH, age, sex, diagnosis (collections) or ischemic time, age, sex (GTEx); VIF = 1 + (m - 1) rho_bar (Wu and Smyth 2012, Nucleic Acids Research 40:e133); z_corrected = z / sqrt(VIF); P = one-sided normal tail",
           f"seed {SEED}; sample bootstraps {NBOOT}; DE-prior null draws {NPERM}"]
    write_tsv(os.path.join(HERE, "hatano_vif.tsv"), hdr, t)
    write_tsv(os.path.join(HERE, "hatano_reference_sets.tsv"), hdr[:2] + ["claim: reference expression sets used for rho_bar"], info)
    write_tsv(os.path.join(HERE, "hatano_rho_bar_blocks.tsv"), hdr[:2] + ["claim: rho_bar of the Hatano SZ gene set split by pair type, primary reference (sign-aligned = sign_g sign_h cor)"], blocks)
    write_tsv(os.path.join(HERE, "hatano_de_prior_null.tsv"), hdr[:2] + ["claim: Hatano score against a DE-prior-matched null (20 quantile bins of de_prior_norm) beside the expression-decile null; new computation, no earlier file holds it", hdr[3], hdr[7]], de)
    write_tsv(os.path.join(HERE, "hatano_vif_simulation_check.tsv"), hdr[:2] + ["claim: direct check of the VIF - SD of the Hatano-style SZ-axis score when 25 samples of a reference set are given a random pH (no association), against the SD for independent genes in Supplementary Table 3", hdr[7]], sim)
