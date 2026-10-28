#!/usr/bin/env python3
"""Step 1: per-gene effects of the Yaqubi et al. 2026 (Nat Commun 17:2848) processing conditions.

Inputs (all open access, from the article page https://www.nature.com/articles/s41467-026-68872-9):
  MOESM3_ESM.xlsx  Supplementary Data 1-7 (DEG tables; per-condition DEG lists)
  MOESM6_ESM.zip   Source Data (raw count matrices behind Supp. Figs 2 and 6; whole-transcriptome log2FC behind Fig. 1D and Supp. Fig. 4C)

Two kinds of effect are produced, both signed "condition minus control" (positive = higher after processing or after death):
  (a) own: OLS on log2(median-of-ratios-normalised count + 1) from the released raw counts (t statistic and log2 difference);
      surgical series: donor fixed effects + condition (paired by donor), 0 h / -80 C is the reference;
      autopsy series: ~ group + age + sex + RIN, as the authors' DESeq2 design (~Age+Sex+RIN+Group), CTRL as reference.
  (b) released: log2FoldChange from the Source Data (24 h 20 C, 36 h 20 C, SPMI, LPMI; whole transcriptome, no p value) for validation.
Usage: python3 01_build_effects.py <dir with MOESM3_ESM.xlsx and unzipped Source Data> <out dir>
"""
import hashlib, os, sys, re, datetime
sys.dont_write_bytecode = True
import numpy as np
import pandas as pd

SRC, OUT = sys.argv[1], sys.argv[2]
AN = os.path.expanduser("<analysis repo>/analyses")
sys.path.insert(0, AN)
import extract
extract.HGNC = extract.load_hgnc()
ENS = extract.HGNC[3]
MIN_MEAN = 10.0   # mean normalised count across the samples of a matrix


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()[:12]


def size_factors(C):
    L = np.log(C.astype(float))
    ok = np.isfinite(L).all(1)
    ref = L[ok].mean(1)
    return np.exp(np.median(L[ok] - ref[:, None], axis=0))


def ols_contrast(Y, X, cvec):
    """Y genes x n, X n x p, cvec p. Returns estimate, t, df."""
    XtXi = np.linalg.pinv(X.T @ X)
    B = XtXi @ X.T @ Y.T                      # p x genes
    R = Y.T - X @ B
    df = X.shape[0] - np.linalg.matrix_rank(X)
    s2 = (R ** 2).sum(0) / df
    est = cvec @ B
    se = np.sqrt(s2 * (cvec @ XtXi @ cvec))
    return est, est / se, df


def map_ids(ids):
    return [ENS.get(i.split(".")[0], ("", ""))[0] for i in ids]


def to_sym_table(ids, cols, meanexpr):
    d = pd.DataFrame(cols)
    d["gene"] = map_ids(ids)
    d["_m"] = meanexpr
    d = d[d.gene != ""].sort_values("_m", ascending=False).drop_duplicates("gene").drop(columns="_m")
    return d.set_index("gene")


res, info = {}, []

# ------------------------------------------------------------------ surgical series
sg = pd.read_pickle(os.path.join(SRC, "surg_counts.pkl")) if os.path.exists(os.path.join(SRC, "surg_counts.pkl")) else \
    pd.read_excel(os.path.join(SRC, "Supplementary Figure 6A_B.xlsx"), sheet_name="All_surgery_samples", engine="openpyxl")
ids = sg["ID"].astype(str).tolist()
C = sg.drop(columns="ID").values.astype(float)
names = sg.columns[1:].tolist()
sf = size_factors(C)
N = C / sf
keep = N.mean(1) >= MIN_MEAN
Y = np.log2(N[keep] + 1)
kid = [i for i, k in zip(ids, keep) if k]
meta = []
for n in names:
    m = re.match(r"(\d+)h_(\d+)C_S(\d+)(_D)?$", n)
    meta.append(dict(sample=n, hours=int(m.group(1)), temp=int(m.group(2)), donor=int(m.group(3)), dup=bool(m.group(4))))
meta = pd.DataFrame(meta)
meta["cond"] = np.where(meta.hours == 0, "0h", meta.hours.astype(str) + "h" + meta.temp.astype(str) + "C")
meta["cond"] = meta.cond.replace({"0h": "0h_-80C"})
conds = ["6h4C", "6h20C", "24h4C", "24h20C", "36h4C", "36h20C"]
donors = sorted(meta.donor.unique())
cols = [f"donor{d}" for d in donors] + conds
X = np.zeros((len(meta), len(cols)))
for i, r in meta.iterrows():
    X[i, donors.index(r.donor)] = 1
    if r.cond in conds:
        X[i, len(donors) + conds.index(r.cond)] = 1
meta.to_csv(os.path.join(OUT, "surgical_sample_design.tsv"), sep="\t", index=False)
mean_expr = N[keep].mean(1)


def cvec(**kw):
    v = np.zeros(len(cols))
    for k, w in kw.items():
        v[cols.index(k)] = w
    return v


surg_contrasts = {
    "S_6h4C_vs_0h": cvec(**{"6h4C": 1}), "S_6h20C_vs_0h": cvec(**{"6h20C": 1}), "S_24h4C_vs_0h": cvec(**{"24h4C": 1}),
    "S_24h20C_vs_0h": cvec(**{"24h20C": 1}), "S_36h4C_vs_0h": cvec(**{"36h4C": 1}), "S_36h20C_vs_0h": cvec(**{"36h20C": 1}),
    "S_20C_minus_4C_6h": cvec(**{"6h20C": 1, "6h4C": -1}), "S_20C_minus_4C_24h": cvec(**{"24h20C": 1, "24h4C": -1}),
    "S_20C_minus_4C_36h": cvec(**{"36h20C": 1, "36h4C": -1}),
    "S_20C_36h_minus_6h": cvec(**{"36h20C": 1, "6h20C": -1}), "S_4C_36h_minus_6h": cvec(**{"36h4C": 1, "6h4C": -1}),
}
for nm, cv in surg_contrasts.items():
    est, t, df = ols_contrast(Y, X, cv)
    d = to_sym_table(kid, {f"lfc_{nm}": est, f"t_{nm}": t}, mean_expr)
    res[nm] = d
    info.append(dict(contrast=nm, series="surgical, donor-paired", n_samples=len(meta), resid_df=df, genes_tested=int(keep.sum()), genes_mapped=len(d),
                     definition=f"{nm}: OLS log2(norm+1) ~ donor + condition, donor fixed effects; est = condition difference (log2); min mean norm count {MIN_MEAN}"))

# ------------------------------------------------------------------ autopsy series
def autopsy(name, xl, count_sheet, meta_sheet, idcol, group_col, label):
    cnt = pd.read_excel(os.path.join(SRC, xl), sheet_name=count_sheet, engine="openpyxl")
    md = pd.read_excel(os.path.join(SRC, xl), sheet_name=meta_sheet, engine="openpyxl")
    gids = cnt.iloc[:, 0].astype(str).tolist()
    cnt = cnt.iloc[:, 1:]
    md = md.set_index(md.columns[0]).loc[cnt.columns]
    C = cnt.values.astype(float)
    N = C / size_factors(C)
    keep = N.mean(1) >= MIN_MEAN
    Y = np.log2(N[keep] + 1)
    grp = (md[group_col].astype(str) != "CTRL").astype(float).values
    sex = md["Sex"].astype(str).str.upper().str[0].eq("M").astype(float).values
    Xa = np.column_stack([np.ones(len(md)), grp, md["Age"].astype(float).values, sex, md["RIN"].astype(float).values])
    est, t, df = ols_contrast(Y, Xa, np.array([0, 1, 0, 0, 0.0]))
    d = to_sym_table([g for g, k in zip(gids, keep) if k], {f"lfc_{name}": est, f"t_{name}": t}, N[keep].mean(1))
    res[name] = d
    info.append(dict(contrast=name, series=label, n_samples=len(md), resid_df=df, genes_tested=int(keep.sum()), genes_mapped=len(d),
                     definition=f"{name}: OLS log2(norm+1) ~ group + age + sex + RIN (authors' DESeq2 design), est = autopsy minus CTRL (log2); groups n = {int(grp.sum())} autopsy, {int((1-grp).sum())} CTRL"))


autopsy("A_LPMI_vs_CTRL", "Supplementary Figure 2F-H.xlsx", "LPMI_CTRL_count", "metadata", None, "Group", "autopsy long PMI (mean 36 h), UK banks, mid temporal gyrus, vs surgical CTRL")
autopsy("A_SPMI_vs_CTRL", "Supplementary Figure 2B-D.xlsx", "SPMI_CTRL_count", "SPMI_CTRL_count_metadata", None, "Group", "autopsy short PMI (mean 6 h), Netherlands Brain Bank, superior temporal gyrus, vs surgical CTRL")

# ------------------------------------------------------------------ released whole-transcriptome log2FC
rel = {}
for sheet, nm in (("24h20c_vs_CTRL", "S_24h20C_vs_0h"), ("36h20c_vs_CTRL", "S_36h20C_vs_0h"), ("SPMI_vs_CTRL", "A_SPMI_vs_CTRL"), ("LPMI_vs_CTRL", "A_LPMI_vs_CTRL")):
    d = pd.read_excel(os.path.join(SRC, "Supplementary Figure 4C.xlsx"), sheet_name=sheet, engine="openpyxl")
    d = d[d.baseMean >= MIN_MEAN]
    t = to_sym_table(d["ID"].astype(str).tolist(), {f"relLFC_{nm}": d["log2FoldChange"].values}, d["baseMean"].values)
    rel[nm] = t
E = pd.concat(list(res.values()) + list(rel.values()), axis=1, join="outer")
E.index.name = "gene"
# validation of own against released log2FC
for r in info:
    nm = r["contrast"]
    if nm in rel:
        j = E[[f"lfc_{nm}", f"relLFC_{nm}"]].dropna()
        from scipy import stats
        r["validation_vs_released_lfc_n"] = len(j)
        r["validation_pearson"] = round(float(np.corrcoef(j.iloc[:, 0], j.iloc[:, 1])[0, 1]), 3)
        r["validation_spearman"] = round(float(stats.spearmanr(j.iloc[:, 0], j.iloc[:, 1])[0]), 3)
today = datetime.date.today().isoformat()
hdr = [f"# script: 01_build_effects.py, run {today}; NEW analysis; no random numbers",
       f"# inputs: Source Data MOESM6_ESM.zip sha256 25591cba3b4a (members: Supplementary Figure 6A_B.xlsx {sha(os.path.join(SRC,'Supplementary Figure 6A_B.xlsx'))}, "
       f"2F-H {sha(os.path.join(SRC,'Supplementary Figure 2F-H.xlsx'))}, 2B-D {sha(os.path.join(SRC,'Supplementary Figure 2B-D.xlsx'))}, 4C {sha(os.path.join(SRC,'Supplementary Figure 4C.xlsx'))})",
       f"# HGNC: data/reference/hgnc_complete_set.txt sha256 {sha(extract.HGNC_FILE)} (Ensembl gene ID without version -> approved symbol; duplicates keep the higher mean count)",
       "# sign: positive = higher in the processed / autopsy sample than in the 0 h -80 C (CTRL) sample; t = OLS t; lfc = log2 difference of log2(norm+1)"]
with open(os.path.join(OUT, "contrast_effects.tsv"), "w") as f:
    f.write("\n".join(hdr) + "\n")
    E.round(5).to_csv(f, sep="\t")
pd.DataFrame(info).to_csv(os.path.join(OUT, "contrast_definitions.tsv"), sep="\t", index=False)
print(pd.DataFrame(info).to_string())
print(E.shape)
