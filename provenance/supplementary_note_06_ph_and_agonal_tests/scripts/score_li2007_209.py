#!/usr/bin/env python3
"""Score the 209 recurrent genes (recurrent_sets.tsv, set study_brain; SZ sign = sign(up - down)) against the Li 2007 agonal-stress gene lists
(li2007_agonal_genes.tsv; s_g = (SZ sign) x t_g, t = low-pH-minus-normal-pH t score, mean over the five regions; positive = the gene moves in SZ the way
it moves with agonal stress / low pH; identical in meaning to -(SZ sign) x (pH-direction statistic) of score_176.py, with the pH-direction statistic = -t).
Null: 10,000 expression-decile-matched random sets from published-brain genes with a SZ sign, statistic over drawn genes that have a Li value (score_176.py).
Rank correlations of the Li statistic with the project control-pH Z and the Pritzker agonal t (gene level).
Variance inflation: rho_bar from hatano_vif.py functions (same five collections + GTEx), VIF = 1 + (m - 1) rho_bar, z_corrected = z / sqrt(VIF).
Run: PYTHONPATH=pyxlsb python3 score_li2007_209.py   (from this folder; pyxlsb is needed only because hatano_vif.py imports it)"""
import collections, gzip, os, random, sys
import numpy as np, pandas as pd
from scipy.stats import norm, spearmanr
HERE = os.path.dirname(os.path.abspath(__file__)); AN = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(AN, "gene_correlation_adjusted_tests"))
import hatano_vif as H
from common import sha, vif, write_tsv, today
SEED, NPERM = 20261012, 10000
SIGN = H.SIGN; REC = H.RECG
PRZ = os.path.join(AN, "axes_pritzker_agonal_2026-10-09/axis_pritzker_agonal_GPL17027.tsv.gz")
CTRL = H.CONF + "/results/combined_controls_pH.tsv"
LI = os.path.join(HERE, "li2007_agonal_genes.tsv"); HG = H.HGNC

# ---- Li statistic by HGNC symbol (approved symbol; else unique previous/alias symbol)
hg = pd.read_csv(HG, sep="\t", usecols=["symbol", "alias_symbol", "prev_symbol"], dtype=str)
appr = set(hg.symbol); alt = collections.defaultdict(set)
for r in hg.itertuples():
    for col in (r.alias_symbol, r.prev_symbol):
        if isinstance(col, str):
            for s in col.split("|"): alt[s].add(r.symbol)
li = pd.read_csv(LI, sep="\t", comment="#"); li = li[li.symbol.notna()]
def to_hgnc(s):
    s = str(s).strip()
    if s in appr: return s, "approved"
    if len(alt.get(s, ())) == 1: return next(iter(alt[s])), "alias_or_previous"
    return None, "unmapped"
m = li.symbol.map(to_hgnc); li["hgnc"] = [a for a, b in m]; li["mapkind"] = [b for a, b in m]
mapinfo = li.mapkind.value_counts().to_dict()
LIT = li.dropna(subset=["hgnc"]).groupby("hgnc").t_mean_5regions.mean().to_dict()
print("Li rows with symbol", len(li), mapinfo, "genes with a value", len(LIT))

# ---- reference axes
ctrl = pd.read_csv(CTRL, sep="\t", comment="#"); ctrl = ctrl[ctrl.factor == "pH"].set_index("gene").Z.to_dict()
with gzip.open(PRZ, "rt") as fh: pz = pd.read_csv(fh, sep="\t", comment="#").groupby("gene").t_agonal.mean().to_dict()

# ---- population and null (score_176.py)
rec = pd.read_csv(REC, sep="\t", comment="#"); rec = rec[rec.scenario == "published, brain"]
dec = {r.gene: r.expr_decile for r in rec.itertuples() if r.gene in SIGN or True}
S = [g for g in SIGN if g in dec and str(dec[g]) not in ("nan", "")]  # the 209 of recurrent_sets.tsv (208 with a SZ sign); deciles from the current recurrence_genes.tsv
print("209-set genes with sign and decile:", len(S))
SA = {r.gene: (1 if r.up > r.down else -1) for r in rec.itertuples() if r.up != r.down}  # all published-brain genes with a SZ sign (score_176.py null pool)
assert all(SA.get(g) == SIGN[g] for g in S), "SZ sign of the 209 differs between recurrent_sets and recurrence_genes"
bg = [g for g in SA if g in dec and str(dec[g]) not in ("nan", "")]
by_dec = collections.defaultdict(list)
for g in bg: by_dec[dec[g]].append(g)
need = collections.Counter(dec[g] for g in S)
def score(genes, ax, mode):
    out = []
    for g in genes:
        if g in ax: out.append(SA[g] * ax[g] if mode == "agonal" else -SA[g] * ax[g])
        elif mode == "agonal0": out.append(0.0)
    return out
def run(name, ax, mode, rng):
    s = score(S, ax, "agonal" if mode == "agonal0" else mode)
    if mode == "agonal0": s = [SIGN[g] * ax.get(g, 0.0) for g in S]
    obs = float(np.mean(s)); conc = sum(x > 0 for x in s)
    null = []
    for _ in range(NPERM):
        d = []
        for k, n in need.items(): d += rng.sample(by_dec[k], n)
        t = ([SA[g] * ax.get(g, 0.0) for g in d] if mode == "agonal0" else score(d, ax, mode))
        if t: null.append(np.mean(t))
    null = np.array(null)
    return dict(axis=name, n_genes=len(S), n_with_value=len(s) if mode != "agonal0" else sum(g in ax for g in S), mean_score=round(obs, 4), concordant=conc, discordant=len(s) - conc,
                null_mean=round(float(null.mean()), 4), null_sd=round(float(null.std()), 4), z=round((obs - null.mean()) / null.std(), 2),
                p_one_sided=round((1 + (null >= obs).sum()) / (1 + len(null)), 4)), s
rng = random.Random(SEED)
rows = []
# validation: the project control-pH row of Supplementary Table 3 must be reproduced in mean score
s3 = pd.read_csv(H.SUPP3, sep="\t", comment="#").set_index("axis")
v, _ = run("project control pH Z (reference)", ctrl, "ph", rng)
print("validation control-pH row:", v["mean_score"], "table", s3.loc["project control pH Z (reference)", "mean_score"], "z", v["z"], "table z", s3.loc["project control pH Z (reference)", "z"])
assert abs(v["mean_score"] - s3.loc["project control pH Z (reference)", "mean_score"]) < 6e-5
rows.append(v)
r1, s1 = run("Li 2007 agonal t (listed genes only)", LIT, "agonal", rng); rows.append(r1)
r2, _ = run("Li 2007 agonal t (genes outside the 2,000 listed = 0)", LIT, "agonal0", rng); rows.append(r2)
LIs = {g: np.sign(t) for g, t in LIT.items()}
r3, _ = run("Li 2007 agonal direction (+1/-1, listed genes only)", LIs, "agonal", rng); rows.append(r3)
sc = pd.DataFrame(rows)

# ---- rank correlations
cors = []
def rc(name, a, b, genes):
    g = [x for x in genes if x in a and x in b and np.isfinite(a[x]) and np.isfinite(b[x])]
    r, p = spearmanr([a[x] for x in g], [b[x] for x in g])
    cors.append(dict(comparison=name, n_genes=len(g), spearman_rho=round(float(r), 3), p=float(p)))
gl = list(LIT)
rc("Li t vs control pH Z (pH-direction: expect negative), all Li genes", LIT, ctrl, gl)
rc("Li t vs Pritzker agonal t (expect positive), all Li genes", LIT, pz, gl)
rc("Li t vs Pritzker agonal t, adjusted for pH (t_agonal_with_pH not used; see script note)", LIT, pz, gl)
rc("Li t vs control pH Z, 209 genes", LIT, ctrl, S)
rc("Li t vs Pritzker agonal t, 209 genes", LIT, pz, S)
rc("control pH Z vs Pritzker agonal t (expect negative), all genes", ctrl, pz, list(ctrl))
rc("SZ-sign-aligned: sign x Li t vs sign x (-control pH Z), 209 genes", {g: SIGN[g] * LIT[g] for g in S if g in LIT}, {g: -SIGN[g] * ctrl[g] for g in S if g in ctrl}, S)
rc("SZ-sign-aligned: sign x Li t vs sign x Pritzker t, 209 genes", {g: SIGN[g] * LIT[g] for g in S if g in LIT}, {g: SIGN[g] * pz[g] for g in S if g in pz}, S)
cors = [c for c in cors if "adjusted for pH" not in c["comparison"]]
cd = pd.DataFrame(cors)
# sign agreement of Li listed genes with the Pritzker agonal t
lg = [g for g in LIT if g in pz]
agree = np.mean([np.sign(LIT[g]) == np.sign(pz[g]) for g in lg])

# ---- variance inflation
print("loading reference expression ...")
R, P = {}, None
for name, gid, cov, reg in H.DS:
    E, cv = H.load_ds(name, gid, cov, reg); R[name] = H.resid_collection(E, cv)
P = pd.concat(R.values(), axis=0, sort=False)
smp = pd.read_csv(os.path.join(H.GT, "samples.tsv"), sep="\t")
G = {k: H.resid_gtex(os.path.join(H.GT, f), smp) for k, f in (("GTEx BA9", "expr_Brain__Frontal_Cortex_BA9.tsv.gz"), ("GTEx Cortex", "expr_Brain__Cortex.tsv.gz"))}
refs = {"five collections pooled (primary)": P, "GTEx BA9": G["GTEx BA9"], "GTEx Cortex": G["GTEx Cortex"], **R}
vrows = []
gs_listed = [g for g in S if g in LIT]
gs_all = list(S)
for ax_name, gs, row in (("Li 2007 agonal t (listed genes only)", gs_listed, r1), ("Li 2007 agonal t (genes outside the 2,000 listed = 0)", gs_all, r2)):
    for rn, Rm in refs.items():
        rb, npairs, mref, C, gref = H.rho_bar(Rm, gs)
        lo, hi = (H.boot_rho(Rm, gs, np.random.default_rng(SEED), 1000) if rn.endswith("(primary)") else (np.nan, np.nan))
        vv = vif(len(gs), rb); zc = row["z"] / np.sqrt(vv)
        vrows.append(dict(axis=ax_name, reference=rn, m_genes=len(gs), genes_in_reference=mref, n_pairs=npairs, rho_bar=round(rb, 4), rho_bar_lo=round(lo, 4), rho_bar_hi=round(hi, 4),
                          VIF=round(vv, 2), z_uncorrected=row["z"], z_corrected=round(zc, 2), p_corrected_one_sided=float(norm.sf(zc))))
vd = pd.DataFrame(vrows)
# the control-pH reference row corrected the same way, for comparison
rb, npairs, mref, C, gref = H.rho_bar(P, [g for g in S if g in ctrl])
vc = vif(len([g for g in S if g in ctrl]), rb)
print("control pH row: rho_bar", round(rb, 4), "VIF", round(vc, 2), "z", v["z"], "->", round(v["z"] / np.sqrt(vc), 2))
hd = [f"claim: 209 recurrent genes (study_brain; {len(SIGN)} with SZ sign) scored against Li 2007 agonal-stress gene lists; null = {NPERM} expression-decile-matched sets, seed {SEED}",
      "script: analyses/li2007_agonal_stress_scoring/score_li2007_209.py", f"source: li2007_agonal_genes.tsv sha256 {sha(LI)}",
      f"source: recurrent_sets.tsv sha256 {sha(H.RSET)}; recurrence_genes.tsv sha256 {sha(REC)}; axis_pritzker_agonal_GPL17027.tsv.gz sha256 {sha(PRZ)}; combined_controls_pH.tsv sha256 {sha(CTRL)}",
      f"Li symbol mapping to HGNC: {mapinfo}; {len(LIT)} genes with a value", f"run: {today()}"]
write_tsv(os.path.join(HERE, "li2007_score_209.tsv"), hd, sc)
write_tsv(os.path.join(HERE, "li2007_rank_correlations.tsv"), hd + [f"sign agreement of Li t with Pritzker agonal t over {len(lg)} shared genes: {agree:.3f}"], cd)
write_tsv(os.path.join(HERE, "li2007_vif.tsv"), hd + [f"control-pH row same reference: rho_bar {rb:.4f}, VIF {vc:.2f}, z {v['z']} -> {v['z']/np.sqrt(vc):.2f}"], vd)
print(sc.to_string()); print(cd.to_string()); print("sign agreement", agree, len(lg)); print(vd.to_string())

# ---- per-gene table for the 209 and listed-membership comparison with the null pool
gt = pd.DataFrame([dict(gene=g, sz_sign=SIGN[g], li_t_mean_5regions=LIT.get(g, np.nan), li_listed=g in LIT,
                        li_direction=("up_in_low_pH" if LIT.get(g, 0) > 0 else "down_in_low_pH") if g in LIT else "",
                        concordant_with_agonal=(SIGN[g] * LIT[g] > 0) if g in LIT else "", control_pH_Z=ctrl.get(g, np.nan), pritzker_agonal_t=pz.get(g, np.nan)) for g in S])
write_tsv(os.path.join(HERE, "li2007_209_genes.tsv"), hd, gt)
pool_listed = np.mean([g in LIT for g in bg])
dn = int(((gt.li_listed) & (gt.li_t_mean_5regions < 0)).sum()); upn = int(((gt.li_listed) & (gt.li_t_mean_5regions > 0)).sum())
print(f"listed among 209: {int(gt.li_listed.sum())} of {len(gt)} ({gt.li_listed.mean():.3f}); listed among the decile-matched null pool: {pool_listed:.3f}; up-in-low-pH {upn}, down-in-low-pH {dn}")
disc = gt[gt.li_listed & (gt.concordant_with_agonal == False)].gene.tolist(); print("discordant listed:", disc)
