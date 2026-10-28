#!/usr/bin/env python3
"""NEW 2026-10-18 (composition_sim). Step 2: composition-only simulation of the bulk pH and agonal axes.
No cell-intrinsic pH effect anywhere (except the labelled positive-control configs). Run:
  OMP_NUM_THREADS=1 nice -n 19 python sim.py <run_dir> [NREP] [NPROC]      seed 20261018 (SeedSequence([seed, config index, rep]))
Model: CLR of NSW control class fractions (7 classes) ~ pH + age + male (OLS). Per synthetic sample: pH drawn from the (bootstrapped)
NSW control pH values, CLR = intercept(mean covariates) + scale*slope*(pH - mean pH) + residual (multivariate normal, residual covariance),
softmax, x (1 - Mural fraction) with Mural (absent from NSW) fixed at the PsychAD control mean fraction of nuclei; RNA share
w_k ~ p_k * UMI-per-nucleus_k (or p_k without RNA weighting); expected CPM = sum_k w_k * class mean CPM; counts ~ NB(mean = CPM*depth/1e6, phi) or Poisson;
log2((count+0.5)/(sum+1)*1e6); per-gene OLS t on pH; Spearman with the axes."""
import sys, os, json, time, numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import rankdata
RUN = sys.argv[1]; NREP = int(sys.argv[2]) if len(sys.argv) > 2 else 200; NPROC = int(sys.argv[3]) if len(sys.argv) > 3 else 8
SEED = 20261018; DEPTH = 25e6; PHI = 0.1
OUT = os.path.join(RUN, "out")
P = np.load(os.path.join(OUT, "prep.npz"), allow_pickle=True)
CLS = list(P["classes"]); K7 = [c for c in CLS if c != "Mural"]; IM = CLS.index("Mural")
I7 = [CLS.index(c) for c in K7]
nsw = pd.read_csv(os.path.join(OUT, "nsw_control_composition.tsv"), sep="\t")
CNT = nsw[K7].values.astype(float); PH = nsw.ph.values; AGE = nsw.age.values; MALE = nsw.male.values.astype(float)
ND = len(nsw)
def clr(c):
    p = (c + 0.5) / (c + 0.5).sum(1, keepdims=True); l = np.log(p); return l - l.mean(1, keepdims=True)
CLR = clr(CNT); PFRAC = CNT / CNT.sum(1, keepdims=True)

def fit_comp(ix):
    """OLS CLR ~ pH + age + male on donors ix. Returns a (intercept at covariate means), b (pH slope), cov (residual covariance), pH values."""
    X = np.column_stack([np.ones(len(ix)), PH[ix], AGE[ix], MALE[ix]])
    if np.linalg.matrix_rank(X) < 4: return None
    B, *_ = np.linalg.lstsq(X, CLR[ix], rcond=None)
    R = CLR[ix] - X @ B
    cov = R.T @ R / (len(ix) - 4)
    a = np.array([1, PH[ix].mean(), AGE[ix].mean(), MALE[ix].mean()]) @ B
    # a is the fitted CLR at mean pH and mean covariates; slope sign: per pH unit
    return dict(a=a, b=B[1], cov=cov, ph=PH[ix], phbar=PH[ix].mean(), sd=np.sqrt(np.diag(cov)))
FULL = fit_comp(np.arange(ND))

def softmax(z):
    z = z - z.max(1, keepdims=True); e = np.exp(z); return e / e.sum(1, keepdims=True)

def props(rng, fit, ph, scale, mode, k_fixed=None, b_override=None):
    n = len(ph); b = fit["b"] if b_override is None else b_override
    w, V = np.linalg.eigh(fit["cov"]); w = np.clip(w, 0, None)
    eps = (rng.standard_normal((n, len(w))) * np.sqrt(w)) @ V.T
    z = fit["a"] + scale * np.outer(ph - fit["phbar"], b) + eps
    if k_fixed is None: return softmax(z)
    f = softmax(fit["a"][None])[0, k_fixed]; o = [j for j in range(7) if j != k_fixed]
    p = np.zeros((n, 7)); p[:, o] = (1 - f) * softmax(z[:, o]); p[:, k_fixed] = f; return p

def to8(p7, m):
    p = np.zeros((p7.shape[0], 8)); p[:, I7] = (1 - m) * p7; p[:, IM] = m; return p

class Prof:
    def __init__(self, name):
        self.name = name; self.cpm = P["cpm_" + name]; self.rna = P["rna_" + name]; self.m = float(P["mural_" + name])
        mix0 = (self.cpm * (to8(PFRAC.mean(0)[None], self.m)[0] * self.rna)[:, None]).sum(0) / (to8(PFRAC.mean(0)[None], self.m)[0] * self.rna).sum()
        self.keep = mix0 >= 1.0
        self.cpm = self.cpm[:, self.keep]; self.genes = P["genes"][self.keep]
        sel = lambda x: x[self.keep]
        self.pH = sel(P["ax_pH"]); self.ag = sel(P["ax_ag"]); self.i383 = sel(P["in383"]); self.i107 = sel(P["in107"])
        self.sets = {}
        for axn, ax in [("pH", self.pH), ("agonal", self.ag)]:
            for sn, msk in [("all", np.ones(len(ax), bool)), ("rec383", self.i383), ("rec107", self.i107)]:
                ok = np.flatnonzero(np.isfinite(ax) & msk)
                self.sets[(axn, sn)] = (ok, rankdata(ax[ok]))
        rng = np.random.default_rng(SEED + 7); ng = len(self.genes)
        self.pc_idx = np.sort(rng.choice(ng, int(round(0.05 * ng)), replace=False))
        self.pc_sign = np.where(rng.random(len(self.pc_idx)) < 0.5, -1.0, 1.0)
PROF = {n: Prof(n) for n in ["HBCC", "MSSM"]}

def mix(prof, p8, rna_w):
    w = p8 * prof.rna[None, :] if rna_w else p8.copy()
    w = w / w.sum(1, keepdims=True)
    return w @ prof.cpm      # expected CPM, N x G (rows sum to ~1e6 over all genes; subset after filter)

def sample_counts(rng, cpm, noise):
    lam = cpm / 1e6 * DEPTH
    if noise == "nb": lam = lam * rng.gamma(1 / PHI, PHI, size=lam.shape)
    return rng.poisson(lam).astype(np.float64)

def tstat(Y, x):
    n = len(x); xc = x - x.mean(); Yc = Y - Y.mean(0)
    r = (xc @ Yc) / (np.sqrt((xc ** 2).sum()) * np.sqrt((Yc ** 2).sum(0)) + 1e-300)
    r = np.clip(r, -0.999999, 0.999999); return r * np.sqrt((n - 2) / (1 - r ** 2)), (xc @ Yc) / (xc ** 2).sum()

def spear(prof, t, dec=False, mu=None):
    out = {}
    for (axn, sn), (ok, rk) in prof.sets.items():
        out[f"rho_{axn}_{sn}"] = np.corrcoef(rankdata(t[ok]), rk)[0, 1]
        out[f"n_{axn}_{sn}"] = len(ok)
        if dec and sn == "all":
            q = np.quantile(mu[ok], np.linspace(0, 1, 11)); d = np.clip(np.searchsorted(q, mu[ok], side="right") - 1, 0, 9)
            for j in range(10):
                s = d == j; out[f"rho_{axn}_decile{j+1:02d}"] = np.corrcoef(rankdata(t[ok][s]), rankdata(prof.__dict__["pH" if axn == "pH" else "ag"][ok][s]))[0, 1]
    return out

def run_one(args):
    cid, cfg, rep = args
    rng = np.random.default_rng(np.random.SeedSequence([SEED, cid, rep]))
    prof = PROF[cfg["profile"]]; var = cfg["variant"]; N = cfg["N"]; scale = cfg["scale"]
    # donor bootstrap (slope uncertainty) unless fixed
    for _ in range(50):
        ix = rng.integers(0, ND, ND) if cfg["boot"] else np.arange(ND)
        fit = fit_comp(ix)
        if fit is not None: break
    if var == "observed":
        p7 = PFRAC[ix]; ph = PH[ix]
    else:
        ph = rng.choice(fit["ph"], N, replace=True)
        b = fit["b"]
        if var == "null_shuffle": b = b[rng.permutation(7)]
        kf = K7.index(cfg["loco"]) if var == "loco" else None
        p7 = props(rng, fit, ph, scale, var, k_fixed=kf, b_override=b)
    cpm = mix(prof, to8(p7, prof.m), cfg["rna"])
    if var == "positive_control":
        beta = np.zeros(cpm.shape[1]); beta[prof.pc_idx] = cfg["pc_beta"] * prof.pc_sign
        cpm = cpm * 2.0 ** (np.outer(ph - ph.mean(), beta))
    C = sample_counts(rng, cpm, cfg["noise"])
    Y = np.log2((C + 0.5) / (C.sum(1, keepdims=True) + 1) * 1e6)
    xph = rng.permutation(ph) if var == "null_pH_perm" else ph
    t, slope = tstat(Y, xph)
    res = spear(prof, t, dec=cfg["dec"], mu=Y.mean(0))
    res["sd_pH"] = float(np.std(ph)); res["n_samples"] = len(ph)
    # composition actually generated: Spearman-free summary, slope of each class fraction (p8) on pH
    pf = to8(p7, prof.m)
    for k, c in enumerate(CLS): res[f"fracslope_{c}"] = np.polyfit(ph, pf[:, k], 1)[0] if np.std(ph) > 0 else np.nan
    if var == "positive_control":
        a = np.abs(t); inj = np.zeros(len(t), bool); inj[prof.pc_idx] = True
        rk = rankdata(a); res["pc_auc_abs_t"] = (rk[inj].sum() - inj.sum() * (inj.sum() + 1) / 2) / (inj.sum() * (~inj).sum())
        res["pc_rho_t_vs_beta"] = np.corrcoef(rankdata(t), rankdata(beta))[0, 1]
        res["pc_sign_agree"] = float(np.mean(np.sign(t[prof.pc_idx]) == prof.pc_sign))
        res["pc_median_abs_t_injected"] = float(np.median(a[inj])); res["pc_median_abs_t_other"] = float(np.median(a[~inj]))
    return cid, rep, res

def configs():
    C = []
    def add(**k):
        d = dict(profile="HBCC", variant="model", N=200, rna=1, scale=1.0, noise="nb", boot=True, loco="", pc_beta=0.0, dec=False, tag="main_grid"); d.update(k); C.append(d)
    for pr in PROF:
        for N in [32, 100, 200, 500]:
            for rna in [1, 0]:
                for sc in [0.0, 0.5, 1.0, 1.5]:
                    add(profile=pr, N=N, rna=rna, scale=sc, dec=(N == 200 and sc == 1.0), tag="main_grid")
        for rna in [1, 0]: add(profile=pr, variant="observed", N=32, rna=rna, tag="observed_proportions")
        for N in [32, 100, 200, 500]:
            add(profile=pr, variant="null_shuffle", N=N, tag="null_slopes_shuffled")
            add(profile=pr, variant="null_pH_perm", N=N, tag="null_pH_permuted")
        add(profile=pr, noise="poisson", tag="noise_poisson")
        add(profile=pr, boot=False, dec=True, tag="slopes_fixed_no_bootstrap")
        add(profile=pr, boot=False, rna=0, tag="slopes_fixed_no_bootstrap")
        for k in K7: add(profile=pr, variant="loco", loco=k, tag="leave_one_class_out")
        for beta in [0.25, 0.5]:
            for sc in [0.0, 1.0]: add(profile=pr, variant="positive_control", pc_beta=beta, scale=sc, tag="positive_control")
    return C
CFG = configs()

def expected_rows():
    """Noise-free expectation: log2 mixture CPM slope on pH over the 32 observed control pH values (no residual noise, fitted slopes)."""
    rows = []
    for pr, prof in PROF.items():
        for rna in [1, 0]:
            for sc in [0.0, 0.5, 1.0, 1.5]:
                for var in ["model", "observed"]:
                    if var == "observed" and sc != 1.0: continue
                    if var == "model":
                        z = FULL["a"] + sc * np.outer(PH - FULL["phbar"], FULL["b"]); p7 = softmax(z)
                    else: p7 = PFRAC
                    Y = np.log2(mix(prof, to8(p7, prof.m), rna))
                    x = PH - PH.mean(); slope = (x @ (Y - Y.mean(0))) / (x @ x)
                    if np.allclose(slope, 0): slope = slope + 1e-12 * np.random.default_rng(1).standard_normal(len(slope))
                    r = spear(prof, slope, dec=True, mu=Y.mean(0)); r.update(profile=pr, rna=rna, scale=sc, variant="expected_" + var); rows.append(r)
    return rows

if __name__ == "__main__":
    t0 = time.time()
    pd.DataFrame([dict(config_id=i, **c) for i, c in enumerate(CFG)]).to_csv(os.path.join(OUT, "sim_configs.tsv"), sep="\t", index=False)
    pd.DataFrame(expected_rows()).to_csv(os.path.join(OUT, "sim_expected_noise_free.tsv"), sep="\t", index=False)
    fit_rows = [dict(class_nsw=c, class_psychad=c, clr_slope_per_pH=FULL["b"][j], resid_sd=FULL["sd"][j], mean_fraction=PFRAC[:, j].mean(), mean_clr=FULL["a"][j]) for j, c in enumerate(K7)]
    pd.DataFrame(fit_rows).to_csv(os.path.join(OUT, "nsw_class_slopes_point.tsv"), sep="\t", index=False)
    tasks = [(i, c, r) for i, c in enumerate(CFG) for r in range(NREP)]
    print("configs", len(CFG), "tasks", len(tasks), "genes", {n: len(p.genes) for n, p in PROF.items()}, flush=True)
    rows = []
    with Pool(NPROC) as pool:
        for n, (cid, rep, res) in enumerate(pool.imap_unordered(run_one, tasks, chunksize=4)):
            for k, v in res.items(): rows.append((cid, rep, k, v))
            if n % 1000 == 0: print(n, len(tasks), f"{time.time()-t0:.0f}s", flush=True)
    pd.DataFrame(rows, columns=["config_id", "rep", "stat", "value"]).to_csv(os.path.join(OUT, "sim_reps_long.tsv.gz"), sep="\t", index=False)
    json.dump({n: dict(n_genes=len(p.genes), sets={f"{a}|{s}": len(v[0]) for (a, s), v in p.sets.items()}, n_positive_control_genes=len(p.pc_idx)) for n, p in PROF.items()}, open(os.path.join(OUT, "sim_gene_counts.json"), "w"), indent=1)
    print("done", f"{time.time()-t0:.0f}s", flush=True)
