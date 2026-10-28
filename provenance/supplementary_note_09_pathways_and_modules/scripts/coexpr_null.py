#!/usr/bin/env python3
"""Supplementary Note 5 Fig. 1a (pathway and gene sharing between lists) against nulls that keep co-expression.

Reuses fairness_lists/a5b_go.py unchanged for: the list set (rule A primary lists, scenario 'published, brain'), the GO BP
collection (terms of 5-200 genes), the DE-prior weights (de_prior_weights + 1, median for genes without a prior), top pathways
(up to 20 terms, hypergeometric p < 0.01, >= 2 hits), the pair statistics (mean Jaccard, shared totals) and the independent-gene
nulls (uniform, DE-prior-weighted). The only change is the added null: each list's hits (within its annotated tested set U_k) are
moved by one circular shift along the GTEx v8 cortex co-expression ordering, within strata
(rotation_lib.py of gene_correlation_adjusted_tests/axes_and_recurrence_dependence_tests/recurrence_overlap_prior, used unchanged).
  S0  one stratum (analogue of the uniform null)
  S3  40 equal-count DE-prior bins over the annotated genes (analogue of the DE-prior null)
  R1  every list shifted on its own;  R2  lists of one cohort cluster shifted together (shared donors kept)
usage: python3 coexpr_null.py MODE NREP WORKERS [SEED] [LISTSET]
  MODE: obs | nulls ;  LISTSET: current | stored62 (lists with included=yes in fairness_lists/a5b_lists.tsv)"""
import sys, os, itertools, json, time, hashlib, multiprocessing
import numpy as np, pandas as pd
from scipy.sparse import csr_matrix
HERE = os.path.dirname(os.path.abspath(__file__))
FL = os.path.abspath(os.path.join(HERE, "..", "fairness_lists"))
ROT = os.path.abspath(os.path.join(HERE, "..", "gene_correlation_adjusted_tests", "axes_and_recurrence_dependence_tests", "recurrence_overlap_prior"))
sys.path.insert(0, FL); sys.path.insert(0, ROT)
try:
    _cwd = os.getcwd(); os.chdir(FL)
    import a5b_go as A
    from a5b_go import C
    R = C.R
    HAVE_REPO = True
except Exception:  # remote runner: no repo, inputs come from the pickle written by `dump`; enrich/top_terms are verbatim copies of a5b_go.py
    HAVE_REPO = False
    from scipy.stats import hypergeom
    class A:
        MIN_HITS = 5; TOP = 20; P_TOP = 0.01
        @staticmethod
        def enrich(S, u, h):
            K = S @ u; k = S @ h; N, n = int(u.sum()), int(h.sum()); p = np.ones(S.shape[0]); m = (K >= 5) & (k >= 2)
            p[m] = hypergeom.sf(k[m] - 1, N, K[m], n); return p, k
        @staticmethod
        def top_terms(p, k):
            idx = np.where((p < A.P_TOP) & (k >= 2))[0]; idx = idx[np.lexsort((-k[idx], p[idx]))][:A.TOP]; return frozenset(idx.tolist())
import rotation_lib as RL
KS = (2, 3, 5, 10)


def build(listset):
    ctx = C.Ctx(); terms, names, members = A.load_library()
    annotated = sorted(set().union(*map(set, members))); gi = {g: i for i, g in enumerate(annotated)}
    rows_, cols_ = [], []
    for ti, m in enumerate(members):
        for g in m: rows_.append(ti); cols_.append(gi[g])
    S = csr_matrix((np.ones(len(rows_)), (rows_, cols_)), shape=(len(terms), len(annotated)))
    size = np.asarray(S.sum(1)).ravel(); keep = (size >= 5) & (size <= 200)
    S = S[np.where(keep)[0]]; terms = [t for t, k in zip(terms, keep) if k]
    w = np.array([ctx.hitn.get(g, np.nan) + 1.0 for g in annotated]); w = np.where(np.isnan(w), np.nanmedian(w), w)
    base = ctx.reps(C.MAIN)
    if listset == "stored62":
        old = pd.read_csv(os.path.join(FL, "a5b_lists.tsv"), sep="\t", comment="#", dtype={"study": str})
        base = [k for k in base if k in set(old.loc[old.included == "yes", "study"])]
    U, H = {}, {}
    for k in base:
        u = np.zeros(len(annotated), bool); u[[gi[g] for g in ctx.M[k] if g in gi]] = True
        h = np.zeros(len(annotated), bool); h[[gi[g] for g in ctx.lists[k] if g in gi and g in ctx.M[k]]] = True
        h &= u
        if h.sum() >= A.MIN_HITS: U[k], H[k] = u, h
    inc = [k for k in base if k in U]
    unit = {}
    for k in inc:
        lab = R.cohort_label(ctx.units, k)
        unit[k] = f"multibank:{k}" if (lab.startswith("per-bank") or lab.startswith("left out")) else lab
    pos, ingtex = RL.gene_positions(annotated)
    B = dict(ctx=ctx, S=S, terms=terms, names=names, annotated=annotated, w=w, U=U, H=H, inc=inc, unit=unit, pos=pos, ingtex=ingtex,
             pairs=list(itertools.combinations(range(len(inc)), 2)))
    B["uf"] = {k: U[k].astype(float) for k in inc}
    B["Kt"] = {k: B["S"] @ B["uf"][k] for k in inc}
    return B


def stats(B, hits):
    """hits: dict k -> bool vector. Returns per-pair arrays and per-gene list counts."""
    inc, S = B["inc"], B["S"]; tops = {}
    for k in inc:
        p, kk = A.enrich(S, B["uf"][k], hits[k].astype(float)); tops[k] = A.top_terms(p, kk)
    P = B["pairs"]; n = len(P)
    jt = np.full(n, np.nan); st = np.zeros(n); jg = np.full(n, np.nan); sg = np.zeros(n)
    for i, (a, b) in enumerate(P):
        ka, kb = inc[a], inc[b]; TA, TB = tops[ka], tops[kb]
        if TA or TB: jt[i] = len(TA & TB) / len(TA | TB)
        st[i] = len(TA & TB)
        uu = B["U"][ka] & B["U"][kb]; ha, hb = hits[ka] & uu, hits[kb] & uu
        un = (ha | hb).sum(); ov = (ha & hb).sum(); sg[i] = ov
        if un: jg[i] = ov / un
    cnt = np.zeros(len(B["annotated"]), int)
    for k in inc: cnt += hits[k]
    ge = np.array([(cnt >= q).sum() for q in KS], float)
    return jt, st, jg, sg, ge


def summarize(jt, st, jg, sg, mask=None):
    m = np.ones(len(jt), bool) if mask is None else mask
    return np.array([np.nanmean(jt[m]), st[m].sum(), np.nanmean(jg[m]), sg[m].sum()])

METRICS = ["mean_pathway_jaccard", "shared_pathways_total", "mean_gene_jaccard", "shared_genes_total"]

_B = None; _CFG = None
def _rep(seed_and_cfg):
    seed, cfg = seed_and_cfg; B = _B; rng = np.random.default_rng(seed); inc = B["inc"]
    if cfg[0] == "indep":
        hits = {}
        for k in inc:
            idx = np.where(B["U"][k])[0]; n = int(B["H"][k].sum())
            if cfg[1] == "uniform": pr = rng.choice(idx, n, replace=False)
            else: pw = B["w"][idx] / B["w"][idx].sum(); pr = rng.choice(idx, n, replace=False, p=pw)
            h = np.zeros(len(B["annotated"]), bool); h[pr] = True; hits[k] = h
        size = 1.0
    elif cfg[0] == "strat":  # control: same strata as rot_*_S3, hits drawn independently (no co-expression) within each stratum
        st = B["strat3"]; hits = {}
        for k in inc:
            h = np.zeros(len(B["annotated"]), bool); Uk = B["U"][k]; Hk = B["H"][k]
            for b in np.unique(st[Hk]):
                idx = np.where(Uk & (st == b))[0]; n = int((Hk & (st == b)).sum()); h[rng.choice(idx, n, replace=False)] = True
            hits[k] = h
        size = 1.0
    else:
        rot = B["rot"][cfg[1:]]; out = rot.draw(rng); hits = {}
        for j, k in enumerate(inc):
            h = np.zeros(len(B["annotated"]), bool); h[out[j]] = True; hits[k] = h
        size = np.mean([hits[k].sum() / B["H"][k].sum() for k in inc])
    jt, st, jg, sg, ge = stats(B, hits)
    return jt.astype(np.float32), st.astype(np.float32), jg.astype(np.float32), sg.astype(np.float32), ge, size


def main():
    mode, nrep, workers = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    seed0 = int(sys.argv[4]) if len(sys.argv) > 4 else 20261013
    listset = sys.argv[5] if len(sys.argv) > 5 else "current"
    out = os.path.join(HERE, "out_" + listset); os.makedirs(out, exist_ok=True)
    global _B
    import pickle
    if mode == "dump":
        B = build(listset); ctx = B.pop("ctx"); B["labels"] = [ctx.label(k) for k in B["inc"]]
        pickle.dump(B, open(os.path.join(HERE, f"base_{listset}.pkl"), "wb")); print("dumped"); return
    if not HAVE_REPO:
        B = pickle.load(open(os.path.join(HERE, f"base_{listset}.pkl"), "rb")); B["ctx"] = type("X", (), {"label": lambda self, k: B["labels"][B["inc"].index(k)]})()
    else:
        B = build(listset)
    inc = B["inc"]; print("lists", len(inc), "pairs", len(B["pairs"]), flush=True)
    hits0 = B["H"]; obs = stats(B, hits0)
    json.dump(dict(inc=inc, labels=[B["ctx"].label(k) for k in inc], unit=[B["unit"][k] for k in inc], observed=summarize(*obs[:4]).tolist(),
                   genes_ge_k=dict(zip(KS, obs[4].tolist())), n_pairs=len(B["pairs"])), open(os.path.join(out, "observed.json"), "w"), indent=1)
    np.savez_compressed(os.path.join(out, "observed_pairs.npz"), jt=obs[0], st=obs[1], jg=obs[2], sg=obs[3], ge=obs[4])
    print("observed", summarize(*obs[:4]), obs[4], flush=True)
    if mode == "obs": return
    # strata and rotators
    prior = B["w"]; ebins = np.zeros(len(B["annotated"]), int)
    meas = np.stack([B["U"][k] for k in inc], 1); hl = [np.where(B["H"][k])[0] for k in inc]
    ukeys = [B["unit"][k] for k in inc]
    B["rot"] = {}
    for sk in ("S0", "S3"):
        strat = RL.strata_labels(sk, ebins, prior)
        B["rot"][("R1", sk)] = RL.Rotator(meas, hl, list(range(len(inc))), strat, B["pos"])
        B["rot"][("R2", sk)] = RL.Rotator(meas, hl, ukeys, strat, B["pos"])
    _B = B
    B["strat3"] = RL.strata_labels("S3", ebins, prior)
    cfgs = [("indep", "uniform"), ("indep", "deprior")] + [("rot", m, s) for m in ("R1", "R2") for s in ("S0", "S3")] + [("strat", "S3")]
    if os.environ.get("CFGS") == "strat": cfgs = [("strat", "S3")]
    for ci, cfg in enumerate(cfgs):
        t0 = time.time(); name = "_".join(cfg[1:])
        seeds = [(seed0 + 1000003 * (ci + (6 if os.environ.get('CFGS') else 0)) + r, cfg) for r in range(nrep)]
        with multiprocessing.get_context("fork").Pool(workers) as pool: res = pool.map(_rep, seeds, chunksize=4)
        np.savez_compressed(os.path.join(out, f"reps_{name}.npz"), jt=np.stack([r[0] for r in res]), st=np.stack([r[1] for r in res]),
                            jg=np.stack([r[2] for r in res]), sg=np.stack([r[3] for r in res]), ge=np.stack([r[4] for r in res]), size=np.array([r[5] for r in res]))
        print(name, "done", round(time.time() - t0), "s", flush=True)
    if not os.environ.get("CFGS"): json.dump(dict(seed0=seed0, nrep=nrep, cfgs=[list(c) for c in cfgs]), open(os.path.join(out, "run.json"), "w"))

if __name__ == "__main__": main()
