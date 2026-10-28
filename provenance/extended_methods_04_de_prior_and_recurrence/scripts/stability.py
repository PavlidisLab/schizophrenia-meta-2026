#!/usr/bin/env python3
"""Stability of the DE prior: (a) bootstrap of the 242-contrast postmortem-brain prior, (b) stabilization curve over random subsets of contrasts.
Prior = de_prior_lists/our_prior/ per-contrast DE calls (cache/calls/<gemma_id>.json.gz, one contrast per experiment; DE = FDR < 0.05 and |log2FC| >= 1,
as build_group_priors.py reads them: up | down); n_de_total of a gene = number of contrasts in which it is DE. Not reimplemented: the calls are read as stored.
Gene key and universe as recurrence_priors.py (HGNC symbol by the de_prior_table rule; universe = measured by >= 1 of the 64 'published, brain' lists and in the prior).
AUROC = Mann-Whitney AUC (ties 1/2) of the 209 (confounds/recurrent_sets.tsv, study_brain) against the other universe genes, the quantity of set_auc.tsv.
usage: stability.py CALLS_DIR INPUTS_DIR OUT_DIR NWORKERS      (OMP/MKL/OPENBLAS threads = 1 set by the caller)"""
import os, sys, gzip, json, csv, hashlib, time, datetime
import numpy as np
from multiprocessing import Pool
from scipy.stats import rankdata
CALLS, INP, OUT, NW = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
SEED = 20261012; NREP = 200; NBOOT = 1000
GRID = {"all experiments": [25, 50, 100, 200, 400, 800, 1600, 3200, 8712], "postmortem brain": [25, 50, 100, 150, 200, 242]}
rd = lambda fn: list(csv.DictReader((l for l in open(fn) if not l.startswith("#")), delimiter="\t"))
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
exps = [r["gemma_id"] for r in rd(f"{INP}/experiments_all.tsv")]; have = set(exps)
def members(gf):
    ids = {r["gemma_id"] for r in rd(f"{INP}/{gf}")}; return [i for i in exps if i in ids]
def used_ids(ids):
    out = []
    for i in ids:
        fp = f"{CALLS}/{i}.json.gz"
        if os.path.exists(fp) and not json.loads(gzip.decompress(open(fp, "rb").read()))["skipped"]: out.append(i)
    return out
pools = {"all experiments": members("group_1_global.tsv"), "postmortem brain": members("group_2b_brain_postmortem_no_sz_banks.tsv")}
pools = {k: used_ids(v) for k, v in pools.items()}
print({k: len(v) for k, v in pools.items()}, flush=True)
assert len(pools["all experiments"]) == 8712 and len(pools["postmortem brain"]) == 242
n2s = {r["ncbi_gene_id"]: r["symbol"] for r in rd(f"{INP}/ncbi_to_symbol.tsv")}
U = {"all experiments": [l.strip() for l in open(f"{INP}/universe_global.txt") if l.strip()], "postmortem brain": [l.strip() for l in open(f"{INP}/universe_2b.txt") if l.strip()]}
S209 = {l.strip() for l in open(f"{INP}/set209.txt") if l.strip()}
usym = set(U["all experiments"]) | set(U["postmortem brain"])
cand = sorted(n for n, s in n2s.items() if s in usym); cidx = {n: j for j, n in enumerate(cand)}
allids = pools["all experiments"]; assert set(pools["postmortem brain"]) <= set(allids)

def load(i):
    g = json.loads(gzip.decompress(open(f"{CALLS}/{i}.json.gz", "rb").read()))["genes"]
    ix, de = [], []
    for n, v in g.items():
        j = cidx.get(n)
        if j is not None: ix.append(j); de.append(1 if (v[1] or v[2]) else 0)
    return np.array(ix, np.int32), np.array(de, np.uint8)
t0 = time.time()
with Pool(NW) as p: res = p.map(load, allids, chunksize=50)
T = np.zeros((len(allids), len(cand)), np.uint8); D = np.zeros_like(T)
for r, (ix, de) in enumerate(res): T[r, ix] = 1; D[r, ix] = de
del res; print("matrix", T.shape, f"{time.time()-t0:.0f}s", flush=True)
row = {i: r for r, i in enumerate(allids)}

# per pool: pick one NCBI column per symbol (largest full count, as de_prior_table keeps the largest), restrict to the pool's universe
P = {}
for name in pools:
    rows = np.array([row[i] for i in pools[name]]); full = D[rows].sum(0); tst = T[rows].sum(0)
    best = {}
    for j, n in enumerate(cand):
        s = n2s[n]
        if tst[j] == 0 or (s in best and full[best[s]] >= full[j]): continue  # columns untested in the pool are not rows of its prior file
        best[s] = j
    genes = U[name]; cols = np.array([best[g] for g in genes])
    Dp, Tp = np.ascontiguousarray(D[rows][:, cols]), np.ascontiguousarray(T[rows][:, cols])
    P[name] = dict(genes=genes, D=Dp, T=Tp, in209=np.array([g in S209 for g in genes]), n=len(rows))
    assert Tp.sum(0).min() >= 1
del T, D
# verification against the shipped prior files (n_de_total, n_tested)
for name, f in (("all experiments", "1_global_our_prior.tsv"), ("postmortem brain", "2b_brain_postmortem_no_sz_banks_our_prior.tsv")):
    sh = {}
    for r in rd(os.path.join(INP, "..", "..", "de_prior_lists", "our_prior", f)) if os.path.exists(os.path.join(INP, "..", "..", "de_prior_lists", "our_prior", f)) else rd(os.path.join(INP, f)):
        g = n2s.get(r["ncbi_gene_id"])
        if g and (g not in sh or int(r["n_de_total"]) > sh[g]): sh[g] = int(r["n_de_total"])
    mine = P[name]["D"].sum(0); ok = sum(int(sh[g]) == int(m) for g, m in zip(P[name]["genes"], mine))
    print("verify", name, ok, "of", len(mine), "genes equal the shipped n_de_total", flush=True)
    assert ok == len(mine)

def rho(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 3: return np.nan
    ra, rb = rankdata(a[m]), rankdata(b[m]); ra -= ra.mean(); rb -= rb.mean()
    d = np.sqrt((ra ** 2).sum() * (rb ** 2).sum()); return float((ra * rb).sum() / d) if d else np.nan
def auroc(x, pos):
    m = np.isfinite(x); x, pos = x[m], pos[m]; n1 = int(pos.sum()); n0 = len(x) - n1
    if n1 == 0 or n0 == 0: return np.nan, n1
    r = rankdata(x); return float((r[pos].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)), n1
def stats(c, t):  # count vector (float, nan where untested), frequency vector
    cnt = np.where(t > 0, c, np.nan); fr = np.where(t > 0, c / np.maximum(t, 1), np.nan); return cnt, fr

def subset_task(a):
    name, n = a; Pn = P[name]; N = Pn["n"]; Dp, Tp, in209 = Pn["D"], Pn["T"], Pn["in209"]
    cf, ff = stats(Dp.sum(0).astype(float), Tp.sum(0).astype(float))
    rng = np.random.default_rng([SEED, 1 if name == "all experiments" else 2, n]); out = []
    reps = 1 if n >= N else NREP
    for k in range(reps):
        perm = rng.permutation(N); A = perm[:n]
        cA, fA = stats(Dp[A].sum(0, dtype=np.int32).astype(float), Tp[A].sum(0, dtype=np.int32).astype(float))
        r = dict(pool=name, n=n, rep=k, rho_full_count=rho(cA, cf), rho_full_freq=rho(fA, ff), rho_disjoint_count=np.nan, rho_disjoint_freq=np.nan)
        if 2 * n <= N:
            B = perm[n:2 * n]; cB, fB = stats(Dp[B].sum(0, dtype=np.int32).astype(float), Tp[B].sum(0, dtype=np.int32).astype(float))
            r["rho_disjoint_count"], r["rho_disjoint_freq"] = rho(cA, cB), rho(fA, fB)
        r["auroc_count"], r["n209_in_universe"] = auroc(cA, in209); r["auroc_freq"], _ = auroc(fA, in209)
        r["n_universe_tested"] = int(np.isfinite(cA).sum()); out.append(r)
    print("done", name, n, flush=True); return out

def boot():
    Pn = P["postmortem brain"]; Dp, Tp, in209, genes = Pn["D"].astype(np.float64), Pn["T"].astype(np.float64), Pn["in209"], Pn["genes"]; N = Pn["n"]
    rng = np.random.default_rng([SEED, 3]); W = np.stack([np.bincount(rng.integers(0, N, N), minlength=N) for _ in range(NBOOT)]).astype(np.float64)
    C, Tt = W @ Dp, W @ Tp  # 1000 x genes
    cf, ff = stats(Dp.sum(0), Tp.sum(0))
    cnt = np.where(Tt > 0, C, np.nan); fr = np.where(Tt > 0, C / np.maximum(Tt, 1), np.nan)
    perm = rng.permutation(NBOOT); perm = np.where(perm == np.arange(NBOOT), np.roll(perm, 1), perm)  # partner draw, never itself
    rows = []
    for b in range(NBOOT):
        a1, _ = auroc(cnt[b], in209); a2, _ = auroc(fr[b], in209)
        rows.append(dict(draw=b, rho_full_count=rho(cnt[b], cf), rho_full_freq=rho(fr[b], ff), rho_pair_count=rho(cnt[b], cnt[perm[b]]), rho_pair_freq=rho(fr[b], fr[perm[b]]),
                         auroc_count=a1, auroc_freq=a2, n_contrasts_distinct=int((W[b] > 0).sum())))
    sd_f = np.nanstd(fr, axis=0, ddof=1); sd_c = np.nanstd(cnt, axis=0, ddof=1); nan_draws = np.isnan(fr).sum(0)
    lo, hi = np.nanpercentile(fr, [2.5, 97.5], axis=0)
    gene_rows = [dict(symbol=g, in_209=int(in209[i]), n_tested=int(Tp.sum(0)[i]), n_de=int(Dp.sum(0)[i]), freq=ff[i], boot_sd_freq=sd_f[i], boot_sd_count=sd_c[i],
                      freq_lo=lo[i], freq_hi=hi[i], draws_untested=int(nan_draws[i])) for i, g in enumerate(genes)]
    return rows, gene_rows

if __name__ == "__main__":
    t0 = time.time(); jobs = [(nm, n) for nm in GRID for n in GRID[nm][::-1]]
    with Pool(NW) as p: sub = p.map(subset_task, jobs, chunksize=1)
    sub = [r for lst in sub for r in lst]
    brows, grows = boot()
    stamp = datetime.date.today().isoformat()
    ins = {f: sha(os.path.join(INP, f)) for f in sorted(os.listdir(INP))}
    mani = hashlib.sha256("\n".join(f"{i} {sha(f'{CALLS}/{i}.json.gz')}" for i in allids).encode()).hexdigest()
    head = [f"# script: analyses/de_prior_stability_2026-10-12/stability.py, run {stamp} on remote; seed {SEED}; {NREP} subsets per n, {NBOOT} bootstrap draws",
            "# source: per-contrast DE calls " + CALLS + f" ({len(allids)} used contrasts of 8,762 experiments; manifest sha256 of the call files {mani})",
            "# inputs sha256: " + "; ".join(f"{k} {v[:12]}" for k, v in ins.items())]
    def wr(fn, rows):
        keys = list(rows[0]); 
        with open(f"{OUT}/{fn}", "w") as fh:
            fh.write("\n".join(head) + "\n" + "\t".join(keys) + "\n")
            for r in rows: fh.write("\t".join(("" if (isinstance(r[k], float) and np.isnan(r[k])) else (f"{r[k]:.6g}" if isinstance(r[k], float) else str(r[k]))) for k in keys) + "\n")
    wr("subset_runs.tsv", sub); wr("bootstrap_draws.tsv", brows); wr("bootstrap_genes.tsv", grows)
    full = {nm: auroc(stats(P[nm]["D"].sum(0).astype(float), P[nm]["T"].sum(0).astype(float))[0], P[nm]["in209"]) for nm in P}
    json.dump(dict(full_auroc_count={k: v[0] for k, v in full.items()}, n209={k: v[1] for k, v in full.items()}, n_universe={k: len(P[k]["genes"]) for k in P}, n_contrasts={k: P[k]["n"] for k in P}, seconds=time.time() - t0), open(f"{OUT}/full_prior_check.json", "w"), indent=1)
    print("full AUROC", full, f"{time.time()-t0:.0f}s", flush=True)
