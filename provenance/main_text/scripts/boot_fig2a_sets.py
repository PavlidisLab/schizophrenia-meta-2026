"""Fig. 2a unit-level interval: observed / expected genes on >= k lists (k = 1..17) under the three analytic nulls of recurrence.py
(uniform, DE prior, expression-matched), with a cohort-cluster bootstrap (the 23 units of ps_core.units_def: 17 clusters plus 6 multi-bank
studies counted singly; a resample draws 23 units with replacement and takes all lists of each drawn unit, duplicated lists are independent
Bernoulli columns). Expected = sum over genes of the exact Poisson-binomial upper tail (recurrence.pb_matrix), RECURRENCE_NK=membership.
usage: python3 boot_fig2a.py B WORKERS -> fig2a_unit_bootstrap.tsv  (percentile 2.5 / 97.5 of observed, expected and ratio over resamples)"""
import sys, os, csv, hashlib, multiprocessing as mp
os.environ.setdefault("RECURRENCE_NK", "membership")
from ps_core import *
B = int(sys.argv[1]); W = int(sys.argv[2]); SEED = 20261013
d = load(); K = len(d.sel)
st = R.study_tier(d.sel, d.M, d.nk, d.lists, d.hitn, d.expr)
assert list(st["genes"]) == list(d.genes), "gene order differs"
NUL = {"uniform": st["P"], "de_prior": st["PP"], "expression": st["PE"]}
H = d.H
KS = list(range(1, 18))
# RA 2026-10-21: cluster-set rule. A multi-bank study is a member of the cohort cluster of each of its banks (R.cluster_set); its list sits in every one of those units.
def units_sets(d):
    m = {}
    for j, k in enumerate(d.sel):
        for c in (R.cluster_set(d.units, k) if k in d.units["multibank"] else [d.cluster[j]]):
            m.setdefault(c, []).append(j)
    return {c: np.array(v) for c, v in m.items()}
ucols = units_sets(d); uniq = sorted(ucols)

def stats(cols):
    o = H[:, cols].sum(1); out = {}
    for name, P in NUL.items():
        D = R.pb_matrix(P[:, cols]); tail = np.cumsum(D[:, ::-1], axis=1)[:, ::-1]
        out[name] = [(int((o >= k).sum()), float(tail[:, k].sum()) if k < tail.shape[1] else 0.0) for k in KS]
    return out

def job(seed):
    r = np.random.default_rng(seed)
    pick = r.integers(0, len(uniq), len(uniq)); cols = np.concatenate([ucols[uniq[i]] for i in pick])  # a spanning study drawn through several units is duplicated, one independent column per unit
    return stats(cols)

if __name__ == "__main__":
    full = stats(np.arange(K))
    seeds = [SEED * 10000 + i for i in range(B)]
    with mp.get_context("fork").Pool(W) as pool: res = pool.map(job, seeds, chunksize=5)
    rows = []
    for name in NUL:
        for ki, k in enumerate(KS):
            o0, e0 = full[name][ki]
            ob = np.array([r_[name][ki][0] for r_ in res], float); ex = np.array([r_[name][ki][1] for r_ in res])
            ok = (ob > 0) & (ex > 0); ra = ob[ok] / ex[ok]
            pc = lambda a, q: float(np.percentile(a, q)) if len(a) else float("nan")
            rows.append(dict(null=name, k=k, observed=o0, expected=round(e0, 3), ratio=round(o0 / e0, 4) if e0 > 0 else "",
                             obs_lo=pc(ob, 2.5), obs_hi=pc(ob, 97.5), exp_lo=round(pc(ex, 2.5), 3), exp_hi=round(pc(ex, 97.5), 3),
                             ratio_lo=round(pc(ra, 2.5), 4), ratio_hi=round(pc(ra, 97.5), 4), n_resamples_ratio_defined=int(ok.sum()), B=B))
    sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
    with open("fig2a_unit_bootstrap.tsv", "w", newline="") as fh:
        fh.write("# figure: fig02 panel a (per-k unit-level interval, replaces the +/-1.96 SD of null draws)\n"
                 f"# script: analyses/extmeth4_alpha_2026-10-13/boot_fig2a.py, run 2026-10-13; scenario 'published, brain', study tier, {K} lists, {len(uniq)} cohort units; RECURRENCE_NK=membership; RECURRENCE_MULTIBANK=sets\n"
                 f"# source: primary_lists.json {sha(os.path.join(A,'primary_lists.json'))}; universes.json {sha(os.path.join(A,'universes.json'))}\n"
                 f"# method: cluster bootstrap B={B} seed base {SEED}; obs_*, exp_*, ratio_* = percentile 2.5 and 97.5 over resamples; ratio only where observed and expected > 0 in the resample\n")
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    print("done", len(rows))
