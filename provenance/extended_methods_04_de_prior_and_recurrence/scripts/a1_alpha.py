"""Step 2-3: expected and observed genes in >=5 and >=10 lists across the DE-prior strength alpha (w = (count+1)^alpha),
profile likelihood, bootstrap intervals over lists and over cohort clusters. Writes alpha_grid.tsv, profile_ll.tsv, alpha_hat.json, ll_list.npz."""
import time, json, multiprocessing as mp
from ps_core import *
d = load(); K = len(d.sel)
GRID = (0, 0.1, 0.25, 0.5, 0.75, 1, 2, 5)
FINE = np.round(np.arange(0, 2.0001, 0.025), 4)
B = 200
rng0 = 20261010

def alpha_hat(llsum):
    i = int(np.argmax(llsum)); 
    if 0 < i < len(FINE) - 1:   # parabolic refinement
        y0, y1, y2 = llsum[i - 1], llsum[i], llsum[i + 1]; den = y0 - 2 * y1 + y2
        return float(FINE[i] + 0.5 * (y0 - y2) / den * (FINE[1] - FINE[0])) if den < 0 else float(FINE[i])
    return float(FINE[i])

# ---- exact grid, both conventions
out = []
Pg = {}
for conv in ("project", "exact"):
    for miss in ("median", "zero"):
        for a in GRID:
            P = probs(d, a, miss, exact_total=(conv == "exact")); e = exp_ge(P, (5, 10)); ll = list_ll(d, P)
            out.append(dict(total=conv, missing_fill=miss, alpha=a, obs5=int((d.obs >= 5).sum()), obs10=int((d.obs >= 10).sum()),
                            exp5=e[5], exp10=e[10], ratio5=(d.obs >= 5).sum() / e[5], ratio10=(d.obs >= 10).sum() / e[10], loglik=ll.sum()))
            if conv == "project" and miss == "median": Pg[a] = P
print("grid done", flush=True)
# ---- fine profile list-level ll
LL = {}
for conv in ("project", "exact"):
    for miss in ("median", "zero"):
        LL[(conv, miss)] = np.array([list_ll(d, probs(d, a, miss, exact_total=(conv == "exact"))) for a in FINE])  # (alpha, K)
np.savez("ll_list.npz", fine=FINE, **{f"{c}_{m}": v for (c, m), v in LL.items()})
prof = []; hat = {}
for (conv, miss), M_ in LL.items():
    s = M_.sum(1); ah = alpha_hat(s); mx = s.max(); inside = FINE[s >= mx - 1.92]
    hat[f"{conv}_{miss}"] = dict(alpha_hat=ah, grid_max=float(FINE[np.argmax(s)]), loglik_max=float(mx), lr_lo=float(inside.min()), lr_hi=float(inside.max()),
                                 ll_at_1=float(s[list(FINE).index(1.0)]), ll_at_0=float(s[0]))
    for a, v in zip(FINE, s): prof.append(dict(total=conv, missing_fill=miss, alpha=a, loglik=v, delta_from_max=v - mx))
import csv
def wtsv(fn, rows):
    with open(fn, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
wtsv("profile_ll.tsv", prof)

# ---- bootstrap of alpha-hat (ll matrix; cheap): over lists and over cohort clusters
uid = units_def(d); uniq = sorted(set(uid)); ucols = {u: np.where(uid == u)[0] for u in uniq}
rng = np.random.default_rng(rng0)
for key in ("project_median",):
    M_ = LL[("project", "median")]
    bl = [alpha_hat(M_[:, rng.integers(0, K, K)].sum(1)) for _ in range(2000)]
    bc = []
    for _ in range(2000):
        pick = rng.integers(0, len(uniq), len(uniq)); cols = np.concatenate([ucols[uniq[i]] for i in pick]); bc.append(alpha_hat(M_[:, cols].sum(1)))
    hat[key]["boot_lists_alpha_hat_95"] = [float(np.percentile(bl, 2.5)), float(np.percentile(bl, 97.5))]
    hat[key]["boot_clusters_alpha_hat_95"] = [float(np.percentile(bc, 2.5)), float(np.percentile(bc, 97.5))]
    hat[key]["boot_clusters_alpha_hat_median"] = float(np.median(bc))
    hat[key]["n_units"] = len(uniq)
ah = hat["project_median"]["alpha_hat"]
Pa = probs(d, ah); ea = exp_ge(Pa, (5, 10))
for lab, a_ in (("alpha_hat", ah), ("lr_lo", hat["project_median"]["lr_lo"]), ("lr_hi", hat["project_median"]["lr_hi"])):
    e = exp_ge(probs(d, a_), (5, 10)); hat["project_median"][f"exp5_at_{lab}"] = e[5]; hat["project_median"][f"exp10_at_{lab}"] = e[10]
out.append(dict(total="project", missing_fill="median", alpha=round(ah, 3), obs5=int((d.obs >= 5).sum()), obs10=int((d.obs >= 10).sum()), exp5=ea[5], exp10=ea[10],
                ratio5=(d.obs >= 5).sum() / ea[5], ratio10=(d.obs >= 10).sum() / ea[10], loglik=list_ll(d, Pa).sum(), note="alpha_hat"))

# The list/cluster bootstrap of the ratios (naive_* columns of the 2026-10-10 alpha_grid.tsv; documented there as biased upward and not used) is not repeated.
keys = sorted({k for r in out for k in r}, key=lambda k: (k not in ("total", "missing_fill", "alpha"), k))
for r in out:
    for k in keys: r.setdefault(k, "")
with open("alpha_grid.tsv", "w", newline="") as fh:
    fh.write(f"# script: a1_alpha.py (copied from mauney_add_2026-10-11/prior_strength_membership/a1_alpha.py, run 2026-10-13 on the 64-list base; ratio bootstrap block absent); inputs: primary_lists.json (64 lists, scenario 'published, brain'), DE prior our_prior/1_global_our_prior.tsv; alpha_hat bootstraps: 2000 resamples of lists and of cohort clusters, seed {rng0}\n")
    w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t"); w.writeheader(); w.writerows(out)
json.dump(hat, open("alpha_hat.json", "w"), indent=1)
print(json.dumps(hat["project_median"], indent=1)); print(json.dumps({k: v for k, v in hat.items() if k != "project_median"}, indent=1))
