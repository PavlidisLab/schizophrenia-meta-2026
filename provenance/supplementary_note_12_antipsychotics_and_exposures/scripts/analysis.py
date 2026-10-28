"""Are the recurrent SZ genes antipsychotic-responsive, and does the drug effect point the way the SZ effect does?

    opt/anaconda3/bin/python3 scripts/analysis.py
Inputs: results/contrasts_long.tsv.gz (build_long.py), ../confounds/recurrent_sets.tsv (study_brain: rule-A recurrent
genes, uniform q < 0.1, with sz_direction), ../recurrence_genes.tsv ('published, brain', q_expr_bh < 0.1: the
expression-matched set; direction = majority of list directions), the DE prior in force in ../recurrence.py (global), LIBD DLPFC mean logCPM
(../gwas_vs_de/results/libd_mean_logcpm.tsv).

Per unit (animals or cultures sharing controls), gene z = mean of its drug contrasts' z (drug vs vehicle).
Responsiveness: mean |z| over the set's measured genes. Direction: mean of sz_sign x z over genes with a SZ direction
(positive = the drug moves genes the way SZ does). Null: NPERM random human gene sets matched gene by gene on DE-prior
quintile x LIBD DLPFC expression decile (genes outside either table form their own bins); for direction the drawn genes
carry the set genes' SZ signs. One draw is evaluated in every unit, so combined tests (mean over a class's units) keep
the between-unit correlation of gene responsiveness. p: responsiveness one-sided (greater), direction two-sided.
Calibration: NCAL random matched sets taken as the 'observed' set give the p distribution under the null.
Leave-one-out: the combined chronic tests with each chronic unit left out in turn.
Writes results/set_tests.tsv, results/gene_effects.tsv, results/calibration.tsv, results/leave_one_out.tsv,
results/analysis_meta.json.
"""
import csv, gzip, hashlib, json, os, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GL = os.path.dirname(HERE)
DATA = os.path.normpath(os.path.join(GL, "..", "data"))  # repo data/ (lists, manifests, reference, raw)
RES = os.path.join(HERE, "results")
NPERM, NCAL, SEED = 5000, 300, 20260926
rng = np.random.default_rng(SEED)
IN = {"long": os.path.join(RES, "contrasts_long.tsv.gz"), "sets": os.path.join(GL, "confounds", "recurrent_sets.tsv"),
      "rec": os.path.join(GL, "recurrence_genes.tsv"), "prior": None,  # set below: the DE prior file of recurrence.py
      "libd": os.path.join(GL, "gwas_vs_de", "results", "libd_mean_logcpm.tsv")}
sys.path.insert(0, GL)
import recurrence as R  # noqa: E402  (de_prior_ranks; DE prior default global since 2026-09-29)
IN["prior"] = R.DE_PRIOR


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def load():
    L = pd.read_csv(IN["long"], sep="\t", comment="#")
    unit_class = L.groupby("unit")["class"].first()
    Z = L.groupby(["unit", "gene"])["z"].mean().unstack(0)          # genes x units
    rs = pd.read_csv(IN["sets"], sep="\t", comment="#")
    s1 = rs[rs.set == "study_brain"].set_index("gene")["sz_direction"]
    rec = pd.read_csv(IN["rec"], sep="\t", comment="#")
    rb = rec[(rec.scenario == "published, brain") & (rec.q_expr_bh < 0.1)].copy()
    rb["sz_direction"] = np.where(rb.up > rb.down, "up", np.where(rb.down > rb.up, "down", "mixed"))
    s2 = rb.set_index("gene")["sz_direction"]
    sets = {"recurrent (study_brain, uniform q < 0.1)": s1, "recurrent (expression-matched q < 0.1)": s2}
    prior = pd.Series(R.de_prior_ranks(), name="DE_Prior_Rank")  # the DE prior in force (normalized rank, 1 = most often DE)
    lb = pd.read_csv(IN["libd"], sep="\t", comment="#")
    expr = lb[lb.region == "DLPFC"].drop_duplicates("symbol").set_index("symbol")["mean_logcpm"]
    return Z, unit_class, sets, prior, expr


def bins_for(genes, prior, expr):
    pq = pd.qcut(prior.reindex(genes), 5, labels=False, duplicates="drop")
    eq = pd.qcut(expr.reindex(genes), 10, labels=False, duplicates="drop")
    return (pq.fillna(-1).astype(int).astype(str) + "_" + eq.fillna(-1).astype(int).astype(str)).values


def main():
    Z, unit_class, sets, prior, expr = load()
    genes = Z.index.values
    b = bins_for(genes, prior, expr)
    by_bin = {k: np.where(b == k)[0] for k in np.unique(b)}
    gidx = {g: i for i, g in enumerate(genes)}
    Zv = Z.values                                           # genes x units, NaN where not measured
    units = list(Z.columns)
    classes = {c: [units.index(u) for u in units if unit_class[u] == c] for c in ("chronic", "acute", "cells")}

    def draw(idx, n):
        """n matched draws for gene indices idx: array n x len(idx)"""
        out = np.empty((n, len(idx)), dtype=int)
        for j, i in enumerate(idx):
            pool = by_bin[b[i]]
            out[:, j] = pool[rng.integers(0, len(pool), n)]
        return out

    def stats(D, sg):
        """D: draws x k gene indices; sg: k SZ signs -> resp[draws x units], dirn[draws x units]"""
        z = Zv[D]                                             # draws x k x units
        resp = np.nanmean(np.abs(z), axis=1)
        ok = (sg != 0)
        dirn = np.nanmean(z[:, ok, :] * sg[ok][None, :, None], axis=1) if ok.any() else np.full(resp.shape, np.nan)
        return resp, dirn

    rows, cal_rows, gene_rows, loo_rows = [], [], [], []
    for sname, sdir in sets.items():
        idx = np.array([gidx[g] for g in sdir.index if g in gidx])
        sg = np.array([{"up": 1, "down": -1}.get(sdir[genes[i]], 0) for i in idx])
        obs_r, obs_d = stats(idx[None, :], sg)
        D = draw(idx, NPERM)
        null_r, null_d = stats(D, sg)
        def p_greater(o, n): return (1 + np.sum(n >= o - 1e-12)) / (len(n) + 1)
        def p_two(o, n):
            m = np.nanmean(n); return (1 + np.sum(np.abs(n - m) >= abs(o - m) - 1e-12)) / (len(n) + 1)
        for u, uname in enumerate(units):
            nm = int(np.sum(~np.isnan(Zv[idx, u])))
            rows.append({"set": sname, "unit": uname, "class": unit_class[uname], "n_set_genes": len(idx), "n_measured": nm,
                         "n_p05": int(np.sum(np.abs(Zv[idx, u]) > 1.959964)),
                         "resp_obs": obs_r[0, u], "resp_null_mean": np.nanmean(null_r[:, u]), "resp_p": p_greater(obs_r[0, u], null_r[:, u]),
                         "dir_obs": obs_d[0, u], "dir_null_mean": np.nanmean(null_d[:, u]), "dir_p": p_two(obs_d[0, u], null_d[:, u]),
                         "n_sig_same_as_sz": int(np.sum((np.abs(Zv[idx, u]) > 1.959964) & (np.sign(Zv[idx, u]) == sg) & (sg != 0))),
                         "n_sig_opposite_sz": int(np.sum((np.abs(Zv[idx, u]) > 1.959964) & (np.sign(Zv[idx, u]) == -sg) & (sg != 0)))})
        for c, cu in classes.items():
            if not cu:
                continue
            o_r, n_r = np.nanmean(obs_r[0, cu]), np.nanmean(null_r[:, cu], axis=1)
            o_d, n_d = np.nanmean(obs_d[0, cu]), np.nanmean(null_d[:, cu], axis=1)
            rows.append({"set": sname, "unit": f"ALL {c} ({len(cu)} units)", "class": c, "n_set_genes": len(idx), "n_measured": "",
                         "n_p05": int(sum(np.sum(np.abs(Zv[idx, u]) > 1.959964) for u in cu)),
                         "resp_obs": o_r, "resp_null_mean": n_r.mean(), "resp_p": p_greater(o_r, n_r),
                         "dir_obs": o_d, "dir_null_mean": np.nanmean(n_d), "dir_p": p_two(o_d, n_d),
                         "n_sig_same_as_sz": int(sum(np.sum((np.abs(Zv[idx, u]) > 1.959964) & (np.sign(Zv[idx, u]) == sg) & (sg != 0)) for u in cu)),
                         "n_sig_opposite_sz": int(sum(np.sum((np.abs(Zv[idx, u]) > 1.959964) & (np.sign(Zv[idx, u]) == -sg) & (sg != 0)) for u in cu))})
            # calibration: random matched sets as the observed set, each against its own null (fewer draws)
            if c == "chronic":
                pr_, pd_ = [], []
                for _ in range(NCAL):
                    fake = draw(idx, 1)[0]
                    fr, fd = stats(fake[None, :], sg)
                    Dn = draw(fake, 500)
                    nr, nd = stats(Dn, sg)
                    pr_.append(p_greater(np.nanmean(fr[0, cu]), np.nanmean(nr[:, cu], axis=1)))
                    ndm = np.nanmean(nd[:, cu], axis=1)
                    pd_.append(p_two(np.nanmean(fd[0, cu]), ndm))
                for name, v in (("responsiveness", pr_), ("direction", pd_)):
                    v = np.array(v)
                    cal_rows.append({"set": sname, "class": c, "test": name, "n": len(v), "frac_p_le_05": np.mean(v <= 0.05),
                                     "frac_p_le_01": np.mean(v <= 0.01), "median_p": np.median(v)})
        # leave one chronic unit out: combined chronic tests without it
        cu = classes["chronic"]
        for drop in cu:
            k = [u for u in cu if u != drop]
            o_r, n_r = np.nanmean(obs_r[0, k]), np.nanmean(null_r[:, k], axis=1)
            o_d, n_d = np.nanmean(obs_d[0, k]), np.nanmean(null_d[:, k], axis=1)
            loo_rows.append({"set": sname, "left_out": units[drop], "resp_obs": o_r, "resp_p": p_greater(o_r, n_r),
                             "dir_obs": o_d, "dir_p": p_two(o_d, n_d)})
        for i, s in zip(idx, sg):
            zz = Zv[i, classes["chronic"]]
            gene_rows.append({"set": sname, "gene": genes[i], "sz_direction": {1: "up", -1: "down", 0: "mixed"}[s],
                              "chronic_units_measured": int(np.sum(~np.isnan(zz))), "chronic_mean_z": np.nanmean(zz) if np.any(~np.isnan(zz)) else np.nan,
                              "chronic_units_up_p05": int(np.sum(zz > 1.959964)), "chronic_units_down_p05": int(np.sum(zz < -1.959964))})
    fmt = lambda df: df.apply(lambda c: c.map(lambda x: f"{x:.4g}" if isinstance(x, float) else x))  # noqa: E731
    head = [f"# medication/scripts/analysis.py, NPERM {NPERM}, seed {SEED}; matched on DE-prior quintile x LIBD DLPFC expression decile",
            "# resp = mean |z| over the set's measured genes (drug vs vehicle); dir = mean sz_sign x z (positive: drug moves genes the way SZ does)"]
    for fn, df in (("set_tests.tsv", pd.DataFrame(rows)), ("gene_effects.tsv", pd.DataFrame(gene_rows)), ("calibration.tsv", pd.DataFrame(cal_rows)),
                   ("leave_one_out.tsv", pd.DataFrame(loo_rows))):
        with open(os.path.join(RES, fn), "w") as fh:
            fh.write("\n".join(head) + "\n")
            fmt(df).to_csv(fh, sep="\t", index=False)
    json.dump({"script": "medication/scripts/analysis.py", "nperm": NPERM, "ncal": NCAL, "seed": SEED,
               "inputs": {k: {"path": os.path.relpath(p, GL), "sha256": sha(p)} for k, p in IN.items()}},
              open(os.path.join(RES, "analysis_meta.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
