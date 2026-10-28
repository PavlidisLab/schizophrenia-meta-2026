#!/usr/bin/env python3
"""COPY of analyses/brain_genes_periphery_2026-10-10/scripts/overlap.py; only the AN path line differs (sets.tsv here holds the frozen 209 as study_133*).
Do brain SZ gene sets turn up in blood, iPSC-derived and other cell-model SZ lists more often than chance? (Paul, 2026-10-01)

    python3 analyses/brain_genes_periphery_2026-10-01/scripts/overlap.py      (after build_sets.py)

Machinery: ../../recurrence.py, imported and unchanged. Study selection (load_units: meta-lane, multi-series, lineage group 3
and non-first-look lists left out; same_data_as groups collapsed to one representative, published > reproduced > gemma, then
FDR, then earliest PMID), lists restricted to their own tested universe M_k with denominator n_k (universes.json), and
per-list inclusion probabilities from recurrence.study_tier:
  uniform        p_gk = |L_k & M_k| / n_k for g in M_k
  de_prior       same list total shared out by (n_de_total + 1) of the global DE prior in force in recurrence.py (ours, 8,712
                 Gemma contrasts)
  expr_libd      same total shared out by the pooled list-membership rate of the gene's LIBD brain expression decile (rates
                 estimated from the scenario's own lists). Brain reference: reported for brain and the iPSC scenarios only
  de_prior_blood same total shared out by (n_de_total + 1) of the blood DE prior (../../de_prior_lists/our_prior/
                 3_blood_our_prior.tsv, 784 Gemma blood contrasts); the tissue-fitting null for blood, reported for blood only

Scenarios (representatives among the studies recurrence.py keeps):
  brain        'published, brain' (recurrence.py's scenario); clusters = recurrence.py's cohort tier (lane clusters,
               multi-bank studies per bank). Reference only: the brain sets were made from these lists
  blood        'published, blood' (recurrence.py's scenario, published lists only)
  ipsc_neural  tissue ipsc, neural cells (neurons, NPCs, neural stem cells, organoids, spheroids, neurospheres), any provenance
  ipsc_all     tissue ipsc, any provenance (adds astrocytes, microglia, glial progenitors, endothelial cells, undifferentiated iPSCs)
  other_cells  LCLs (tissue lcl) and tissue 'other' studies of fibroblasts or olfactory cells/epithelium, any provenance
               (Oldmeadow 2014, blood lymphocytes and brain exon splicing, left out)
  nonbrain     blood + ipsc_all + other_cells

Re-use rule (brief item 4): non-brain representatives are clustered, and a cluster counts once, by union-find over lineage
edges between non-brain nodes (lineage.json; every edge kind but reviewed / resource paper / included in), plus each
representative's same_data_as group. Edges through meta-lane nodes and lineage group 3 nodes (no new data) do not chain.
Donor edges (shared subjects, lines, collection) chain through any other non-brain node; data-reuse edges (public dataset,
pooled study, same dataset, same samples) count only when the re-using study is itself a modelled representative, since a
reanalysis that pools two datasets does not make their donors shared. One edge whose evidence shows re-use of results only
is left out (NOT_DONORS). P(g in cluster) = 1 - prod over member lists
(1 - p_gk); a cluster holds g if any member list does; its direction is the majority of the member lists holding g.

Set statistics, per scenario x set x null: hits = sum over set genes of the clusters holding the gene; expected = sum of
P(g in cluster); exact Poisson-binomial tails over the (gene, cluster) pairs (genes treated as independent), and a 95%
interval for observed / expected by resampling clusters (2,000 draws), which keeps any co-listing of set genes within a
cluster. Direction: each hit whose cluster direction and brain direction are both up or down is concordant or discordant;
two-sided exact binomial test against 0.5.

Outputs (results/): set_overlap.tsv, direction.tsv, ph_contrast.tsv, gene_hits.tsv, units.tsv
"""
import collections
import csv
import datetime
import hashlib
import json
import os
import re
import sys

import numpy as np
from scipy.stats import binomtest

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.environ["OUTDIR"]; AN = os.environ["ANALYSES"]   # refreeze copy: output folder and analyses dir from the environment
ROOT = os.path.dirname(AN)
sys.path.insert(0, AN)
import recurrence as R  # noqa: E402

RES = os.path.join(OUT, "results")
SETS = os.environ["SETSTSV"]
BLOOD_PRIOR = os.path.join(AN, "de_prior_lists", "our_prior", "3_blood_our_prior.tsv")
NONCHAIN_KINDS = {"reviewed", "resource paper", "included in"}
REUSE_KINDS = {"public dataset", "pooled study", "same samples (reanalysis)", "same dataset"}
# lineage edges whose evidence says the later study used the earlier one's results, not its donors or lines
NOT_DONORS = {("29263384", "34010744"): "Tiihonen 2021 downloaded Hoffman 2017's neuron DEG table for comparison (edge "
              "evidence); its own lines are Finnish twins (Tiihonen 2019), not Hoffman's childhood-onset SZ lines"}
NEURAL = re.compile(r"neuron|\bnpcs?\b|neural progenitor|neural stem|\bnscs?\b|organoid|spheroid|neurosphere|neuron-committed|\bnccs?\b|\bcins?\b|interneuron", re.I)
OTHER_CELLS = re.compile(r"fibroblast|olfactory|lymphoblast|\blcls?\b", re.I)
OTHER_LEFT_OUT = {"24507884": "Oldmeadow 2014: blood lymphocytes and DLPFC, exon splicing (not a cell model)"}
NBOOT = 2000
SEED = 20261001
SCEN_ORDER = ["brain", "blood", "ipsc_neural", "ipsc_all", "other_cells", "nonbrain"]
NULLS = {"brain": ["uniform", "de_prior", "expr_libd"], "blood": ["uniform", "de_prior", "de_prior_blood"],
         "ipsc_neural": ["uniform", "de_prior", "expr_libd"], "ipsc_all": ["uniform", "de_prior", "expr_libd"],
         "other_cells": ["uniform", "de_prior"], "nonbrain": ["uniform", "de_prior"]}


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]


def node_key(n):
    return n["pmid"] or "doi_" + n["id"].split(":", 1)[-1].replace("/", "_")


def read_sets():
    with open(SETS, encoding="utf-8") as fh:
        rows = list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))
    sets = collections.OrderedDict()
    for r in rows:
        sets.setdefault(r["set"], {})[r["gene"]] = r["brain_direction"]
    return sets


def nonbrain_clusters(nodes, edges, units, modelled):
    """pmid -> cluster root over non-brain lineage edges and same_data_as groups (module docstring). `modelled`: the
    modelled non-brain representatives. Donor edges (shared subjects, lines or collection) chain through any non-brain
    node. Data-reuse edges (REUSE_KINDS: B analysed A's data) link A and B only when B is a modelled representative, whose
    list may then hold A's donors; a non-modelled node that re-used other data (Ji 2026, Gasso 2017, English 2015,
    Tiihonen 2019, which re-analysed Hoffman 2017's series as a replication set) does not bridge its sources and its own
    donors."""
    by = {node_key(n): n for n in nodes}
    id2k = {n["id"]: node_key(n) for n in nodes}
    par = {}

    def find(x):
        par.setdefault(x, x)
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    def ok(k):
        n = by.get(k)
        return n is not None and n.get("tissue") not in ("brain", None) and n.get("lane") != "meta" and n.get("category") != "3"
    used, skipped = [], []
    for e in edges:
        a, b = id2k.get(e["from"], e["from"]), id2k.get(e["to"], e["to"])
        if e.get("kind") in NONCHAIN_KINDS or not (ok(a) and ok(b)):
            continue
        if (a, b) in NOT_DONORS:
            skipped.append((a, b, e.get("kind"), NOT_DONORS[(a, b)]))
            continue
        if e.get("kind") in REUSE_KINDS and b not in modelled:
            skipped.append((a, b, e.get("kind"), "data re-use by a node that is not a modelled representative"))
            continue
        par[find(a)] = find(b)
        used.append((a, b, e.get("kind")))
    for g in units["groups"]:
        for x, y in zip(g, g[1:]):
            par[find(x)] = find(y)
    return find, used, skipped


def tails(ps, obs):
    """Exact Poisson-binomial P(X >= obs) and P(X <= obs) for independent Bernoulli(ps)."""
    ps = ps[ps > 0]
    m = obs + 2  # 0..obs exact, last = more than obs
    d = np.zeros(m)
    d[0] = 1.0
    for q in ps:
        last = d[m - 1] + d[m - 2] * q
        d[1:m - 1] = d[1:m - 1] * (1 - q) + d[0:m - 2] * q
        d[0] *= (1 - q)
        d[m - 1] = last
    return float(min(1.0, d[obs] + d[m - 1])), float(min(1.0, d[:obs + 1].sum()))


def main():
    os.makedirs(RES, exist_ok=True)
    pl = json.load(open(os.path.join(AN, "primary_lists.json")))
    lin = json.load(open(R.LINEAGE))
    nodes, edges = lin["nodes"], lin["edges"]
    tissue_n = {node_key(n): (n.get("tissue_n") or "") for n in nodes}
    units = R.load_units(pl, nodes)
    uni = json.load(open(os.path.join(AN, "universes.json")))["lists"]
    lists, M, nk, utype, nouni = R.load_lists(units, uni)
    studies = units["studies"]
    hit_glob = R.de_prior_weights()
    saved = R.DE_PRIOR
    R.DE_PRIOR = BLOOD_PRIOR  # same columns as the global file; de_prior_table reads R.DE_PRIOR
    hit_blood = R.de_prior_weights()
    R.DE_PRIOR = saved
    expr = R.expression_levels()
    sets = read_sets()

    def is_neural(k):
        return bool(NEURAL.search(tissue_n.get(k, "")))

    def other_ok(k):
        return k not in OTHER_LEFT_OUT and bool(OTHER_CELLS.search(tissue_n.get(k, "")))
    filt = {
        "brain": dict(R.SCENARIOS)["published, brain"],
        "blood": dict(R.SCENARIOS)["published, blood"],
        "ipsc_neural": lambda s: s["tissue"] == "ipsc",
        "ipsc_all": lambda s: s["tissue"] == "ipsc",
        "other_cells": lambda s: s["tissue"] in ("lcl", "other"),
    }
    extra = {"ipsc_neural": is_neural, "other_cells": other_ok}
    usable = R.usable(units, lists)
    reps = {}
    for sc, f in filt.items():
        okf = (lambda e: (lambda k: usable(k) and e(k)))(extra[sc]) if sc in extra else usable
        reps[sc] = R.representatives(units, f, okf)
    reps["nonbrain"] = {**reps["blood"], **reps["ipsc_all"], **reps["other_cells"]}
    find, used_edges, skipped_edges = nonbrain_clusters(nodes, edges, units, {k for k in reps["nonbrain"] if k in M})

    unit_rows, set_rows, dir_rows, gene_rows, contrast_rows = [], [], [], [], []
    rng = np.random.default_rng(SEED)
    for sc in SCEN_ORDER:
        sel = sorted((k for k in reps[sc] if k in M), key=lambda k: (0, int(k)) if k.isdigit() else (1, 0))
        st = R.study_tier(sel, M, nk, lists, hit_glob, expr)
        genes, gi = st["genes"], st["gi"]
        mats = {"uniform": st["P"], "de_prior": st["PP"], "expr_libd": st["PE"]}
        if "de_prior_blood" in NULLS[sc]:
            mats["de_prior_blood"] = R.study_tier(sel, M, nk, lists, hit_blood, expr)["PP"]
        jx = {k: j for j, k in enumerate(sel)}
        # clusters: name -> member indices into sel (+ per-bank lists for brain multi-bank studies)
        members = collections.OrderedDict()
        if sc == "brain":
            cm, mb_out = R.cohort_members(sel, units, st, M, nk)
            # cohort_members gives (name, p unif, p prior, meas, genes, p expr) per member; rebuild with our own matrices
            for c in sorted(cm):
                members[c] = [m[0] for m in cm[c]]
            cname = {m: c for c, ms in members.items() for m in ms}
        else:
            mb_out = []
            for k in sel:
                c = find(k)
                members.setdefault(c, []).append(k)
            cname = {k: c for c, ms in members.items() for k in ms}

        def member_cols(name, mat):
            """inclusion-probability column, list genes and direction map of one cluster member."""
            if ":" in name and not name.startswith("doi_"):  # multi-bank per-bank list (brain)
                k, b = name.split(":", 1)
                bl = units["bank_lists"][k][b]
                Lb = {}
                for n in bl:
                    for g, v in R.read_list(n).items():
                        Lb[g] = R.merge_dir(Lb.get(g), v)
                Lb = {g: v for g, v in Lb.items() if g in M[k]}
                idx = np.array([gi[g] for g in M[k]])
                w = {"uniform": None, "de_prior": st["wall"], "expr_libd": st["wexp"]}
                col = np.zeros(len(genes))
                tot = len(Lb) / nk[k] * len(idx)
                col[idx] = (len(Lb) / nk[k]) if w.get(mat) is None else R.prior_probs(tot, w[mat][idx])
                return col, Lb
            k = name
            return mats[mat][:, jx[k]], {g: v for g, v in lists[k].items() if g in M[k]}
        cl = list(members)

        def mlabel(m):
            k, _, b = m.partition(":") if not m.startswith("doi_") else (m, "", "")
            return studies[k]["label"] + (f" [{b}]" if b else "")

        def clabel(c):
            return mlabel(members[c][0]) if c in members else c

        def L_dir(m, g):
            return member_cols(m, "uniform")[1].get(g, "")
        # observed / direction per cluster (null independent)
        obsC = np.zeros((len(genes), len(cl)), dtype=bool)
        dirC = np.zeros((len(genes), len(cl)), dtype=int)  # +1 up, -1 down, 0 mixed / none
        measC = np.zeros((len(genes), len(cl)), dtype=bool)
        for ci, c in enumerate(cl):
            up, dn = np.zeros(len(genes)), np.zeros(len(genes))
            for m in members[c]:
                col, L = member_cols(m, "uniform")
                measC[:, ci] |= col > 0
                for g, v in L.items():
                    obsC[gi[g], ci] = True
                    if v == "up":
                        up[gi[g]] += 1
                    elif v == "down":
                        dn[gi[g]] += 1
            dirC[:, ci] = np.sign(up - dn).astype(int)
        PC = {}
        for nm in NULLS[sc]:
            Q = np.ones((len(genes), len(cl)))
            for ci, c in enumerate(cl):
                for m in members[c]:
                    Q[:, ci] *= 1 - member_cols(m, nm)[0]
            PC[nm] = 1 - Q
        # units
        for k in sorted(reps[sc], key=lambda k: (0, int(k)) if k.isdigit() else (1, 0)):
            c = cname.get(k) or next((cname[m] for m in cname if m.startswith(k + ":")), "")
            unit_rows.append({"scenario": sc, "representative": k, "label": studies[k]["label"], "tissue": studies[k]["tissue"],
                              "cells_or_tissue": tissue_n.get(k, "")[:160], "provenance": studies[k]["provenance"],
                              "threshold_kind": studies[k]["threshold_kind"],
                              "modelled": "yes" if k in M else f"no ({nouni.get(k) or utype[k]})",
                              "list_genes_in_universe": len(set(lists[k]) & M[k]) if k in M else "",
                              "universe_genes": len(M[k]) if k in M else "", "n_k": round(nk[k]) if k in M else "",
                              "cluster": (c if sc == "brain" else clabel(c)) if k in M else "",
                              "cluster_members": "; ".join(studies[m.split(":")[0]]["label"] + (":" + m.split(":")[1] if ":" in m and not m.startswith("doi_") else "") for m in members.get(c, [])) if k in M else "",
                              "same_data_as_group": "; ".join(studies[m]["label"] for m in reps[sc][k])})
        for k in mb_out:
            unit_rows.append({"scenario": sc, "representative": k, "label": studies[k]["label"], "tissue": studies[k]["tissue"],
                              "modelled": "study tier only (multi-bank, no per-bank lists): left out of clusters"})
        n_cl = len(cl)
        # per-set statistics
        for sname, sg in sets.items():
            idx = np.array([gi[g] for g in sg if g in gi], dtype=int)
            meas_any = measC[idx].any(axis=1) if len(idx) else np.array([], bool)
            obs_hits = int(obsC[idx].sum()) if len(idx) else 0
            obs_genes = int(obsC[idx].any(axis=1).sum()) if len(idx) else 0
            per_cl_obs = obsC[idx].sum(axis=0) if len(idx) else np.zeros(n_cl)
            boot_idx = rng.integers(0, n_cl, size=(NBOOT, n_cl)) if n_cl else None
            for nm in NULLS[sc]:
                P = PC[nm][idx] if len(idx) else np.zeros((0, n_cl))
                exp_hits = float(P.sum())
                exp_genes = float((1 - np.prod(1 - P, axis=1)).sum()) if len(idx) else 0.0
                pu, pd = tails(P.ravel(), obs_hits) if len(idx) else (1.0, 1.0)
                per_cl_exp = P.sum(axis=0)
                if n_cl and exp_hits > 0:
                    bo, be = per_cl_obs[boot_idx].sum(axis=1), per_cl_exp[boot_idx].sum(axis=1)
                    rr = bo[be > 0] / be[be > 0]
                    lo, hi = np.percentile(rr, [2.5, 97.5])
                else:
                    lo = hi = float("nan")
                set_rows.append({"scenario": sc, "set": sname, "null": nm, "n_set": len(sg), "n_measured": int(meas_any.sum()),
                                 "n_clusters": n_cl, "n_lists": len(sel), "obs_hits": obs_hits, "exp_hits": round(exp_hits, 2),
                                 "obs_exp": round(obs_hits / exp_hits, 3) if exp_hits else "",
                                 "ci95_lo": round(float(lo), 3) if exp_hits else "", "ci95_hi": round(float(hi), 3) if exp_hits else "",
                                 "p_enriched": f"{pu:.3g}", "p_depleted": f"{pd:.3g}",
                                 "obs_genes": obs_genes, "exp_genes": round(exp_genes, 2)})
            # direction (null independent)
            bd = np.array([{"up": 1, "down": -1}.get(sg[g], 0) for g in sg if g in gi], dtype=int)
            if len(idx):
                hit = obsC[idx] & (dirC[idx] != 0) & (bd[:, None] != 0)
                conc = int((hit & (dirC[idx] == bd[:, None])).sum())
                disc = int((hit & (dirC[idx] == -bd[:, None])).sum())
                # cluster-level vote: per cluster, concordant minus discordant hits
                net = ((hit & (dirC[idx] == bd[:, None])).sum(axis=0) - (hit & (dirC[idx] == -bd[:, None])).sum(axis=0))
            else:
                conc = disc = 0
                net = np.zeros(n_cl)
            p = binomtest(conc, conc + disc, 0.5).pvalue if conc + disc else float("nan")
            dir_rows.append({"scenario": sc, "set": sname, "hits_with_direction": conc + disc, "concordant": conc, "discordant": disc,
                             "share_concordant": round(conc / (conc + disc), 3) if conc + disc else "",
                             "binom_p": f"{p:.3g}" if conc + disc else "",
                             "clusters_net_concordant": int((net > 0).sum()), "clusters_net_discordant": int((net < 0).sum()),
                             "clusters_tied_or_none": int((net == 0).sum())})
        # pH contrast: genes not pH-associated (|Z| < 1) vs pH-associated (|Z| >= 1)
        for base in ("study_133", "lists5plus"):
            a, b = sets[f"{base}_ph_lt1"], sets[f"{base}_ph_ge1"]
            ia = np.array([gi[g] for g in a if g in gi], dtype=int)
            ib = np.array([gi[g] for g in b if g in gi], dtype=int)
            for nm in NULLS[sc]:
                oa, ob = int(obsC[ia].sum()), int(obsC[ib].sum())
                ea, eb = float(PC[nm][ia].sum()), float(PC[nm][ib].sum())
                p = binomtest(oa, oa + ob, ea / (ea + eb)).pvalue if oa + ob and ea + eb else float("nan")
                ca, cb = obsC[ia].sum(axis=0), obsC[ib].sum(axis=0)
                xa, xb = PC[nm][ia].sum(axis=0), PC[nm][ib].sum(axis=0)
                bi = rng.integers(0, n_cl, size=(NBOOT, n_cl))
                ra = ca[bi].sum(axis=1) / np.maximum(xa[bi].sum(axis=1), 1e-12)
                rb = cb[bi].sum(axis=1) / np.maximum(xb[bi].sum(axis=1), 1e-12)
                okb = (xa[bi].sum(axis=1) > 0) & (xb[bi].sum(axis=1) > 0) & (rb > 0)
                lo, hi = np.percentile(ra[okb] / rb[okb], [2.5, 97.5]) if okb.any() else (float("nan"),) * 2
                contrast_rows.append({"scenario": sc, "genes": base, "null": nm, "n_ph_lt1": len(a), "n_ph_ge1": len(b),
                                      "obs_lt1": oa, "exp_lt1": round(ea, 2), "obs_exp_lt1": round(oa / ea, 3) if ea else "",
                                      "obs_ge1": ob, "exp_ge1": round(eb, 2), "obs_exp_ge1": round(ob / eb, 3) if eb else "",
                                      "ratio_lt1_over_ge1": round((oa / ea) / (ob / eb), 3) if ea and eb and ob else "",
                                      "ci95_lo": round(float(lo), 3), "ci95_hi": round(float(hi), 3),
                                      "p_conditional_binomial": f"{p:.3g}"})
        # per-gene rows for the named small sets and the pH / careful sets
        keep = set().union(*(sets[s] for s in ("interneuron_14", "ieg_4", "careful_25", "study_133", "ph_top200")))
        for g in sorted(keep):
            if g not in gi:
                gene_rows.append({"scenario": sc, "gene": g, "n_clusters_measured": 0})
                continue
            i = gi[g]
            # the member lists holding the gene, grouped by cluster
            hitc = [" + ".join(mlabel(m) + f" ({L_dir(m, g)})" for m in members[cl[ci]] if g in member_cols(m, "uniform")[1])
                    for ci in range(n_cl) if obsC[i, ci]]
            gene_rows.append({"scenario": sc, "gene": g, "n_clusters_measured": int(measC[i].sum()),
                              "n_clusters_observed": int(obsC[i].sum()),
                              "expected_uniform": round(float(PC["uniform"][i].sum()), 3),
                              "expected_de_prior": round(float(PC["de_prior"][i].sum()), 3),
                              "clusters_up": int(((dirC[i] == 1) & obsC[i]).sum()), "clusters_down": int(((dirC[i] == -1) & obsC[i]).sum()),
                              "lists_holding_gene": "; ".join(hitc)})
        print(sc, "lists", len(sel), "clusters", n_cl, "genes", len(genes))

    today = datetime.date.today().isoformat()
    head = [f"# script: analyses/brain_genes_periphery_2026-10-01/scripts/overlap.py, run {today}",
            f"# input: sets.tsv sha256 {sha(SETS)}; analyses/primary_lists.json {sha(os.path.join(AN, 'primary_lists.json'))}; "
            f"analyses/universes.json {sha(os.path.join(AN, 'universes.json'))}; corpus/lineage.json {sha(R.LINEAGE)}",
            f"# input: {os.path.relpath(R.DE_PRIOR, ROOT)} {sha(R.DE_PRIOR)} (de_prior); {os.path.relpath(BLOOD_PRIOR, ROOT)} "
            f"{sha(BLOOD_PRIOR)} (de_prior_blood); {os.path.relpath(R.EXPR_FILE, ROOT)} {sha(R.EXPR_FILE)} (expr_libd)",
            "# machinery: analyses/recurrence.py (load_units, load_lists, representatives, study_tier, cohort_members), unchanged; "
            "clusters: lineage edges between non-brain nodes (not through meta-lane or group 3 nodes) + same_data_as groups; "
            "brain: recurrence.py cohort tier",
            "# expr_libd uses LIBD brain expression deciles (reported for brain and iPSC); de_prior_blood is the blood DE prior "
            "(reported for blood)"]

    def write(fn, rows, extra=()):
        fields = list(dict.fromkeys(k for r in rows for k in r))
        with open(os.path.join(RES, fn), "w", encoding="utf-8") as fh:
            fh.write("\n".join(head + list(extra)) + "\n")
            w = csv.DictWriter(fh, fields, delimiter="\t", restval="")
            w.writeheader()
            w.writerows(rows)
    write("set_overlap.tsv", set_rows, ["# obs_hits: (gene, cluster) pairs where the cluster holds the gene; exp_hits: sum of "
                                        "P(g in cluster); p_enriched = P(X >= obs), p_depleted = P(X <= obs), exact "
                                        "Poisson-binomial over pairs (genes independent); ci95: obs/exp over 2,000 cluster "
                                        "resamples; obs_genes / exp_genes: set genes in at least one cluster"])
    write("direction.tsv", dir_rows, ["# hits whose cluster direction (majority of member lists) and brain direction are both "
                                      "up/down; ph_* sets: brain direction = the low-pH prediction; binom_p two-sided vs 0.5"])
    write("ph_contrast.tsv", contrast_rows, ["# genes of the 133 / 1,151 with control pH |Z| < 1 (lt1) vs >= 1 (ge1); "
                                             "p_conditional_binomial: obs_lt1 ~ Binomial(obs_lt1 + obs_ge1, exp_lt1 / "
                                             "(exp_lt1 + exp_ge1)); ci95 of ratio_lt1_over_ge1 over 2,000 cluster resamples"])
    write("gene_hits.tsv", gene_rows, ["# genes of interneuron_14, ieg_4, careful_25, study_133 and ph_top200; clusters_up / "
                                       "_down: clusters holding the gene, by majority direction"])
    write("units.tsv", unit_rows, ["# one row per representative and scenario; cluster = re-use cluster (lineage edges between "
                                   "non-brain nodes, same_data_as groups; brain: recurrence.py cohort clusters)",
                                   "# lineage edges used for non-brain clusters: " + "; ".join(
                                       f"{studies[a]['label'] if a in studies else a} - {studies[b]['label'] if b in studies else b} ({k})"
                                       for a, b, k in used_edges),
                                   "# lineage edges not used (data re-use by a non-modelled node, or evidence of result re-use "
                                   "only): " + "; ".join(
                                       f"{studies[a]['label'] if a in studies else a} - {studies[b]['label'] if b in studies else b} ({k}; {why})"
                                       for a, b, k, why in skipped_edges)])


if __name__ == "__main__":
    main()
