#!/usr/bin/env python3
"""Analysis 1: agreement among Gemma's uniform reanalyses vs among the published lists of the same datasets.

Datasets: studies whose rule A primary list is published (primary_lists.json) and that have a Gemma result set for
the same GEO series and, for multi-subset series, the same region or cell type (DATASETS below; one entry per series).
Left out: Iwamoto 2005 (published list covers mitochondria-related probe sets only), Topol 2016 (miRNA lists).

For each pair of datasets (i, j), over the common universe U = M_i & M_j & G_i & G_j (M = the published list's tested
gene set from universes.json, G = genes with a contrast p in the Gemma table):
  published      the published lists restricted to U (sizes n_i, n_j);
  gemma_sized    Gemma's n_i and n_j genes with the smallest contrast p within U (ties by |log2 fold change|), so each
                 Gemma list has exactly the size of the published list of the same dataset;
  published_topN / gemma_topN   N = min(n_i, n_j, 100): each published list's N genes by p (or fdr; h1.top_genes,
                 pairs with an unranked published list are left out) and Gemma's top N by p;
  gemma_fdr      Gemma's own lists (BH q < 0.05, lists/*gemma*.tsv) restricted to U, where both are non-empty.
  overlap, expected |a||b|/|U|, fold, log2 fold ((obs + 0.5) / (exp + 0.5)), concordance of direction over shared genes;
  rho_gemma      Spearman correlation of Gemma log2 fold changes over all genes in both Gemma tables.
Pairs of datasets in one same_data_as group (shared cell lines) are reported but left out of the summaries.
Dataset-level bootstrap (datasets resampled with replacement) of gemma minus published.

    python3 a1_uniform_vs_published.py   -> a1_datasets.tsv, a1_pairs.tsv, a1_summary.tsv, a1_boot.tsv,
                                            a1_uniform_vs_published_meta.json
"""
import itertools
import math
import os

import numpy as np

import fl_common as C

# (study pmid, Gemma list names (several = combined by smallest p per gene), contrast of the Gemma subset matches the
#  published list's, note)
DATASETS = [
    ("19255580", ["19255580__gemma_GSE17612_all_fdr05"], "yes", "Maycox 2009 CCHPC BA10 (GSE17612); a member of Glatt 2005's group"),
    ("18778695", ["18778695__gemma_GSE21138_all_fdr05"], "no", "published list: short illness duration subgroup (8 vs 8); Gemma: all 30 vs 29"),
    ("21538462", ["21538462__gemma_GSE21935_all_fdr05"], "yes", "BA22"),
    ("24702465", ["24702465__gemma_GSE37981_all_fdr05"], "yes", "layer 3 pyramidal neurons, BA42"),
    ("25560755", ["25560755__gemma_GSE93987_dorsolateral_prefrontal_cortex_layer_3_fdr05",
                  "25560755__gemma_GSE93987_dorsolateral_prefrontal_cortex_layer_5_fdr05"], "yes",
     "published list pools layer 3 and layer 5; Gemma's two layer result sets combined by smallest p per gene"),
    ("31123247", ["31123247__gemma_GSE53987_dorsolateral_prefrontal_cortex_fdr05"], "yes", "DLPFC subset"),
    ("31288836", ["31288836__gemma_GSE107638_all_fdr05"], "no", "published list: NeuN+ nuclei; Gemma subset 'all' (43 vs 53 samples)"),
    ("32152472", ["32152472__gemma_GSE138082_dentate_gyrus_of_hippocampal_formation_fdr05"], "yes", "dentate gyrus"),
    ("36302706", ["36302706__gemma_GSE202537_caudate_nucleus_fdr05"], "yes", "caudate"),
    ("21490598", ["21490598__gemma_GSE25673_all_fdr05"], "yes", "hiPSC neurons"),
    ("24686136", ["24686136__gemma_GSE40102_all_fdr05"], "yes", "hiPSC NPCs"),
    ("28736215", ["28736215__gemma_GSE86906_all_fdr05"], "no", "published list: genes shared by four per-patient comparisons"),
    ("29263384", ["29263384__gemma_GSE106589_neural_progenitor_cell_fdr05"], "yes", "hiPSC NPCs"),
    ("30664768", ["30664768__gemma_GSE118313_all_fdr05"], "yes", "hiPSC cINs, pilot cohort"),
    ("33677214", ["33677214__gemma_GSE155067_cin_migrated_fdr05"], "yes", "hiPSC migrated cINs"),
    ("35701597", ["35701597__gemma_GSE184102_all_fdr05"], "yes", "hiPSC glutamatergic neurons"),
    ("31019265", ["31019265__gemma_GSE125999_all_fdr05"], "yes", "hiPSC cINs, expanded cohort"),
    ("32186681", ["32186681__gemma_GSE133534_all_fdr05"], "yes", "cerebral organoids"),
    ("33288841", ["33288841__gemma_GSE182875_all_fdr05"], "yes", "hiPSC cortical neurons"),
    ("34936134", ["34936134__gemma_GSE191248_all_fdr05"], "yes", "hiPSC astrocytes"),
    ("37117634", ["37117634__gemma_GSE213232_all_fdr05"], "yes", "iPSC microglia-like cells"),
    ("39307399", ["39307399__gemma_GSE282524_glutamatergic_neuron_derives_from_cell_l_fdr05"], "yes", "NGN2 glutamatergic neurons"),
    ("22393424", ["22393424__gemma_GSE27383_all_fdr05"], "yes", "PBMC"),
    ("25658856", ["25658856__gemma_GSE62333_all_fdr05"], "yes", "fibroblasts"),
    ("27727244", ["27727244__gemma_GSE73129_all_fdr05"], "no", "published list: microarray cohort 15 vs 15; Gemma 55 vs 60 samples"),
]
LEFT_OUT = [("15563509", "Iwamoto 2005: published list covers mitochondria-related probe sets only (targeted universe)"),
            ("27117414", "Topol 2016: miRNA lists (universe type exclude)")]


def gemma_table(ctx, names):
    """symbol -> (effect, p): smallest p across the named Gemma result sets."""
    out = {}
    for n in names:
        for g, (e, p) in ctx.table(n).items():
            if p is None:
                continue
            if g not in out or p < out[g][1]:
                out[g] = (e, p)
    return out


def gemma_fdr_genes(names):
    import recurrence as R
    s = {}
    for n in names:
        for g, d in R.read_list(n).items():
            s[g] = d
    return s


def top_by_p(T, U, n):
    gs = sorted((g for g in U if g in T), key=lambda g: (T[g][1], -(abs(T[g][0]) if T[g][0] is not None else 0), g))
    return gs[:n]


def score(a, b, U, da, db):
    """overlap metrics for sets a, b over universe U; da, db: gene -> 'up'/'down'."""
    ov = a & b
    exp_ = len(a) * len(b) / len(U) if U else 0
    same = sum(1 for g in ov if da.get(g) in ("up", "down") and da.get(g) == db.get(g))
    opp = sum(1 for g in ov if {da.get(g), db.get(g)} == {"up", "down"})
    return {"size_a": len(a), "size_b": len(b), "overlap": len(ov), "expected": exp_,
            "fold": len(ov) / exp_ if exp_ else None, "log2_fold": math.log2((len(ov) + 0.5) / (exp_ + 0.5)),
            "n_dir": same + opp, "n_same": same, "concordance": same / (same + opp) if same + opp else None,
            "p_hyper": C.H.hyper_sf(len(ov), len(U), len(a), len(b))}


def sign_dir(e):
    return "up" if e is not None and e > 0 else ("down" if e is not None and e < 0 else "")


def main():
    ctx = C.Ctx()
    ds = []
    for k, names, cmatch, note in DATASETS:
        s = ctx.studies[k]
        T = gemma_table(ctx, names)
        G = set(T)
        ds.append({"study": k, "label": s["label"], "tissue": s["tissue"], "published_list": s["primary"][0],
                   "published_threshold": s["threshold_kind"], "published_universe": ctx.utype[k],
                   "n_published_genes": len(set(ctx.lists[k]) & ctx.M[k]), "gemma_lists": ";".join(names),
                   "n_gemma_genes_tested": len(G), "n_gemma_fdr_genes": len(set(gemma_fdr_genes(names)) & G),
                   "region_published": ctx.region(s["primary"][0]), "region_gemma": ctx.region(names[0]),
                   "contrast_matches": cmatch, "group": ctx.units["group_of"][k], "note": note, "_T": T})
    pairs = []
    for x, y in itertools.combinations(ds, 2):
        i, j = x["study"], y["study"]
        U = ctx.M[i] & ctx.M[j] & set(x["_T"]) & set(y["_T"])
        if len(U) < 100:
            continue
        a = set(ctx.lists[i]) & U
        b = set(ctx.lists[j]) & U
        da = {g: ctx.stats[i][g][0] for g in a}
        db = {g: ctx.stats[j][g][0] for g in b}
        Ta, Tb = x["_T"], y["_T"]
        ga, gb = set(top_by_p(Ta, U, len(a))), set(top_by_p(Tb, U, len(b)))
        gda = {g: sign_dir(Ta[g][0]) for g in U}
        gdb = {g: sign_dir(Tb[g][0]) for g in U}
        base = {"study_a": i, "label_a": x["label"], "tissue_a": x["tissue"], "study_b": j, "label_b": y["label"],
                "tissue_b": y["tissue"], "n_universe": len(U),
                "tissue_pair": "-".join(sorted((x["tissue"], y["tissue"]))),
                "same_group": "yes" if x["group"] == y["group"] else "no",
                "contrast_matches_both": "yes" if x["contrast_matches"] == y["contrast_matches"] == "yes" else "no",
                "region_pair": f"{x['region_published']}/{y['region_published']}"}
        rows = {"published": score(a, b, U, da, db), "gemma_sized": score(ga, gb, U, gda, gdb)}
        n = min(len(a), len(b), C.N_TOP)
        ta, _ = C.H.top_genes(ctx.stats[i], a, n)
        tb, _ = C.H.top_genes(ctx.stats[j], b, n)
        if n and ta is not None and tb is not None:
            rows["published_topN"] = score(ta, tb, U, da, db)
            rows["gemma_topN"] = score(set(top_by_p(Ta, U, n)), set(top_by_p(Tb, U, n)), U, gda, gdb)
        fa = set(gemma_fdr_genes(x["gemma_lists"].split(";"))) & U
        fb = set(gemma_fdr_genes(y["gemma_lists"].split(";"))) & U
        if fa and fb:
            rows["gemma_fdr"] = score(fa, fb, U, gda, gdb)
        sh = [g for g in Ta if g in Tb and Ta[g][0] is not None and Tb[g][0] is not None]
        rho = C.spearman([Ta[g][0] for g in sh], [Tb[g][0] for g in sh]) if len(sh) >= 3 else None
        for meth, m in rows.items():
            pairs.append(dict(base, method=meth, top_n=n if "topN" in meth else "", rho_gemma=rho, n_rho_gemma=len(sh), **m))
    # summaries
    summ, bt = [], []
    classes = [("all pairs, different groups", lambda r: r["same_group"] == "no"),
               ("brain-brain, different groups", lambda r: r["same_group"] == "no" and r["tissue_pair"] == "brain-brain"),
               ("ipsc-ipsc, different groups", lambda r: r["same_group"] == "no" and r["tissue_pair"] == "ipsc-ipsc"),
               ("all pairs, different groups, contrasts match", lambda r: r["same_group"] == "no" and r["contrast_matches_both"] == "yes"),
               ("brain-brain, different groups, contrasts match", lambda r: r["same_group"] == "no" and r["tissue_pair"] == "brain-brain"
                and r["contrast_matches_both"] == "yes"),
               ("same group (shared lines)", lambda r: r["same_group"] == "yes")]
    for cname, f in classes:
        for meth in ("published", "gemma_sized", "published_topN", "gemma_topN", "gemma_fdr"):
            rows = [r for r in pairs if f(r) and r["method"] == meth]
            if not rows:
                continue
            s = C.summarize(rows, extra_median=("rho_gemma",))
            summ.append(dict({"pair_class": cname, "method": meth}, **s))
        for ma, mb in (("gemma_sized", "published"), ("gemma_topN", "published_topN")):
            rows = [r for r in pairs if f(r) and r["method"] in (ma, mb)]
            keys = {(r["study_a"], r["study_b"]) for r in rows if r["method"] == ma} & {(r["study_a"], r["study_b"]) for r in rows if r["method"] == mb}
            rows = [r for r in rows if (r["study_a"], r["study_b"]) in keys]
            if not rows:
                continue
            la = {(r["study_a"], r["study_b"]): r["log2_fold"] for r in rows if r["method"] == ma}
            lb = {(r["study_a"], r["study_b"]): r["log2_fold"] for r in rows if r["method"] == mb}
            d = np.array([la[k] - lb[k] for k in keys])
            b = C.boot_diff(rows, lambda r, ma=ma: r["method"] == ma, lambda r, mb=mb: r["method"] == mb)
            for k, v in b.items():
                if k == "pooled_top_fold":  # a1 rows carry their own sets in overlap / expected; no separate top-N columns
                    continue
                bt.append({"pair_class": cname, "comparison": f"{ma} minus {mb}", "statistic": k, "estimate": v[0],
                           "boot_lo95": v[1], "boot_hi95": v[2], "share_boot_gt0": v[3], "n_pairs": len(keys),
                           "median_paired_log2_diff": float(np.median(d)), "n_gemma_higher": int((d > 0).sum()),
                           "n_published_higher": int((d < 0).sum()), "n_equal": int((d == 0).sum())})
    for x in ds:
        x.pop("_T")
    head = [f"# a1_uniform_vs_published.py, run {C.TODAY}; rule A primary lists (published side); Gemma result sets from "
            "manifests/batch6_gemma.json (gemma/<GSE>/*_sz_vs_control.tsv)",
            "# left out: " + "; ".join(t for _, t in LEFT_OUT)]
    outs = []
    for fn, rows, extra in (("a1_datasets.tsv", ds, ["# one row per dataset: published primary list and the Gemma result set(s) matched to it"]),
                            ("a1_pairs.tsv", pairs, ["# one row per dataset pair and method (see script docstring); U = common universe"]),
                            ("a1_summary.tsv", summ, ["# summaries per pair class and method (fl_common.summarize); median_rho_gemma = "
                                                      "Spearman of Gemma log2 fold changes over all shared genes"]),
                            ("a1_boot.tsv", bt, [f"# dataset-level bootstrap ({C.NBOOT}, seed {C.SEED}) of gemma minus published, same pairs; "
                                                 "median_paired_log2_diff and n_*_higher are over pairs, no resampling"])):
        p = os.path.join(C.HERE, fn)
        C.write_tsv(p, rows, head + extra)
        outs.append(p)
    C.write_meta("a1_uniform_vs_published.py", outs, extra_inputs=["../data/manifests/batch6_gemma.json"],
                 params={"n_top": C.N_TOP, "min_universe": 100, "datasets": [d[0] for d in DATASETS]})
    print(len(ds), "datasets", len({(r['study_a'], r['study_b']) for r in pairs}), "pairs")


if __name__ == "__main__":
    main()
