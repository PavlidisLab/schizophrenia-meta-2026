#!/usr/bin/env python3
# copy of analyses/de_prior_lists/our_prior/build_group_priors.py (RA 2026-10-23); the default HGNC path was replaced by the repository-relative data/reference/hgnc_complete_set.txt, nothing else changed
"""Per-group DE priors from the per-contrast calls of ../test_rnaseq_1000/build_test_prior.py (cache/calls/, one contrast per experiment).

    python3 build_group_priors.py        (needs numpy, pandas, scipy)

Membership comes from ../group_*.tsv, restricted to the experiments in experiments_all.tsv (public, not troubled, not single-cell).
A gene is in a group's file if at least one used contrast of the group tested it (n_tested >= 1), including genes with n_de = 0.

Columns follow the external files (../de_prior_ranks/<n>_<group>_de_prior.tsv):
  n_contrasts_used   contrasts of the group used, INCLUDING contrasts with no DE gene (the external file leaves those out)
  n_de_total, n_de_up, n_de_down   contrasts in which the gene is DE (FDR < 0.05, |log2FC| >= 1), and up / down (a gene with elements going
                     both ways counts in both, so up + down can exceed total)
  mf_score           sum over contrasts of hit / (N_in * N_out), N_in = genes DE in the contrast, N_out = genes tested and not DE in it
                     (the external build's N_out also counts untested genes as not DE, because its hit matrix codes them 0)
  full_rank, full_norm   rank of mf_score over the group's genes, all contrasts at once; rank 1 = lowest score, ties get the average rank
                     (R rank() default; the external full_rank has average ranks such as 3.5; its de_prior_rank ties come from the minimum over held-out runs), full_norm = (rank - min) / (max - min)
  de_prior_rank, de_prior_norm   the minimum, over held-out contrasts, of the rank after that contrast is left out and the scores recomputed
                     (Crow et al. 2019, as the external readme describes); same tie rule and rescaling; sorted by de_prior_rank descending
  n_tested, frac_de (= n_de_total / n_tested), n_de_4fold (contrasts with |log2FC| >= 2), locus_group (HGNC)
"""
import csv
import gzip
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.stats import rankdata, spearmanr

HERE = Path(__file__).resolve().parent
DL = HERE.parent
CACHE = Path(os.environ.get("CACHE", DL / "test_rnaseq_1000" / "cache"))
OUT = Path(os.environ.get("OUT", HERE))
HGNC = Path(os.environ.get("HGNC_FILE", "data/reference/hgnc_complete_set.txt"))
GROUPS = [("1_global", "group_1_global"), ("1_global_no_microarray", "group_1_global_no_microarray"),
          ("2_brain_postmortem", "group_2_brain_postmortem"), ("2_brain_postmortem_no_microarray", "group_2_brain_postmortem_no_microarray"),
          ("2b_brain_postmortem_no_sz_banks", "group_2b_brain_postmortem_no_sz_banks"),
          ("2b_brain_postmortem_no_sz_banks_no_microarray", "group_2b_brain_postmortem_no_sz_banks_no_microarray"),
          ("3_blood", "group_3_blood"), ("3_blood_no_microarray", "group_3_blood_no_microarray"),
          ("4_global_minus_brain_blood", "group_4_global_minus_brain_blood"),
          ("4_global_minus_brain_blood_no_microarray", "group_4_global_minus_brain_blood_no_microarray")]
NAMED = ["SMPD4", "WBP1", "TNIP2", "SNRPN", "SNURF", "XIST"]


def rd(fn):
    with open(fn, encoding="utf-8") as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def main():
    exps = rd(HERE / "experiments_all.tsv")
    fetched = {e["gemma_id"] for e in exps}
    calls = {}
    for e in exps:
        fp = CACHE / "calls" / f"{e['gemma_id']}.json.gz"
        if fp.exists():
            calls[e["gemma_id"]] = json.loads(gzip.decompress(fp.read_bytes()))
    lg, sym_h = {}, {}
    for r in rd(HGNC):
        if r["entrez_id"]:
            lg[r["entrez_id"]], sym_h[r["entrez_id"]] = r["locus_group"], r["symbol"]
    stats = {}
    for stem, gf in GROUPS:
        ids = [r["gemma_id"] for r in rd(DL / f"{gf}.tsv")]
        inlist = [i for i in ids if i in fetched]
        used = [i for i in inlist if i in calls and not calls[i]["skipped"]]
        st = {"experiments_in_group_file": len(ids), "not_in_experiments_all (private)": len(ids) - len(inlist), "experiments": len(inlist),
              "not_fetched": sum(1 for i in inlist if i not in calls), "skipped": {}, "contrasts_used": len(used)}
        for i in inlist:
            if i in calls and calls[i]["skipped"]:
                st["skipped"][calls[i]["skipped"]] = st["skipped"].get(calls[i]["skipped"], 0) + 1
        gi, gene_ids = {}, []
        for i in used:
            for nid in calls[i]["genes"]:
                if nid not in gi:
                    gi[nid] = len(gene_ids); gene_ids.append(nid)
        G = len(gene_ids)
        ntest, nde, nup, ndown, nde4, score = (np.zeros(G) for _ in range(6))
        sym = {}
        hits = []  # per contrast: (gene index array of DE genes, weight)
        zero_de = 0
        for i in used:
            gs = calls[i]["genes"]
            idx = np.fromiter((gi[n] for n in gs), int, len(gs))
            arr = np.array([v[1:] for v in gs.values()], int).reshape(-1, 4)
            up, down, up4, down4 = arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]
            de = (up | down).astype(bool)
            ntest[idx] += 1; nup[idx] += up; ndown[idx] += down; nde[idx] += de; nde4[idx] += (up4 | down4)
            for n, v in gs.items():
                sym.setdefault(n, v[0])
            n_in, n_out = int(de.sum()), len(gs) - int(de.sum())
            if n_in == 0:
                zero_de += 1
                continue
            w = 1.0 / (n_in * n_out) if n_out else 0.0
            score[idx[de]] += w
            hits.append((idx[de], w))
        st["contrasts_zero_de"] = zero_de
        st["genes"] = G
        full = rankdata(score) if G else score
        best = full.copy()
        for hidx, w in hits:  # leave the contrast out, recompute the scores, re-rank, keep each gene's minimum
            s2 = score.copy(); s2[hidx] -= w
            best = np.minimum(best, rankdata(s2))
        norm = lambda r: (r - r.min()) / (r.max() - r.min()) if G and r.max() > r.min() else np.zeros(G)
        order = np.lexsort((-nde, -best))  # rank descending, then n_de_total
        cols = ["ncbi_gene_id", "gene_symbol", "n_contrasts_used", "n_de_total", "n_de_up", "n_de_down", "mf_score", "full_rank", "full_norm",
                "de_prior_rank", "de_prior_norm", "n_tested", "frac_de", "n_de_4fold", "locus_group"]
        fn = OUT / f"{stem}_our_prior.tsv"
        with open(fn, "w", newline="", encoding="utf-8") as fh:
            w_ = csv.writer(fh, delimiter="\t"); w_.writerow(cols)
            fnorm, bnorm = norm(full), norm(best)
            for j in order:
                n = gene_ids[j]
                w_.writerow([n, sym_h.get(n) or sym[n], len(used), int(nde[j]), int(nup[j]), int(ndown[j]), f"{score[j]:.6g}", f"{full[j]:g}", f"{fnorm[j]:.6g}",
                             f"{best[j]:g}", f"{bnorm[j]:.6g}", int(ntest[j]), f"{nde[j] / ntest[j]:.6g}", int(nde4[j]), lg.get(n, "not in HGNC")])
        # comparison with the matching external file
        ext_fn = DL / "de_prior_ranks" / f"{stem}_de_prior.tsv"
        if ext_fn.exists() and G:
            e = {r["ncbi_gene_id"]: r for r in rd(ext_fn)}
            shared = [j for j in range(G) if gene_ids[j] in e]
            ex = lambda k: np.array([float(e[gene_ids[j]][k]) for j in shared])
            sp = lambda a, b: float(spearmanr(a, b)[0]) if len(shared) > 2 else None
            st["external_genes"] = len(e)
            st["shared_genes"] = len(shared)
            st["ours_outside_external"] = G - len(shared)
            st["external_not_in_ours"] = len(e) - len(shared)
            st["spearman_n_de_total"] = sp(nde[shared], ex("n_de_total"))
            st["spearman_n_de_4fold_vs_external_n_de_total"] = sp(nde4[shared], ex("n_de_total"))
            st["spearman_frac_de_vs_external_n_de_total"] = sp((nde / np.maximum(ntest, 1))[shared], ex("n_de_total"))
            st["spearman_de_prior_rank"] = sp(best[shared], ex("de_prior_rank"))
            st["spearman_full_rank"] = sp(full[shared], ex("full_rank"))
            st["external_n_contrasts_used"] = int(float(next(iter(e.values()))["n_contrasts_used"]))
        st["locus_groups"] = {}
        for j in range(G):
            lgp = lg.get(gene_ids[j], "not in HGNC"); st["locus_groups"][lgp] = st["locus_groups"].get(lgp, 0) + 1
        st["named"] = {s: None for s in NAMED}
        for j in range(G):
            s = sym_h.get(gene_ids[j]) or sym[gene_ids[j]]
            if s in st["named"]:
                st["named"][s] = {"n_tested": int(ntest[j]), "n_de": int(nde[j]), "n_de_4fold": int(nde4[j]), "de_prior_rank": float(best[j]),
                                  "external_n_de_total": int(float(e[gene_ids[j]]["n_de_total"])) if ext_fn.exists() and gene_ids[j] in e else None}
        stats[stem] = st
        print(stem, st["contrasts_used"], G, flush=True)
    json.dump(stats, open(OUT / "group_stats.json", "w"), indent=1)


if __name__ == "__main__":
    main()
