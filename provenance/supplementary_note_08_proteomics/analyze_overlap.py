#!/usr/bin/env python3
"""Step 2 (2026-10-14): protein hits among the 209 RNA recurrent genes: recompute the side study's overlap, add a
matched-random-gene null (expression decile x DE-prior quintile), direction vs the RNA call, sensitivity to the universe.
Read-only on the project repository. Run: PYTHONDONTWRITEBYTECODE=1 python3 -I analyze_overlap.py"""
import sys, os, csv, json, hashlib
import numpy as np, pandas as pd
from scipy.stats import hypergeom, binomtest
REPO = os.path.expanduser("repo")
OUT = os.environ["OUTDIR"]; W = OUT + "/work"
sys.path.insert(0, REPO + "/analyses"); sys.path.insert(0, REPO + "/analyses/recurrence_sources"); sys.path.insert(0, REPO + "/side-studies/proteomics")
import analyse as A, recurrence as R, extract
extract.HGNC = extract.load_hgnc()
rng = np.random.default_rng(20261014)
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
REC = os.environ.get("REC_OVERRIDE", "core80_work_2026-10-16/P5/M03/inputs_new/recurrence_genes_final.tsv"); TAG = os.environ.get("TAG", ""); EXPECT_N = int(os.environ.get("EXPECT_N", "0"))
# PATCH: set definition from the final recurrent_sets.tsv (env SETNAME); without SETNAME the original q_bh < 0.1 rule on REC is used (old 176-gene file)
SETNAME = os.environ.get("SETNAME")
SETS_TSV = os.environ.get("SETS_TSV", "core80_work_2026-10-16/P5/M03/inputs_new/sets_final.tsv")
SETGENES = set(pd.read_csv(SETS_TSV, sep="\t", comment="#").query("set == @SETNAME").gene) if SETNAME else None
def in_set(df):
    return df.gene.isin(SETGENES) if SETNAME else df.q_bh < 0.1
rg = pd.DataFrame([r for r in csv.DictReader((l for l in open(REC) if not l.startswith("#")), delimiter="\t") if r["scenario"] == "published, brain"])
for c in ["q_bh", "de_prior_hitlists", "up", "down", "expr_decile"]:
    rg[c] = pd.to_numeric(rg[c], errors="coerce")
rg["in209"] = in_set(rg)
rg["rna_dir"] = np.where(rg.up > rg.down, "up", np.where(rg.down > rg.up, "down", "tie"))
med = rg.de_prior_hitlists.median(); rg["prior"] = rg.de_prior_hitlists.fillna(med)
rg["pq"] = pd.qcut(rg.prior.rank(method="first"), 5, labels=False)
rg["stratum"] = rg.expr_decile.astype(int).astype(str) + "_" + rg.pq.astype(str)
S209 = set(rg.gene[rg.in209]); RNAG = set(rg.gene); strat = dict(zip(rg.gene, rg.stratum)); rdir = dict(zip(rg.gene, rg.rna_dir))
assert (EXPECT_N == 0 and (SETNAME is not None)) or len(S209) == EXPECT_N, len(S209)  # PATCH: old assertion len == 209; size is now the size of the named set
NSET = len(S209)
# units
M = {pm: R.read_set("../side-studies/proteomics/universe_sets/%s.txt" % pm) if False else {l.strip() for l in open(REPO + "/side-studies/proteomics/universe_sets/%s.txt" % pm) if l.strip() and not l.startswith("#")} for pm in ["22248022","25048004","27236410","31055969","31642882","32639965","35354897","41611717","42748139"]}
L = {pm: A.read_list(pm, A.UNITS[pm][1]) for pm in A.UNITS}
vm = pd.read_csv(W + "/protein_table__Vera-Montecinos_2021.tsv", sep="\t")
M["34576238"] = set(vm.symbol); L["34576238"] = R.read_list("34576238__dataset_s2_fdr01")
LAB = {pm: A.UNITS[pm][0] for pm in A.UNITS}; LAB["34576238"] = "Vera-Montecinos 2021"
BANK = {pm: A.UNITS[pm][4] for pm in A.UNITS}; BANK["34576238"] = "barcelona"
order = ["22248022","25048004","27236410","31055969","31642882","32639965","35354897","41611717","42748139","34576238"]
def pooled_dist(terms):
    d = np.array([1.0])
    for N, K, n in terms:
        k = np.arange(0, min(K, n) + 1); d = np.convolve(d, hypergeom.pmf(k, N, K, n))
    return d
def tail(d, o): return float(d[o:].sum())
def qint(d, lo=.025, hi=.975):
    c = np.cumsum(d); return int(np.searchsorted(c, lo)), int(np.searchsorted(c, hi))
def unit_terms(pm, mode):
    """hypergeometric terms for one unit. mode: 'plain' (project: universe = the study's tested set), 'matched' (strata), 'allgenes' (universe = all RNA genes), 'allgenes_matched'."""
    U = M[pm]; Ls = set(L[pm]) & U
    if mode == "plain":
        return [(len(U), len(S209 & U), len(Ls))]
    if mode == "allgenes":
        return [(len(RNAG), NSET, len(Ls & RNAG))]  # PATCH: 209 -> set size
    Ls = Ls & RNAG; U = (U if mode == "matched" else RNAG) & RNAG
    out = []
    for s in sorted({strat[g] for g in Ls}):
        Ns = sum(1 for g in U if strat[g] == s); Ks = sum(1 for g in U if strat[g] == s and g in S209); m = sum(1 for g in Ls if strat[g] == s)
        out.append((Ns, Ks, m))
    return out
def hits(pm):
    return sorted((set(L[pm]) & M[pm]) & S209)
def dirs(pm): return {g: v for g, v in L[pm].items() if v in ("up", "down")}
def agree_stats(pm, nperm=20000):
    d = dirs(pm); H = [g for g in hits(pm) if g in d and rdir[g] in ("up", "down")]
    a = sum(d[g] == rdir[g] for g in H)
    genes = list(d); labs = np.array([d[g] == "up" for g in genes]); idx = {g: i for i, g in enumerate(genes)}
    hi = np.array([idx[g] for g in H]); rup = np.array([rdir[g] == "up" for g in H])
    sims = np.empty(nperm, dtype=int)
    for i in range(nperm):
        p = rng.permutation(labs); sims[i] = int((p[hi] == rup).sum()) if len(H) else 0
    return len(H), a, sims
rows = []; sim_pool = {m: np.zeros(20000, dtype=int) for m in ["all"]}
res = {}
for pm in order:
    r = {"unit": LAB[pm], "pmid": pm, "bank": BANK[pm], "universe_tested": len(M[pm]), "209_in_universe": len(S209 & M[pm]), "list_proteins_in_universe": len(set(L[pm]) & M[pm]),
         "hits_in_209": len(hits(pm))}
    for mode in ["plain", "matched", "allgenes", "allgenes_matched"]:
        T = unit_terms(pm, mode); d = pooled_dist(T); E = float((d * np.arange(len(d))).sum())
        lo_, hi_ = qint(d); r[f"null95_{mode}"] = f"{lo_} to {hi_}"
        r[f"expected_{mode}"] = round(E, 3); r[f"p_{mode}"] = tail(d, r["hits_in_209"]) if mode != "allgenes" and mode != "allgenes_matched" else tail(d, len(set(L[pm]) & RNAG & S209))
    n, a, sims = agree_stats(pm)
    r.update({"hits_with_direction": n, "same_direction_as_rna": a, "expected_same_direction_shuffle": round(float(sims.mean()), 2),
              "p_direction_shuffle": float((np.sum(sims >= a) + 1) / (len(sims) + 1)),
              "p_direction_vs_half": float(binomtest(a, n, .5).pvalue) if n else np.nan})
    rows.append(r); res[pm] = (n, a, sims)
df = pd.DataFrame(rows)
# pooled sets
def pooled(sel, label):
    out = {"unit": label, "units_n": len(sel)}
    out["hits_in_209"] = int(sum(hits(pm).__len__() for pm in sel))
    for mode in ["plain", "matched", "allgenes", "allgenes_matched"]:
        T = [t for pm in sel for t in unit_terms(pm, mode)]; d = pooled_dist(T); E = float((d * np.arange(len(d))).sum())
        o = out["hits_in_209"] if mode in ("plain", "matched") else sum(len(set(L[pm]) & RNAG & S209) for pm in sel)
        lo, hi = qint(d)
        out[f"expected_{mode}"] = round(E, 2); out[f"null95_{mode}"] = f"{lo} to {hi}"; out[f"p_{mode}"] = tail(d, o)
        out[f"ratio_{mode}"] = round(o / E, 2) if E else np.nan
    n = sum(res[pm][0] for pm in sel); a = sum(res[pm][1] for pm in sel); sims = sum(res[pm][2] for pm in sel)
    out.update({"hits_with_direction": n, "same_direction_as_rna": a, "expected_same_direction_shuffle": round(float(sims.mean()), 2),
                "null95_direction": "%d to %d" % tuple(np.quantile(sims, [.025, .975])), "p_direction_shuffle": float((np.sum(sims >= a) + 1) / (len(sims) + 1)),
                "p_direction_vs_half": float(binomtest(a, n, .5).pvalue) if n else np.nan})
    return out
P = []
nine = order[:9]
P.append(pooled(nine, "pooled, 9 units (the side study's 9 units with a tested set)"))
P.append(pooled(order, "pooled, 10 units (9 + Vera-Montecinos 2021)"))
P.append(pooled([p for p in nine if p != "42748139"], "9 units without Testa 2026 (RNA used as a prior)"))
P.append(pooled([p for p in nine if BANK[p] != "pitt"], "9 units without the Pittsburgh group (MacDonald 2020, Nucifora 2019, Testa 2026)"))
P.append(pooled([p for p in nine if p != "32639965"], "9 units without Vidal-Domenech (same donors as Pinacho)"))
P.append(pooled([p for p in nine if p not in ("42748139",) and BANK[p] != "pitt"], "9 units without Testa and the Pittsburgh group"))
P.append(pooled(["27236410","35354897","25048004","41611717","42748139"], "one unit per donor group (Pinacho, Alnafisah, Focking PSD, Koopmans, Testa; Krishnamurthy, Nucifora, MacDonald, Vidal left out)"))
P.append(pooled(["27236410","35354897","25048004","41611717","31642882"], "one unit per donor group, MacDonald 2020 for Pittsburgh (no Testa)"))
for pm in nine:
    P.append(pooled([p for p in nine if p != pm], "leave out " + LAB[pm]))
pdf = pd.DataFrame(P)
def hdr(extra=""):
    return ["# built 2026-10-14 by scripts/analyze_overlap.py; RNA set = recurrence_genes.tsv (sha256 %s) scenario 'published, brain', q_bh < 0.1 = %d genes; RNA direction = majority of the gene's up/down lists" % (sha(REC)[:12], len(S209)),
            "# plain = project null (list proteins drawn from the study's tested set, hypergeometric, convolved over units); matched = each list protein drawn from the same expression-decile x DE-prior-quintile stratum of the tested set (exact stratified hypergeometric; expression decile and de_prior_hitlists from recurrence_genes.tsv, 231 genes without a prior count given the median); allgenes = universe is the 15,229 genes in the RNA scenario; allgenes_matched = stratified over all of them",
            "# direction null: the unit's protein up/down labels shuffled among its own directed list proteins (20,000x, seed 20261014), count of hits whose label equals the RNA call; ties in the RNA call excluded" + extra]
def w(fn, d, h):
    with open((f"{W}/{TAG}_{fn}" if TAG else f"{OUT}/{fn}"), "w") as f:
        for l in h: f.write(l + "\n")
        d.to_csv(f, sep="\t", index=False)
w("overlap_209_by_unit.tsv", df, hdr()); w("overlap_209_pooled_and_sensitivity.tsv", pdf, hdr())
# gene-level hits table
gl = []
for pm in order:
    d = dirs(pm)
    for g in hits(pm):
        gl.append((LAB[pm], g, L[pm][g], rdir[g], int(rg.set_index("gene").loc[g, "up"]), int(rg.set_index("gene").loc[g, "down"]), strat[g]))
gdf = pd.DataFrame(gl, columns=["unit", "gene", "protein_direction", "rna_direction", "rna_lists_up", "rna_lists_down", "stratum_expr_prior"])
w("overlap_209_hit_genes.tsv", gdf, hdr()[:1])
meta = {"recurrence_genes_sha256": sha(REC), "S209": len(S209), "units": {LAB[p]: {"universe": len(M[p]), "list": len(L[p])} for p in order}}
json.dump(meta, open(W + "/overlap_meta.json", "w"), indent=1)
pd.set_option("display.width", 250, "display.max_columns", 40)
print(df[["unit","universe_tested","209_in_universe","list_proteins_in_universe","hits_in_209","expected_plain","p_plain","expected_matched","p_matched","expected_allgenes","expected_allgenes_matched","hits_with_direction","same_direction_as_rna","expected_same_direction_shuffle","p_direction_shuffle"]].to_string())
print(pdf[["unit","hits_in_209","expected_plain","p_plain","expected_matched","null95_matched","p_matched","ratio_matched","expected_allgenes","expected_allgenes_matched","hits_with_direction","same_direction_as_rna","expected_same_direction_shuffle","null95_direction","p_direction_shuffle"]].to_string())
