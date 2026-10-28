#!/usr/bin/env python3
"""Task 2: concordance of the LIV-PM statistic (Liharska ST1 moderated t; positive = higher in postmortem) with the recurrent SZ direction,
recomputed with leave-bank-out recurrence runs. Same definitions as liharska_living_vs_pm_2026-10-18/scripts/s1_main.py (383 row):
  SZ sign = sign(up - down) over lists; concordant = sign(t) == SZ sign; genes with sign 0 or no ST1 symbol are left out.
  Matched null = 10,000 sets from published-brain genes with a SZ sign outside the set, same expression decile x DE-prior quintile strata.
  Interval = percentile bootstrap over genes (2,000), as in s1_main.
Variants per run (runs/<name>/recurrence_genes.tsv, scenario 'published, brain'):
  (a) SIGNREDERIVED : set = the original 383 genes (full run); sign from the lists remaining in the run; genes with no remaining list or a tie dropped.
  (b) SETREBUILT    : set = genes of the run with q_bh < 0.1; sign from the run's lists.
Runs: full (validation: must give 296/374), no_bank9 (all nine bank-sharing studies out), drop_<pmid> (one study out). Seed 20261018.   python3 -I s2_recompute.py"""
import sys, os, csv, hashlib, random, collections, math
sys.dont_write_bytecode = True
import numpy as np, pandas as pd
OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.expanduser("<analysis repo>"); AN = f"{REPO}/analyses"; sys.path.insert(0, AN)
LV = os.path.expanduser("~/analysis-delivery/liharska_living_vs_pm_2026-10-18")
SEED, NSETS, NBOOT = 20261018, 10000, 2000
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
def tsv(p): return list(csv.DictReader((l for l in open(p, encoding="utf8") if not l.startswith("#")), delimiter="\t"))
import extract
extract.HGNC = extract.load_hgnc(); ENS2SYM = {e: s for e, (s, _) in extract.HGNC[3].items()}
l = pd.read_csv(f"{LV}/data/liv_pm_de_ST1.tsv", sep="\t", comment="#"); l["sym"] = l.gene.map(ENS2SYM)
l = l.dropna(subset=["sym"]).sort_values("P.Value").drop_duplicates("sym").set_index("sym"); T = l.t
NAMES = {"full": "none (validation)", "no_bank9": "all nine bank-sharing studies", "drop_24702465": "Pietersen 2014", "drop_24628518": "Pietersen 2014 (PV neurons)",
         "drop_34440415": "Lindholm Carlstrom 2021", "drop_38448582": "Ling 2024", "drop_38781388": "Ruzicka 2024", "drop_31288836": "Mendizabal 2019",
         "drop_19255580": "Maycox 2009", "drop_17822540": "Saetre 2007"}
def load(run):
    rows = [r for r in tsv(f"{OUT}/runs/{run}/recurrence_genes.tsv") if r["scenario"] == "published, brain"]
    d = {}
    for r in rows:
        dd = int(r["up"]) - int(r["down"])
        prior = float(r["de_prior_hitlists"]) if r["de_prior_hitlists"] != "" else np.nan
        d[r["gene"]] = dict(sign=(1 if dd > 0 else -1) if dd else 0, dec=r["expr_decile"], prior=prior, q=float(r["q_bh"]))
    med = float(np.nanmedian([v["prior"] for v in d.values()]))
    for v in d.values():
        if np.isnan(v["prior"]): v["prior"] = med
    return d
def analyse(rd, setgenes, seed):
    """rd: gene -> dict(sign, dec, prior); setgenes: list. Returns stats dict."""
    S = [g for g in setgenes if g in rd and rd[g]["sign"] != 0 and g in T.index]
    n = len(S); sg = np.array([rd[g]["sign"] for g in S]); tv = np.array([T[g] for g in S])
    conc = int(((np.sign(tv) == sg)).sum()); frac = conc / n
    rng = np.random.default_rng(seed); bs = [ (np.sign(tv[i]) == sg[i]).mean() for i in (rng.integers(0, n, n) for _ in range(NBOOT))]
    inset = set(setgenes)
    pool = [g for g, v in rd.items() if v["sign"] != 0 and g not in inset and g in T.index]
    qc = np.quantile([rd[g]["prior"] for g in pool], [0.2, 0.4, 0.6, 0.8])
    st = lambda g: (rd[g]["dec"], int(np.searchsorted(qc, rd[g]["prior"], side="right")))
    by = collections.defaultdict(list)
    for g in pool: by[st(g)].append(g)
    need = collections.Counter(st(g) for g in S if g in rd)
    pz = {g: (np.sign(T[g]) == rd[g]["sign"]) for g in pool}
    r = random.Random(seed); nulls = []
    for _ in range(NSETS):
        row = []
        for s, k in need.items():
            c = by[s]
            if not c: continue
            row += r.sample(c, k) if len(c) >= k else r.choices(c, k=k)
        nulls.append(np.mean([pz[g] for g in row]))
    nulls = np.array(nulls)
    return dict(n_set_genes=len(setgenes), n_kept=n, n_concordant=conc, frac=round(frac, 4), ci_lo=round(float(np.percentile(bs, 2.5)), 4), ci_hi=round(float(np.percentile(bs, 97.5)), 4),
                null_mean=round(float(nulls.mean()), 4), null_lo=round(float(np.percentile(nulls, 2.5)), 4), null_hi=round(float(np.percentile(nulls, 97.5)), 4),
                frac_minus_null=round(frac - float(nulls.mean()), 4), pool_genes=len(pool), strata_without_pool_gene=sum(1 for s in need if not by[s]))
full = load("full"); S383 = [g for g, v in full.items() if v["q"] < 0.1]; assert len(S383) == 383
rows = []
for run in NAMES:
    rd = load(run)
    for variant, sset in (("SIGNREDERIVED", S383), ("SETREBUILT", [g for g, v in rd.items() if v["q"] < 0.1])):
        if run == "full" and variant == "SETREBUILT": pass
        res = analyse(rd, sset, SEED)
        rows.append(dict(run=run, studies_removed=NAMES[run], variant=variant, **res)); print(rows[-1], flush=True)
        if run == "full" and variant == "SIGNREDERIVED":
            assert (res["n_kept"], res["n_concordant"]) == (374, 296), res
# overlap of rebuilt set with the original 383
extra = []
for run in NAMES:
    rd = load(run); R = {g for g, v in rd.items() if v["q"] < 0.1}
    extra.append((run, len(R), len(R & set(S383))))
hdr = ["# script: scripts/s2_recompute.py, run 2026-10-18; seed 20261018; 10,000 matched sets; bootstrap over genes 2,000 (percentile 95%)",
       f"# inputs: runs/*/recurrence_genes.tsv (recurrence_fork.py sha 044e2c702342 on the private copy of B80/ws on the remote environment; include file = primary_pmids.txt minus the studies named; validation run 'full' equals B80/out/recurrence_genes.tsv, sha dda5eb621c15, on all data rows); Liharska ST1 data/liv_pm_de_ST1.tsv in liharska_living_vs_pm_2026-10-18 (sha {sha(LV+'/data/liv_pm_de_ST1.tsv')})",
       "# concordant = sign(LIV-PM t) equals SZ sign (sign(up - down)); positive t = higher in postmortem. SIGNREDERIVED: original 383 genes, sign from remaining lists (genes with no remaining list or a tie dropped = n_set_genes - n_kept). SETREBUILT: genes with q_bh < 0.1 in the run.",
       "# the 'full' / SIGNREDERIVED row reproduces the original 296 of 374 (0.7914)"]
with open(f"{OUT}/concordance_leave_bank_out.tsv", "w") as fh:
    fh.write("\n".join(hdr) + "\n"); w = csv.DictWriter(fh, list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
with open(f"{OUT}/rebuilt_set_sizes.tsv", "w") as fh:
    fh.write("# script: s2_recompute.py; genes with q_bh<0.1 (published, brain) per run and overlap with the original 383\nrun\tstudies_removed\tn_q10_genes\tn_in_original_383\n" + "\n".join(f"{a}\t{NAMES[a]}\t{b}\t{c}" for a, b, c in extra) + "\n")
