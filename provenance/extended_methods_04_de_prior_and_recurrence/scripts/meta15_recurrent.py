"""P10 / META15: how often the core-80 recurrent sets (383 study tier, 107 in >=10 lists) appear in the 15 brain meta-analysis lists.
Adapted from figure_data/mab_numbers_2026-10-12/task_c_meta_agreement.py (lists, k-way counts, permutation) and h3.py expected_deprior
(DE-prior null: each list's n genes shared out over its universe with weight n_de_total + 1, median for genes missing from the prior,
inclusion probability capped at 1 by recurrence.prior_probs; here the list is DRAWN as independent Bernoulli(p_g), so list size varies
around n; mean size = n).  Uniform null: n genes without replacement from the same universe.
All 15 units are run jointly (15 independent draws per replicate); variants drop two units (see units_rule table).
Output paths are in the private folder only."""
import csv, hashlib, os, sys, re, json, datetime
import numpy as np

OUT = os.path.dirname(os.path.abspath(__file__)) + "/out"
os.makedirs(OUT, exist_ok=True)
NPERM, SEED = 5000, 20261016
REPO = "sztr-rebuild-infra/repo"
AN = REPO + "/analyses"
WS = "base78_2026-10-16/B80/ws"
B80 = "base78_2026-10-16/B80/out"
MAIL = "recurrence_on_80_studies"
UNITS, INV = AN + "/h3/h3_units.tsv", AN + "/h3/h3_inventory.tsv"
SETS = B80 + "/recurrent_sets.tsv"
DEP = AN + "/de_prior_lists/our_prior/1_global_our_prior.tsv"
HGNC = REPO + "/data/reference/hgnc_complete_set.txt"
PRIMARY, UNIJ = WS + "/analyses/primary_lists.json", WS + "/analyses/universes.json"
BASE80 = MAIL + "/base80_studies.tsv"
GENEDIFF = MAIL + "/set_gene_differences.tsv"
SHAS = {}

def sha(p):
    if p not in SHAS:
        SHAS[p] = hashlib.sha256(open(p, "rb").read()).hexdigest()
    return SHAS[p]

def read_tsv(p):
    return list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t"))

def prior_probs(total, w):  # copy of recurrence.prior_probs
    p = np.zeros(len(w)); free = np.ones(len(w), bool); rem = total
    while True:
        p[free] = rem * w[free] / w[free].sum()
        over = free & (p > 1)
        if not over.any():
            return p
        p[over] = 1.0; free &= ~over; rem = total - p[~free].sum()

def de_prior_weights():  # copy of recurrence.de_prior_table (counts only)
    ent, approved = {}, set()
    for r in csv.DictReader(open(HGNC, encoding="utf-8"), delimiter="\t"):
        if r["status"] == "Approved":
            approved.add(r["symbol"])
            if r["entrez_id"]: ent[r["entrez_id"]] = r["symbol"]
    out = {}
    for r in read_tsv(DEP):
        g = ent.get((r.get("ncbi_gene_id") or "").strip())
        if not g:
            s = (r.get("gene_symbol") or "").strip(); g = s if s in approved else None
        if g:
            c = int(float(r["n_de_total"]))
            if g not in out or c > out[g]: out[g] = c
    return out

# ---- core-80 brain universe (union of tested sets of the 72 non-empty base lists) for units with no own universe
p = json.load(open(PRIMARY))["studies"]; uj = json.load(open(UNIJ))["lists"]
b80 = read_tsv(BASE80)
assert len(b80) == 80 and sum(r["empty_list_dropped"] == "True" for r in b80) == 8
UB = set(); ub_src = []
for r in b80:
    if r["empty_list_dropped"] == "True": continue
    for l in p[r["pmid"]]["primary"]:
        f = WS + "/data/" + uj[l]["gene_set"]
        UB |= {x.strip() for x in open(f) if x.strip() and not x.startswith("#")}
        ub_src.append(f)
assert len(UB) > 30000

units = [r for r in read_tsv(UNITS) if r["tissue_set"] == "brain"]
inv = {r["h3_unit"]: r for r in read_tsv(INV) if r["h3_tissue_set"] == "brain" and r["in_pairs"] == "yes" and r["h3_class"].startswith("cross-cohort")}
assert len(units) == 15
S383 = {r["gene"] for r in read_tsv(SETS) if r["set"] == "study_brain"}
S107 = {r["gene"] for r in read_tsv(SETS) if r["set"] == "study_brain_ge10"}
assert len(S383) == 383 and len(S107) == 107 and S107 <= S383 | S107
W = de_prior_weights()

LABEL = {u["unit"]: u["label"] for u in units}
ids = [u["unit"] for u in units]
info = []
srcfiles = []
for u in units:
    r = inv[u["unit"]]
    genes = set()
    for lid in r["h3_list"].split(";"):
        f = REPO + "/data/lists/" + lid + ".tsv"; srcfiles.append(f)
        genes |= {x["symbol"] for x in read_tsv(f) if x["symbol"]}
    assert len(genes) == int(u["list_genes"]), (u["label"], len(genes), u["list_genes"])
    m = re.search(r"\((universe_sets/[^)]+)\)", r["universe_rule"])
    if m and "_OR_" not in m.group(1) or (m and os.path.exists(WS + "/data/" + m.group(1))):
        f = WS + "/data/" + m.group(1); srcfiles.append(f)
        U = {x.strip() for x in open(f) if x.strip() and not x.startswith("#")}
        usrc = "own tested universe: " + m.group(1)
    else:
        U = set(UB); usrc = "FALLBACK core-80 brain union universe (no own universe: %s)" % r["universe_rule"][:60]
    L = genes & U
    info.append(dict(unit=u["unit"], label=u["label"], tech=u["technology"], genes=genes, L=L, U=U, usrc=usrc, own=usrc.startswith("own"),
                     ncoh=int(r["n_input_cohorts"]), ndat=int(r["n_input_datasets"]), role=r["h3_list_role"]))

# ---- units rule table
A_DROP = {"30545856", "40796563"}      # <4 input cohorts
B_DROP = {"40796563", "doi_10.64898_2026.09.18.752621"}   # single-nucleus units
def crit(i):
    return dict(ge4_cohorts="yes" if i["ncoh"] >= 4 else "no", bulk_tissue="no" if i["tech"] == "snrnaseq" else "yes",
                list_role_first_look="yes" if i["role"] == "first_look" else "no",
                list_genes_ge10="yes" if len(i["genes"]) >= 10 else "no", peer_reviewed="no" if "bioRxiv" in i["label"] else "yes")
rule_rows = []
for i in info:
    c = crit(i)
    rule_rows.append(dict(unit=i["unit"], label=i["label"], technology=i["tech"], n_input_datasets=i["ndat"], n_input_cohorts=i["ncoh"],
                          h3_list_role=i["role"], list_genes=len(i["genes"]), own_universe="yes" if i["own"] else "no", **c,
                          dropped_in_variant_A="yes" if i["unit"] in A_DROP else "no", dropped_in_variant_B="yes" if i["unit"] in B_DROP else "no"))
for k in ("ge4_cohorts", "bulk_tissue", "list_role_first_look", "list_genes_ge10", "peer_reviewed"):
    print("rule", k, "keeps", sum(r[k] == "yes" for r in rule_rows))

# ---- global gene index
allg = sorted(set().union(*[i["U"] for i in info]) | UB | S383 | S107 | set().union(*[i["genes"] for i in info]))
gidx = {g: k for k, g in enumerate(allg)}
G = len(allg)
m383 = np.zeros(G, bool); m107 = np.zeros(G, bool)
for g in S383: m383[gidx[g]] = True
for g in S107: m107[gidx[g]] = True
for i in info:
    i["idx"] = np.array(sorted(gidx[g] for g in i["U"]))
    i["n"] = len(i["L"])
    w = np.array([W.get(g, np.nan) + 1.0 for g in sorted(i["U"], key=lambda g: gidx[g])])
    w = np.where(np.isnan(w), np.nanmedian(w), w)
    i["p"] = prior_probs(float(i["n"]), w)
    i["nS383"] = int(m383[i["idx"]].sum()); i["nS107"] = int(m107[i["idx"]].sum())
    i["obs383"] = sum(g in S383 for g in i["L"]); i["obs107"] = sum(g in S107 for g in i["L"])
    i["exp_u383"] = i["n"] * i["nS383"] / len(i["U"]); i["exp_u107"] = i["n"] * i["nS107"] / len(i["U"])
    i["exp_p383"] = float(i["p"][m383[i["idx"]]].sum()); i["exp_p107"] = float(i["p"][m107[i["idx"]]].sum())

# ---- observed counts
cnt = {}
for i in info:
    for g in i["L"]: cnt[g] = cnt.get(g, 0) + 1
cnt_unrestricted = {}
for i in info:
    for g in i["genes"]: cnt_unrestricted[g] = cnt_unrestricted.get(g, 0) + 1
def variant_idx(v):
    if v == "all15": return list(range(15))
    drop = A_DROP if v == "drop_A_13" else B_DROP
    return [k for k, i in enumerate(info) if i["unit"] not in drop]
VARS = ["all15", "drop_A_13", "drop_B_13"]
def obs_counts(v):
    c = {}
    for k in variant_idx(v):
        for g in info[k]["L"]: c[g] = c.get(g, 0) + 1
    return c
OBSK = {v: obs_counts(v) for v in VARS}

# ---- simulation (joint draws, 15 lists per replicate)
rng = np.random.default_rng(SEED)
res = {}
for null in ("uniform", "deprior"):
    hits = np.zeros((NPERM, 15, 2), int)
    ge = {v: np.zeros((NPERM, 16, 3), int) for v in VARS}   # [perm, k, (all, in383, in107)]
    for b in range(NPERM):
        c = np.zeros(G, np.int16)
        percell = []
        for k, i in enumerate(info):
            if null == "uniform":
                sel = i["idx"][np.argpartition(rng.random(len(i["idx"])), i["n"] - 1)[:i["n"]]] if i["n"] else i["idx"][:0]
            else:
                sel = i["idx"][rng.random(len(i["idx"])) < i["p"]]
            hits[b, k, 0] = m383[sel].sum(); hits[b, k, 1] = m107[sel].sum()
            ind = np.zeros(G, np.int16); ind[sel] = 1
            percell.append(ind)
        for v in VARS:
            c = np.sum([percell[k] for k in variant_idx(v)], axis=0)
            for t in range(1, 16):
                m = c >= t
                ge[v][b, t, 0] = m.sum(); ge[v][b, t, 1] = (m & m383).sum(); ge[v][b, t, 2] = (m & m107).sum()
    res[null] = dict(hits=hits, ge=ge)
    print("done", null, flush=True)

def summ(x):
    return x.mean(), np.percentile(x, 2.5), np.percentile(x, 97.5)
def f(x, d=2): return f"{x:.{d}f}"
def pge(sim, o): return (1 + (sim >= o).sum()) / (len(sim) + 1)
def pse(sim, o): return (1 + (sim <= o).sum()) / (len(sim) + 1)

per_unit, pooled, byk, genes = [], [], [], []
for si, (sname, Sset, key) in enumerate((("study_brain_383", S383, "383"), ("study_brain_ge10_107", S107, "107"))):
    for k, i in enumerate(info):
        o = i["obs" + key]
        row = dict(meta_analysis=i["label"], unit=i["unit"], technology=i["tech"], set=sname, list_genes=len(i["genes"]), list_genes_in_universe=i["n"],
                   universe_used=i["usrc"], universe_n=len(i["U"]), set_genes_in_universe=i["nS" + key], observed_hits=o,
                   share_of_list_in_set=f(o / i["n"], 4) if i["n"] else "")
        for null in ("uniform", "deprior"):
            s = res[null]["hits"][:, k, si]; m, lo, hi = summ(s)
            row[f"exp_{null}"] = f(i[("exp_u" if null == "uniform" else "exp_p") + key], 3)
            row[f"sim_mean_{null}"] = f(m, 3); row[f"lo95_{null}"] = f(lo, 1); row[f"hi95_{null}"] = f(hi, 1)
            row[f"fold_{null}"] = f(o / i[("exp_u" if null == "uniform" else "exp_p") + key], 2) if i[("exp_u" if null == "uniform" else "exp_p") + key] > 0 else ""
            row[f"P_ge_obs_{null}"] = f"{pge(s, o):.4f}"
        per_unit.append(row)
    for v in VARS:
        vi = variant_idx(v)
        o = sum(info[k]["obs" + key] for k in vi)
        row = dict(variant=v, n_meta_analyses=len(vi), dropped=";".join(info[k]["label"] for k in range(15) if k not in vi) or "none", set=sname,
                   observed_pooled_hits=o, sum_list_genes_in_universe=sum(info[k]["n"] for k in vi))
        for null in ("uniform", "deprior"):
            s = res[null]["hits"][:, vi, si].sum(1); m, lo, hi = summ(s)
            e = sum(info[k][("exp_u" if null == "uniform" else "exp_p") + key] for k in vi)
            row[f"exp_{null}"] = f(e, 2); row[f"sim_mean_{null}"] = f(m, 2); row[f"lo95_{null}"] = f(lo, 1); row[f"hi95_{null}"] = f(hi, 1)
            row[f"fold_{null}"] = f(o / e, 2); row[f"P_ge_obs_{null}"] = f"{pge(s, o):.4f}"
        pooled.append(row)
for v in VARS:
    ok = OBSK[v]
    for t in range(1, 16):
        gs = [g for g, c in ok.items() if c >= t]
        row = dict(variant=v, k_min=t, genes_in_k_or_more=len(gs), in_383=sum(g in S383 for g in gs), in_107=sum(g in S107 for g in gs))
        for null in ("uniform", "deprior"):
            for j, nm in enumerate(("all", "in383", "in107")):
                s = res[null]["ge"][v][:, t, j]; m, lo, hi = summ(s)
                o = (len(gs), row["in_383"], row["in_107"])[j]
                row[f"{nm}_exp_{null}"] = f(m, 2); row[f"{nm}_lo95_{null}"] = f(lo, 1); row[f"{nm}_hi95_{null}"] = f(hi, 1); row[f"{nm}_P_ge_obs_{null}"] = f"{pge(s, o):.4f}"
        byk.append(row)
allc = OBSK["all15"]; A = OBSK["drop_A_13"]; B = OBSK["drop_B_13"]
for g in sorted({g for g in allc if allc[g] >= 3} | {g for g in A if A[g] >= 3} | {g for g in B if B[g] >= 3}, key=lambda g: (-allc.get(g, 0), g)):
    genes.append(dict(gene=g, k_all15=allc.get(g, 0), k_drop_A_13=A.get(g, 0), k_drop_B_13=B.get(g, 0), in_383="yes" if g in S383 else "no", in_107="yes" if g in S107 else "no",
                      meta_analyses=";".join(i["label"] for i in info if g in i["L"])))
ge5_unrestricted = sum(1 for v in cnt_unrestricted.values() if v >= 5)

# ---- write
TODAY = "2026-10-16"
srcs = [UNITS, INV, SETS, MAIL + "/set_gene_differences.tsv", BASE80, DEP, HGNC, PRIMARY, UNIJ] + sorted(set(srcfiles)) + sorted(set(ub_src))
def write(name, rows, claim, extra=()):
    head = [f"# figure/claim: {claim}", "# script: meta15_recurrent.py (adapted from figure_data/mab_numbers_2026-10-12/task_c_meta_agreement.py; DE prior as in h3.py expected_deprior / recurrence.prior_probs)",
            "# sources (sha256): " + "; ".join(f"{os.path.relpath(s, '$HOME')} {sha(s)}" for s in srcs[:9]) + "; plus per-list files listed in sources_sha256.tsv",
            f"# date: {TODAY}", f"# parameters: NPERM {NPERM}, seed {SEED}, base = 80 studies / 72 contributing lists; sets from core-80 recurrent_sets.tsv (383 study tier q<0.1; 107 in >=10 lists)"] + [f"# {e}" for e in extra]
    with open(OUT + "/" + name, "w", newline="") as fh:
        fh.write("\n".join(head) + "\n")
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
nullnote = ["uniform null: n genes (the list's in-universe size) drawn without replacement from the unit's universe; DE-prior null: each gene included independently with the h3/recurrence probability (weight n_de_total + 1, median for genes missing from the prior, capped at 1, sum = n); exp_* = analytic expectation, sim_mean = simulated mean, lo95/hi95 = 2.5/97.5 percentiles of the simulation, P_ge_obs = (1 + #sim >= observed)/(NPERM + 1)",
            "universe: the unit's own tested universe (universes.json gene set) where it has one; FALLBACK units use the union of the tested sets of the 72 non-empty base-80 lists (%d genes)" % len(UB),
            "DE prior file: analyses/de_prior_lists/our_prior/1_global_our_prior.tsv (default 'ours' of recurrence.py, same file the B80 run used)"]
write("meta15_per_meta_analysis_hits.tsv", per_unit, "META15: per brain meta-analysis, genes of the 383 / 107 core-80 recurrent sets on its list, observed vs expected under the uniform and DE-prior nulls", nullnote)
write("meta15_pooled_hits.tsv", pooled, "META15: hits pooled over the meta-analyses (sum of per-list hits; a gene on several lists counts once per list), observed vs expected under both nulls; all 15 and two 13-unit variants", nullnote + ["variant drop_A_13 = without Gandal 2018 (RNA-seq) and Nakatsuka 2025 (rule: fewer than 4 input cohorts); drop_B_13 = without Nakatsuka 2025 and Endresz 2026 (rule: single-nucleus units). No rule is declared the one; see meta15_units_rule.tsv"])
write("meta15_genes_by_k.tsv", byk, "META15: genes on k or more of the meta-analysis lists (each meta-analysis once; up and down pooled), how many are in the 383 / 107, and the joint-draw expectation under both nulls", nullnote + ["all15 k>=5 is the Note 5 '147 genes in 5 or more of the 15 lists' (earlier version, 209-era base: 147, of which 68 in the earlier 375-gene set)", f"observed counts use each list restricted to its unit's universe; without the restriction genes in >=5 lists = {ge5_unrestricted}"])
write("meta15_genes_k_ge3.tsv", genes, "META15: genes on 3 or more of the meta-analysis lists with their membership of the 383 / 107", nullnote)
write("meta15_units_rule.tsv", rule_rows, "META15: which of the 15 units a rule from the unit definitions (h3_units / h3_inventory) would drop to give 13", ["criteria columns are read from h3_units.tsv / h3_inventory.tsv; all 15 are cross-cohort, pool 2 or more cohorts, publish a list and are in_pairs = yes, so those criteria drop none"])
with open(OUT + "/sources_sha256.tsv", "w") as fh:
    fh.write("# P10 META15 input files, sha256 (full)\nsha256\tpath\n")
    for s in srcs: fh.write(f"{sha(s)}\t{s}\n")
print("universe fallback n", len(UB), "ge5 unrestricted", ge5_unrestricted)
for r in pooled: print(r)
for r in byk:
    if r["k_min"] in (5, 10): print(r["variant"], r["k_min"], r["genes_in_k_or_more"], r["in_383"], r["in_107"], r["all_exp_uniform"], r["all_exp_deprior"], r["in383_exp_uniform"], r["in383_exp_deprior"], r["in107_exp_deprior"])
