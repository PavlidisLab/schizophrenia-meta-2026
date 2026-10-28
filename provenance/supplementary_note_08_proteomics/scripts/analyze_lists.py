#!/usr/bin/env python3
"""Step 4 (2026-10-14): direction agreement of every directed proteomic list with the top of each axis (Extended Fig. 4 test, both axes),
and protein vs RNA direction on the 209 genes against each axis. Read-only on the repository.
Run: PYTHONDONTWRITEBYTECODE=1 python3 -I analyze_lists.py"""
import sys, os, csv, hashlib, collections
import numpy as np, pandas as pd
from scipy.stats import binomtest
REPO = os.path.expanduser("sztr-rebuild-infra/repo")
OUT = os.environ["OUTDIR"]; W = OUT + "/work"
sys.path.insert(0, REPO + "/analyses"); sys.path.insert(0, REPO + "/analyses/recurrence_sources"); sys.path.insert(0, REPO + "/side-studies/proteomics")
import analyse as A, recurrence as R, extract
extract.HGNC = extract.load_hgnc()
rng = np.random.default_rng(20261014)
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
ax = pd.read_csv(W + "/axis_table_hgnc.tsv", sep="\t").set_index("sym")
cols = ["NSW", "Stanley", "Pitt", "BA10", "Pritzker"]
def stouffer(excl=None, minc=2):
    num = np.zeros(len(ax)); den = np.zeros(len(ax)); k = np.zeros(len(ax))
    for c in cols:
        if c == excl: continue
        t = ax["t_pH_" + c].values; n = ax["n_pH_" + c].values; m = ~np.isnan(t) & ~np.isnan(n)
        num += np.where(m, np.sqrt(np.nan_to_num(n)) * np.nan_to_num(t), 0); den += np.where(m, np.nan_to_num(n), 0); k += m
    return pd.Series(np.where((k >= minc) & (den > 0), num / np.sqrt(np.where(den > 0, den, 1)), np.nan), index=ax.index)
AX = pd.DataFrame({"pH": -ax.pH_axis, "agonal": ax.agonal_axis})
for c in ["NSW", "Stanley", "Pitt", "Pritzker"]: AX["pH_LOCO_" + c] = -stouffer(c)
BANK_COLL = {"stanley": "Stanley", "pitt": "Pitt", "nsw": "NSW", "aus": "NSW", "uci": "Pritzker"}
# units: project units + extras
units = {pm: {"label": u[0], "lists": u[1], "bank": u[4] or ("solo:" + pm), "class": u[2]} for pm, u in A.UNITS.items()}
extra = {"26232077": ("Saia-Cereda 2015", ["table2"], "solo:26232077", "MS proteome"), "19405953": ("Martins-de-Souza 2009 (STG)", ["additional_file2_table2"], "mannheim", "gel"),
         "37076706": ("Khan 2023", ["tables1_sn", "tables2_hippocampus", "tables3_pfc"], "solo:37076706", "gel"), "28214564": ("Hirayama-Kurogi 2017", ["table_s2"], "solo:28214564", "panel"),
         "34576238": ("Vera-Montecinos 2021", ["dataset_s2_fdr01"], "barcelona", "MS proteome")}
for pm, (l, ls, b, c) in extra.items(): units[pm] = {"label": l, "lists": ls, "bank": b, "class": c}
def udirs(pm):
    d = collections.defaultdict(list)
    for lid in units[pm]["lists"]:
        if pm in A.UNITS: lst = A.read_list(pm, [lid])
        else: lst = R.read_list(f"{pm}__{lid}")
        for g, v in lst.items():
            if v in ("up", "down"): d[g].append(v)
    out = {}
    for g, v in d.items():
        c = collections.Counter(v)
        if c["up"] != c["down"]: out[g] = "up" if c["up"] > c["down"] else "down"
    return out
UD = {pm: udirs(pm) for pm in units}
def gcall(pms):
    d = collections.defaultdict(list)
    for pm in pms:
        for g, v in UD[pm].items(): d[g].append(v)
    out = {}
    for g, v in d.items():
        c = collections.Counter(v)
        if c["up"] != c["down"]: out[g] = "up" if c["up"] > c["down"] else "down"
    return out
N_PERM = 5000
def topn(series, n):
    s = series.dropna(); s = s.reindex(s.abs().sort_values(ascending=False).index)[:n]
    return {g: ("up" if v > 0 else "down") for g, v in s.items()}
variants = [("pH", 200), ("pH", 500), ("agonal", 200), ("agonal", 500)]
def predicted(variant, bank):
    name, n = variant
    if name == "pH" and bank in BANK_COLL: return topn(AX["pH_LOCO_" + BANK_COLL[bank]], n), "pH_leaving_out_" + BANK_COLL[bank]
    return topn(AX[name], n), name
def run(pm_list, label):
    groups = collections.defaultdict(list)
    for pm in pm_list: groups[units[pm]["bank"]].append(pm)
    grp_rows = []
    tot = collections.defaultdict(lambda: [0, 0, np.zeros(N_PERM)])
    for b, pms in sorted(groups.items()):
        call = gcall(pms); genes = list(call); labs = np.array([call[g] == "up" for g in genes]); gi = {g: i for i, g in enumerate(genes)}
        for variant in variants:
            pred, nm = predicted(variant, b)
            H = [g for g in genes if g in pred]
            conc = sum(call[g] == pred[g] for g in H)
            hi = np.array([gi[g] for g in H], dtype=int); pu = np.array([pred[g] == "up" for g in H])
            sims = np.empty(N_PERM)
            for i in range(N_PERM):
                sims[i] = (rng.permutation(labs)[hi] == pu).sum() if len(H) else 0
            grp_rows.append({"set": label, "axis": nm, "top_n": variant[1], "group": b, "units": ";".join(units[p]["label"] for p in pms), "genes_called": len(call), "hits_in_top_n": len(H), "concordant": conc,
                             "share": round(conc / len(H), 3) if H else np.nan, "expected_share": round(float(sims.mean()) / len(H), 3) if H else np.nan})
            key = (variant[0] + ("_matched" if variant[0] == "pH" else ""), variant[1])
            t = tot[key]; t[0] += len(H); t[1] += conc; t[2] += sims
    G = pd.DataFrame(grp_rows); pool = []
    for key, (h, c, sims) in tot.items():
        if h == 0: continue
        sub = G[(G.top_n == key[1]) & (G.axis.str.startswith("pH") if key[0].startswith("pH") else G.axis == "agonal")]
        sub = sub[sub.hits_in_top_n > 0]
        full = c / h - sims.mean() / h
        jk = []
        for g in sub.group:
            s2 = sub[sub.group != g]
            if s2.hits_in_top_n.sum() == 0: continue
            jk.append(s2.concordant.sum() / s2.hits_in_top_n.sum() - (s2.expected_share * s2.hits_in_top_n).sum() / s2.hits_in_top_n.sum())
        jk = np.array(jk); se = np.sqrt((len(jk) - 1) / len(jk) * ((jk - jk.mean()) ** 2).sum()) if len(jk) > 1 else np.nan
        pool.append({"set": label, "axis": key[0], "top_n": key[1], "donor_groups_with_hits": len(sub), "hits": h, "concordant": c, "share": round(c / h, 3), "expected_share": round(sims.mean() / h, 3),
                     "excess": round(full, 3), "excess_ci_lo_jackknife_groups": round(full - 1.96 * se, 3), "excess_ci_hi_jackknife_groups": round(full + 1.96 * se, 3),
                     "p_one_sided_shuffle": (np.sum(sims >= c) + 1) / (N_PERM + 1), "p_binomial_vs_half_two_sided": binomtest(c, h, .5).pvalue})
    return G, pd.DataFrame(pool)
ALL = list(units); ORIG = list(A.UNITS)
cfg = [("all units (31 + 5 added)", ALL), ("the 31 units of the side study (13 donor groups, as Extended Fig. 4)", ORIG),
       ("all units without the Pittsburgh group", [p for p in ALL if units[p]["bank"] != "pitt"]),
       ("all units, Pittsburgh group without Testa 2026", [p for p in ALL if p != "42748139"])]
Gs, Ps = [], []
for lab, pl in cfg:
    G_, P_ = run(pl, lab); Gs.append(G_); Ps.append(P_); print(lab, flush=True)
G = pd.concat(Gs); PL = pd.concat(Ps)
H2 = ["built 2026-10-14 by scripts/analyze_lists.py; axes: axis_gene_table.tsv (sha256 %s); top_n = genes with the largest |axis value|; predicted SZ direction = the one the axis predicts (positive axis value = up); pH axis = the tool's pH axis, recomputed leaving out the collection that shares a bank with the group (Stanley, Pitt/Pittsburgh, NSW, Pritzker for the Irvine study) where one does; 'pH' rows without a match use all 5 collections" % sha(W + "/axis_gene_table.tsv")[:12],
      "unit of test = donor group (units sharing a bank pooled; a gene's direction = majority over the group's units, ties dropped); null: shuffle the group's up/down labels among its own called genes (%d shuffles, seed 20261014), so the group's balance of up and down and the set of proteins it reports are held fixed" % N_PERM,
      "the excess is concordant share minus expected share; its interval is a jackknife over donor groups (+/-1.96 se); one group is one bank, so groups are not independent of the RNA collections"]
open(OUT + "/axis_direction_lists_by_group.tsv", "w").write("".join("# " + x + "\n" for x in H2) + G.to_csv(sep="\t", index=False))
open(OUT + "/axis_direction_lists_pooled.tsv", "w").write("".join("# " + x + "\n" for x in H2) + PL.to_csv(sep="\t", index=False))
pd.set_option("display.width", 250, "display.max_columns", 30, "display.max_rows", 200)
print(PL.to_string())
# ---- (e) the 209: protein vs RNA direction against each axis (sign of axis value)
REC = "core80_work_2026-10-16/P5/M03/inputs_new/recurrence_genes_final.tsv"
# PATCH: set definition from the final recurrent_sets.tsv (env SETNAME); without SETNAME the original q_bh < 0.1 rule on REC is used (old 176-gene file)
SETNAME = os.environ.get("SETNAME")
SETS_TSV = os.environ.get("SETS_TSV", "core80_work_2026-10-16/P5/M03/inputs_new/sets_final.tsv")
SETGENES = set(pd.read_csv(SETS_TSV, sep="\t", comment="#").query("set == @SETNAME").gene) if SETNAME else None
def in_set(df):
    return df.gene.isin(SETGENES) if SETNAME else df.q_bh < 0.1
rg = pd.DataFrame([r for r in csv.DictReader((l for l in open(REC) if not l.startswith("#")), delimiter="\t") if r["scenario"] == "published, brain"])
for c in ["q_bh", "up", "down"]: rg[c] = pd.to_numeric(rg[c])
rg = rg[in_set(rg)].set_index("gene"); rdir = {g: ("up" if u > d else "down" if d > u else "tie") for g, u, d in zip(rg.index, rg.up, rg.down)}
E_rows = []
def sgn(v): return "up" if v > 0 else "down"
sets = {"all 209": [g for g in rg.index if g in AX.index]}
for axn in ["pH", "agonal"]:
    v = AX[axn]; gg = [g for g in rg.index if g in AX.index and not np.isnan(v[g]) and rdir[g] in ("up", "down")]
    c = sum(rdir[g] == sgn(v[g]) for g in gg)
    E_rows.append({"unit": "RNA consensus, all 209 with an axis value", "axis": axn, "genes": len(gg), "agree_with_axis_sign": c, "share": round(c / len(gg), 3), "p_vs_half": binomtest(c, len(gg), .5).pvalue})
for pm in units:
    for axn in ["pH", "agonal"]:
        bank = units[pm]["bank"]; v = AX["pH_LOCO_" + BANK_COLL[bank]] if (axn == "pH" and bank in BANK_COLL) else AX[axn]
        gg = [g for g in UD[pm] if g in rg.index and g in AX.index and not np.isnan(v[g]) and rdir[g] in ("up", "down")]
        if len(gg) < 3: continue
        cp = sum(UD[pm][g] == sgn(v[g]) for g in gg); cr = sum(rdir[g] == sgn(v[g]) for g in gg); cpr = sum(UD[pm][g] == rdir[g] for g in gg)
        E_rows.append({"unit": units[pm]["label"], "axis": axn if not (axn == "pH" and bank in BANK_COLL) else "pH leaving out " + BANK_COLL[bank], "genes": len(gg), "protein_agree_with_axis_sign": cp, "rna_agree_with_axis_sign_same_genes": cr,
                       "protein_same_direction_as_rna": cpr, "share_protein": round(cp / len(gg), 3), "share_rna": round(cr / len(gg), 3), "p_protein_vs_half": binomtest(cp, len(gg), .5).pvalue})
ED = pd.DataFrame(E_rows)
open(OUT + "/axis_direction_on_209.tsv", "w").write("# built 2026-10-14 by scripts/analyze_lists.py; genes = the unit's directed proteins that are among the 209 recurrent genes, have an axis value and an RNA call (up or down); agreement = direction equals the sign of the axis value (positive = predicted up); same genes scored for the RNA consensus (majority of the gene's lists)\n" + ED.to_csv(sep="\t", index=False))
print(ED.to_string())
