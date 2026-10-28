#!/usr/bin/env python3
"""Step 3 (2026-10-14): score proteomic case-control results against the pH axis and the agonal (GTEx slow-minus-fast death) axis.
(a) effect-size correlation over all tested proteins, per unit and pooled by donor group (units with a full table on disk);
(b) the same over the 209 recurrent genes; (c) comparison with the RNA effect (CommonMind DxSCZ t) over the same genes;
(d) direction agreement of every directed protein list with the genes at the top of each axis (the manuscript's Extended Fig. 4 test, for both axes, with the pH axis leaving out the collection a study's bank shares);
(e) protein vs RNA direction on the 209 against each axis.
Orientation: positive axis value = predicted SZ-up. lowpH = -pH_axis (pH_axis > 0 means higher at higher pH); slow = agonal_axis (> 0 higher after slow death).
Read-only on the repository. Run: PYTHONDONTWRITEBYTECODE=1 python3 -I analyze_axes.py"""
import sys, os, csv, json, hashlib, glob, collections
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata, binomtest
REPO = os.path.expanduser("sztr-rebuild-infra/repo")
OUT = os.environ["OUTDIR"]; W = OUT + "/work"
sys.path.insert(0, REPO + "/analyses"); sys.path.insert(0, REPO + "/analyses/recurrence_sources"); sys.path.insert(0, REPO + "/side-studies/proteomics")
import analyse as A, recurrence as R, extract
extract.HGNC = extract.load_hgnc()
rng = np.random.default_rng(20261014)
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
# ---- axes
ax = pd.read_csv(W + "/axis_table_hgnc.tsv", sep="\t").set_index("sym")
cols = ["NSW", "Stanley", "Pitt", "BA10", "Pritzker"]
def stouffer(excl=None, minc=2):
    num = np.zeros(len(ax)); den = np.zeros(len(ax)); k = np.zeros(len(ax))
    for c in cols:
        if c == excl: continue
        t = ax["t_pH_" + c].values; n = ax["n_pH_" + c].values; m = ~np.isnan(t) & ~np.isnan(n)
        num += np.where(m, np.sqrt(np.nan_to_num(n)) * np.nan_to_num(t), 0); den += np.where(m, np.nan_to_num(n), 0); k += m
    z = np.where((k >= minc) & (den > 0), num / np.sqrt(np.where(den > 0, den, 1)), np.nan); return pd.Series(z, index=ax.index)
AX = pd.DataFrame({"pH": -ax.pH_axis, "agonal": ax.agonal_axis})
for c in ["NSW", "Stanley", "Pitt", "Pritzker"]:
    AX["pH_LOCO_" + c] = -stouffer(c)
BANK_COLL = {"stanley": "Stanley", "pitt": "Pitt", "nsw": "NSW", "aus": "NSW", "uci": "Pritzker"}   # collection to leave out for a bank
def axvec(name, bank=None):
    if name == "pH_matched" : return AX["pH_LOCO_" + BANK_COLL[bank]] if bank in BANK_COLL else AX["pH"]
    return AX[name]
# ---- RNA comparator
cmc = pd.read_csv(REPO + "/side-studies/proteomics/data/cmc_mssm_penn_pitt_dxscz.tsv", sep="\t")
cmc["sym"] = [extract.symbol_to_hgnc(s, official=True)[0] if isinstance(s, str) else "" for s in cmc.symbol]
cmc = cmc[cmc.sym != ""].dropna(subset=["t"]).drop_duplicates("sym").set_index("sym")
# ---- 209
REC = "core80_work_2026-10-16/P5/M03/inputs_new/recurrence_genes_final.tsv"
# PATCH: set definition from the final recurrent_sets.tsv (env SETNAME); without SETNAME the original q_bh < 0.1 rule on REC is used (old 176-gene file)
SETNAME = os.environ.get("SETNAME")
SETS_TSV = os.environ.get("SETS_TSV", "core80_work_2026-10-16/P5/M03/inputs_new/sets_final.tsv")
SETGENES = set(pd.read_csv(SETS_TSV, sep="\t", comment="#").query("set == @SETNAME").gene) if SETNAME else None
def in_set(df):
    return df.gene.isin(SETGENES) if SETNAME else df.q_bh < 0.1
rg = pd.DataFrame([r for r in csv.DictReader((l for l in open(REC) if not l.startswith("#")), delimiter="\t") if r["scenario"] == "published, brain"])
for c in ["q_bh", "up", "down"]: rg[c] = pd.to_numeric(rg[c])
rg = rg[in_set(rg)].set_index("gene"); rna_net = (rg.up - rg.down) / (rg.up + rg.down)
rna_dir = {g: ("up" if u > d else "down" if d > u else "tie") for g, u, d in zip(rg.index, rg.up, rg.down)}
S209 = set(rg.index)
# ---- full-table units
UF = {  # label: (file label, bank/group)
 "Pinacho 2016": "upv", "Vidal-Domenech 2020": "upv", "Alnafisah 2022": "maryland", "Nucifora 2019": "pitt", "MacDonald 2020 homogenate": "pitt",
 "MacDonald 2020 synaptosome": "pitt", "Focking 2015 PSD": "stanley", "Testa 2026": "pitt", "Koopmans 2026": "koopmans", "Vera-Montecinos 2021": "barcelona"}
T = {k: pd.read_csv(f"{W}/protein_table__{k.replace(' ', '_')}.tsv", sep="\t").set_index("symbol") for k in UF}
# Lorca 2024 (extracellular vesicles; pooled preparations, spectral counts): exploratory only
def lorca():
    out = {}
    d1 = W + "/zips/lorca/Supplementary Dataset 1.xlsx"; d2 = W + "/zips/lorca/Supplementary Dataset 2.xlsx"
    def sheet(f, s):
        x = pd.read_excel(f, sheet_name=s, header=None); i = x.index[x.iloc[:, 0] == "GS"][0]
        y = x.iloc[i + 1:, [0, 3]].copy(); y.columns = ["gs", "sc"]; y["sc"] = pd.to_numeric(y.sc, errors="coerce"); y = y.dropna()
        y["sym"] = [extract.symbol_to_hgnc(g, official=True)[0] if isinstance(g, str) else "" for g in y.gs]; return y[y.sym != ""].groupby("sym").sc.sum()
    for reg in ["PFC", "HC", "CAU"]:
        C = pd.concat([sheet(d1, f"{reg}_C{i}") for i in (1, 2, 3)], axis=1).fillna(0); S = pd.concat([sheet(d2, f"{reg}_SZ{i}") for i in (1, 2, 3)], axis=1).fillna(0)
        idx = C.index.union(S.index); C = C.reindex(idx).fillna(0); S = S.reindex(idx).fillna(0)
        keep = ((C > 0).sum(axis=1) >= 2) | ((S > 0).sum(axis=1) >= 2)
        e = np.log2((S.mean(axis=1) + 0.5) / (C.mean(axis=1) + 0.5))[keep]
        out["Lorca 2024 EVs, " + reg] = pd.DataFrame({"effect": e, "p": np.nan})
    return out
LO = lorca()
for k, v in LO.items(): T[k] = v; UF[k] = "lorca"
# ---- (a)(b)(c) per unit
def boot_sp(x, y, B=2000):
    n = len(x); out = np.empty(B)
    for b in range(B):
        i = rng.integers(0, n, n); out[b] = spearmanr(x[i], y[i])[0]
    return np.quantile(out, [.025, .975])
rows = []; bootdraws = {}
for k, tab in T.items():
    bank = UF[k]
    for axname in ["pH", "pH_matched", "agonal"]:
        v = axvec(axname, bank)
        if axname == "pH_matched" and bank not in BANK_COLL: continue
        m = tab.join(v.rename("axis"), how="inner").dropna(subset=["effect", "axis"])
        if len(m) < 20: continue
        rho, p = spearmanr(m.effect, m.axis); lo, hi = boot_sp(m.effect.values, m.axis.values)
        r = {"unit": k, "group": bank, "axis": axname if axname != "pH_matched" else "pH leaving out " + BANK_COLL[bank], "set": "all tested proteins", "n": len(m), "rho": rho, "ci_lo": lo, "ci_hi": hi, "p_nominal": p}
        rows.append(r)
        # 209 subset
        s2 = m[m.index.isin(S209)]
        if len(s2) >= 8:
            rho2, p2 = spearmanr(s2.effect, s2.axis); lo2, hi2 = boot_sp(s2.effect.values, s2.axis.values, 1000)
            rows.append({"unit": k, "group": bank, "axis": r["axis"], "set": "209 recurrent genes", "n": len(s2), "rho": rho2, "ci_lo": lo2, "ci_hi": hi2, "p_nominal": p2})
        # RNA comparator over the same genes (CommonMind DxSCZ t)
        mm = m.join(cmc.t.rename("rna_t"), how="inner")
        if len(mm) >= 20 and k not in LO:
            ra, _ = spearmanr(mm.rna_t, mm.axis); lo3, hi3 = boot_sp(mm.rna_t.values, mm.axis.values, 1000)
            rows.append({"unit": k, "group": bank, "axis": r["axis"], "set": "RNA (CommonMind t) on the same proteins", "n": len(mm), "rho": ra, "ci_lo": lo3, "ci_hi": hi3, "p_nominal": spearmanr(mm.rna_t, mm.axis)[1]})
            if axname == "pH":
                rpp, ppp = spearmanr(mm.effect, mm.rna_t); lo4, hi4 = boot_sp(mm.effect.values, mm.rna_t.values, 1000)
                rows.append({"unit": k, "group": bank, "axis": "(protein effect vs RNA t)", "set": "same proteins", "n": len(mm), "rho": rpp, "ci_lo": lo4, "ci_hi": hi4, "p_nominal": ppp})
    print(k, flush=True)
A1 = pd.DataFrame(rows)
# pooled by donor group: mean over groups of the mean rho over units in a group (all tested proteins)
def pooled(sel_units, axis_label, setname="all tested proteins", B=1000):
    sub = A1[(A1.axis == axis_label) & (A1.set == setname) & A1.unit.isin(sel_units)]
    if sub.empty: return None
    est = sub.groupby("group").rho.mean().mean()
    # bootstrap genes within each unit
    ests = []
    for b in range(B):
        g = collections.defaultdict(list)
        for k in sub.unit:
            bank = UF[k]; an = "pH" if axis_label == "pH" else ("agonal" if axis_label == "agonal" else None)
            v = axvec("pH_matched" if axis_label.startswith("pH leaving") else an, bank)
            m = T[k].join(v.rename("axis"), how="inner").dropna(subset=["effect", "axis"])
            i = rng.integers(0, len(m), len(m)); g[bank].append(spearmanr(m.effect.values[i], m.axis.values[i])[0])
        ests.append(np.mean([np.mean(x) for x in g.values()]))
    lo, hi = np.quantile(ests, [.025, .975]); return est, lo, hi, sub.group.nunique(), len(sub)
FT = [k for k in UF if k not in LO]
P = []
for axl in ["pH", "agonal"]:
    for lab, sel in [("all 10 full-table units", FT), ("without Testa 2026", [k for k in FT if k != "Testa 2026"]), ("without the Pittsburgh group", [k for k in FT if UF[k] != "pitt"]),
                     ("one-line summary: 6 donor groups, MacDonald 2020 synaptosome + homogenate averaged within group (as all units)", FT)]:
        if lab.startswith("one-line"): continue
        r = pooled(sel, axl)
        if r: P.append({"axis": axl, "units": lab, "pooled_rho_mean_over_donor_groups": r[0], "ci_lo": r[1], "ci_hi": r[2], "donor_groups": r[3], "units_n": r[4]})
# leave-one-group-out
for axl in ["pH", "agonal"]:
    for gname in sorted(set(UF[k] for k in FT)):
        r = pooled([k for k in FT if UF[k] != gname], axl, B=300)
        if r: P.append({"axis": axl, "units": "leave out group " + gname, "pooled_rho_mean_over_donor_groups": r[0], "ci_lo": r[1], "ci_hi": r[2], "donor_groups": r[3], "units_n": r[4]})
PD = pd.DataFrame(P)
def hdr(h): return "\n".join("# " + x for x in h) + "\n"
H1 = ["built 2026-10-14 by scripts/analyze_axes.py; axes from axis_gene_table.tsv (sha256 %s, built 2026-10-14 by the axis tool); orientation: positive axis value = predicted higher in SZ (pH: minus the tool's pH_axis; agonal: slow-minus-fast death t)" % sha(W + "/axis_gene_table.tsv")[:12],
      "rho = Spearman between the protein's case-control effect (signed log2 fold change or the study's effect measure, see protein_tables_meta) and the axis value; positive = the study's protein effects lie along the axis in the direction the pH / slow-death account predicts; ci = 95% bootstrap over proteins (2,000 resamples; 1,000 for subsets); p_nominal treats proteins as independent",
      "pH leaving out X: the Stouffer pH axis recomputed without collection X (the collection sharing a brain bank with the study), genes with data in at least 2 remaining collections",
      "RNA row: CommonMind MSSM-Penn-Pitt DxSCZ t (Hoffman 2021 per-cohort) for the proteins of the unit, with the same axis; a Pittsburgh-bank study's donors overlap CommonMind Pittsburgh donors (donor_overlap.tsv). Lorca 2024 rows: log2 ratio of mean spectral counts in 3 pooled extracellular-vesicle preparations per group, exploratory"]
open(OUT + "/axis_effect_correlation_by_unit.tsv", "w").write(hdr(H1) + A1.to_csv(sep="\t", index=False))
open(OUT + "/axis_effect_correlation_pooled.tsv", "w").write(hdr(H1[:3] + ["pooled = mean over donor groups of the mean rho over a group's units (all tested proteins); ci = bootstrap over proteins within each unit (1,000 resamples; 300 for leave-one-group-out), which does not include between-study or between-donor variation"]) + PD.to_csv(sep="\t", index=False))
pd.set_option("display.width", 220, "display.max_columns", 30, "display.max_rows", 200)
print(A1.round(3).to_string()); print(PD.round(3).to_string())
A1.to_pickle(W + "/A1.pkl")
