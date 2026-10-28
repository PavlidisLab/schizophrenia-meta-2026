#!/usr/bin/env python3
"""Per-gene table for the 383 recurrent genes: Liharska living-versus-postmortem statistic (original and with RNA-quality surrogate variables), SZ direction and z, pH axis v2 and agonal axis (RA, 2026-10-19). Reuses the loaders and rules of liharska_living_vs_pm_2026-10-18/scripts/s1_main.py."""
import sys, os, csv, hashlib, datetime
sys.dont_write_bytecode = True
import numpy as np, pandas as pd
H = os.path.expanduser("~"); REPO = f"{H}/Dev/schizophrenia-transcriptome-replication"; AN = f"{REPO}/analyses"
sys.path.insert(0, AN)
LV = f"{H}/analysis-delivery/liharska_living_vs_pm_2026-10-18"
B80 = f"{H}/out"
REC = f"{B80}/recurrence_genes.tsv"
ST1 = f"{LV}/data/liv_pm_de_ST1.tsv"; CT4 = f"{LV}/data/collado_torres/TableS4.csv"
PHV2 = f"{H}/axis_v2_controls.tsv.gz"
GTEX = f"{AN}/gtex_ischemia_2026-10-08/hardy_t_by_region.tsv.gz"
GS = f"{H}/fig07/fig07_a2_a4_gene_scatter.tsv"
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
import extract
extract.HGNC = extract.load_hgnc()
ENS2SYM = {e: s for e, (s, _) in extract.HGNC[3].items()}
_c = {}
def hg(s):
    if s not in _c: _c[s] = extract.symbol_to_hgnc(s)[0] or None
    return _c[s]
def best(df, sym, p):
    df = df.assign(_s=sym).dropna(subset=["_s", p]); return df.sort_values(p).drop_duplicates("_s").set_index("_s")
rec = [r for r in csv.DictReader((l for l in open(REC) if not l.startswith("#")), delimiter="\t") if r["scenario"] == "published, brain"]
S383 = [r["gene"] for r in rec if float(r["q_bh"]) < 0.1]; assert len(S383) == 383
sign = {r["gene"]: (1 if int(r["up"]) - int(r["down"]) > 0 else -1 if int(r["up"]) - int(r["down"]) < 0 else 0) for r in rec}
l = pd.read_csv(ST1, sep="\t", comment="#"); l["sym"] = l.gene.map(ENS2SYM); l1 = best(l, l.sym, "P.Value")
c = pd.read_csv(CT4, low_memory=False); c = c[c.isExprs].copy(); c["sym"] = c.ensemblID.map(ENS2SYM)
c["sym"] = c.sym.where(c.sym.notna(), c.Symbol.map(lambda s: hg(str(s)) if isinstance(s, str) else None))
cw = best(c, c.sym, "with_qSVA_P.Value").dropna(subset=["with_qSVA_t"]); cn = best(c, c.sym, "without_qSVA_P.Value").dropna(subset=["without_qSVA_t"])
ph = pd.read_csv(PHV2, sep="\t", comment="#", index_col=0).z_stouffer; ph.index = [hg(s) or s for s in ph.index]; ph = ph.groupby(level=0).first()
Hd = pd.read_csv(GTEX, sep="\t", index_col=0)
ag = Hd[["Brain - Frontal Cortex (BA9)", "Brain - Anterior cingulate cortex (BA24)", "Brain - Cortex"]].mean(axis=1, skipna=False).dropna(); ag.index = [hg(s) or s for s in ag.index]; ag = ag.groupby(level=0).first()
gs = pd.read_csv(GS, sep="\t", comment="#").set_index("gene")
rows = []
for g in S383:
    r = dict(gene=g, sz_direction={1: "up", -1: "down", 0: "none"}[sign[g]], z_sz_raw=gs.z_sz_raw.get(g, np.nan))
    r["liv_pm_t"] = l1.t.get(g, np.nan); r["liv_pm_adj_p"] = l1["adj.P.Val"].get(g, np.nan)
    r["liv_pm_t_without_qsva"] = cn.without_qSVA_t.get(g, np.nan)
    r["liv_pm_t_with_qsva"] = cw.with_qSVA_t.get(g, np.nan); r["liv_pm_adj_p_with_qsva"] = cw["with_qSVA_adj.P.Val"].get(g, np.nan)
    r["in_374_with_direction"] = "TRUE" if (sign[g] != 0 and not np.isnan(r["liv_pm_t"])) else "FALSE"
    r["liv_pm_in_sz_direction"] = ("" if r["in_374_with_direction"] == "FALSE" else "TRUE" if np.sign(r["liv_pm_t"]) == sign[g] else "FALSE")
    r["liv_pm_with_qsva_in_sz_direction"] = ("" if (sign[g] == 0 or np.isnan(r["liv_pm_t_with_qsva"])) else "TRUE" if np.sign(r["liv_pm_t_with_qsva"]) == sign[g] else "FALSE")
    r["ph_axis_v2_z"] = ph.get(g, np.nan); r["agonal_axis_t"] = ag.get(g, np.nan)
    rows.append(r)
d = pd.DataFrame(rows)
n374 = (d.in_374_with_direction == "TRUE").sum(); c374 = (d.liv_pm_in_sz_direction == "TRUE").sum()
from scipy import stats
rp = stats.spearmanr(d.liv_pm_t, d.ph_axis_v2_z, nan_policy="omit")[0]; ra = stats.spearmanr(d.liv_pm_t, d.agonal_axis_t, nan_policy="omit")[0]
hd = ["# figure: Fig. 4 panel, freshly resected against postmortem (Liharska et al. 2025, Mol Psychiatry, doi 10.1038/s41380-025-03163-1), per-gene scatter over the 383 recurrent genes",
 "# script: fig4_liharska_per_gene_2026-10-19/build_liharska_per_gene.py, run %s; loaders and rules reused from liharska_living_vs_pm_2026-10-18/scripts/s1_main.py (best row per symbol by smallest P; Ensembl ids mapped through HGNC)" % datetime.date.today(),
 f"# input: liharska_living_vs_pm_2026-10-18/data/liv_pm_de_ST1.tsv sha256 {sha(ST1)[:12]} (Liharska Supplementary Table 1, MOESM8: per-gene moderated t of postmortem against living prefrontal cortex); data/collado_torres/TableS4.csv sha256 {sha(CT4)[:12]} (Collado-Torres et al. 2023, living_brain_reanalysis, Zenodo 10.5281/zenodo.10010510: the same contrast without and with quality surrogate variables, qSVA)",
 f"# input: recurrence_on_80_studies/out/recurrence_genes.tsv sha256 {sha(REC)[:12]} (383 genes; sz_direction = sign of n_up - n_down lists); figure_data_core80_2026-10-17/fig07/fig07_a2_a4_gene_scatter.tsv sha256 {sha(GS)[:12]} (z_sz_raw: Stouffer z of the schizophrenia effect, same z used elsewhere in Fig. 4); ph_axis_from_ten_dataset_regions/axis_v2_controls.tsv.gz sha256 {sha(PHV2)[:12]}; analyses/gtex_ischemia_2026-10-08/hardy_t_by_region.tsv.gz sha256 {sha(GTEX)[:12]}; data/reference/hgnc_complete_set.txt via extract.py",
 "# sign conventions: liv_pm_t positive = higher in postmortem, negative = higher in freshly resected (living) tissue; liv_pm_t_with_qsva is the same contrast after adjusting for quality surrogate variables from a degradation experiment (Collado-Torres 2023) and removes most of the signal by construction; ph_axis_v2_z positive = higher at higher brain pH (control-only, ten dataset-regions); agonal_axis_t positive = higher after slow death (mean over BA9, BA24, Cortex of GTEx slow-minus-fast natural death t); z_sz_raw positive = higher in schizophrenia",
 "# in_374_with_direction = gene has an SZ direction (up or down) and a Liharska statistic; liv_pm_in_sz_direction = sign(liv_pm_t) equals the SZ direction; blank where a gene has no direction or no statistic; genes without a value are blank (not measured in that source)",
 f"# checks: genes in the 374 = {n374}; in the SZ direction = {c374} ({c374/n374:.3f}); Spearman of liv_pm_t with ph_axis_v2_z = {rp:.3f} and with agonal_axis_t = {ra:.3f} over the genes with both (the delivered Pearson values, -0.51 and +0.66, come from liharska_living_vs_pm_2026-10-18/fig_livpm_data.tsv)"]
with open("liharska_per_gene_383.tsv", "w") as f:
    f.write("\n".join(hd) + "\n"); d.round(5).to_csv(f, sep="\t", index=False)
print(n374, c374, round(c374/n374, 4), round(rp, 3), round(ra, 3)); print(d.notna().mean().round(2).to_dict())
