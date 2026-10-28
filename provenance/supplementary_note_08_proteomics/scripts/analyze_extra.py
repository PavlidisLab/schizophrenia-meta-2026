#!/usr/bin/env python3
"""Step 5 (2026-10-14): (i) agreement of protein effects between studies and with RNA (ceiling for the axis correlations), (ii) pooled RNA comparator,
(iii) donor overlap with the Pittsburgh transcriptome studies by case number (MacDonald 2020 Table S1), (iv) inventory of omics-scale proteomic SZ brain studies and what is local.
Read-only on the repository. Run: PYTHONDONTWRITEBYTECODE=1 python3 -I analyze_extra.py"""
import sys, os, csv, glob, json, hashlib, collections, itertools
import numpy as np, pandas as pd
from scipy.stats import spearmanr
REPO = os.path.expanduser("sztr-rebuild-infra/repo")
OUT = os.environ["OUTDIR"]; W = OUT + "/work"
rng = np.random.default_rng(20261014)
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
UFG = {"Pinacho 2016": "upv", "Vidal-Domenech 2020": "upv", "Alnafisah 2022": "maryland", "Nucifora 2019": "pitt", "MacDonald 2020 homogenate": "pitt", "MacDonald 2020 synaptosome": "pitt",
       "Focking 2015 PSD": "stanley", "Testa 2026": "pitt", "Koopmans 2026": "koopmans", "Vera-Montecinos 2021": "barcelona"}
T = {k: pd.read_csv(f"{W}/protein_table__{k.replace(' ', '_')}.tsv", sep="\t").set_index("symbol").effect for k in UFG}
# (i) pairwise effect agreement
rows = []
for a, b in itertools.combinations(T, 2):
    m = pd.concat([T[a], T[b]], axis=1, join="inner").dropna()
    if len(m) < 30: continue
    x, y = m.iloc[:, 0].values, m.iloc[:, 1].values; rho = spearmanr(x, y)[0]
    bs = [spearmanr(x[i], y[i])[0] for i in (rng.integers(0, len(x), len(x)) for _ in range(500))]
    rows.append({"unit_a": a, "unit_b": b, "same_donor_group": UFG[a] == UFG[b], "shared_proteins": len(m), "rho": round(rho, 3), "ci_lo": round(np.quantile(bs, .025), 3), "ci_hi": round(np.quantile(bs, .975), 3)})
R1 = pd.DataFrame(rows)
cmc = pd.read_csv(REPO + "/side-studies/proteomics/data/cmc_mssm_penn_pitt_dxscz.tsv", sep="\t").dropna(subset=["t"]).drop_duplicates("symbol").set_index("symbol").t
# (ii) pooled RNA comparator: from the by-unit table
A1 = pd.read_csv(OUT + "/axis_effect_correlation_by_unit.tsv", sep="\t", comment="#")
cmp_rows = []
for axl in ["pH", "agonal"]:
    for setn, lab in [("all tested proteins", "protein effect vs axis"), ("RNA (CommonMind t) on the same proteins", "CommonMind RNA t vs axis, same proteins")]:
        s = A1[(A1.axis == axl) & (A1.set == setn) & ~A1.unit.str.startswith("Lorca")]
        cmp_rows.append({"axis": axl, "quantity": lab, "units": len(s), "donor_groups": s.group.nunique(), "pooled_rho_mean_over_donor_groups": round(s.groupby("group").rho.mean().mean(), 3),
                         "unit_rho_min": round(s.rho.min(), 3), "unit_rho_max": round(s.rho.max(), 3)})
s = A1[(A1.axis == "(protein effect vs RNA t)")]
cmp_rows.append({"axis": "-", "quantity": "protein effect vs CommonMind RNA t, same proteins", "units": len(s), "donor_groups": s.group.nunique(), "pooled_rho_mean_over_donor_groups": round(s.groupby("group").rho.mean().mean(), 3), "unit_rho_min": round(s.rho.min(), 3), "unit_rho_max": round(s.rho.max(), 3)})
C1 = pd.DataFrame(cmp_rows)
# (iii) donor overlap, MacDonald 2020 Table S1
cw = list(csv.DictReader((l for l in open(REPO + "/analyses/donor_overlap/pitt_crosswalk.tsv") if not l.startswith("#")), delimiter="\t"))
by = collections.defaultdict(set)
for x in cw: by[x["dataset_or_study"]].add(str(x["biobank_id"]).strip())
f = REPO + "/data/raw/31642882/preprint/table_s1_subjects_gsheet_1HEHXGVT9y10GdcmckSbDwU9OttfC_DPuc7SXiEC91PI.xlsx"
x = pd.read_excel(f, header=None); hdr = x.iloc[0].tolist()
left = x.iloc[1:, :13].copy(); left.columns = hdr[:13]; right = x.iloc[1:, 14:27].copy(); right.columns = hdr[14:27]
ids = {str(int(v)) for v in list(left["HU"].dropna()) + list(right["HU"].dropna())}
anyp = set().union(*by.values())
D = [{"proteomic_study": "MacDonald 2020 (LC-SRM, auditory cortex; Table S1, 48 pairs)", "donors_listed": len(ids), "transcriptome_study": "any Pittsburgh study with printed case numbers", "study_donors": len(anyp), "donors_in_both": len(ids & anyp), "basis": "exact case-number (HU number) match, " + os.path.relpath(f, os.path.expanduser("~")) }]
for s, v in sorted(by.items(), key=lambda t: -len(ids & t[1])):
    if ids & v: D.append({"proteomic_study": "MacDonald 2020", "donors_listed": len(ids), "transcriptome_study": s.split(" (")[0] if "CommonMind" not in s else "Fromer 2016 / CommonMind Pittsburgh DLPFC arm", "study_donors": len(v), "donors_in_both": len(ids & v), "basis": "exact case-number match"})
_INV = os.path.expanduser("review_claims_work/noncoding_modalities_inventory_2026-10-14/inventory.tsv")  # PATCH (P5): file absent on remote; set-independent, feeds studies_local_tables.tsv only
inv = [l for l in csv.DictReader((l for l in open(_INV) if not l.startswith("#")), delimiter="\t")] if os.path.exists(_INV) else []
DI = pd.DataFrame(D)
# (iv) inventory of omics-scale proteomic studies
units = pd.read_csv(REPO + "/side-studies/proteomics/units.tsv", sep="\t", comment="#").set_index("unit")
tab = {"27236410": "full table, per-protein log2 ratio (1,208 genes)", "32639965": "full table, per-protein log2 ratio (1,370 genes)", "35354897": "full table, limma logFC (487 genes)", "31055969": "full table, abundance ratio (868 genes)",
       "31642882": "full table, fold change (homogenate 402; synaptosome 155)", "25048004": "full table, fold change (698 genes)", "42748139": "full table, isoform effect (5,763 genes)", "41611717": "full table, log2 fold change by layer (5,243 genes)",
       "34576238": "per-subject intensities, 1,474 proteins (Supplementary dataset 1); effect computed here (1,449 genes)", "38255234": "pooled spectral counts, 3 preparations per group and region (Supplementary Datasets 1 and 2); exploratory effect"}
INV = []
seen = set()
for r in inv:
    if not r["group"].startswith("proteomics"): continue
    pm = r["pmid"]
    ls = sorted(glob.glob(f"{REPO}/data/lists/{pm}__*.tsv"))
    nl = 0; nd = 0
    for l in ls:
        rr = list(csv.DictReader((q for q in open(l) if not q.startswith("#")), delimiter="\t")); nl += len({q["symbol"] for q in rr if q["symbol"]}); nd += len({q["symbol"] for q in rr if q["symbol"] and q["direction"] in ("up", "down")})
    kind = tab.get(pm, "list of significant proteins only" if ls else ("no protein table on disk" if pm not in ("38255234",) else tab["38255234"]))
    if pm == "16624532": kind = "none: SELDI peaks not identified as proteins"
    if pm == "40385394": kind = "no table on disk (non-corpus study; text not extracted)"
    full = pm in tab
    INV.append({"study": r["study"], "pmid": pm, "region": r["region"], "bank_inventory": r["collection_study_linked"][:60], "evidence_level_donor_overlap_inventory": r["evidence_level"], "in_side_study_31_units": pm in units.index.astype(str).tolist() or pm in ("19405953",) and False,
                "local_table": kind, "proteins_in_lists_with_symbol": nl, "proteins_with_direction": nd, "scored_here": ("effect-size correlation, 209 overlap and top-N direction" if full and pm != "38255234" else "exploratory effect-size correlation" if pm == "38255234" else "top-N direction" if nd >= 5 else "not scored (no usable table or fewer than 5 directed proteins)")})
INVD = pd.DataFrame(INV)
if len(INVD): INVD["in_side_study_31_units"] = INVD.pmid.isin([str(i) for i in units.index])
def hdr(h): return "".join("# " + x + "\n" for x in h)
open(OUT + "/protein_effect_agreement_between_studies.tsv", "w").write(hdr(["built 2026-10-14 by scripts/analyze_extra.py; Spearman between two studies' protein effects over the proteins both report (gene symbols, HGNC); ci = bootstrap over proteins (500); same_donor_group = same bank / donor group as listed in the report"]) + R1.to_csv(sep="\t", index=False))
open(OUT + "/axis_comparison_with_rna_pooled.tsv", "w").write(hdr(["built 2026-10-14 by scripts/analyze_extra.py from axis_effect_correlation_by_unit.tsv; pooled = mean over donor groups of the mean rho over a group's units; Lorca 2024 (exploratory) left out"]) + C1.to_csv(sep="\t", index=False))
open(OUT + "/donor_overlap.tsv", "w").write(hdr(["built 2026-10-14 by scripts/analyze_extra.py; case numbers: MacDonald 2020 Table S1 ('HU' column), compared with analyses/donor_overlap/pitt_crosswalk.tsv (exact identifier rows; sha256 %s)" % sha(REPO + "/analyses/donor_overlap/pitt_crosswalk.tsv")[:12],
        "other proteomic studies: no per-subject table on disk, so donor overlap is not measured; see report.md for what each paper states"]) + DI.to_csv(sep="\t", index=False))
if os.path.exists(_INV): open(OUT + "/studies_local_tables.tsv", "w").write(hdr(["built 2026-10-14 by scripts/analyze_extra.py; rows = the inventory's proteomics_corpus and proteomics_noncorpus groups (noncoding_modalities_inventory_2026-10-14/inventory.tsv, sha256 %s); lists counted from data/lists/<pmid>__*.tsv" % sha(os.path.expanduser("review_claims_work/noncoding_modalities_inventory_2026-10-14/inventory.tsv"))[:12]]) + INVD.to_csv(sep="\t", index=False))
pd.set_option("display.width", 250, "display.max_rows", 200, "display.max_columns", 20)
print(R1.to_string()); print(C1.to_string()); print(DI.to_string()); print(INVD.drop(columns=["bank_inventory"]).to_string())
