#!/usr/bin/env python3
"""Fig. 3 panels a and b of the sample-level reanalysis on the 81 base.

Panel a is built with the logic of fig9_a() of analyses/figure_exports/export_fig08_fig09_fig10.py (copied; only the
dataset set, the pair relation table and the paths are inputs). Panel b (PsychAD, independent of the bulk base) is rebuilt
with fig9_b() unchanged and delivered as is.

Run from the run folder:  python assemble_a.py
Reads : in/stored/*  (stored 64-base inputs and delivered fig09 files), in/single_cell/*,
        out/validate28_b/b_pairs.tsv (rerun of b_covariate_series.R, levels L1 and L2),
        out/new_pairs/pairs_rho_from_b_de.tsv, in/derived_stage2/*_meta.tsv,
        repo/analyses/subject_relations/subject_relations*.tsv
Writes: deliver/ (see README)
"""
import csv, hashlib, itertools, json, os, sys
import numpy as np
import pandas as pd

RUN = os.environ["RUN"]
OUT = os.path.join(RUN, "deliver")
os.makedirs(OUT, exist_ok=True)
SEED = 20261015
rng = np.random.default_rng(SEED)
VAL = []     # validation rows


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def pdr(path, **kw):
    return pd.read_csv(path, sep="\t", comment="#", **kw)


def fmt(v):
    if isinstance(v, (float, np.floating)):
        return "" if np.isnan(v) else repr(round(float(v), 6))
    if isinstance(v, (np.integer,)):
        return str(int(v))
    return v


def write_tsv(path, header, fields, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        for l in header:
            fh.write("# " + l + "\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(fields)
        n = 0
        for r in rows:
            w.writerow([("" if r.get(f) is None else fmt(r.get(f))) for f in fields])
            n += 1
    return n


def body(path):
    """file body without '#' lines (header row and data)"""
    return [l for l in open(path, encoding="utf-8") if not l.startswith("#")]


IN = lambda *p: os.path.join(RUN, "in", *p)
ST = lambda f: IN("stored", f)
AN = os.path.join(RUN, "repo", "analyses")

# ------------------------------------------------------------------ inputs
b_val = pd.read_csv(os.path.join(RUN, "out", "validate28_b", "b_pairs.tsv"), sep="\t")
b_old = pd.read_csv(ST("b_pairs.tsv"), sep="\t")
DU = pdr(ST("datasets_used.tsv"))
NP = pd.read_csv(os.path.join(RUN, "out", "new_pairs", "pairs_rho_from_b_de.tsv"), sep="\t")
SR = pdr(os.path.join(AN, "subject_relations", "subject_relations.tsv"), dtype=str)
UN = pdr(os.path.join(AN, "subject_relations", "subject_relations_units.tsv"), dtype=str)
CP = pdr(IN("single_cell", "sc_cross_pairs.tsv"))

# ---- V: rerun of b_covariate_series.R (unchanged) against the stored b_pairs.tsv (levels L1, L2)
k = ["dataset_a", "dataset_b", "level"]
m = b_old[b_old.level.isin(["L1", "L2"])].merge(b_val, on=k, suffixes=("_s", "_v"))
mx = max((m[c + "_s"] - m[c + "_v"]).abs().max() for c in ["n_common", "rho_t", "rho_logfc", "r_logfc", "top100_overlap", "fdr_overlap", "pi1"])
VAL.append(("V01", "fairness_sample_level/scripts/b_covariate_series.R unchanged, SERIES_LEVELS=L1,L2: all metrics of b_pairs.tsv",
            f"{int((b_old.level.isin(['L1','L2'])).sum())} stored rows (sha {sha(ST('b_pairs.tsv'))[:12]})", f"{len(m)} rows matched; max abs difference over n_common, rho_t, rho_logfc, r_logfc, top100_overlap, fdr_overlap, pi1 = {mx:g}",
            "PASS" if mx == 0 and len(m) == int(b_old.level.isin(['L1', 'L2']).sum()) else "FAIL"))

# ---- units
NEW = ["GSE62191", "GSE87610_L3", "GSE87610_L5", "GSE78936_BA9", "GSE78936_BA11", "GSE78936_BA24", "GSE12654", "GSE12649", "GSE120340"]
NEWBULK = [u for u in NEW if not u.startswith("GSE87610")]
meta = {}
for f in ("GSE87610_meta.tsv", "GSE62191_meta.tsv", "GSE78936_meta.tsv", "GSE12649_meta.tsv", "extra_microarray_meta.tsv"):
    for r in pd.read_csv(IN("derived_stage2", f), sep="\t", dtype=str).itertuples():
        meta[r.dataset] = r
nsz = {d: int(a) for d, a in zip(DU.dataset, DU.n_sz)}
nct = {d: int(a) for d, a in zip(DU.dataset, DU.n_ctl)}
for u in NEW:
    nsz[u], nct[u] = int(meta[u].n_sz), int(meta[u].n_ctl)
comp = dict(zip(DU.dataset, DU.compartment))
comp.update({u: "bulk" for u in NEWBULK})
comp.update({"GSE87610_L3": "LCM layer 3 pyramidal neurons", "GSE87610_L5": "LCM layer 5 pyramidal neurons"})
REGB = [d for d in DU.dataset if comp[d] == "bulk"]
assert len(DU) == 28 and len(REGB) == 22
OUTSIDE = ["GSE92538_GPL570", "GSE92538_GPL96", "GSE208338"]   # Hagenauer 2018 and Worf 2025: Table 6 (excluded), not in base81.tsv
DERIVED7 = ["GSE35977", "GSE35974", "GSE80655_DLPFC", "GSE80655_AnCg", "GSE80655_nAcc", "GSE12679_neuron", "GSE12679_endothelial"]
SETS = {
    "S28_previous": REGB,                                                   # 22 bulk of the 28 dataset-regions
    "S34_base81": [d for d in REGB if d not in OUTSIDE] + NEWBULK,           # datasets of studies in the 81 base
    "S29_28plus9": REGB + NEWBULK,                                           # 28 dataset-regions + 9 derived units, as proposed
}
ALLUNITS = {"S28_previous": list(DU.dataset), "S34_base81": [d for d in DU.dataset if d not in OUTSIDE] + NEW, "S29_28plus9": list(DU.dataset) + NEW}

# ---- relations: stored table; GSE120340 by the same bank rule (copied from subject_relations.py bank_relation)
BANK_SAME = {frozenset(("cmc", "pitt")), frozenset(("cmc", "mssm"))}


def bank_relation(ba, ca, bb, cb):
    if ba is None or bb is None:
        return "bank not stated", "subject overlap unknown, bank not stated"
    if ba == bb == "stanley":
        if "unknown" in (ca, cb):
            return "Stanley, sub-collection not stated", "same collection, subject overlap unknown"
        if ca == cb:
            return f"Stanley, same sub-collection ({ca})", "same collection, subject overlap unknown"
        return "Stanley, different sub-collections (Array vs Neuropathology Consortium)", "subject-separate"
    if ba == bb:
        return f"same bank ({ba})", "same collection, subject overlap unknown"
    if frozenset((ba, bb)) in BANK_SAME:
        return f"CommonMind with a source bank ({ba}, {bb})", "same collection, subject overlap unknown"
    if frozenset((ba, bb)) == frozenset(("hbcc", "libd")):
        return "HBCC with LIBD (donor_overlap: fingerprint overlap)", "same collection, subject overlap unknown"
    return f"different banks ({ba}, {bb})", "subject-separate"


UNIT = UN[UN.kind == "dataset"].set_index("unit")
REL = {}
for r in SR.itertuples():
    REL[frozenset((r.dataset_a, r.dataset_b))] = dict(group=r.relation_group, evidence=r.evidence, bank_a=r.bank_a, bank_b=r.bank_b, src="subject_relations.tsv",
                                                    order=(r.dataset_a, r.dataset_b))
bank_of = lambda d: None if (d not in UNIT.index or UNIT.loc[d, "bank"] in ("not stated",)) else UNIT.loc[d, "bank"]
coll_of = lambda d: UNIT.loc[d, "stanley_collection"] if d in UNIT.index and isinstance(UNIT.loc[d, "stanley_collection"], str) else ""
added = []
for d in list(DU.dataset) + NEW:
    if d == "GSE120340":
        continue
    key = frozenset(("GSE120340", d))
    if key in REL:
        continue
    br, rel = bank_relation("stanley", "unknown", bank_of(d), coll_of(d) or None)
    REL[key] = dict(group=rel, evidence="collection", bank_a="stanley", bank_b=bank_of(d) or "not stated", src="added here: subject_relations.py bank_relation rule, GSE120340 (Abdolmaleky 2019, cohort cluster stanley, sub-collection not stated; no age or sex recorded, so no fingerprint test)",
                    order=("GSE120340", d))
    added.append((d, rel))
pd.DataFrame(added, columns=["partner", "relation_group"]).to_csv(os.path.join(OUT, "tables_relations_GSE120340_added.tsv"), sep="\t", index=False)


def b2_separate(key):
    """B2 reading of the bank rule: shared persons only within one bank; a fingerprint match between two different banks does not count"""
    r = REL[key]
    if r["group"] == "subject-separate":
        return True
    if r["evidence"] == "fingerprint" and r["bank_a"] != r["bank_b"] and "not stated" not in (r["bank_a"], r["bank_b"]) \
            and frozenset((r["bank_a"], r["bank_b"])) not in BANK_SAME and r["group"].startswith("partial subject overlap"):
        return True
    return False


# ---- rho lookups
def rho_reg(a, b, lv):
    x = b_val[(b_val.level == lv) & (((b_val.dataset_a == a) & (b_val.dataset_b == b)) | ((b_val.dataset_a == b) & (b_val.dataset_b == a)))]
    assert len(x) == 1, (a, b, lv, len(x))
    return float(x.rho_t.iloc[0])


NPI = {}
for r in NP.itertuples():
    NPI[(r.variant, frozenset((r.dataset_a, r.dataset_b)))] = r.rho_t


def rho_pair(a, b, mode):
    """mode 'L2': registry pairs at L2, pairs with a derived unit: partner at L2 (derived own limma fit); 'L1': every partner at L1"""
    if a in NEW and b in NEW:
        return NPI[("new_vs_new", frozenset((a, b)))]
    if a in NEW or b in NEW:
        return NPI[(f"new_vs_reg_{mode}", frozenset((a, b)))]
    return rho_reg(a, b, mode)


# registry pair order of the stored file (b_pairs, level L2)
base_order = b_old[(b_old.level == "L2")][["dataset_a", "dataset_b"]].values.tolist()
ORD = {frozenset(x): tuple(x) for x in base_order}


def pair_list(units):
    ps = []
    for a, b in itertools.combinations(units, 2):
        key = frozenset((a, b))
        if key in ORD:
            ps.append(ORD[key])
        elif (a in NEW) != (b in NEW):
            ps.append((a, b) if a in NEW else (b, a))
        else:
            ps.append((a, b))
    # registry pairs first in the order of the stored file, then pairs with derived units
    pos = {tuple(p): i for i, p in enumerate(base_order)}
    ps.sort(key=lambda p: (0, pos[tuple(p)]) if tuple(p) in pos else (1, p[0], p[1]))
    return ps


# ------------------------------------------------------------------ panel a for one set
def sn_rows():
    CPn = CP.copy()
    CPn["n_min"] = CPn[["n_sz_a", "n_ctl_a", "n_sz_b", "n_ctl_b"]].min(axis=1)
    E = CPn[(CPn.group_a == "endresz") & (CPn.group_b == "endresz") & (CPn.relation == "subject-separate")]
    rows = [{"data_type": "single-nucleus cohort pair (Endresz 2026, one pipeline)", "class_or_bulk": c, "unit_a": a, "unit_b": b, "relation": "subject-separate", "n_per_group": n, "rho_t": r}
            for c, a, b, n, r in zip(E.gemma_class, E.source_a, E.source_b, E.n_min, E.rho_t)]
    return E, rows


E_SN, SN_ROWS = sn_rows()


def panel_a(setname, rule, mode, header_extra):
    units = SETS[setname]
    rows, nb = [], 0
    pairs = pair_list(units)
    sep = lambda p: (REL[frozenset(p)]["group"] == "subject-separate") if rule == "R1" else b2_separate(frozenset(p))
    keep = [p for p in pairs if sep(p)]
    for a, b in keep:
        n = min(min(nsz[a], nct[a]), min(nsz[b], nct[b]))
        if a in NEW or b in NEW:
            dt = "bulk dataset pair (derived-tier limma fit with a derived unit, level L1 as recorded)"
        else:
            dt = "bulk dataset pair (Gemma, level L2)" if mode == "L2" else "bulk dataset pair (Gemma, level L1)"
        rows.append({"data_type": dt, "class_or_bulk": "bulk tissue", "unit_a": a, "unit_b": b, "relation": "subject-separate", "n_per_group": n, "rho_t": rho_pair(a, b, mode)})
    BP = pd.DataFrame(rows)
    BP["n_min"] = BP.n_per_group
    ee = E_SN.groupby(["source_a", "source_b"]).agg(rho_t=("rho_t", "median"), n_min=("n_min", "first")).reset_index()
    bins = []
    for lo, hi, lab in ((0, 10, "0-10"), (11, 30, "11-30"), (31, 1000, "31+")):
        bb = BP[(BP.n_min >= lo) & (BP.n_min <= hi)]
        e2 = ee[(ee.n_min >= lo) & (ee.n_min <= hi)]
        for dt, d in (("bulk dataset pairs", bb), ("single-nucleus cohort pairs (median over classes)", e2)):
            bins.append({"donors_per_group_bin": lab, "data_type": dt, "n_pairs": len(d), "median_rho_t": float(d.rho_t.median()) if len(d) else None,
                         "q25": float(d.rho_t.quantile(.25)) if len(d) else None, "q75": float(d.rho_t.quantile(.75)) if len(d) else None})
    cls = E_SN.groupby("gemma_class").rho_t.median()
    crow = [{"class": c, "n_pair_rows": int((E_SN.gemma_class == c).sum()), "median_rho_t": float(v)} for c, v in cls.items()]
    return BP, rows, bins, crow, keep


def write_a(dirpath, setname, rule, mode, label):
    BP, rows, bins, crow, keep = panel_a(setname, rule, mode, label)
    hdr = [f"figure: fig09 panel a  ({label})",
           "script: sample_level_81_2026-10-15/code/assemble_a.py (panel logic of analyses/figure_exports/export_fig08_fig09_fig10.py fig9_a); run 2026-10-15",
           f"dataset set: {setname} ({len(SETS[setname])} bulk datasets); pair rule: {rule} ({'relation_group == subject-separate of subject_relations.tsv, the rule behind the stored 192 pairs' if rule == 'R1' else 'B2 reading: fingerprint matches between two different banks do not count as shared persons'}); covariate level of bulk partners: {mode}",
           f"input: out/validate28_b/b_pairs.tsv (rerun of fairness_sample_level/scripts/b_covariate_series.R, unchanged) sha256 {sha(os.path.join(RUN,'out','validate28_b','b_pairs.tsv'))[:12]}",
           f"input: out/new_pairs/pairs_rho_from_b_de.tsv sha256 {sha(os.path.join(RUN,'out','new_pairs','pairs_rho_from_b_de.tsv'))[:12]}",
           f"input: single_cell/results/sc_cross_pairs.tsv sha256 {sha(IN('single_cell','sc_cross_pairs.tsv'))[:12]}",
           f"input: subject_relations/subject_relations.tsv sha256 {sha(os.path.join(AN,'subject_relations','subject_relations.tsv'))[:12]}",
           f"input: fairness_sample_level/results/datasets_used.tsv sha256 {sha(ST('datasets_used.tsv'))[:12]}; derived-unit sizes from derived_lists_2026-10-14/work/stage2/*meta.tsv"]
    fields = ["data_type", "class_or_bulk", "unit_a", "unit_b", "relation", "n_per_group", "rho_t"]
    n_b = len(rows)
    ext = [f"one row per pair of datasets with separate subjects: Spearman rho of the SZ vs control t statistics over shared genes against n_per_group = the smaller of (SZ, control) group sizes over the two datasets. Bulk: {n_b} pairs; single nucleus: Endresz 2026 cohort pairs (7 cohorts), one row per Gemma cell class (14) and pair, relation subject-separate ({len(SN_ROWS)} rows), unchanged from the stored file"]
    write_tsv(os.path.join(dirpath, "fig09_a_agreement_pairs.tsv"), hdr + ext, fields, rows + SN_ROWS)
    write_tsv(os.path.join(dirpath, "fig09_a_agreement_by_n_bin.tsv"), hdr + ["medians of fig09_a_agreement_pairs.tsv within bins of the smaller group size (bins 0-10, 11-30, 31+); single nucleus = median over classes within each cohort pair first"],
              ["donors_per_group_bin", "data_type", "n_pairs", "median_rho_t", "q25", "q75"], bins)
    write_tsv(os.path.join(dirpath, "fig09_a_single_nucleus_by_class.tsv"), hdr[:2] + [hdr[5]] + ["median over the subject-separate Endresz cohort pairs, per Gemma cell class (14); unchanged from the stored file"],
              ["class", "n_pair_rows", "median_rho_t"], crow)
    return BP, bins, keep


# ------------------------------------------------------------------ validation of panel a on the previous base (S28, R1, L2) against the delivered files
tmp = os.path.join(RUN, "out", "validate_panel_a")
BP28, bins28, keep28 = write_a(tmp, "S28_previous", "R1", "L2", "validation against the delivered 64-base files")
for f in ("fig09_a_agreement_pairs.tsv", "fig09_a_agreement_by_n_bin.tsv", "fig09_a_single_nucleus_by_class.tsv"):
    same = body(os.path.join(tmp, f)) == body(ST(f))
    VAL.append((f"V0{2 + ['fig09_a_agreement_pairs.tsv','fig09_a_agreement_by_n_bin.tsv','fig09_a_single_nucleus_by_class.tsv'].index(f)}",
                f"{f}: rebuilt from the rerun with the stored dataset set (22 bulk of the 28 regions), rule R1, level L2; body (without '#' lines) against the delivered 64-base file",
                f"sha256 {sha(ST(f))[:12]}, {len(body(ST(f)))-1} rows", f"{len(body(os.path.join(tmp, f)))-1} rows; body identical = {same}", "PASS" if same else "FAIL"))
VAL.append(("V05", "number of bulk subject-separate pairs on the previous base (S28, rule R1)", "192 (fig09_a legend; delivered file comment)", f"{len(BP28)}", "PASS" if len(BP28) == 192 else "FAIL"))

# ------------------------------------------------------------------ panel b: unchanged
import re
pbs = pdr(IN("single_cell", "sc_split_half_summary.tsv"))
x = pbs[pbs.n_per_group == 20].copy()
cols = ["cell_type", "contrast", "comparison", "n_per_group", "reps", "rho_t", "rho_t_q25", "rho_t_q75", "rho_logfc", "top100_sign_in_other", "top200_fold", "fdr_overlap", "frac_reps_fdr_any"]
x = x[cols]
cls8 = x[x.cell_type != "total"]
med = cls8.groupby(["comparison", "contrast"]).rho_t.agg(["median", "count"]).reset_index()
tmpb = os.path.join(RUN, "out", "validate_panel_b")
write_tsv(os.path.join(tmpb, "fig09_b_psychad_20_per_group.tsv"), ["x"], cols, x.to_dict("records"))
write_tsv(os.path.join(tmpb, "fig09_b_medians_over_8_classes.tsv"), ["x"], ["comparison", "contrast", "median_over_classes", "n_classes"],
          [{"comparison": a, "contrast": b, "median_over_classes": mm, "n_classes": n} for a, b, mm, n in zip(med.comparison, med.contrast, med["median"], med["count"])])
for i, f in enumerate(("fig09_b_psychad_20_per_group.tsv", "fig09_b_medians_over_8_classes.tsv")):
    same = body(os.path.join(tmpb, f)) == body(ST(f))
    VAL.append((f"V0{6 + i}", f"{f}: rebuilt from single_cell/results/sc_split_half_summary.tsv (sha {sha(IN('single_cell','sc_split_half_summary.tsv'))[:12]}); body against the delivered file",
                f"sha256 {sha(ST(f))[:12]}", f"body identical = {same}", "PASS" if same else "FAIL"))

# panel b files into every delivery folder that has panel a (the content does not depend on the bulk base)
HDRB = ["figure: fig09 panel b  (unchanged by the base)",
        "script: sample_level_81_2026-10-15/code/assemble_a.py (panel logic of analyses/figure_exports/export_fig08_fig09_fig10.py fig9_b, unchanged); run 2026-10-15",
        f"input: single_cell/results/sc_split_half_summary.tsv sha256 {sha(IN('single_cell','sc_split_half_summary.tsv'))[:12]}",
        "PsychAD pseudobulk (Mount Sinai and HBCC banks), 20 donors per group per fit, medians over replicates: not part of the bulk dataset set, so the same for the previous base and the 81 base; body identical to the delivered 64-base file"]
for folder in ("base81", "previous_base_28_regions"):
    write_tsv(os.path.join(OUT, folder, "fig09", "fig09_b_psychad_20_per_group.tsv"), HDRB, cols, x.to_dict("records"))
    write_tsv(os.path.join(OUT, folder, "fig09", "fig09_b_medians_over_8_classes.tsv"), HDRB, ["comparison", "contrast", "median_over_classes", "n_classes"],
              [{"comparison": a, "contrast": b, "median_over_classes": mm, "n_classes": n} for a, b, mm, n in zip(med.comparison, med.contrast, med["median"], med["count"])])

# ------------------------------------------------------------------ deliveries
summary = []
combos = [("base81", "S34_base81", "R1", "L2", "81 base: datasets of the studies in the 81 base"),
          ("previous_base_28_regions", "S28_previous", "R1", "L2", "previous base: the 28 dataset-regions (22 bulk), as delivered"),
          ("sensitivity_S29_28plus9", "S29_28plus9", "R1", "L2", "28 dataset-regions plus the 9 derived units"),
          ("sensitivity_S34_B2", "S34_base81", "B2", "L2", "81 base, B2 reading of the bank rule"),
          ("sensitivity_S29_B2", "S29_28plus9", "B2", "L2", "28 dataset-regions plus 9 derived units, B2 reading"),
          ("sensitivity_S28_B2", "S28_previous", "B2", "L2", "previous base, B2 reading"),
          ("sensitivity_S34_allL1", "S34_base81", "R1", "L1", "81 base, every dataset at level L1 (age and sex)"),
          ("sensitivity_S28_allL1", "S28_previous", "R1", "L1", "previous base, every dataset at level L1")]
res_rows = []
for folder, sname, rule, mode, label in combos:
    BP, bins, keep = write_a(os.path.join(OUT, folder, "fig09"), sname, rule, mode, label)
    r = BP.rho_t
    bs = np.median(rng.choice(r.values, size=(10000, len(r)), replace=True), axis=1)
    res_rows.append(dict(folder=folder, dataset_set=sname, rule=rule, level=mode, n_bulk_datasets=len(SETS[sname]), n_pairs_total=len(list(itertools.combinations(SETS[sname], 2))),
                         n_pairs_no_shared_donors=len(BP), median_rho_t=float(r.median()), q25=float(r.quantile(.25)), q75=float(r.quantile(.75)),
                         median_ci_lo_pairs_bootstrap=float(np.percentile(bs, 2.5)), median_ci_hi_pairs_bootstrap=float(np.percentile(bs, 97.5)),
                         **{f"bin_{b['donors_per_group_bin']}_{'bulk' if b['data_type'].startswith('bulk') else 'sn'}": f"{b['n_pairs']}|{b['median_rho_t']:.4f}" if b['median_rho_t'] is not None else f"{b['n_pairs']}|" for b in bins}))
pd.DataFrame(res_rows).to_csv(os.path.join(OUT, "tables_panel_a_summary.tsv"), sep="\t", index=False)

# all pairs with every rule, for the S29 universe (bulk) and all units
allp = []
U = SETS["S29_28plus9"]
for a, b in pair_list(U):
    key = frozenset((a, b))
    r = REL[key]
    n = min(min(nsz[a], nct[a]), min(nsz[b], nct[b]))
    allp.append(dict(unit_a=a, unit_b=b, relation_group=r["group"], evidence=r["evidence"], bank_a=r["bank_a"], bank_b=r["bank_b"], relation_source=r["src"],
                     separate_R1=r["group"] == "subject-separate", separate_B2=b2_separate(key), n_per_group=n,
                     rho_t_L2=rho_pair(a, b, "L2"), rho_t_allL1=rho_pair(a, b, "L1"),
                     in_S28_previous=(a in SETS["S28_previous"] and b in SETS["S28_previous"]), in_S34_base81=(a in SETS["S34_base81"] and b in SETS["S34_base81"])))
pd.DataFrame(allp).to_csv(os.path.join(OUT, "tables_bulk_pairs_all_rules.tsv"), sep="\t", index=False)

# relation group counts per set
cnt = []
for sname in ("S28_previous", "S34_base81", "S29_28plus9"):
    for scope, units in (("bulk datasets", SETS[sname]), ("all units incl. LCM and sorted", [u for u in ALLUNITS[sname]])):
        for a, b in itertools.combinations(units, 2):
            pass
        g = pd.Series([REL[frozenset(p)]["group"] for p in itertools.combinations(units, 2)]).value_counts()
        b2 = sum(b2_separate(frozenset(p)) for p in itertools.combinations(units, 2))
        for k_, v in g.items():
            cnt.append(dict(dataset_set=sname, scope=scope, relation_group=k_, n_pairs=int(v)))
        cnt.append(dict(dataset_set=sname, scope=scope, relation_group="TOTAL pairs", n_pairs=int(g.sum())))
        cnt.append(dict(dataset_set=sname, scope=scope, relation_group="separate under B2 reading", n_pairs=int(b2)))
pd.DataFrame(cnt).to_csv(os.path.join(OUT, "tables_pair_relation_counts.tsv"), sep="\t", index=False)

# rounding check of the new pairs' registry side
chk = []
for lv in ("L1", "L2"):
    r = NP[NP.variant == "reg_rounded_" + lv]
    dif = []
    for t in r.itertuples():
        dif.append(abs(t.rho_t - rho_reg(t.dataset_a, t.dataset_b, lv)))
    chk.append((lv, len(dif), max(dif)))
VAL.append(("V08", "rho_t of registry pairs recomputed from the 5-significant-digit b_de tables (the route used for pairs with derived units) against the full-precision values of the rerun",
            "full-precision in-memory fits", "; ".join(f"{lv}: {n} pairs, max abs difference {d:.2e}" for lv, n, d in chk), "PASS" if max(c[2] for c in chk) < 1e-5 else "CHECK"))
pd.DataFrame(VAL, columns=["id", "check", "expected", "found", "result"]).to_csv(os.path.join(OUT, "VALIDATION_a_b.tsv"), sep="\t", index=False)
print(pd.DataFrame(VAL, columns=["id", "check", "expected", "found", "result"]).to_string())
print(pd.DataFrame(res_rows).iloc[:, :12].to_string())
