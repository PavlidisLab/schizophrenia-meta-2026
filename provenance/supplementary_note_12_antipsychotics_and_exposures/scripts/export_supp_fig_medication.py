#!/usr/bin/env python3
"""Per-panel tidy exports for the medication Supplementary Figure, panels a-i.

Reads existing analysis outputs only (no model fits, no network, no data/raw). Rerun with

    python3 analyses/figure_data_exports/export_supp_fig_medication.py

Writes figure_data/supp_fig_medication/*.tsv, _meta.json and _status_currency.txt.
Inputs in analyses/antipsychotic_medication_analyses/ and analyses/ph_and_confound_analyses/ are read only.
"""
import csv, hashlib, json, os

HERE = os.path.dirname(os.path.abspath(__file__))
AN = os.path.dirname(HERE)
REPO = os.path.dirname(AN)
RUN = "2026-10-09"
SCRIPT = "analyses/figure_data_exports/export_supp_fig_medication.py"
FIG = "supp_fig_medication"
OUT = os.path.expanduser(f"figure_data/{FIG}")
os.makedirs(OUT, exist_ok=True)
MR = os.path.join(AN, "medication", "results")
CR = os.path.join(AN, "confounds", "results")
STRICT = "recurrent (study_brain, uniform q < 0.1)"


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def rel(p):
    return os.path.relpath(p, REPO)


def rd(p):
    with open(p, encoding="utf-8") as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


META = {"script": SCRIPT, "run_id": RUN, "figure": FIG, "inputs": {}, "panels": {}}
S = []


def write(name, panel, fields, rows, inputs, comments=()):
    for p in inputs:
        META["inputs"][rel(p)] = sha(p)
    lines = [f"# figure: {FIG} panel {panel}", f"# script: {SCRIPT}, run {RUN}"]
    lines += [f"# input: {rel(p)} sha256 {sha(p)[:12]}" for p in inputs]
    lines += [f"# {c}" for c in comments]
    with open(os.path.join(OUT, name), "w", encoding="utf-8", newline="") as fh:
        for l in lines:
            fh.write(l + "\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(fields)
        for r in rows:
            w.writerow(["" if r.get(f) is None else r.get(f) for f in fields])
    META["panels"][name] = {"panel": panel, "rows": len(rows)}
    print(name, len(rows))


# a: per-unit and pooled tests, both recurrent sets
p = os.path.join(MR, "set_tests.tsv")
st_rows = rd(p)
for r in st_rows:
    r["dir_minus_null"] = round(float(r["dir_obs"]) - float(r["dir_null_mean"]), 4)
write("a_unit_tests.tsv", "a",
      ["set", "unit", "class", "n_set_genes", "n_measured", "n_p05", "resp_obs", "resp_null_mean", "resp_p",
       "dir_obs", "dir_null_mean", "dir_minus_null", "dir_p", "n_sig_same_as_sz", "n_sig_opposite_sz"],
      st_rows, [p],
      ["resp = mean |z| (drug vs vehicle) over the set's measured genes; dir = mean SZ sign x z (positive: the drug moves genes the way SZ does)",
       "nulls matched on DE-prior quintile x LIBD DLPFC expression decile, 5000 draws, seed 20260926; dir_minus_null computed in the export script",
       "sets: recurrent (study_brain, uniform q < 0.1) = the 176 genes; recurrent (expression-matched q < 0.1) = 99 genes (those present in the drug data)"])

# per-unit sign summary (chronic units), computed here
for label, pred in (("strict set (176 genes)", lambda s: s == STRICT),
                    ("expression-matched set", lambda s: s.startswith("recurrent (expression"))):
    ch = [r for r in st_rows if pred(r["set"]) and r["class"] == "chronic" and not r["unit"].startswith("ALL")]
    above = sum(1 for r in ch if float(r["dir_minus_null"]) > 0)
    sig = [r["unit"] for r in ch if float(r["dir_p"]) < 0.05]
    S.append(f"- chronic units, {label}: direction score above its matched-null mean in {above} of {len(ch)}; units with direction p < 0.05: {', '.join(sig)}")

# b: leave-one-out
p = os.path.join(MR, "leave_one_out.tsv")
write("b_leave_one_out.tsv", "b", ["set", "left_out", "resp_obs", "resp_p", "dir_obs", "dir_p"], rd(p), [p],
      ["chronic units pooled with one unit left out; same matched nulls"])

# c: gene effects (chronic units)
p = os.path.join(MR, "gene_effects.tsv")
ge = rd(p)
write("c_gene_effects.tsv", "c",
      ["set", "gene", "sz_direction", "chronic_units_measured", "chronic_mean_z", "chronic_units_up_p05", "chronic_units_down_p05"],
      ge, [p], ["mean drug z over the chronic units that measured the gene; both recurrent sets"])
for d, sgn in (("down", -1), ("up", 1)):
    s = [r for r in ge if r["set"] == STRICT and r["sz_direction"] == d and r["chronic_mean_z"] != ""]
    S.append(f"- strict set, SZ-{d} genes with mean chronic drug z in the SZ direction: {sum(1 for r in s if float(r['chronic_mean_z']) * sgn > 0)} of {len(s)}")

# d: calibration
p = os.path.join(MR, "calibration.tsv")
write("d_calibration.tsv", "d", ["set", "class", "test", "n", "frac_p_le_05", "frac_p_le_01", "median_p"], rd(p), [p],
      ["random matched sets taken as the observed set (300 each, chronic units), scored with the same test"])

# e: unit design
p = os.path.join(MR, "units_table.tsv")
seen = {}
for r in rd(p):
    u = seen.setdefault(r["unit"], {"unit": r["unit"], "class": r["class"], "geo": r["geo"], "species": r["species"],
                                    "region": r["region"], "drugs": [], "n_drug_n_baseline": [], "duration": r["duration"], "dose": r["dose"]})
    u["drugs"].append(r["drug"])
    u["n_drug_n_baseline"].append(f"{r['n_drug']}/{r['n_baseline']}")
rows = []
for u in seen.values():
    u["drugs"] = "; ".join(u["drugs"])
    u["n_drug_n_baseline"] = "; ".join(u["n_drug_n_baseline"])
    rows.append(u)
write("e_units.tsv", "e", ["unit", "class", "geo", "species", "region", "drugs", "n_drug_n_baseline", "duration", "dose"], rows, [p],
      ["one row per unit; drugs and n (drug/vehicle) listed per contrast in the same order"])

# f: Fromer 2016 rhesus
p = os.path.join(MR, "primate_fromer2016.tsv")
write("f_fromer2016_rhesus.tsv", "f",
      ["drug", "n_control", "n_treated", "genes_affected", "enrichment_p", "n_concordant_with_sz", "concordance_p_enrichment", "concordance_p_depletion"],
      rd(p), [p], ["Fromer et al. 2016 (PMID 27668389) Supplementary Table 2A as read in primate_fromer2016.tsv; no gene list published"])

# g: primate named-gene tally and changed genes
p1 = os.path.join(MR, "primate_tally.tsv")
p2 = os.path.join(MR, "primate_genes_level.tsv")
rows = [{"kind": "tally", "key": r["measure"], "value": r["value"]} for r in rd(p1)]
for r in rd(p2):
    if r["changed"] == "True":
        rows.append({"kind": "changed_gene", "key": r["symbol"], "monkey_change": r["effects"],
                     "relation_to_sz_in_papers": r["relation_to_sz_in_papers"], "recurrent_set_direction": r["rec"],
                     "libd_dlpfc_direction": r["libd"]})
write("g_primate_named_genes.tsv", "g",
      ["kind", "key", "value", "monkey_change", "relation_to_sz_in_papers", "recurrent_set_direction", "libd_dlpfc_direction"], rows, [p1, p2],
      ["named-gene results from papers on antipsychotic-exposed macaques (model extraction, quote-checked, hand-checked); not a random sample of genes",
       "changed_gene rows: genes with a reported change; monkey_change holds the 'effects' field as written"])

# h: dose within SZ
p = os.path.join(CR, "medication_dose.tsv")
write("h_dose_within_sz.tsv", "h",
      ["dataset", "factor", "set", "n_sz", "n_set_measured", "genome_rate_p05", "n_set_p05", "expected_p05", "fisher_p",
       "n_p05_same_sign_as_sz", "score_dose_follows_sz", "sample_perm_p_dose_follows_sz"], rd(p), [p],
      ["GSE21138, chlorpromazine-equivalent dose within SZ cases (the one dataset with dose); the gene-list column of the source is omitted"])

# i: list-based check against published antipsychotic gene lists
p = os.path.join(CR, "geneset_enrichment.tsv")
rows = [r for r in rd(p) if r["factor"] == "antipsychotic" and r["recurrent_set"] == "study_brain"]
rows.sort(key=lambda r: float(r["p_hypergeom"]))
write("i_published_antipsychotic_lists.tsv", "i",
      ["gene_set", "universe", "n_set_in_universe", "n_recurrent", "n_overlap", "expected", "fold", "p_hypergeom",
       "n_with_direction", "n_direction_agree", "q_bh_within_universe_and_set"], rows, [p],
      ["hypergeometric overlap of published antipsychotic gene lists with the 176-gene recurrent set; no DE-prior matching; all rows for factor antipsychotic, sorted by p"])

with open(os.path.join(OUT, "_meta.json"), "w") as fh:
    json.dump(META, fh, indent=1)

# currency
rt = os.path.join(AN, "recurrence.tsv")
S.append(f"- recurrence.tsv sha256 {sha(rt)[:12]} (expected d11761227abd): {'match' if sha(rt).startswith('d11761227abd') else 'MISMATCH'}")
rs = os.path.join(AN, "confounds", "recurrent_sets.tsv")
S.append(f"- confounds/recurrent_sets.tsv sha256 {sha(rs)[:12]} (recorded in the medication report: 18d6535ad65f): {'match' if sha(rs).startswith('18d6535ad65f') else 'MISMATCH'}")
cl = os.path.join(MR, "contrasts_long.tsv.gz")
S.append(f"- medication/results/contrasts_long.tsv.gz sha256 {sha(cl)[:12]} (recorded: 05618fa1f172): {'match' if sha(cl).startswith('05618fa1f172') else 'MISMATCH'}")
with open(os.path.join(OUT, "_status_currency.txt"), "w") as fh:
    fh.write("\n".join(S) + "\n")
print("\n".join(S))
