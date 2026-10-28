"""02_stats.py: per-list, per-study and per-year counts and the year trend. Reads out/rows_classified.tsv.gz and out/list_set.tsv.
Writes out/lists_summary.tsv, out/studies_summary.tsv, out/by_year.tsv, out/by_stratum.tsv, out/special_classes.tsv, out/ensembl_id_drift.tsv, out/trend.tsv"""
import sys, os, csv, gzip, collections, math, random, datetime, json
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hlib as lib
import numpy as np
from scipy import stats
OUT = lib.WD + "/out"
H = lib.HGNC()
today = datetime.date.today()
def hdr(fh, what):
    fh.write("# %s; written %s by scripts/11_stats.py\n" % (what, today))
    fh.write("# input: out/rows_classified.tsv.gz sha256 %s (from work/sztr/data/lists/*.tsv (GitHub commit 746d6740) via scripts/10_classify.py); out/list_set.tsv sha256 %s\n" % (lib.sha(OUT + "/rows_classified.tsv.gz"), lib.sha(OUT + "/list_set.tsv")))
    fh.write("# hgnc: hgnc_complete_set.txt sha256 %s; withdrawn.txt sha256 %s; accessed %s\n" % (H.sha, H.sha_w, open(lib.HGNC_DIR + "/access_date.txt").read().strip()))
rows = list(csv.DictReader((l for l in gzip.open(OUT + "/rows_classified.tsv.gz", "rt") if not l.startswith("#")), delimiter="\t"))
meta = {r["list"]: r for r in csv.DictReader((l for l in open(OUT + "/list_set.tsv") if not l.startswith("#")), delimiter="\t")}
def wilson(k, n, z=1.96):
    if n == 0: return ("", "")
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (round(max(0, c - h), 4), round(min(1, c + h), 4))
CL = ["current", "previous", "alias", "ambiguous", "withdrawn_merged", "withdrawn", "excel_date", "excel_serial", "multi_token", "not_found", "descriptive", "empty"]
# distinct (list, token) with class
tok = {}
for r in rows:
    key = (r["list"], r["gene_reported"].strip().rstrip("*").strip())
    if key not in tok:
        tok[key] = r
bylist = collections.defaultdict(list)
for (n, t), r in tok.items():
    bylist[n].append(r)
def counts(rs):
    c = collections.Counter(r["cls"] for r in rs)
    den = sum(v for k, v in c.items() if k not in ("descriptive", "empty"))
    nar = sum(c[k] for k in lib.NARROW); brd = sum(c[k] for k in lib.BROAD)
    return c, den, nar, brd
# --- per list
with open(OUT + "/lists_summary.tsv", "w") as fh:
    hdr(fh, "per-list counts of distinct printed gene names by class against current HGNC")
    fh.write("# n_names = distinct gene_reported strings; denominator = names that are not descriptive text or empty; narrow = previous+alias+ambiguous+withdrawn+excel dates; broad = narrow + not_found + multi_token; CI = Wilson 95%\n")
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["list", "label", "year", "stratum", "in64", "n_names", "denominator"] + CL + ["narrow", "narrow_frac", "narrow_lo", "narrow_hi", "broad", "broad_frac", "broad_lo", "broad_hi"])
    for n in sorted(bylist):
        rs = bylist[n]; c, den, nar, brd = counts(rs); m = meta[n]
        lo, hi = wilson(nar, den); lo2, hi2 = wilson(brd, den)
        w.writerow([n, m["label"], m["year"], rs[0]["stratum"], m["in64"], len(rs), den] + [c[k] for k in CL] + [nar, round(nar / den, 4) if den else "", lo, hi, brd, round(brd / den, 4) if den else "", lo2, hi2])
# --- per study (in64), pooled distinct names across the study's primary lists
stud = collections.defaultdict(lambda: collections.defaultdict(dict))
for (n, t), r in tok.items():
    m = meta[n]
    if m["in64"] != "yes": continue
    stud[m["study_pmid"]]["names"].setdefault(t, r)
    stud[m["study_pmid"]]["meta"] = m
    stud[m["study_pmid"]]["strata"] = stud[m["study_pmid"]].get("strata", set()) | {r["stratum"]}
S = []
with open(OUT + "/studies_summary.tsv", "w") as fh:
    hdr(fh, "per-study counts (distinct printed names pooled over the study's primary lists), 64 studies of the 746d6740 selection")
    fh.write("# study set = recurrence.scenario_selection('published, brain') on the 746d6740 copy. strata: id types of the study's lists\n")
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["pmid", "label", "year", "strata", "n_names", "denominator"] + CL + ["narrow", "narrow_frac", "narrow_lo", "narrow_hi", "broad", "broad_frac", "broad_lo", "broad_hi"])
    for p, d in sorted(stud.items(), key=lambda x: (int(x[1]["meta"]["year"] or 0), x[0])):
        rs = list(d["names"].values()); c, den, nar, brd = counts(rs); m = d["meta"]
        lo, hi = wilson(nar, den); lo2, hi2 = wilson(brd, den)
        w.writerow([p, m["label"], m["year"], "+".join(sorted(d["strata"])), len(rs), den] + [c[k] for k in CL] + [nar, round(nar / den, 4) if den else "", lo, hi, brd, round(brd / den, 4) if den else "", lo2, hi2])
        S.append(dict(p=p, year=int(m["year"]), den=den, nar=nar, brd=brd, strata=d["strata"]))
# --- per year (pooled names over studies; names counted once per study)
def agg(keyf, rows_by_study, fn):
    pass
with open(OUT + "/by_year.tsv", "w") as fh:
    hdr(fh, "per publication year and 5-year bin, names pooled over the 64 studies (a name counted once per study); Wilson 95% intervals treat names as independent")
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["group", "n_studies", "denominator"] + CL + ["narrow", "narrow_frac", "narrow_lo", "narrow_hi", "broad", "broad_frac", "broad_lo", "broad_hi"])
    groups = collections.defaultdict(list)
    for p, d in stud.items():
        y = int(d["meta"]["year"]); b = "%d-%d" % (y - (y - 2000) % 5, y - (y - 2000) % 5 + 4)
        for g in (str(y), "bin " + b):
            groups[g].append(d)
    for g in sorted(groups, key=lambda s: (s.startswith("bin"), s)):
        rs = [r for d in groups[g] for r in d["names"].values()]; c, den, nar, brd = counts(rs)
        lo, hi = wilson(nar, den); lo2, hi2 = wilson(brd, den)
        w.writerow([g, len(groups[g]), den] + [c[k] for k in CL] + [nar, round(nar / den, 4) if den else "", lo, hi, brd, round(brd / den, 4) if den else "", lo2, hi2])
# --- per stratum (64 studies, names pooled by (study,list-stratum))
with open(OUT + "/by_stratum.tsv", "w") as fh:
    hdr(fh, "per identifier stratum of the list (symbol = no identifier column; ensembl; probe; accession; entrez); 94 lists = 64-study lists plus recurrence-selection lists")
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["set", "stratum", "n_lists", "denominator"] + CL + ["narrow", "narrow_frac", "narrow_lo", "narrow_hi", "broad", "broad_frac", "broad_lo", "broad_hi"])
    for setn, f in (("64 studies", lambda n: meta[n]["in64"] == "yes"), ("all lists", lambda n: True)):
        for st in ("symbol", "accession", "probe", "ensembl", "entrez"):
            ns = [n for n in bylist if f(n) and bylist[n][0]["stratum"] == st]
            rs = [r for n in ns for r in bylist[n]]
            if not rs: continue
            c, den, nar, brd = counts(rs); lo, hi = wilson(nar, den); lo2, hi2 = wilson(brd, den)
            w.writerow([setn, st, len(ns), den] + [c[k] for k in CL] + [nar, round(nar / den, 4) if den else "", lo, hi, brd, round(brd / den, 4) if den else "", lo2, hi2])
# --- special classes
with open(OUT + "/special_classes.tsv", "w") as fh:
    hdr(fh, "named problem classes among distinct (list, name) pairs, 64-study lists; counts of names and of lists affected")
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["problem_class", "definition", "names", "lists_affected", "studies_affected", "detail"])
    rs62 = [(k, r) for k, r in tok.items() if meta[k[0]]["in64"] == "yes"]
    def emit(name, defn, pred, detail=lambda sel: ""):
        sel = [(k, r) for k, r in rs62 if pred(r)]
        w.writerow([name, defn, len(sel), len({k[0] for k, r in sel}), len({meta[k[0]]["study_pmid"] for k, r in sel}), detail(sel)])
    emit("excel_date", "name matches d-Mon or Mon-d, or a five-digit 4xxxx serial", lambda r: r["cls"] in ("excel_date", "excel_serial"),
         lambda sel: "recoverable_unique %d; recoverable_ambiguous %d; unrecoverable/serial %d; names: %s" % (
             sum(1 for k, r in sel if "recoverable_unique" in r["cls_flags"]), sum(1 for k, r in sel if "recoverable_ambiguous" in r["cls_flags"]),
             sum(1 for k, r in sel if r["cls"] == "excel_serial" or "unrecoverable" in r["cls_flags"]), ",".join(sorted({k[1] for k, r in sel}))))
    emit("ambiguous_shared_name", "previous or alias symbol shared by more than one current gene", lambda r: r["cls"] == "ambiguous",
         lambda sel: "with exactly one previous-symbol candidate (resolvable by previous-over-alias precedence) %d; names: %s" % (sum(1 for k, r in sel if r["cls_flags"].startswith("prev:1")), ",".join(sorted({k[1] for k, r in sel})[:40])))
    emit("current_but_also_old_name_of_another_gene", "approved symbol that is also a previous or alias symbol of a different gene", lambda r: r["cls"] == "current" and "also_prev_or_alias" in r["cls_flags"],
         lambda sel: "names: %s" % ",".join(sorted({k[1] for k, r in sel})[:40]))
    emit("multi_gene_token", "name contains / , or ;", lambda r: r["cls"] == "multi_token", lambda sel: ",".join(sorted({k[1] for k, r in sel})[:30]))
    emit("withdrawn", "name only in HGNC withdrawn.txt", lambda r: r["cls"] in ("withdrawn", "withdrawn_merged"), lambda sel: "merged into a current gene %d; entry withdrawn %d" % (sum(1 for k, r in sel if r["cls"] == "withdrawn_merged"), sum(1 for k, r in sel if r["cls"] == "withdrawn")))
    emit("loc_ncbi_only", "LOC<number> name not in HGNC", lambda r: r["cls"] == "not_found" and "loc_ncbi" in r["cls_flags"])
    emit("not_found_other", "name not in HGNC approved, previous, alias or withdrawn", lambda r: r["cls"] == "not_found" and "loc_ncbi" not in r["cls_flags"])
    emit("descriptive_text", "name contains a space or is longer than 40 characters (a protein name or annotation, not a symbol)", lambda r: r["cls"] == "descriptive")
    emit("readthrough_gene", "resolved gene has HGNC locus_type readthrough or 'readthrough' in its name", lambda r: "readthrough" in r["gene_flags"],
         lambda sel: "classes: %s; names: %s" % (dict(collections.Counter(r["cls"] for k, r in sel)), ",".join(sorted({k[1] for k, r in sel}))))
    emit("antisense_gene", "resolved gene symbol ends in -AS<n> or its name contains 'antisense'", lambda r: "antisense" in r["gene_flags"],
         lambda sel: "classes: %s" % dict(collections.Counter(r["cls"] for k, r in sel)))
    # case differs
    emit("case_differs_from_approved", "name matches an approved symbol only when case is ignored", lambda r: r["cls"] == "current" and "case_differs" in r["cls_flags"], lambda sel: ",".join(sorted({k[1] for k, r in sel})[:30]))
    # dated: previous names retired before vs after publication year
    def yr(k):
        return int(meta[k[0]]["year"])
    def changed_after(r, y):
        ids = r["hgnc_ids"].split(",")
        if len(ids) != 1 or not ids[0]: return None
        d = H.info[ids[0]]["date_symbol_changed"].strip('"')
        return int(d[:4]) > y if d else None
    for lab, flag in (("previous_or_alias_symbol_changed_after_publication_year", True), ("previous_or_alias_symbol_changed_in_or_before_publication_year", False)):
        sel = [(k, r) for k, r in rs62 if r["cls"] in ("previous", "alias") and changed_after(r, yr(k)) is flag]
        w.writerow([lab, "HGNC date_symbol_changed (latest symbol change of the gene) vs the study's year; a gene whose symbol changed twice shows only its latest date", len(sel), len({k[0] for k, r in sel}), len({meta[k[0]]["study_pmid"] for k, r in sel}), ""])
# --- Ensembl and entrez id drift
with open(OUT + "/ensembl_id_drift.tsv", "w") as fh:
    hdr(fh, "ENSG ids in the ensembl-stratum lists (64 studies): present in current HGNC or not, and whether the printed symbol resolves to the same gene as the id")
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["list", "label", "year", "distinct_ensg", "in_hgnc", "not_in_hgnc", "not_in_hgnc_frac", "not_in_hgnc_lo", "not_in_hgnc_hi", "printed_symbol_resolves", "agree_with_id", "disagree_with_id"])
    seen = collections.defaultdict(dict)
    for r in rows:
        if r["ens_cls"] and meta[r["list"]]["in64"] == "yes":
            seen[r["list"]].setdefault(r["id_reported"].split(".")[0], r)
    tot = [0, 0, 0, 0]
    for n in sorted(seen):
        rs = list(seen[n].values()); a = sum(1 for r in rs if r["ens_cls"] == "in_hgnc"); b = len(rs) - a
        ag = sum(1 for r in rs if r["ens_symbol_agree"] == "yes"); dg = sum(1 for r in rs if r["ens_symbol_agree"] == "no")
        lo, hi = wilson(b, len(rs))
        w.writerow([n, meta[n]["label"], meta[n]["year"], len(rs), a, b, round(b / len(rs), 4), lo, hi, ag + dg, ag, dg])
        tot = [tot[0] + len(rs), tot[1] + b, tot[2] + ag, tot[3] + dg]
    lo, hi = wilson(tot[1], tot[0])
    w.writerow(["TOTAL", "", "", tot[0], tot[0] - tot[1], tot[1], round(tot[1] / tot[0], 4), lo, hi, tot[2] + tot[3], tot[2], tot[3]])
# --- trend (Spearman, study level) with bootstrap CI over studies
def spear(X, Y):
    r, p = stats.spearmanr(X, Y); return r, p
res = []
def trend(label, ss, key):
    ss = [s for s in ss if s["den"] > 0]
    X = np.array([s["year"] for s in ss]); Y = np.array([s[key] / s["den"] for s in ss])
    r, p = spear(X, Y)
    rng = np.random.default_rng(20261014); bs = []
    for _ in range(10000):
        i = rng.integers(0, len(ss), len(ss))
        if len(set(X[i])) < 2 or len(set(Y[i])) < 2: continue
        bs.append(stats.spearmanr(X[i], Y[i])[0])
    lo, hi = np.percentile(bs, [2.5, 97.5])
    # slope (fraction per 10 years) by least squares, as effect size
    sl = np.polyfit(X, Y, 1)[0] * 10
    res.append([label, key, len(ss), round(r, 3), round(lo, 3), round(hi, 3), "%.3g" % p, round(sl, 4)])
for key in ("nar", "brd"):
    trend("all 64 studies", S, key)
    trend("studies with >= 20 names", [s for s in S if s["den"] >= 20], key)
    trend("symbol-only studies (all lists have no identifier column)", [s for s in S if s["strata"] == {"symbol"}], key)
    trend("studies with any ensembl, probe, accession or entrez list", [s for s in S if s["strata"] - {"symbol"}], key)
with open(OUT + "/trend.tsv", "w") as fh:
    hdr(fh, "Spearman correlation between publication year and the fraction of distinct printed names affected, study level; CI = percentile bootstrap over studies (10000 resamples, seed 20261014); slope = least-squares change in fraction per 10 years")
    fh.write("# nar = narrow definition (previous, alias, ambiguous, withdrawn, excel dates); brd = broad (narrow + not_found + multi_token)\n")
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["subset", "definition", "n_studies", "spearman_rho", "ci_lo", "ci_hi", "p", "slope_per_10y"])
    w.writerows(res)
for r in res: print(r)
