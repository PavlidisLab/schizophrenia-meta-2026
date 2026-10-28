#!/usr/bin/env python3
"""Pairwise reproducibility and recurrence of SZ, sex and age effects across the same datasets.

    python3 scripts/pairs.py [<matched workdir>]      (run from controls_age_sex/)
    python3 scripts/pairs.py --summary-only           (summary.tsv only, from pairs.tsv, pairs_matched.tsv, subsample/rep_metrics.tsv.gz)

Inputs (read only):
  data/datasets.tsv                 registry (scripts/build_covariates.py): collection cluster(s), compartment
  de/<dataset>.tsv                  full-data limma tables (scripts/de_contrasts.R full)
  results/de_n.tsv                  n per contrast
  results/sex_check_samples.tsv     the samples used per dataset (and their sex-label check)
  <workdir>/rep_XX/<dataset>.bin    matched-subsample t and p (scripts/de_contrasts.R matched <workdir>); optional
  ../subject_relations/subject_relations.tsv   subject relation of every dataset pair (the project's one table)
  ../../data/reference/hgnc_complete_set.txt  chromosome (HGNC 'location'); autosomal = location starting 1-22
  DE prior                                   recurrence.de_prior_weights (global)
Outputs:
  pairs.tsv                 one row per dataset pair x contrast (full data)
  pairs_matched.tsv         one row per dataset pair x contrast (matched subsamples: mean and SD over repeats)
  subsample/rep_metrics.tsv.gz   per repeat x pair x contrast metrics
  summary.tsv               per contrast x analysis x tier x subject relation (column donor_class = relation_group of
                            ../subject_relations/subject_relations.tsv): distribution of pair metrics
  results/recurrence.tsv    recurrence of genes across collection clusters vs uniform and DE-prior nulls
  results/recurrence_genes_<contrast>.tsv   per gene: units in top 200 / at q < 0.05, direction, Stouffer z
  results/sex_sanity.tsv    X/Y marker genes in the full-data sex contrast, per dataset
  results/known_sets.tsv    recurrent age / autosomal-sex genes vs Mistry & Pavlidis 2010 top-10 sets (../confounds/genesets)

Pair metrics, over the genes measured in both datasets (restricted to autosomal genes for '_auto' contrasts):
  rho_t = Spearman correlation of the moderated t; rho_logfc = the same for log2 fold change (full data only);
  top-N (N = 100, 200, 500): each dataset's N genes with the smallest p (ties by larger |t|); overlap, expected N^2 / |U|,
  fold = overlap / expected, hypergeometric p (h1.hyper_sf); dir_N = fraction of the overlapping genes with the same sign;
  pi1_ab = 1 - pi0 of B's p values over A's top 200 genes, pi0 = Storey's estimator at lambda = 0.5 (capped at 1);
  pi1 = mean of pi1_ab and pi1_ba. Full data also: pi1q_* from R qvalue with defaults (h1.run_pi1).
Subject relation of a pair: relation_group of ../subject_relations/subject_relations.tsv (subject_relations.py; tracker,
  Vocabulary). The earlier classifier here (same series / same collection / pH-and-PMI fingerprint, columns donor_relation
  and fingerprint_* of pairs.tsv) is removed (2026-09-25); those columns in pairs.tsv files written before then are not read.
"""
import csv, datetime, gzip, hashlib, json, math, os, re, sys
import numpy as np
from scipy.stats import rankdata

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
GL = os.path.join(ROOT, "..")
DATA = os.path.normpath(os.path.join(GL, "..", "data"))  # repo data/ (lists, manifests, reference, raw)
sys.path.insert(0, GL)
import recurrence  # noqa: E402
import h1  # noqa: E402
sys.path.insert(0, os.path.join(GL, "subject_relations"))
import relations as SR  # noqa: E402

SEED = 20260925
NS = (100, 200, 500)
NPI = 200
FULL = [("sz_auto", "sz", "auto"), ("sz_all", "sz", "all"), ("sex_auto", "sex", "auto"), ("sex_all", "sex", "all"),
        ("sex_ctl_auto", "sex_ctl", "auto"), ("sex_ctl_all", "sex_ctl", "all"), ("age_auto", "age", "auto"), ("age_all", "age", "all"),
        ("age_ctl_auto", "age_ctl", "auto"), ("age2_auto", "age2", "auto")]
MCON = ["sz_m", "sex_m", "age2_m", "age_m"]  # order in the .bin files
MATCHED = [("sz_m_auto", "sz_m", "auto"), ("sz_m_all", "sz_m", "all"), ("sex_m_auto", "sex_m", "auto"), ("sex_m_all", "sex_m", "all"),
           ("age2_m_auto", "age2_m", "auto"), ("age_m_auto", "age_m", "auto")]
MARKERS = {"XIST": -1, "RPS4Y1": 1, "DDX3Y": 1, "KDM5D": 1, "UTY": 1, "EIF1AY": 1}  # expected sign of male - female


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def read_tsv(p):
    with open(p) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def chrom_classes():
    out = {}
    for r in read_tsv(os.path.join(DATA, "reference", "hgnc_complete_set.txt")):
        loc = r.get("location") or ""
        m = re.match(r"^(\d+)", loc)
        out[r["symbol"]] = "auto" if m and 1 <= int(m.group(1)) <= 22 else ("X" if loc.startswith("X") else "Y" if loc.startswith("Y") else "other")
    return out


def load_de(ds):
    rows = read_tsv(os.path.join(ROOT, "de", f"{ds}.tsv"))
    genes = [r["gene"] for r in rows]
    cols = {}
    for k in rows[0]:
        if k.startswith(("t_", "p_", "logfc_", "q_")):
            cols[k] = np.array([float(r[k]) if r[k] not in ("", "NA") else np.nan for r in rows])
    return genes, cols


def topn(t, p, n):
    o = np.lexsort((-np.abs(t), p))
    return o[:n]


def storey_pi1(p):
    if len(p) == 0:
        return np.nan
    return 1.0 - min(1.0, float(np.mean(p > 0.5)) / 0.5)


def pair_metrics(tA, pA, tB, pB, full=False, lA=None, lB=None):
    U = len(tA)
    out = {"n_genes_shared": U}
    rA, rB = rankdata(tA), rankdata(tB)
    out["rho_t"] = float(np.corrcoef(rA, rB)[0, 1])
    if full and lA is not None:
        out["rho_logfc"] = float(np.corrcoef(rankdata(lA), rankdata(lB))[0, 1])
    for n in NS:
        a, b = topn(tA, pA, n), topn(tB, pB, n)
        ov = np.intersect1d(a, b)
        e = n * n / U
        out[f"overlap{n}"] = len(ov)
        out[f"exp{n}"] = round(e, 3)
        out[f"fold{n}"] = len(ov) / e
        out[f"dir{n}"] = float(np.mean(np.sign(tA[ov]) == np.sign(tB[ov]))) if len(ov) else np.nan
        if full:
            out[f"hyper_p{n}"] = h1.hyper_sf(len(ov), U, n, n)
    a, b = topn(tA, pA, NPI), topn(tB, pB, NPI)
    out["pi1_ab"], out["pi1_ba"] = storey_pi1(pB[a]), storey_pi1(pA[b])
    out["pi1"] = (out["pi1_ab"] + out["pi1_ba"]) / 2
    out["_topA"], out["_topB"] = a, b
    return out


# ---------------------------------------------------------------- recurrence
def recur(lists, meas, genes, hitn):
    """lists/meas: per unit, sets of gene indices; returns obs counts, uniform and DE-prior tails."""
    G, K = len(genes), len(lists)
    wall = np.array([hitn.get(g, np.nan) + 1.0 for g in genes])
    wall = np.where(np.isnan(wall), np.nanmedian(wall), wall)
    P = np.zeros((G, K)); PP = np.zeros((G, K)); obs = np.zeros(G, dtype=int)
    for j, (L, M) in enumerate(zip(lists, meas)):
        idx = np.array(sorted(M))
        P[idx, j] = len(L) / len(M)
        PP[idx, j] = recurrence.prior_probs(float(len(L)), wall[idx])
        obs[list(L)] += 1
    pu, tail = recurrence.upper_tail(recurrence.pb_matrix(P), obs)
    pp, tailP = recurrence.upper_tail(recurrence.pb_matrix(PP), obs)
    return obs, pu, tail, pp, tailP


SUBJECT_RELATIONS = SR.GROUPS
_REL = None


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def subject_relation(r):
    """Subject relation of a pair row: relation_group of ../subject_relations/subject_relations.tsv (tracker, Vocabulary)."""
    global _REL
    if _REL is None:
        _REL = SR.pairs()
    return _REL[(r["dataset_a"], r["dataset_b"])]["relation_group"]


def summarize(rows, mrows, reprow, reps, today):
    """summary.tsv from the full-data and matched pair rows (called by main, or alone with --summary-only)."""
    summ = []
    # one group per tier x subject relation; relations are never pooled (tracker, Vocabulary, 2026-09-25)
    relmap = {(r["dataset_a"], r["dataset_b"]): subject_relation(r) for r in rows}
    groups = [(tier, cls, (lambda r, t=tier, c=cls: (t == "all brain" or r["tier"] == "bulk") and relmap[(r["dataset_a"], r["dataset_b"])] == c))
              for cls in SUBJECT_RELATIONS for tier in ("bulk", "all brain")]
    for analysis, data, specs, suffix, ref in (("full data", rows, FULL, "", "sz_auto"), ("matched", mrows, MATCHED, "_mean", "sz_m_auto")):
        if not data:
            continue
        for tier, dcls, f in groups:
            sel = [r for r in data if f(r)]
            refv = {(r["dataset_a"], r["dataset_b"]): r[f"rho_t{suffix}"] for r in sel if r["contrast"] == ref}
            reff = {(r["dataset_a"], r["dataset_b"]): r[f"fold200{suffix}"] for r in sel if r["contrast"] == ref}
            for key, _, _ in specs:
                s = [r for r in sel if r["contrast"] == key]
                if not s:
                    continue
                row = {"analysis": analysis, "contrast": key, "tier": tier, "donor_class": dcls, "n_pairs": len(s),
                       "n_datasets": len({r["dataset_a"] for r in s} | {r["dataset_b"] for r in s})}
                for k in ("rho_t", "fold100", "fold200", "fold500", "dir200", "pi1"):
                    v = np.array([r[k + suffix] for r in s], float)
                    v = v[~np.isnan(v)]
                    row[f"{k}_median"] = float(np.median(v)) if len(v) else np.nan
                    row[f"{k}_q25"] = float(np.percentile(v, 25)) if len(v) else np.nan
                    row[f"{k}_q75"] = float(np.percentile(v, 75)) if len(v) else np.nan
                row["frac_rho_t_gt0"] = float(np.mean([r["rho_t" + suffix] > 0 for r in s]))
                d = [r["rho_t" + suffix] - refv[(r["dataset_a"], r["dataset_b"])] for r in s if (r["dataset_a"], r["dataset_b"]) in refv]
                row["rho_t_minus_sz_median"] = float(np.median(d)) if d else np.nan
                row["n_pairs_rho_t_above_sz"] = int(sum(x > 0 for x in d))
                d2 = [r["fold200" + suffix] - reff[(r["dataset_a"], r["dataset_b"])] for r in s if (r["dataset_a"], r["dataset_b"]) in reff]
                row["fold200_minus_sz_median"] = float(np.median(d2)) if d2 else np.nan
                row["n_pairs_fold200_above_sz"] = int(sum(x > 0 for x in d2))
                if analysis == "matched" and reprow:
                    pairset = {(r["dataset_a"], r["dataset_b"]) for r in s}
                    med = []
                    for rep in reps:
                        v = [x["rho_t"] for x in reprow if x["rep"] == rep and x["contrast"] == key and (x["dataset_a"], x["dataset_b"]) in pairset]
                        med.append(np.median(v))
                    row["rho_t_median_rep_q025"] = float(np.percentile(med, 2.5)); row["rho_t_median_rep_q975"] = float(np.percentile(med, 97.5))
                summ.append(row)
    write(os.path.join(ROOT, "summary.tsv"), summ, [
        f"# scripts/pairs.py run {today}; distribution over dataset pairs of each metric (median, 25th and 75th percentiles)",
        "# matched rows summarise per-pair means over repeats; rho_t_median_rep_q025/q975 = 2.5th and 97.5th percentiles over repeats of the median over pairs",
        "# rho_t_minus_sz_median / n_pairs_rho_t_above_sz: per pair, contrast minus the SZ contrast of the same analysis (sz_auto or sz_m_auto), same pairs; "
        "pairs are not independent (a dataset enters many pairs)",
        f"# donor_class = subject relation of the pair: relation_group of ../subject_relations/subject_relations.tsv ({SR.built()[0]}; sha256 "
        f"{SR.built()[1][:12]}). Relations are not pooled."])



def summary_only():
    """Rebuild summary.tsv without the fingerprint, pi1q (R) or recurrence steps. Full-data pair metrics are recomputed from de/*.tsv
    with pair_metrics (pairs.tsv holds them rounded to 4 significant digits); matched per-pair means come from
    subsample/rep_metrics.tsv.gz (full precision); pair relation from ../subject_relations/subject_relations.tsv (subject_relation)."""
    today = datetime.date.today().isoformat()
    MET = ("rho_t", "fold100", "fold200", "fold500", "dir200", "pi1")
    reg = read_tsv(os.path.join(ROOT, "data", "datasets.tsv"))
    names = [r["dataset"] for r in reg]
    chrom = chrom_classes()
    DE = {d: load_de(d) for d in names}
    gidx = {d: {g: i for i, g in enumerate(DE[d][0])} for d in names}
    prow = {}
    for r in read_tsv(os.path.join(ROOT, "pairs.tsv")):
        prow.setdefault((r["dataset_a"], r["dataset_b"]), r)
    keep = ("dataset_a", "dataset_b", "tier")
    rows = []
    for (a, b), pr in prow.items():
        common = [g for g in DE[a][0] if g in gidx[b]]
        ia = np.array([gidx[a][g] for g in common]); ib = np.array([gidx[b][g] for g in common])
        au = np.array([chrom.get(g) == "auto" for g in common])
        sh = {"all": (ia, ib), "auto": (ia[au], ib[au])}
        for key, con, gs in FULL:
            ta, tb = DE[a][1].get(f"t_{con}"), DE[b][1].get(f"t_{con}")
            if ta is None or tb is None:
                continue
            xa, xb = sh[gs]
            m = pair_metrics(ta[xa], DE[a][1][f"p_{con}"][xa], tb[xb], DE[b][1][f"p_{con}"][xb], full=True,
                             lA=DE[a][1][f"logfc_{con}"][xa], lB=DE[b][1][f"logfc_{con}"][xb])
            row = {k: pr[k] for k in keep}
            row["contrast"] = key
            row.update({k: m[k] for k in MET})
            rows.append(row)
    with gzip.open(os.path.join(ROOT, "subsample", "rep_metrics.tsv.gz"), "rt") as fh:
        reprow = list(csv.DictReader(fh, delimiter="\t"))
    acc = {}
    for r in reprow:
        for k in MET:
            r[k] = _num(r[k])
        acc.setdefault((r["dataset_a"], r["dataset_b"], r["contrast"]), []).append(r)
    mrows = []
    for (a, b, key), lst in acc.items():
        row = {k: prow[(a, b)][k] for k in keep}
        row["contrast"] = key
        for k in MET:
            row[f"{k}_mean"] = float(np.nanmean(np.array([x[k] for x in lst], float)))
        mrows.append(row)
    reps = sorted({r["rep"] for r in reprow})
    summarize(rows, mrows, reprow, reps, today)
    print("wrote summary.tsv")


def main():
    work = sys.argv[1] if len(sys.argv) > 1 else None
    today = datetime.date.today().isoformat()
    reg = read_tsv(os.path.join(ROOT, "data", "datasets.tsv"))
    R = {r["dataset"]: r for r in reg}
    names = [r["dataset"] for r in reg]
    chrom = chrom_classes()
    den = {(r["dataset"], r["contrast"]): r for r in read_tsv(os.path.join(ROOT, "results", "de_n.tsv"))}
    used = {}
    for r in read_tsv(os.path.join(ROOT, "results", "sex_check_samples.tsv")):
        used.setdefault(r["dataset"], []).append(r["gsm"])
    DE = {d: load_de(d) for d in names}
    gidx = {d: {g: i for i, g in enumerate(DE[d][0])} for d in names}
    auto = {d: np.array([chrom.get(g) == "auto" for g in DE[d][0]]) for d in names}

    # subject relation per pair: ../subject_relations/subject_relations.tsv (the earlier classifier is removed)
    pairs = [(a, b) for i, a in enumerate(names) for b in names[i + 1:]]
    SRP = SR.pairs()
    rel = {}
    for a, b in pairs:
        tier = "bulk" if R[a]["compartment"] == "bulk" and R[b]["compartment"] == "bulk" else "cell-enriched"
        rel[(a, b)] = {"relation": SRP[(a, b)]["relation_group"], "tier": tier, "k": SRP[(a, b)]["k_shared"],
                       "evidence": SRP[(a, b)]["evidence"]}

    # shared-gene index per pair and gene set
    shared = {}
    for a, b in pairs:
        common = [g for g in DE[a][0] if g in gidx[b]]
        ia = np.array([gidx[a][g] for g in common]); ib = np.array([gidx[b][g] for g in common])
        au = np.array([chrom.get(g) == "auto" for g in common])
        shared[(a, b)] = {"all": (ia, ib, np.array(common)), "auto": (ia[au], ib[au], np.array(common)[au])}

    base_cols = lambda a, b: {"dataset_a": a, "dataset_b": b, "gse_a": R[a]["gse"], "gse_b": R[b]["gse"], "cluster_a": R[a]["cluster"],
                              "cluster_b": R[b]["cluster"], "compartment_a": R[a]["compartment"], "compartment_b": R[b]["compartment"],
                              "tier": rel[(a, b)]["tier"], "donor_relation": rel[(a, b)]["relation"],
                              "subjects_shared": rel[(a, b)]["k"], "subject_evidence": rel[(a, b)]["evidence"]}
    # ---------------- full data
    rows, pijobs = [], {}
    for a, b in pairs:
        for key, con, gs in FULL:
            ta, tb = DE[a][1].get(f"t_{con}"), DE[b][1].get(f"t_{con}")
            if ta is None or tb is None:
                continue
            ia, ib, gn = shared[(a, b)][gs]
            m = pair_metrics(ta[ia], DE[a][1][f"p_{con}"][ia], tb[ib], DE[b][1][f"p_{con}"][ib], full=True,
                             lA=DE[a][1][f"logfc_{con}"][ia], lB=DE[b][1][f"logfc_{con}"][ib])
            j = f"{a}|{b}|{key}"
            pijobs[j + "|ab"] = list(DE[b][1][f"p_{con}"][ib][m["_topA"]])
            pijobs[j + "|ba"] = list(DE[a][1][f"p_{con}"][ia][m["_topB"]])
            del m["_topA"], m["_topB"]
            row = base_cols(a, b)
            row.update({"contrast": key, "model_contrast": con, "genes": gs, "n_a": den[(a, con)]["n"], "n_b": den[(b, con)]["n"]})
            row.update(m)
            rows.append(row)
    piq = h1.run_pi1(pijobs)
    for row in rows:
        j = f"{row['dataset_a']}|{row['dataset_b']}|{row['contrast']}"
        row["pi1q_ab"] = piq.get(j + "|ab", (None,))[0]
        row["pi1q_ba"] = piq.get(j + "|ba", (None,))[0]
    write(os.path.join(ROOT, "pairs.tsv"), rows, [
        f"# scripts/pairs.py run {today}; full-data limma tables de/*.tsv; one row per dataset pair x contrast",
        "# contrasts: sz = SZ vs control (~ sz + age + sex); sex = male vs female, all subjects (~ male + sz + age); sex_ctl = controls only; "
        "age = per year, all subjects (~ age + sz + sex); age_ctl = controls only; age2 = older vs younger half; _auto = autosomal genes only (HGNC location 1-22), _all = all genes",
        f"# metrics over genes measured in both datasets; top-N by p; fold = overlap / (N^2 / n_genes_shared); pi1 = Storey lambda 0.5 over A's top {NPI} in B's p (and reverse); pi1q = R qvalue defaults",
        f"# input sha256: data/datasets.tsv {sha(os.path.join(ROOT, 'data', 'datasets.tsv'))[:12]}; results/de_n.tsv {sha(os.path.join(ROOT, 'results', 'de_n.tsv'))[:12]}"])

    # ---------------- matched subsamples
    mrows, reprow = [], []
    mdp = os.path.join(ROOT, "results", "matched_design.tsv")
    mdesign = {r["dataset"]: r for r in read_tsv(mdp)} if os.path.exists(mdp) else {}
    mnames = [d for d in names if d in mdesign and mdesign[d]["used"] == "TRUE"]
    reps = sorted(x for x in os.listdir(work) if x.startswith("rep_")) if work and os.path.isdir(work) else []
    if reps:
        mg = {d: [l.strip() for l in open(os.path.join(work, f"{d}_genes.txt"))] for d in mnames}
        for d in mnames:
            assert mg[d] == DE[d][0], d
        acc = {}
        for rep in reps:
            M = {}
            for d in mnames:
                x = np.fromfile(os.path.join(work, rep, f"{d}.bin"), dtype="<f4").astype(float)
                G = len(mg[d])
                M[d] = {c: (x[(2 * k) * G:(2 * k + 1) * G], x[(2 * k + 1) * G:(2 * k + 2) * G]) for k, c in enumerate(MCON)}
            for a, b in pairs:
                if a not in M or b not in M:
                    continue
                for key, con, gs in MATCHED:
                    ia, ib, _ = shared[(a, b)][gs]
                    (ta, pa), (tb, pb) = M[a][con], M[b][con]
                    m = pair_metrics(ta[ia], pa[ia], tb[ib], pb[ib])
                    rr = {"rep": rep, "dataset_a": a, "dataset_b": b, "contrast": key, "tier": rel[(a, b)]["tier"],
                          "donor_relation": rel[(a, b)]["relation"]}
                    rr.update({k: m[k] for k in ("rho_t", "fold100", "fold200", "fold500", "dir200", "pi1")})
                    reprow.append(rr)
                    acc.setdefault((a, b, key), []).append(rr)
            print("matched", rep)
        for (a, b, key), lst in acc.items():
            row = base_cols(a, b)
            row.update({"contrast": key, "reps": len(lst), "n_matched_a": mdesign[a]["n_matched"], "n_matched_b": mdesign[b]["n_matched"],
                        "a_b_a": f"{mdesign[a]['a']}/{mdesign[a]['b']}", "a_b_b": f"{mdesign[b]['a']}/{mdesign[b]['b']}"})
            for k in ("rho_t", "fold100", "fold200", "fold500", "dir200", "pi1"):
                v = np.array([x[k] for x in lst], float)
                row[f"{k}_mean"] = float(np.nanmean(v)); row[f"{k}_sd"] = float(np.nanstd(v, ddof=1))
            mrows.append(row)
        write(os.path.join(ROOT, "pairs_matched.tsv"), mrows, [
            f"# scripts/pairs.py run {today}; matched subsamples ({len(reps)} repeats; scripts/de_contrasts.R matched); mean and SD over repeats per pair",
            "# per dataset every contrast uses the same group sizes a/b (results/matched_design.tsv); sz_m ~ sz + age + sex; sex_m ~ male + sz + age; "
            "age2_m older vs younger half ~ old + sz + sex; age_m continuous age on a + b random subjects ~ age + sz + sex"])
        os.makedirs(os.path.join(ROOT, "subsample"), exist_ok=True)
        with gzip.open(os.path.join(ROOT, "subsample", "rep_metrics.tsv.gz"), "wt") as fh:
            w = csv.DictWriter(fh, fieldnames=list(reprow[0].keys()), delimiter="\t")
            w.writeheader(); w.writerows(reprow)

    # ---------------- summary
    summarize(rows, mrows, reprow, reps, today)

    # ---------------- recurrence over collection clusters (bulk), one representative per cluster
    hitn = recurrence.de_prior_weights()
    bulk = [d for d in names if R[d]["compartment"] == "bulk"]
    rep_of, rep_m = {}, {}
    for d in sorted(bulk, key=lambda d: -int(den[(d, "sz")]["n"])):
        rep_of.setdefault(R[d]["cluster"], d)
        if d in mnames:
            rep_m.setdefault(R[d]["cluster"], d)
    units, units_m = sorted(rep_of.values()), sorted(rep_m.values())
    rrows, gene_tables = [], {}
    specs = [("sz_auto", "sz", "auto"), ("sex_auto", "sex", "auto"), ("sex_all", "sex", "all"), ("sex_ctl_auto", "sex_ctl", "auto"),
             ("age_auto", "age", "auto"), ("age_ctl_auto", "age_ctl", "auto"), ("age2_auto", "age2", "auto")]
    for key, con, gs in specs:
        U = sorted(set().union(*({g for g, a in zip(DE[d][0], auto[d]) if gs == "all" or a} for d in units)))
        ui = {g: i for i, g in enumerate(U)}
        for kind in ("q05", "top200"):
            lists, meas, signs = [], [], []
            for d in units:
                g, c = DE[d]
                keep = np.array([gs == "all" or a for a in auto[d]])
                idx = np.where(keep)[0]
                if kind == "q05":
                    sel = idx[c[f"q_{con}"][idx] < 0.05]
                else:
                    sel = idx[topn(c[f"t_{con}"][idx], c[f"p_{con}"][idx], 200)]
                lists.append({ui[g[i]] for i in sel}); meas.append({ui[g[i]] for i in idx})
                signs.append({ui[g[i]]: np.sign(c[f"t_{con}"][i]) for i in sel})
            obs, pu, tail, pp, tailP = recur(lists, meas, U, hitn)
            same = np.zeros(len(U), bool)
            for gi in np.where(obs >= 3)[0]:
                sg = {s[gi] for s in signs if gi in s}
                same[gi] = len(sg) == 1
            row = {"analysis": "full data", "contrast": key, "list": kind, "units": len(units), "unit_datasets": ",".join(units),
                   "median_list_size": int(np.median([len(L) for L in lists])), "list_sizes": ",".join(str(len(L)) for L in lists),
                   "genes_in_universe": len(U)}
            for k in (3, 5):
                row[f"obs_ge{k}"] = int((obs >= k).sum())
                row[f"exp_uniform_ge{k}"] = round(float(tail[:, k].sum()), 2) if k < tail.shape[1] else 0.0
                row[f"exp_prior_ge{k}"] = round(float(tailP[:, k].sum()), 2) if k < tailP.shape[1] else 0.0
            row["obs_ge3_same_sign"] = int(same.sum())
            row["genes_q10_uniform"] = int((recurrence.bh(pu) < 0.1).sum())
            row["genes_q10_prior"] = int((recurrence.bh(pp) < 0.1).sum())
            rrows.append(row)
            if kind == "top200":
                gene_tables[key] = (U, obs, pu, pp, units, con, gs)
    # matched recurrence at top 200
    if reps:
        for key, con, gs in [("sz_m_auto", "sz_m", "auto"), ("sex_m_auto", "sex_m", "auto"), ("sex_m_all", "sex_m", "all"),
                             ("age2_m_auto", "age2_m", "auto"), ("age_m_auto", "age_m", "auto")]:
            U = sorted(set().union(*({g for g, a in zip(DE[d][0], auto[d]) if gs == "all" or a} for d in units_m)))
            ui = {g: i for i, g in enumerate(U)}
            res = []
            for rep in reps:
                lists, meas = [], []
                for d in units_m:
                    x = np.fromfile(os.path.join(work, rep, f"{d}.bin"), dtype="<f4").astype(float)
                    G = len(DE[d][0]); k = MCON.index(con)
                    t, p = x[(2 * k) * G:(2 * k + 1) * G], x[(2 * k + 1) * G:(2 * k + 2) * G]
                    idx = np.where(np.array([gs == "all" or a for a in auto[d]]))[0]
                    sel = idx[topn(t[idx], p[idx], 200)]
                    lists.append({ui[DE[d][0][i]] for i in sel}); meas.append({ui[DE[d][0][i]] for i in idx})
                obs, pu, tail, pp, tailP = recur(lists, meas, U, hitn)
                res.append([(obs >= 3).sum(), (obs >= 5).sum(), tail[:, 3].sum(), tail[:, 5].sum(), tailP[:, 3].sum(), tailP[:, 5].sum()])
            res = np.array(res, float)
            rrows.append({"analysis": "matched", "contrast": key, "list": "top200", "units": len(units_m), "unit_datasets": ",".join(units_m),
                          "median_list_size": 200, "list_sizes": "200", "genes_in_universe": len(U),
                          "obs_ge3": round(res[:, 0].mean(), 2), "exp_uniform_ge3": round(res[:, 2].mean(), 2), "exp_prior_ge3": round(res[:, 4].mean(), 2),
                          "obs_ge5": round(res[:, 1].mean(), 2), "exp_uniform_ge5": round(res[:, 3].mean(), 2), "exp_prior_ge5": round(res[:, 5].mean(), 2),
                          "obs_ge3_rep_range": f"{int(res[:, 0].min())}-{int(res[:, 0].max())}", "obs_ge5_rep_range": f"{int(res[:, 1].min())}-{int(res[:, 1].max())}"})
    write(os.path.join(ROOT, "results", "recurrence.tsv"), rrows, [
        f"# scripts/pairs.py run {today}; units = one bulk dataset per collection cluster, the one with the largest SZ contrast n (matched: largest among datasets in the matched analysis)",
        "# list q05 = genes at BH q < 0.05 in the unit; top200 = the unit's 200 genes with the smallest p; universe = union of the units' measured genes (autosomal for _auto)",
        "# exp_uniform = expected genes in >= k lists when each unit draws its list uniformly from its measured genes (Poisson-binomial, recurrence.pb_matrix); "
        f"exp_prior = the same with draws in proportion to the DE-prior count + 1 ({recurrence.DE_PRIOR_SOURCES[recurrence.DE_PRIOR_SOURCE]['label']}; recurrence.prior_probs; median weight for genes not in the prior)",
        "# obs_ge3_same_sign = genes in >= 3 lists whose sign (t) agrees in every list; matched rows = mean over repeats"])

    # per-gene tables for top-200 recurrence
    for key, (U, obs, pu, pp, us, con, gs) in gene_tables.items():
        out = []
        qa, qb = recurrence.bh(pu), recurrence.bh(pp)
        for gi in np.argsort(-obs, kind="stable"):
            if obs[gi] < 2:
                break
            g = U[gi]
            ts = [DE[d][1][f"t_{con}"][gidx[d][g]] for d in us if g in gidx[d]]
            qs = [DE[d][1][f"q_{con}"][gidx[d][g]] for d in us if g in gidx[d]]
            out.append({"gene": g, "chrom": chrom.get(g, ""), "units_top200": int(obs[gi]), "units_measured": len(ts),
                        "units_q05": int(sum(q < 0.05 for q in qs)), "n_pos": int(sum(t > 0 for t in ts)), "n_neg": int(sum(t < 0 for t in ts)),
                        "stouffer_t": round(float(np.sum(ts) / math.sqrt(len(ts))), 3), "p_uniform": pu[gi], "q_uniform": qa[gi],
                        "p_prior": pp[gi], "q_prior": qb[gi], "de_prior_hitlists": hitn.get(g, "")})
        out.sort(key=lambda r: (-r["units_top200"], -abs(r["stouffer_t"])))
        write(os.path.join(ROOT, "results", f"recurrence_genes_{key}.tsv"), out, [
            f"# scripts/pairs.py run {today}; genes in the top 200 of >= 2 units ({', '.join(us)}); stouffer_t = sum of t over units measuring the gene / sqrt(units); "
            "t sign: sz = higher in SZ, sex = higher in males, age = higher with age"])

    # X/Y sanity
    san = []
    for d in names:
        g, c = DE[d]
        for mk, s in MARKERS.items():
            if mk in gidx[d] and "t_sex" in c:
                i = gidx[d][mk]
                san.append({"dataset": d, "gene": mk, "measured": "yes", "logfc_sex": c["logfc_sex"][i], "t_sex": c["t_sex"][i], "q_sex": c["q_sex"][i],
                            "expected_sign": "+" if s > 0 else "-", "pass": bool(c["q_sex"][i] < 0.05 and np.sign(c["t_sex"][i]) == s)})
            else:
                san.append({"dataset": d, "gene": mk, "measured": "no", "logfc_sex": "", "t_sex": "", "q_sex": "", "expected_sign": "+" if s > 0 else "-", "pass": ""})
    write(os.path.join(ROOT, "results", "sex_sanity.tsv"), san, [
        f"# scripts/pairs.py run {today}; full-data sex contrast (male vs female, all subjects); pass = q < 0.05 and expected sign (XIST lower in males, Y genes higher); "
        "'measured' = gene-level value present after the expression filter (gene level keeps the probe with the highest mean)"])

    # known sets (Mistry & Pavlidis 2010 top-10 per direction)
    known = []
    for setname, key in (("mistry2010_top10_age", "age_auto"), ("mistry2010_top10_sex", "sex_auto"), ("mistry2010_top10_sex", "sex_all")):
        ks = read_tsv(os.path.join(GL, "confounds", "genesets", f"{setname}.tsv"))
        rg = {r["gene"]: r for r in read_tsv(os.path.join(ROOT, "results", f"recurrence_genes_{key}.tsv"))}
        for r in ks:
            g = r["gene"]; x = rg.get(g)
            known.append({"set": setname, "contrast": key, "gene": g, "chrom": chrom.get(g, ""), "known_direction": r["direction"],
                          "units_top200": x["units_top200"] if x else "<2", "stouffer_t": x["stouffer_t"] if x else "",
                          "direction_agrees": (("pos" if float(x["stouffer_t"]) > 0 else "neg") == r["direction"]) if x else ""})
    write(os.path.join(ROOT, "results", "known_sets.tsv"), known, [
        f"# scripts/pairs.py run {today}; known sets from ../confounds/genesets (Mistry & Pavlidis 2010 Table 7, top 10 per direction); "
        "units_top200 from results/recurrence_genes_<contrast>.tsv ('<2' = in the top 200 of fewer than 2 units or not measured); sex direction pos = higher in males"])
    print("done")


def write(path, rows, header):
    if not rows:
        return
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w") as fh:
        fh.writelines(h + "\n" for h in header)
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4g}" if isinstance(v, (float, np.floating)) and not isinstance(v, bool) else v) for k, v in r.items()})


if __name__ == "__main__":
    if "--summary-only" in sys.argv[1:]:
        summary_only()
    else:
        main()
