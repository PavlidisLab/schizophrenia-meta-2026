#!/usr/bin/env python3
"""Per-panel tidy exports for manuscript Figures 8 (pH and agonal axes), 9 (single nucleus) and 10 (cell-type proportions).

Request: MAB_TO_RA_2026_10_08_FRESH_FIGURES_PANEL_DATA.md, Priority 2 (Fig 8, 9, 10).
Population: the 64 published brain studies (recurrence.py scenario 'published, brain', study tier); recurrent set = the 176 study-tier genes at
q < 0.1 (uniform null); recurrence.tsv sha256 d11761227abd.

Reads existing analysis outputs only (no model fits, no permutation reruns, no network). The only recomputation is light:
  8a  per-gene log2 change at 24 h from the Dachet tar (the gene-level table behind the scatter is not saved by dachet_live_tissue_2026-10-08)
  8c  Spearman of control-pH Z with each class specificity, to confirm the 2026-10-03 values against today's combined_controls_pH.tsv
  8f  per-subject table of GSE92538 pivoted from the saved GEO sample characteristics
  9a/9b  medians of saved per-pair / per-class tables
Writes to figure_data/fig08, fig09, fig10, each TSV with '#' provenance lines, plus _meta.json and STATUS.md per
folder. One function per panel; rerun with

    python3 analyses/figure_data_exports/export_fig08_fig09_fig10.py
"""
import csv
import gzip
import hashlib
import io
import json
import os
import re
import statistics
import sys
import tarfile
import zipfile
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
AN = os.path.dirname(HERE)            # analyses/
REPO = os.path.dirname(AN)
RUN = "2026-10-08"
SCRIPT = "analyses/figure_data_exports/export_fig08_fig09_fig10.py"
OUT = os.path.expanduser("figure_data")
POP = ("64 published brain studies (recurrence.py scenario 'published, brain', study tier); recurrent set = the 176 "
       "study-tier genes at q < 0.1 (uniform null); recurrence.tsv sha256 d11761227abd")
DACHET_TAR = "dachet/GSE153669_RAW.tar"      # sha256 cf6f11a4099f, as in dachet_live_tissue_2026-10-08

_sha_cache = {}


def sha(path):
    path = os.path.abspath(path)
    if path not in _sha_cache:
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        _sha_cache[path] = h.hexdigest()
    return _sha_cache[path]


def rel(path):
    p = os.path.abspath(path)
    return os.path.relpath(p, REPO) if p.startswith(REPO) else p.replace(os.path.expanduser("~"), "~")


def rd(path):
    """Rows of a TSV ('.gz' allowed), '#' lines skipped, values kept as the strings in the file."""
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", encoding="utf-8") as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def pdr(path, **kw):
    return pd.read_csv(path, sep="\t", comment="#", **kw)


def A(*p):
    return os.path.join(AN, *p)


def fmt(v):
    if isinstance(v, (float, np.floating)):
        return "" if np.isnan(v) else repr(round(float(v), 6))
    if isinstance(v, (np.integer,)):
        return str(int(v))
    return v


class Figure:
    def __init__(self, fig):
        self.fig = fig
        self.dir = os.path.join(OUT, fig)
        os.makedirs(self.dir, exist_ok=True)
        self.meta = {"script": SCRIPT, "run_id": RUN, "figure": fig, "population": POP, "inputs": {}, "panels": {}}
        self.status = []

    def note(self, panel, text):
        self.status.append(f"- **{panel}**: {text}")

    def write(self, name, panel, fields, rows, inputs, extra_comments=()):
        for p in inputs:
            self.meta["inputs"][rel(p)] = sha(p)
        lines = [f"# figure: {self.fig} panel {panel}", f"# script: {SCRIPT}, run {RUN}", f"# population: {POP}"]
        lines += [f"# input: {rel(p)} sha256 {sha(p)[:12]}" for p in inputs]
        lines += [f"# {c}" for c in extra_comments]
        path = os.path.join(self.dir, name)
        with open(path, "w", encoding="utf-8", newline="") as fh:
            for l in lines:
                fh.write(l + "\n")
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(fields)
            n = 0
            for r in rows:
                w.writerow([("" if r.get(f) is None else fmt(r.get(f))) for f in fields])
                n += 1
        self.meta["panels"][name] = {"panel": panel, "rows": n}
        print(f"{self.fig}/{name}: {n} rows")
        return n

    def df(self, name, panel, d, inputs, extra_comments=()):
        return self.write(name, panel, list(d.columns), d.to_dict("records"), inputs, extra_comments)

    def finish(self):
        with open(os.path.join(self.dir, "_meta.json"), "w") as fh:
            json.dump(self.meta, fh, indent=1)
        with open(os.path.join(self.dir, "STATUS.md"), "w") as fh:
            fh.write(f"# {self.fig} status ({RUN})\n\n" + "\n".join(self.status) + "\n")


# ------------------------------------------------------------------ shared inputs
P_CTLPH = A("confounds", "results", "combined_controls_pH.tsv")
P_STUDY_CLAIMS = A("study_claims_2026-10-08", "study_claims.tsv")
P_DACHET_DIR = A("dachet_live_tissue_2026-10-08", "results")
P_AX_REMOTE = A("axes_pritzker_agonal_2026-10-09")
P_CLAIMS_AX = A("claims_vs_ph_axis_2026-10-08")
P_HAG = A("hagihara_replication_2026-10-08", "results")
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def ctl_ph(factor):
    """gene -> (Z, q) from combined_controls_pH.tsv (Z > 0: higher with higher pH / PMI in controls)."""
    out = {}
    for r in rd(P_CTLPH):
        if r["factor"] == factor and r["Z"]:
            out[r["gene"]] = (float(r["Z"]), float(r["q"]) if r["q"] else None)
    return out


def xl_pairs(data):
    """gene -> value (columns A, B) of the first sheet of an xlsx given as bytes (same reading as dachet_live_tissue_2026-10-08/build_report.py)."""
    z = zipfile.ZipFile(io.BytesIO(data))
    ss = [''.join(t.text or '' for t in si.iter(NS + 't')) for si in ET.parse(z.open('xl/sharedStrings.xml')).getroot()] if 'xl/sharedStrings.xml' in z.namelist() else []
    rels = {r.get('Id'): r.get('Target') for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
    s = ET.fromstring(z.read('xl/workbook.xml')).find(NS + 'sheets')[0]
    tgt = rels[s.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')].lstrip('/')
    tgt = tgt if tgt.startswith('xl/') else 'xl/' + tgt
    o = {}
    for r in ET.parse(z.open(tgt)).getroot().iter(NS + 'row'):
        row = {}
        for c in r.findall(NS + 'c'):
            v = c.find(NS + 'v')
            if v is not None:
                row[re.match(r'[A-Z]+', c.get('r')).group()] = ss[int(v.text)] if c.get('t') == 's' else v.text
        if 'A' in row and 'B' in row:
            try:
                o[row['A']] = float(row['B'])
            except ValueError:
                pass
    return o


# ================================================================== Figure 8
def fig8_a(F):
    """Dachet fresh cortex: gene-level pH association vs 24 h change; correlations by time point; direction of recurrent genes."""
    # gene-level table: recomputed from the tar (the gene table is not saved by build_report.py)
    assert sha(DACHET_TAR)[:12] == "cf6f11a4099f", "unexpected Dachet tar"
    E = {}
    with tarfile.open(DACHET_TAR) as tf:
        for m in tf.getmembers():
            t = re.search(r'_(\d+)H\.xlsx$', m.name)
            if t:
                E[int(t.group(1))] = xl_pairs(tf.extractfile(m).read())
    hours = [0, 1, 2, 4, 8, 12, 24]
    genes = sorted(set.intersection(*[set(E[h]) for h in hours]))
    ph, pmi = ctl_ph("pH"), ctl_ph("PMI")
    rows = []
    for g in genes:
        row = {"gene": g, "ph_assoc": ph[g][0] if g in ph else None, "q_ph": ph[g][1] if g in ph else None,
               "pmi_assoc": pmi[g][0] if g in pmi else None}
        for h in hours[1:]:
            row[f"change_{h}h"] = float(np.log2(E[h][g] + 1) - np.log2(E[0][g] + 1))
        rows.append(row)
    d = pd.DataFrame(rows).rename(columns={"change_24h": "change_24h"})
    k = d.dropna(subset=["ph_assoc"])
    rho, _ = spearmanr(k.ph_assoc, k.change_24h)
    F.write("fig08_a_dachet_gene_scatter.tsv", "a (left)",
            ["gene", "ph_assoc", "q_ph", "pmi_assoc", "change_1h", "change_2h", "change_4h", "change_8h", "change_12h", "change_24h"],
            d.to_dict("records"), [P_CTLPH],
            [f"Dachet et al. 2021 GSE153669_RAW.tar sha256 cf6f11a4099f (one patient, 0 h to 24 h at room temperature; not an analysis output, recomputed here as in dachet_live_tissue_2026-10-08/build_report.py)",
             "ph_assoc = control-only Stouffer Z of the gene with brain pH (Z > 0: higher at higher pH; five collections); pmi_assoc likewise for PMI; blank = gene not in the control table",
             "changeNh = log2(x_N + 1) - log2(x_0 + 1) of the unlogged Dachet values; 15,655 genes are in all seven time points",
             f"Spearman rho(ph_assoc, change_24h) = {rho:+.3f} over {len(k)} genes (legend: -0.50, 13,027)"])
    F.note("a (left)", f"{len(d)} genes in all seven time points; {len(k)} with a control pH Z; Spearman rho of ph_assoc with change_24h = {rho:+.3f} (legend: -0.50 over 13,027 genes: match). "
                       "The gene-level table was recomputed here from the tar in dachet (sha256 cf6f11a4099f); dachet_live_tissue_2026-10-08 saves only the correlations.")
    p = os.path.join(P_DACHET_DIR, "alignment_by_time.tsv")
    r = rd(p)
    F.write("fig08_a_alignment_by_time.tsv", "a (right)", ["hours", "axis", "genes", "rho"], r, [p],
            ["Spearman rho of the change from 0 h with the control pH Z, the control PMI Z and Gandal 2018 SCZ / BD / ASD t statistics; Z > 0 = higher with higher pH / PMI in controls"])
    p2 = os.path.join(P_DACHET_DIR, "recurrent_direction.tsv")
    F.write("fig08_a_recurrent_direction.tsv", "a (text)", ["set", "genes", "same_direction", "share", "expected", "p_binomial"], rd(p2), [p2],
            ["share of recurrent genes whose SZ direction equals their direction in the dying tissue at 24 h (not drawn in the Fig 8a legend; the manuscript text quotes 62% of 172, p 0.0028)"])
    rr = {x["set"]: x for x in rd(p2)}["study tier (q < 0.1)"]
    F.note("a (text)", f"study tier: {rr['same_direction']} of {rr['genes']} genes same direction ({rr['share']}, expected {rr['expected']}, p {rr['p_binomial']}). The manuscript says 172 study-tier genes; the file has 172 genes with a direction and a Dachet measurement, out of the 176 study-tier genes.")


def fig8_b(F):
    """Axis projection: per-dataset alignment of the SZ t with each axis, pair agreement before and after projecting an axis out."""
    # per-dataset alignment, remote run 2026-10-08/09, four axes
    pf = os.path.join(P_AX_REMOTE, "results_remote", "alignment_by_dataset.tsv")
    f = pdr(pf)
    pc = os.path.join(P_CLAIMS_AX, "alignment_two_axes.tsv")
    c = pdr(pc).rename(columns={"rho_pH_L0": "claims_rho_pH_L0", "rho_agonal_L0": "claims_rho_gtex_L0", "rho_pH_L4": "claims_rho_pH_L4", "rho_agonal_L4": "claims_rho_gtex_L4"})
    pbd = os.path.join(P_CLAIMS_AX, "alignment_by_dataset.tsv")
    b = pdr(pbd).rename(columns={"rho_axis_L0": "claims_by_dataset_rho_pH_L0", "rho_axis_L4": "claims_by_dataset_rho_pH_L4"})[["dataset", "axis_collections", "claims_by_dataset_rho_pH_L0", "claims_by_dataset_rho_pH_L4"]]
    m = f.rename(columns={"rho_pH_L0": "rho_ph_axis_L0", "rho_pH_L4": "rho_ph_axis_L4", "rho_GTEx_L0": "rho_gtex_agonal_L0", "rho_GTEx_L4": "rho_gtex_agonal_L4",
                          "rho_Pritz_L0": "rho_pritzker_agonal_L0", "rho_Pritz_L4": "rho_pritzker_agonal_L4", "rho_Dai_L0": "rho_dai_L0", "rho_Dai_L4": "rho_dai_L4"})
    m = m.merge(b, on="dataset", how="left").merge(c, on="dataset", how="left")
    cols = ["dataset", "axis_collections", "in_ALL", "in_NOPRITZ", "in_NODAI", "rho_ph_axis_L0", "rho_ph_axis_L4", "rho_gtex_agonal_L0", "rho_gtex_agonal_L4",
            "rho_pritzker_agonal_L0", "rho_pritzker_agonal_L4", "rho_dai_L0", "rho_dai_L4",
            "claims_by_dataset_rho_pH_L0", "claims_by_dataset_rho_pH_L4", "claims_rho_pH_L0", "claims_rho_pH_L4", "claims_rho_gtex_L0", "claims_rho_gtex_L4"]
    F.write("fig08_b_dataset_axis_alignment.tsv", "b (left)", cols, m.to_dict("records"), [pf, pbd, pc],
            ["rho_*: Spearman rho of the dataset's SZ t with the axis over shared genes; L0 = unadjusted DE, L4 = all measured covariates in the model. 28 datasets (GSE92538_GPL570 is the 29th row of nothing; the remote file has 28)",
             "ph_axis = control-only pH Z, Stouffer over the collections sharing no donors with the target (leave-collection-out); expected sign negative. gtex_agonal = GTEx v8 slow vs fast natural death t (Hardy scale), mean of BA9, BA24, Cortex; "
             "pritzker_agonal = per-gene t of the agonal factor in GSE92538 GPL17027 (104 subjects); dai = Dai 2021 fever/infection +1/-1 symbol sets; expected sign positive for all three",
             "Columns without 'claims_': axes_pritzker_agonal_2026-10-09/results_remote (remote-regenerated inputs, run 2026-10-08/09, README says identical to the laptop run to four decimals). "
             "Columns 'claims_by_dataset_*': claims_vs_ph_axis_2026-10-08/alignment_by_dataset.tsv (pH only); 'claims_*': alignment_two_axes.tsv of the same folder (pH and GTEx). The two claims_ pH columns are NOT the same (see STATUS.md)"])
    # reconcile medians
    med = lambda s: float(np.nanmedian(pd.to_numeric(s)))                       # noqa: E731
    F.note("b (left)", "median rho by source (L0 / L4): remote pH %+.3f / %+.3f, remote GTEx %+.3f / %+.3f, Pritzker %+.3f / %+.3f; claims_vs_ph_axis alignment_by_dataset pH %+.3f / %+.3f; "
           "alignment_two_axes pH %+.3f / %+.3f, agonal %+.3f / %+.3f. The claims_vs_ph_axis README quotes pH -0.20 / -0.10 and agonal +0.24 / +0.19."
           % (med(m.rho_ph_axis_L0), med(m.rho_ph_axis_L4), med(m.rho_gtex_agonal_L0), med(m.rho_gtex_agonal_L4), med(m.rho_pritzker_agonal_L0), med(m.rho_pritzker_agonal_L4),
              med(m.claims_by_dataset_rho_pH_L0), med(m.claims_by_dataset_rho_pH_L4), med(m.claims_rho_pH_L0), med(m.claims_rho_pH_L4), med(m.claims_rho_gtex_L0), med(m.claims_rho_gtex_L4)))
    # pairs: claims_vs_ph_axis (GTEx axis) medians over 192 both-bulk pairs
    pp = os.path.join(P_CLAIMS_AX, "pairs_two_axes.tsv")
    d = pdr(pp)
    bb = d[d.both_bulk.astype(str) == "True"]
    rows = []
    for lev in ("L0", "L4"):
        for lab, col in (("none (raw)", f"raw_{lev}"), ("pH axis removed", f"partial_pH_{lev}"), ("GTEx agonal axis removed", f"partial_agonal_{lev}"), ("pH and GTEx agonal removed", f"partial_both_{lev}")):
            rows.append({"axis_version": "GTEx Hardy agonal axis + control pH axis (claims_vs_ph_axis_2026-10-08, laptop)", "pair_class": "both-bulk donor-separate", "level": lev,
                         "axis_removed": lab, "n_pairs": len(bb), "median_rho": float(bb[col].median())})
            rows.append({"axis_version": "GTEx Hardy agonal axis + control pH axis (claims_vs_ph_axis_2026-10-08, laptop)", "pair_class": "all donor-separate", "level": lev,
                         "axis_removed": lab, "n_pairs": len(d), "median_rho": float(d[col].median())})
    F.write("fig08_b_pair_agreement_gtex_axis.tsv", "b (right)", ["axis_version", "pair_class", "level", "axis_removed", "n_pairs", "median_rho"], rows, [pp],
            ["median over donor-separate pairs of the Spearman rho of the two datasets' SZ t, raw and after the rank-residual projection of the axis out of both profiles; L0 unadjusted DE, L4 all measured covariates",
             "computed here as medians of pairs_two_axes.tsv (310 lines = 308 pairs; 192 both-bulk)"])
    g = {(r["axis_removed"], r["level"], r["pair_class"]): r["median_rho"] for r in rows}
    F.note("b (right)", "GTEx-axis medians over 192 both-bulk pairs recomputed from pairs_two_axes.tsv: L0 raw %.3f -> pH removed %.3f, GTEx removed %.3f, both %.3f; L4 raw %.3f -> pH %.3f, GTEx %.3f, both %.3f "
           "(legend: 0.118 to 0.072 [0.076 to 0.048 with covariates in the model], 0.104 with the pH axis removed)."
           % tuple(g[(a, l, "both-bulk donor-separate")] for l, a in (("L0", "none (raw)"), ("L0", "pH axis removed"), ("L0", "GTEx agonal axis removed"), ("L0", "pH and GTEx agonal removed"),
                                                                       ("L4", "none (raw)"), ("L4", "pH axis removed"), ("L4", "GTEx agonal axis removed"), ("L4", "pH and GTEx agonal removed"))))
    # remote same-gene single-axis table
    ps = os.path.join(P_AX_REMOTE, "single_axis_same_genes_remote.tsv")
    s = pdr(ps)
    s.insert(0, "axis_version", s.axis_removed.map({"pH only": "control pH axis", "GTEx": "GTEx Hardy agonal axis", "Pritz": "Pritzker agonal axis (GSE92538)", "Dai": "Dai 2021 fever/infection symbols"}))
    s["set_description"] = s.set.map({"ALL": "all 28 datasets (axis source's own donors included for Pritzker)", "NOPRITZ": "23 datasets: minus GSE92538_*, GSE80655_*",
                                      "NOPRITZ2": "22 datasets: also minus GSE224683"})
    F.df("fig08_b_pair_agreement_single_axis_same_genes.tsv", "b (right)", s, [ps],
         ["axes_pritzker_agonal_2026-10-09/single_axis_same_genes_remote.tsv (remote inputs, run 2026-10-08/09): both-bulk donor-separate pairs, each axis on its own gene set; share_removed_of_median negative = agreement rises after projection",
          "ALL/GTEx L0 0.1177 -> 0.0717 and L4 0.0755 -> 0.0477 reproduce the laptop numbers of claims_vs_ph_axis"])
    # pairs, remote, all combinations
    pf2 = os.path.join(P_AX_REMOTE, "results_remote", "pairs_partial.tsv")
    pairs = pdr(pf2)
    F.df("fig08_b_pairs_remote.tsv", "b (right)", pairs, [pf2],
         ["one row per donor-separate pair (308; both_bulk TRUE for 192): median-ready values for every axis combination at L0 and L4 (remote inputs); genes_Lk = genes in common"])
    # parse summary.txt: axis correlations and per-target-set blocks
    st = os.path.join(P_AX_REMOTE, "results_remote", "summary.txt")
    corr_rows, set_rows, pair_rows = [], [], []
    tset, nds = None, None
    for ln in open(st):
        ln = ln.rstrip("\n")
        m1 = re.search(r"\s*(Pritzker agonal vs GTEx Hardy|Pritzker agonal vs pH axis \(all five collections\)|GTEx Hardy vs pH axis|Dai .*?) ([+-][\d.]+) \(n genes (\d+)\)", ln)
        if m1:
            corr_rows.append({"pair_of_axes": m1.group(1), "spearman": m1.group(2), "n_genes": m1.group(3)})
            continue
        m2 = re.match(r"== target set (\w+): (\d+) datasets", ln)
        if m2:
            tset, nds = m2.group(1), m2.group(2)
            continue
        m3 = re.match(r"\s*rho\(SZ t, (\w+) axis\) (L\d): median ([+-][\d.]+), range ([+-][\d.]+) to ([+-][\d.]+), expected sign in (\d+) of (\d+)", ln)
        if m3 and tset:
            set_rows.append({"target_set": tset, "n_datasets": nds, "axis": m3.group(1), "level": m3.group(2), "median_rho": m3.group(3), "range_lo": m3.group(4), "range_hi": m3.group(5),
                             "n_expected_sign": m3.group(6), "n_of": m3.group(7)})
            continue
        m4 = re.match(r"\s*(both-bulk pairs|all donor-separate pairs) \((\d+)\) (L\d): (.*)", ln)
        if m4 and tset:
            for kv in re.finditer(r"([\w+]+) ([+-]?[\d.]+)", m4.group(4)):
                pair_rows.append({"target_set": tset, "pair_class": m4.group(1), "n_pairs": m4.group(2), "level": m4.group(3), "axes_removed": kv.group(1), "median_rho": kv.group(2)})
    F.write("fig08_b_axis_correlations.tsv", "b (context)", ["pair_of_axes", "spearman", "n_genes"], corr_rows, [st],
            ["parsed from summary.txt of axes_pritzker_agonal_2026-10-09/results_remote; Spearman correlation between axes over genes"])
    F.write("fig08_b_alignment_medians_by_target_set.tsv", "b (context)", ["target_set", "n_datasets", "axis", "level", "median_rho", "range_lo", "range_hi", "n_expected_sign", "n_of"], set_rows, [st],
            ["parsed from summary.txt: median rho of the SZ t with each axis over the datasets of the target set; 'expected sign' = negative for pH, positive for the others"])
    F.write("fig08_b_pair_medians_by_target_set.tsv", "b (context)", ["target_set", "pair_class", "n_pairs", "level", "axes_removed", "median_rho"], pair_rows, [st],
            ["parsed from summary.txt: on the gene set common to all four axes (differs from single_axis_same_genes by gene restriction), median pair agreement with the listed axes projected out ('none' = raw)"])
    F.note("b", "TWO AGONAL AXES, labelled in every file: 'GTEx Hardy' (gtex_ischemia_2026-10-08: GTEx v8 slow vs fast natural death, donors separate from every array dataset) and 'Pritzker agonal' "
                "(GSE92538 agonal-factor t, 104 subjects; shares donors with GSE80655 and part of GSE224683, so only the NOPRITZ / NOPRITZ2 target sets are independent of it). "
                "Numbers are from remote, run 2026-10-08/09 (axes_pritzker_agonal_2026-10-09/README.md); the README says the laptop copy is identical to four decimals in every summary line.")
    x = {(r["set"], r["axis_removed"], r["level"]): r for r in s.to_dict("records")}
    F.note("b", "remote same-gene medians (both-bulk): ALL GTEx L0 %.4f -> %.4f; NOPRITZ (109 pairs) GTEx L0 %.4f -> %.4f, Pritz L0 %.4f -> %.4f, pH L0 %.4f -> %.4f; at L4 removing Pritz / pH raises the agreement in NOPRITZ "
                "(%.4f -> %.4f / %.4f)." % (x[("ALL", "GTEx", "L0")]["median_raw_same_genes"], x[("ALL", "GTEx", "L0")]["median_after_removal"],
                                            x[("NOPRITZ", "GTEx", "L0")]["median_raw_same_genes"], x[("NOPRITZ", "GTEx", "L0")]["median_after_removal"],
                                            x[("NOPRITZ", "Pritz", "L0")]["median_raw_same_genes"], x[("NOPRITZ", "Pritz", "L0")]["median_after_removal"],
                                            x[("NOPRITZ", "pH only", "L0")]["median_raw_same_genes"], x[("NOPRITZ", "pH only", "L0")]["median_after_removal"],
                                            x[("NOPRITZ", "Pritz", "L4")]["median_raw_same_genes"], x[("NOPRITZ", "Pritz", "L4")]["median_after_removal"], x[("NOPRITZ", "pH only", "L4")]["median_after_removal"]))


def fig8_c(F):
    """Signed control-pH Z against neuronal and astrocyte specificity; class structure; hypoxia reference-set tests."""
    ps = A("interneuron_ph_2026-10-01", "results", "specificity_MSSM.tsv.gz")
    z = {r["gene"]: r for r in rd(P_CTLPH) if r["factor"] == "pH" and r["Z"]}
    spec = {r["gene"]: r for r in rd(ps)}
    classes = ["EN", "coarse_IN", "Astro", "Oligo", "OPC", "Micro", "Vascular"]
    rows = []
    for g, r in z.items():
        if g in spec:
            row = {"gene": g, "z_ph": r["Z"], "q_ph": r["q"], "n_datasets": r["n_datasets"]}
            for c in classes:
                row[f"spec_{c}"] = spec[g][f"spec_{c}"]
            rows.append(row)
    F.write("fig08_c_gene_ph_specificity.tsv", "c (left)", ["gene", "z_ph", "q_ph", "n_datasets"] + [f"spec_{c}" for c in classes], rows, [P_CTLPH, ps],
            ["z_ph = control-only Stouffer Z with brain pH over the five collections (Z > 0: RNA higher at higher pH, so neuronal genes at low pH are lower); spec_* = EWCE specificity (mean CPM in class / sum over classes) from PsychAD MSSM control donors "
             "(interneuron_ph_2026-10-01/scripts/specificity.py); EN excitatory, coarse_IN inhibitory",
             "the legend's 'neuron_spec' = spec_EN (rho +0.53), 'astro_spec' = spec_Astro (rho -0.35); the running median is not saved (compute from this table)"])
    # recompute class rho against the 2026-10-03 table
    old = {r["statistic"]: r for r in rd(A("confounds", "ph_pHindex_2026-10-03", "results", "results.tsv")) if r["test"] == "T1_cells"}
    zz = np.array([float(r["z_ph"]) for r in rows])
    crow = []
    for c in classes:
        sp = np.array([float(r[f"spec_{c}"]) for r in rows])
        rho = spearmanr(sp, zz)[0]
        o = old.get(f"Z_pH ~ spec_{'coarse_IN' if c == 'coarse_IN' else c} (signed)")
        crow.append({"class": c, "n_genes": len(rows), "spearman_z_ph_vs_spec": rho, "value_in_results_2026-10-03": o["value"] if o else "", "residual_after_expression_fit_2026-10-03": (re.search(r"residual ([+-][\d.]+)", o["detail"]).group(1) if o else "")})
    p_old = A("confounds", "ph_pHindex_2026-10-03", "results", "results.tsv")
    F.write("fig08_c_class_direction.tsv", "c (right)", ["class", "n_genes", "spearman_z_ph_vs_spec", "value_in_results_2026-10-03", "residual_after_expression_fit_2026-10-03"], crow, [P_CTLPH, ps, p_old],
            ["direction structure by cell class: Spearman rho of signed control-pH Z with the class's specificity over the genes tested (pH x specificity); recomputed here from the table fig08_c_gene_ph_specificity.tsv and set beside the 2026-10-03 values of ph_pHindex_2026-10-03/results/results.tsv"])
    F.note("c (left/right)", "Rho recomputed on today's combined_controls_pH.tsv: " + ", ".join(f"{r['class']} {r['spearman_z_ph_vs_spec']:+.3f} (2026-10-03: {r['value_in_results_2026-10-03']})" for r in crow)
           + f"; {len(rows)} genes (2026-10-03: 13398). ph_pHindex_2026-10-03 was built on combined_controls_pH.tsv sha 8e029621de28, today's is 3b969e3ac151; the pH Z column is identical for all 15,187 genes (max abs difference 0), so the 2026-10-03 results stand and nothing was rerun.")
    # reference-set tests
    srcs = [("pr_auc_lower_full.tsv", "4-set family (what_are_the_pH_genes.py; end = RNA higher at low pH)"), ("pr_auc_external.tsv", "6-set family incl. MSigDB Hallmark hypoxia"),
            ("pr_auc_postmortem.tsv", "postmortem pH signatures (replication cross-cohort)")]
    out, ins = [], []
    for fn, fam in srcs:
        p = A("confounds", "ph_pHindex_2026-10-03", "results", fn)
        ins.append(p)
        for r in rd(p):
            out.append({"source_file": fn, "test_family": fam, "set_id": r["set_id"], "set_name": r["set_name"], "set_size": r["set_size"], "n_in_ranking": r["set_size_after_filter"],
                        "pr_auc": r["pr_auc"], "p": r["pvalue"], "q": r["qvalue"], "mf_pr_auc": r["mf_pr_auc"], "mf_corrected_p": r["mf_corrected_pvalue"], "mf_corrected_q": r["mf_corrected_qvalue"]})
    F.write("fig08_c_reference_set_tests.tsv", "c (right)", ["source_file", "test_family", "set_id", "set_name", "set_size", "n_in_ranking", "pr_auc", "p", "q", "mf_pr_auc", "mf_corrected_p", "mf_corrected_q"], out, ins,
            ["PR-AUC enrichment (Ballouz 2017; ranking of 15,187 genes by control-pH Z, interesting end = RNA higher at low pH; permutations; BH q within the family of the source file; mf_* = multifunctionality-corrected). Run 2026-10-03"])
    q = {(r["source_file"], r["set_id"]): r["q"] for r in out}
    F.note("c (right)", "GO:0001666 q = %s in pr_auc_lower_full.tsv (4-set family) but %s in pr_auc_external.tsv (6-set family); Hallmark hypoxia q = %s (external only). The Fig 8 legend and report text give GO q = 0.0016 and Hallmark q = 0.0007, "
                        "the report's figure caption says 'q 0.0007 each': the two q values for GO come from different BH families. Both sets are flagged multifunctionality-driven (MF-corrected p %s and %s)."
        % (q[("pr_auc_lower_full.tsv", "GO_0001666_response_to_hypoxia")], q[("pr_auc_external.tsv", "GO_0001666_response_to_hypoxia")], q[("pr_auc_external.tsv", "msigdb_hallmark_hypoxia_up")],
           next(r["mf_corrected_p"] for r in out if r["source_file"] == "pr_auc_external.tsv" and r["set_id"] == "GO_0001666_response_to_hypoxia"),
           next(r["mf_corrected_p"] for r in out if r["source_file"] == "pr_auc_external.tsv" and r["set_id"] == "msigdb_hallmark_hypoxia_up")))


def fig8_d(F):
    """Interneuron marker scores by model, pH-matched random genes; single-gene papers by pH handling."""
    d = A("interneuron_ph_2026-10-01", "results")
    p = os.path.join(d, "subject_level_summary.tsv")
    rows = rd(p)
    for r in rows:
        r["is_ph_matched_random"] = "yes" if r["score"].endswith("__random") else "no"
    F.write("fig08_d_interneuron_score_by_model.tsv", "d (left)",
            ["set", "score", "is_ph_matched_random", "clusters", "d_base", "d_pH", "d_pH_hyp", "rho_ctl", "rho_sz", "slope_ctl", "slope_sz", "share_pH_median", "clusters_d_base_negative"], rows, [p],
            ["d_* = case-minus-control difference of the mean row-z score of the marker group (SD units), median over collection clusters: d_base adjusted for age and sex, d_pH + pH, d_pH_hyp + pH and S_hypoxia score; "
             "'__random' rows = 200 random gene sets matched to the group on pH Z, age Z, hypoxia loading and expression. Run 2026-10-08 (subject_level.py); 18 dataset-regions that record pH, 1,034 subjects",
             "the legend's rows: CLEAN cortex (Stanley, Pritzker DLPFC; 2 clusters) IN_markers and IN_markers__random; the other sets are the five-cluster and leave-same-data-out variants"])
    byset = {}
    for r in rows:
        byset.setdefault(r["set"], {})[r["score"]] = r
    cl = byset["CLEAN cortex (Stanley, Pritzker DLPFC)"]
    F.note("d (left)", "CLEAN cortex IN_markers d_base/d_pH/d_pH_hyp = %s / %s / %s (legend -0.53, -0.39, -0.32); IN_markers__random %s / %s / %s (legend -0.38, -0.20, -0.09): match." % (
        cl["IN_markers"]["d_base"], cl["IN_markers"]["d_pH"], cl["IN_markers"]["d_pH_hyp"], cl["IN_markers__random"]["d_base"], cl["IN_markers__random"]["d_pH"], cl["IN_markers__random"]["d_pH_hyp"]))
    p2 = os.path.join(d, "subject_level_by_dataset.tsv")
    F.write("fig08_d_interneuron_score_by_dataset.tsv", "d (left, detail)",
            ["dataset", "cluster", "score", "n_ctl", "n_sz", "rho_ctl", "rho_sz", "d_base", "d_pH", "d_pH_hyp", "slope_ctl", "slope_sz"], rd(p2), [p2],
            ["per dataset-region and score group, same quantities as fig08_d_interneuron_score_by_model.tsv"])
    p3 = os.path.join(d, "taskB_rates.tsv")
    F.write("fig08_d_single_gene_papers_by_ph_handling.tsv", "d (right)",
            ["gene_group", "tier", "n_calls", "n_papers", "n_down", "n_up", "n_ns", "n_mixed", "rate_down_calls", "n_bank_clusters", "rate_down_bank_clustered", "n_calls_pitt", "rate_down_pitt", "rate_down_non_pitt"], rd(p3), [p3],
            ["H4 single-gene brain calls (paper x gene x region x molecule) split by how the paper handled pH: A in the gene model or pairs matched on pH, B groups matched, C reported not used, D open methods read and pH not reported, E no open methods text; "
             "rate_down_* = share of usable calls that report a decrease. File header says run 2026-10-01; the folder README says the H4 tally was rebuilt on the 2026-10-03 h4_calls.tsv"])
    F.note("d (right)", "taskB_rates.tsv has 52 rows (gene groups x tiers). Its provenance lines record analyses/single_gene_study_agreement/h4_calls.tsv sha256 e7ebd1a44b2c; I did not check it against the current h4_calls.tsv.")


def fig8_e(F):
    """GTEx: Hardy scale (slow vs fast natural death) and ischemic time, per brain region, against the control pH association and the Dachet change."""
    p = A("gtex_ischemia_2026-10-08", "results_by_region.tsv")
    d = pdr(p)
    hd = d[d.statistic == "slow vs fast natural death"].set_index("region")
    isch = d[d.statistic == "ischemic time"].set_index("region")
    rows = []
    for reg in hd.index:
        rows.append({"region": reg, "samples": hd.loc[reg, "samples"], "genes": hd.loc[reg, "genes"],
                     "rho_hardy_vs_ph": hd.loc[reg, "rho_pH Z"], "n_genes_ph": hd.loc[reg, "n_pH Z"], "rho_hardy_vs_pmi": hd.loc[reg, "rho_PMI Z"],
                     "rho_hardy_vs_dachet": hd.loc[reg, "rho_Dachet 24 h change"], "n_genes_dachet": hd.loc[reg, "n_Dachet 24 h change"],
                     "rho_ischemic_vs_ph": isch.loc[reg, "rho_pH Z"], "rho_ischemic_vs_pmi": isch.loc[reg, "rho_PMI Z"], "rho_ischemic_vs_dachet": isch.loc[reg, "rho_Dachet 24 h change"]})
    F.write("fig08_e_gtex_by_region.tsv", "e",
            ["region", "samples", "genes", "rho_hardy_vs_ph", "n_genes_ph", "rho_hardy_vs_pmi", "rho_hardy_vs_dachet", "n_genes_dachet", "rho_ischemic_vs_ph", "rho_ischemic_vs_pmi", "rho_ischemic_vs_dachet"], rows, [p],
            ["GTEx v8 brain: per gene (median TPM >= 1 in the region) OLS of log2(TPM + 1) on ischemic time, age, sex and Hardy scale; 'hardy' = slow death (4) vs fast natural death (2) t statistic; Spearman rho with the control pH Z, PMI Z and the Dachet 24 h change. "
             "13 regions; run 2026-10-08 (gtex_ischemic_time.py)"])
    v = [float(r["rho_hardy_vs_ph"]) for r in rows]
    w = [float(r["rho_hardy_vs_dachet"]) for r in rows]
    F.note("e", f"{len(rows)} regions; Hardy vs pH negative in {sum(x < 0 for x in v)} of {len(v)}, median {statistics.median(v):+.3f} (legend: negative in all, median -0.28: match); Hardy vs Dachet positive in {sum(x > 0 for x in w)} of {len(w)}, median {statistics.median(w):+.3f} (manuscript text: 12 of 13, +0.30: match).")


def fig8_f(F):
    """GSE92538 (Pritzker): agonal factor, pH and diagnosis; pH effect by dataset with and without agonal state."""
    pg = os.path.join(P_HAG, "geo_sample_characteristics.tsv")
    L = pd.read_csv(pg, sep="\t", dtype=str, comment="#")
    L = L[L.series == "GSE92538"]
    W = L.pivot_table(index=["file", "gsm"], columns="key", values="value", aggfunc="first").reset_index()
    num = lambda s: pd.to_numeric(s.replace({"NA": np.nan, "N/A": np.nan}), errors="coerce")      # noqa: E731
    D = pd.DataFrame({"platform": W.file.str.extract(r"(GPL\d+)")[0], "gsm": W.gsm, "subject": W["subject id"], "diagnosis": W["diagnosis"], "age": num(W["age"]), "sex": W["gender"],
                      "pmi_h": num(W["post-mortem interval"]), "ph_cerebellum": num(W["tissue ph (cerebellum)"]), "agonal_factor": num(W["agonal factor"]), "suicide": num(W["suicide (1=yes)"])})
    D = D.sort_values("platform", ascending=False).drop_duplicates("subject")          # GPL17027 preferred, as in gse92538_ph_agonal.py
    D = D.sort_values(["diagnosis", "subject"])
    D["complete_case_for_models"] = D[["ph_cerebellum", "age", "pmi_h", "agonal_factor"]].notna().all(axis=1).map({True: "yes", False: "no"})
    F.df("fig08_f_gse92538_subjects.tsv", "f", D[["subject", "platform", "gsm", "diagnosis", "age", "sex", "pmi_h", "ph_cerebellum", "agonal_factor", "suicide", "complete_case_for_models"]], [pg],
         ["one row per subject ('subject id'; the series carry technical replicates; GPL17027 sample preferred over GPL10526), pivoted from the saved GEO sample characteristics; pH is cerebellar tissue pH, agonal factor is the 0-4 score, expression is DLPFC; complete_case_for_models = pH, age, PMI and agonal factor all present (the rows the models of gse92538_ph_models.tsv use)"])
    cc = D[D.complete_case_for_models == "yes"]
    cnt = cc.groupby("diagnosis").agonal_factor.agg(lambda s: f"{int((s > 0).sum())}/{len(s)}").to_dict()
    F.note("f", f"{len(D)} unique subjects, {len(cc)} complete cases (legend: 157 subjects: match to the complete cases). Agonal factor > 0 by diagnosis among complete cases (n>0 / n): {cnt}. SUMMARY_AGONAL_AND_PH_BY_DATASET.md section 4 gives 19 of 71 controls, 6/40 MDD, 6/24 bipolar, 5/22 SZ.")
    pm = os.path.join(P_HAG, "gse92538_ph_models.tsv")
    F.write("fig08_f_ph_models.tsv", "f", ["subset", "formula", "term", "n", "coef", "se", "p"], rd(pm), [pm],
            ["OLS (HC3 errors) of tissue pH on diagnosis (Control reference) and covariates, per platform and for both platforms with one row per subject; rows 'agonal' = pH per agonal-score point. Run 2026-10-08/09 on remote"])
    r = [x for x in rd(pm) if x["subset"].startswith("both") and x["term"] == "Schizophrenia"]
    ag = [x for x in rd(pm) if x["subset"].startswith("both") and x["term"] == "agonal"]
    F.note("f", "both platforms, SZ pH term: " + "; ".join(f"{x['formula'].split('+', 1)[1].strip() if '+' in x['formula'] else x['formula']}: {float(x['coef']):+.3f} (SE {float(x['se']):.3f}, p {float(x['p']):.2f}, n {x['n']})" for x in r)
           + f"; agonal per point {float(ag[0]['coef']):+.3f} (SE {float(ag[0]['se']):.3f}, p {float(ag[0]['p']):.1g}). These tables are the saved results of 2026-10-08; SUMMARY_AGONAL_AND_PH_BY_DATASET.md is dated 2026-10-09 (remote). No file marks them as final.")
    pe = os.path.join(P_HAG, "ph_effect_by_dataset.tsv")
    F.write("fig08_f_ph_effect_by_dataset.tsv", "f", ["dataset", "model", "extra_terms", "n_sz", "n_ctl", "coef_sz", "se", "p", "agonal_recorded", "death_proxy_recorded"], rd(pe), [pe],
            ["schizophrenia minus control brain pH per dataset-region; A age + sex, B + PMI, C + PMI + death proxy, D + PMI + agonal factor, E + PMI + agonal + suicide (only GSE92538 records the agonal factor). GSE80655, GSE202537 and GSE53987 regions repeat the same donors"])
    pp = os.path.join(P_HAG, "ph_effect_pooled.tsv")
    F.write("fig08_f_ph_effect_pooled.tsv", "f", ["model", "n_sz", "n_ctl", "n_studies", "coef_sz", "se", "p"], rd(pp), [pp], ["pooled over 11 studies with study fixed effects; row 1 and 2 use one row per study and subject"])
    pm2 = os.path.join(P_HAG, "ph_by_manner_of_death.tsv")
    F.write("fig08_f_ph_by_manner_of_death.tsv", "f", ["dataset", "group", "manner", "n", "mean_ph", "sd_ph", "mean_age", "mean_pmi"], rd(pm2), [pm2], ["pH by diagnosis group and manner of death (GSE208338, GSE202537)"])
    pn = os.path.join(P_HAG, "ph_natural_death_only.tsv")
    F.write("fig08_f_ph_natural_death_only.tsv", "f", ["dataset", "model", "n_sz", "n_ctl", "coef_sz", "se", "p"], rd(pn), [pn], ["SZ minus control pH among natural deaths only, age + sex"])
    pe2 = os.path.join(P_HAG, "gse92538_expression_summary.tsv")
    F.write("fig08_f_expression_summary.tsv", "f", list(rd(pe2)[0].keys()), rd(pe2), [pe2],
            ["GSE92538 expression models M0-M3 per platform and diagnosis: genes at FDR < 0.05, correlation of the diagnosis t across models, agonal-term t over the Dai 2021 fever symbols"])
    F.note("f", "SUMMARY_AGONAL_AND_PH_BY_DATASET.md section 4 says 'closest, GPL17027 sample-level adjusting for PMI and agonal: +0.013 (p 0.88)'; the saved tables give +0.022 (p 0.83, n 92) for pH ~ dx + age + sex + pmi + agonal on GPL17027 (and -0.013, p 0.92 with agonal > 0) "
                "and +0.028 (p 0.79; 9 SZ / 45 control) for GSE92538_GPL96 model D in ph_effect_by_dataset.tsv. The +0.013 (p 0.88) is not in any exported table.")


def fig8_g(F):
    """Per-study agonal score of the hit lists (no panel yet): median +0.72, 46 of 55."""
    rows = rd(P_STUDY_CLAIMS)
    keep = ["label", "pmid", "year", "cluster", "lane", "ph_status", "n_signed_genes", "ph_score", "agonal_score", "agonal_share", "agonal_share_z", "agonal_share_p", "dachet_share", "dachet_share_z", "dachet_share_p", "themes"]
    F.write("fig08_g_hitlist_agonal_score_by_study.tsv", "g (no panel yet)", keep, rows, [P_STUDY_CLAIMS],
            ["agonal_score = mean over the study's listed genes with a direction of (SZ sign) x GTEx agonal t (slow vs fast natural death, cortical mean); positive = the list moves with the agonal state the way it moves in SZ. "
             "agonal_share = share of listed genes whose SZ direction equals the sign of the agonal t; dachet_share = share whose SZ direction equals the direction in dying fresh cortex at 24 h; *_z, *_p against 200 sign shuffles within the list. "
             "Source script analyses/study_claims_2026-10-08/build_study_claims.py (output built 2026-10-08 18:15 from recurrence_units.tsv / primary_lists.json e282904dcacb)"])
    v = [float(r["agonal_score"]) for r in rows if r["agonal_score"]]
    sh = [float(r["agonal_share"]) for r in rows if r["agonal_share"]]
    dsh = [float(r["dachet_share"]) for r in rows if r["dachet_share"]]
    n10 = sum(1 for r in rows if r["n_signed_genes"] and int(r["n_signed_genes"]) >= 10)
    F.note("g", f"'median +0.72, positive in 46 of 55' IS in a file: study_claims.tsv column agonal_score: {len(v)} studies with a score, median {statistics.median(v):.3f}, positive in {sum(x > 0 for x in v)}. Quoted in report/report.md line 538 and report/numbers.tsv line 219 "
                f"(source recorded there: analyses/study_claims_2026-10-08/study_claims.tsv), manuscript.md line 98 (Fig 8 paragraph; the HTML comment there says no panel exists). Agonal share median {statistics.median(sh):.3f} "
                f"(above 50% in {sum(x > 0.5 for x in sh)} of {len(sh)}; individually p < 0.05 in {sum(float(r['agonal_share_p']) < 0.05 for r in rows if r['agonal_share_p'])}); Dachet share median {statistics.median(dsh):.3f} (above 50% in {sum(x > 0.5 for x in dsh)} of {len(dsh)}, "
                f"p < 0.05 in {sum(float(r['dachet_share_p']) < 0.05 for r in rows if r['dachet_share_p'])}). Report says '55 studies list at least 10 genes with a direction'; the file has {n10} studies with n_signed_genes >= 10 and {len(v)} with an agonal score "
                "(the score is blank for studies without GTEx-measured listed genes or too few genes; the filter is in build_study_claims.py).")


# ================================================================== Figure 9
def fig9_a(F):
    """Agreement (Spearman rho of t) between cohorts with separate subjects against donors per group; single-nucleus classes and bulk."""
    pb = A("fairness_sample_level", "results", "b_pairs.tsv")
    pd_ = A("fairness_sample_level", "results", "datasets_used.tsv")
    pc = A("single_cell", "results", "sc_cross_pairs.tsv")
    BP = pdr(pb)
    DU = pdr(pd_)
    nmin = dict(zip(DU.dataset, DU[["n_sz", "n_ctl"]].min(axis=1)))
    BP = BP[(BP.level == "L2") & (BP.both_bulk.astype(str) == "True") & (BP.relation == "subject-separate")].copy()
    BP["n_min"] = [min(nmin.get(a, np.nan), nmin.get(b, np.nan)) for a, b in zip(BP.dataset_a, BP.dataset_b)]
    CP = pdr(pc)
    CP["n_min"] = CP[["n_sz_a", "n_ctl_a", "n_sz_b", "n_ctl_b"]].min(axis=1)
    E = CP[(CP.group_a == "endresz") & (CP.group_b == "endresz") & (CP.relation == "subject-separate")]
    rows = [{"data_type": "bulk dataset pair (Gemma, level L2)", "class_or_bulk": "bulk tissue", "unit_a": a, "unit_b": b, "relation": "subject-separate", "n_per_group": n, "rho_t": r}
            for a, b, n, r in zip(BP.dataset_a, BP.dataset_b, BP.n_min, BP.rho_t)]
    rows += [{"data_type": "single-nucleus cohort pair (Endresz 2026, one pipeline)", "class_or_bulk": c, "unit_a": a, "unit_b": b, "relation": "subject-separate", "n_per_group": n, "rho_t": r}
             for c, a, b, n, r in zip(E.gemma_class, E.source_a, E.source_b, E.n_min, E.rho_t)]
    F.write("fig09_a_agreement_pairs.tsv", "a", ["data_type", "class_or_bulk", "unit_a", "unit_b", "relation", "n_per_group", "rho_t"], rows, [pb, pd_, pc],
            ["one row per pair of datasets with separate subjects: Spearman rho of the SZ vs control t statistics over shared genes against n_per_group = the smaller of (SZ, control) group sizes over the two datasets. "
             "Bulk: Gemma datasets, level L2 (sz + age + sex + PMI), both bulk, 192 pairs; single nucleus: Endresz 2026 cohort pairs (7 cohorts), one row per Gemma cell class (14) and pair, relation subject-separate (266 rows)"])
    cls = E.groupby("gemma_class").rho_t.median()
    bins = []
    for lo, hi, lab in ((0, 10, "0-10"), (11, 30, "11-30"), (31, 1000, "31+")):
        bb = BP[(BP.n_min >= lo) & (BP.n_min <= hi)]
        ee = E.groupby(["source_a", "source_b"]).agg(rho_t=("rho_t", "median"), n_min=("n_min", "first")).reset_index()
        ee = ee[(ee.n_min >= lo) & (ee.n_min <= hi)]
        for dt, d in (("bulk dataset pairs", bb), ("single-nucleus cohort pairs (median over classes)", ee)):
            bins.append({"donors_per_group_bin": lab, "data_type": dt, "n_pairs": len(d), "median_rho_t": float(d.rho_t.median()) if len(d) else None,
                         "q25": float(d.rho_t.quantile(.25)) if len(d) else None, "q75": float(d.rho_t.quantile(.75)) if len(d) else None})
    F.write("fig09_a_agreement_by_n_bin.tsv", "a", ["donors_per_group_bin", "data_type", "n_pairs", "median_rho_t", "q25", "q75"], bins, [pb, pd_, pc],
            ["medians of fig09_a_agreement_pairs.tsv within bins of the smaller group size (same bins as sc_vs_bulk.py: 0-10, 11-30, 31+); single nucleus = median over classes within each cohort pair first"])
    crow = [{"class": c, "n_pair_rows": int((E.gemma_class == c).sum()), "median_rho_t": float(v)} for c, v in cls.items()]
    F.write("fig09_a_single_nucleus_by_class.tsv", "a", ["class", "n_pair_rows", "median_rho_t"], crow, [pc], ["median over the subject-separate Endresz cohort pairs, per Gemma cell class (14)"])
    ps = A("single_cell", "results", "sc_vs_bulk.tsv")
    F.note("a", f"single nucleus median rho(t) over {len(E)} class x pair rows {E.rho_t.median():.3f}, per-class medians {cls.min():.3f} to {cls.max():.3f} over {len(cls)} classes (legend 0.07, range 0.02 to 0.16, 14 classes: match); bulk {len(BP)} pairs median {BP.rho_t.median():.3f} (legend 0.11: match). "
                "Bins: " + "; ".join(f"{b['donors_per_group_bin']} {b['data_type'].split(' (')[0]} {b['median_rho_t']:.2f} (n={b['n_pairs']})" for b in bins if b['median_rho_t'] is not None)
                + ". Legend: 0.16 against 0.14 at 11 to 30 donors per group (match). report.md line 388 gives -0.01 against 0.05 at 0 to 10; the file value is -0.015 (rounds to -0.02 at two decimals; -0.01 truncated). Bin boundaries are those of sc_vs_bulk.py; sc_vs_bulk.tsv (saved) holds the same bin medians.")


def fig9_b(F):
    """PsychAD HBCC against Mount Sinai at 20 donors per group, and split halves."""
    p = A("single_cell", "results", "sc_split_half_summary.tsv")
    d = pdr(p)
    x = d[d.n_per_group == 20].copy()
    x["cell_class"] = x.cell_type
    cols = ["cell_type", "contrast", "comparison", "n_per_group", "reps", "rho_t", "rho_t_q25", "rho_t_q75", "rho_logfc", "top100_sign_in_other", "top200_fold", "fdr_overlap", "frac_reps_fdr_any"]
    F.df("fig09_b_psychad_20_per_group.tsv", "b", x[cols], [p],
         ["PsychAD pseudobulk, 20 donors per group per fit, medians over replicates (20 per cell): HBCC_vs_MSSM = one half from HBCC against one from Mount Sinai (main = with the dementia-term SZ donors, nodementia = without), split_HBCC / split_MSSM = split halves within one bank; "
          "contrasts sz, age, male (autosomal genes); cell_type 'total' = whole-tissue pseudobulk, the other eight are PsychAD classes"])
    cls8 = x[x.cell_type != "total"]
    med = cls8.groupby(["comparison", "contrast"]).rho_t.agg(["median", "count"]).reset_index()
    F.write("fig09_b_medians_over_8_classes.tsv", "b", ["comparison", "contrast", "median_over_classes", "n_classes"],
            [{"comparison": a, "contrast": b, "median_over_classes": m, "n_classes": n} for a, b, m, n in zip(med.comparison, med.contrast, med["median"], med["count"])], [p],
            ["median over the 8 PsychAD cell classes (whole-tissue 'total' excluded) of the per-class medians at 20 donors per group"])
    g = {(a, b): m for a, b, m in zip(med.comparison, med.contrast, med["median"])}
    F.note("b", "medians over the 8 classes at 20 per group: HBCC vs MSSM SZ %.3f (legend 0.04), age %.3f (legend 0.17); HBCC split halves SZ %.3f (legend 0.19); MSSM split halves SZ %.3f (legend -0.01); without dementia-term SZ donors HBCC vs MSSM SZ %.3f (single_cell_report.md: 0.12)."
           % (g[("HBCC_vs_MSSM_main", "sz")], g[("HBCC_vs_MSSM_main", "age")], g[("split_HBCC_main", "sz")], g[("split_MSSM_main", "sz")], g[("HBCC_vs_MSSM_nodementia", "sz")]))
    F.note("b", "source rows come from sc_split_half_reps.tsv.gz (remote run, 2026-09-27); the summary was rebuilt by sc_vs_bulk.py on 2026-10-07 (single_cell/sc_vs_bulk_meta.json). The older analyses/psychad_pseudobulk/results (2026-09-25, a_hbcc_vs_mssm.tsv etc.) are a different, earlier analysis and were not used.")


def fig9_c(F):
    """Bulk control-pH association against each subclass's single-nucleus SZ effect (Endresz 2026, Stouffer over cohorts)."""
    p1 = A("confounds", "results", "ph_singlecell_by_class.tsv")
    p2 = A("confounds", "results", "ph_singlecell_by_class_genomewide.tsv")
    a = pdr(p1).rename(columns={"cell_type": "subclass", "n_genes": "n_genes_recurrent_universe", "r": "r_ph_vs_sn_effect_5plus_universe", "median_n_cohorts": "median_n_cohorts"})
    b = pdr(p2).rename(columns={"cell_type": "subclass", "n_genes": "n_genes_genomewide", "r": "r_ph_vs_sn_effect_genomewide"})
    m = a.merge(b, on="subclass")
    F.df("fig09_c_ph_vs_sn_effect_by_subclass.tsv", "c", m, [p1, p2],
         ["Pearson r (see ph_singlecell.R) between each gene's bulk control-only pH Z (five collections) and its single-nucleus SZ effect, Stouffer-combined (weights sqrt(N)) over the Endresz 2026 cohorts holding the subclass, per Endresz subclass; "
          "restricted universe = the genes in 5+ lists (study_brain_5plus, 1,641 genes) that have a control-only pH Z (Z_pH_controls in recurrent_genes_pH.tsv; 1,583 genes; was 1,431 of 1,478); genome-wide = 15,187 genes measured in >= 4 of 5 bulk pH datasets. Rerun 2026-10-11 on the 515-study corpus"])
    r1 = m.r_ph_vs_sn_effect_5plus_universe
    r2 = m.r_ph_vs_sn_effect_genomewide
    F.note("c", f"{len(m)} subclasses; restricted universe (5+ lists genes with a control pH Z): negative in {int((r1 < 0).sum())}, median {r1.median():.3f}, range {r1.min():.3f} to {r1.max():.3f}; genome-wide: negative in {int((r2 < 0).sum())}, median {r2.median():.3f}. "
                "Per-subclass n_genes is the number of those genes measured in that subclass. Tables rebuilt 2026-10-11 by confounds/ph_singlecell_2026-10-02/scripts/ph_singlecell.R on the 5+ lists set of the 515-study corpus (1,641 genes, 1,583 with a control pH Z; before 1,478 and 1,431). r is Pearson (R cor default).")


def fig9_d(F):
    """Noise-corrected SD of the true SZ log2 fold changes by data type (ashr)."""
    pe = A("effect_sizes", "effects_by_study.tsv")
    pg = A("effect_sizes", "data", "by_group.tsv")
    ps = A("effect_sizes", "summary_by_type.tsv")
    E = pdr(pe)
    H = E[E.headline_C.astype(str) == "True"].copy()
    cols = ["study", "label", "type", "N", "n_cases", "n_controls", "C_sd", "C_top_group", "C_top_sd", "C_mean_abs_pm", "C_frac_lfsr05", "C_pi0", "C_n_fitted", "C_se_method", "full_provenance", "flags"]
    F.df("fig09_d_effect_sd_by_study.tsv", "d (large points)", H[cols], [pe],
         ["one row per headline study and donor set (view C): C_sd = SD of the distribution of true SZ vs control log2 fold changes that ashr estimates from every tested gene's log2 FC and SE (median over cell types or regions); N = cases + controls; "
          "log scale; C_sd = 0 sits on the dashed line. effect_sizes/ built 2026-10-07 20:27 from primary_lists.json e282904dcacb (current)"])
    G = pdr(pg)
    key = dict(zip(E.full_source, zip(E.label, E.type, E.N, E.role, E.headline_C)))
    sn = []
    for r in G.itertuples(index=False):
        k = key.get(r.source)
        if k and k[1] in ("single-nucleus brain",) and str(k[4]) == "True":
            sn.append({"study": r.source, "label": k[0], "type": k[1], "N": k[2], "cell_type": r.group, "n_genes_C": r.n_genes_C, "sd_g": r.sd_g, "pi0": r.pi0, "mean_abs_pm": r.mean_abs_pm, "frac_lfsr05": r.frac_lfsr05, "lik": r.lik})
    F.write("fig09_d_effect_sd_single_nucleus_cell_types.tsv", "d (small points)", ["study", "label", "type", "N", "cell_type", "n_genes_C", "sd_g", "pi0", "mean_abs_pm", "frac_lfsr05", "lik"], sn, [pg, pe],
            ["per cell type SD of true effects (ashr sd_g) for the three headline single-nucleus brain studies"])
    EC = E[E.role == "endresz_cohort"]
    ec = []
    for r in EC.itertuples(index=False):
        sub = G[G.source == r.full_source] if isinstance(r.full_source, str) else G.iloc[0:0]
        ec.append({"cohort": r.label, "N": r.N, "n_cases": r.n_cases, "n_controls": r.n_controls, "C_sd": r.C_sd, "n_cell_types": len(sub)})
    F.write("fig09_d_endresz_cohorts.tsv", "d (panel b of the effect-size figure)", ["cohort", "N", "n_cases", "n_controls", "C_sd", "n_cell_types"], ec, [pe, pg],
            ["the seven Endresz 2026 snRNA-seq cohorts in one pipeline (not headline: MSSM1 and MSSM2 share Mount Sinai donors); effect_sizes Fig 1b / Fig 4 open circles"])
    S = pdr(ps)
    F.df("fig09_d_summary_by_type.tsv", "d", S, [ps], ["medians and bootstrap 95% CI by type, as in effect_sizes/summary_by_type.tsv"])
    med = H.groupby("type").C_sd.agg(["count", "median"])
    F.note("d", "headline_C studies by type (n, median SD): " + "; ".join(f"{t} {int(c)} / {m:.3f}" for t, (c, m) in med.iterrows()) + ". CAPTIONS.md: brain array 0.035 (14), RNA-seq 0.048 (13), single-nucleus 0.021 (3), Endresz cohorts 0.070 (7): match. "
                f"Single-nucleus small points: {len(sn)} cell-type rows from {len({r['study'] for r in sn})} studies. This panel is the effect_sizes Figure 1 (and 4); the legend 'large points sized by N, small points single-nucleus cell types' fits Fig 1a.")


# ================================================================== Figure 10
def fig10(F):
    """Cell-type proportions: what exists as a table."""
    # a. Endresz 2026 abundance (crumblr, per supertype, per dataset), from the bioRxiv Supplemental Table 4
    z = os.path.join(REPO, "data", "raw", "doi_10.64898_2026.09.18.752621", "biorxiv", "media-4.zip")
    import openpyxl
    wb = openpyxl.load_workbook(io.BytesIO(zipfile.ZipFile(z).read("Supplemental_table_4.xlsx")), read_only=True)
    rows = []
    for ws in wb:
        it = ws.iter_rows(values_only=True)
        hdr = next(it)
        for r in it:
            if r[0] is None:
                continue
            ct = str(r[0])
            rows.append({"dataset": ws.title, "supertype": ct, "subclass": re.sub(r"_\d+$", "", ct), "estimate": r[1], "SE": r[2], "p": r[3], "FDR": r[4]})
    F.write("fig10_a_endresz_abundance_by_dataset.tsv", "a / b", ["dataset", "supertype", "subclass", "estimate", "SE", "p", "FDR"], rows, [z],
            ["Endresz 2026 (bioRxiv 10.64898/2026.09.18.752621) Supplemental Table 4: donor-level proportion of each supertype, SCZ vs control, crumblr adjusted for age, sex and PMI, per dataset and fixed-effect 'Meta-analysis'; "
             "estimate = SCZ coefficient (log-scale), positive = more abundant in SCZ. Donors younger than 70 years only (HBCC older than 20) in every row: the paper's inclusion rule, so no 'with the oldest donors' version exists in this table. "
             "Source file is a data/raw supplement, not an analysis output of this repository"])
    meta = [r for r in rows if r["dataset"] == "Meta-analysis"]
    sst = sorted([r for r in meta if r["supertype"].startswith("Sst_")], key=lambda r: r["p"])[:5]
    F.note("a", f"Endresz Supplemental Table 4: {len(rows)} rows over {len({r['dataset'] for r in rows})} sheets (7 cohorts + Xenium + Meta-analysis), {len(meta)} supertypes in the meta-analysis. Lowest-p Sst supertypes in the meta-analysis: "
                + ", ".join(f"{r['supertype']} {float(r['estimate']):+.2f} (FDR {float(r['FDR']):.3g})" for r in sst) + ". Text of the paper: Sst_2, Sst_22, Sst_25, Sst_3 FDR < 0.05 (beta -0.19 to -0.27); no Pvalb or non-neuronal supertype at FDR < 0.10.")
    # a. Lee 2026 PsychAD compositional variation, SCZ
    pl = os.path.join(REPO, "data", "raw", "42778699", "publisher", "moesm3", "2024-07-15352A-s3", "PsychAD_SupplementaryTable3.csv")
    L = pd.read_csv(pl)
    Ls = L[L.coef == "SCZ"].rename(columns={"AnnoLevel": "annotation_level", "assay": "cell_type", "std.error": "std_error"})
    F.df("fig10_a_lee_psychad_composition_scz.tsv", "a", Ls[["annotation_level", "cell_type", "coef", "estimate", "std_error", "FDR", "log10FDR"]], [pl],
         ["Lee et al. 2026 (PsychAD; Nature) Supplementary Table 3, 'Compositional variation analysis across 8 cross-disorder traits', rows for SCZ only (all 10 traits in the file; 27 subclasses and 65 subtypes for SCZ). "
          "estimate = SCZ coefficient on the cell-type proportion model; positive = more abundant in SCZ. The paper's age handling is not in this table (see STATUS.md)"])
    pv = Ls[Ls.cell_type.str.contains("PVALB|SST")][["cell_type", "estimate", "FDR"]]
    F.note("a", "Lee Supp Table 3 SCZ rows: " + "; ".join(f"{r.cell_type} {r.estimate:+.3f} (FDR {r.FDR:.3g})" for r in pv.itertuples()) + ". The only PVALB-class rows with FDR < 0.05 are IN_PVALB_CHC, IN_PVALB_HSPA1A (subtype); IN_PVALB (subclass) +0.032, FDR 0.40. "
                "The manuscript statement 'Lee report an increase in PVALB cells' is read from the paper's text/Fig; which Lee result carries it was not traced here.")
    # c. PsychAD donor-level composition (repo analysis, 2026-09-25)
    pc = A("psychad_cellxgene", "results", "composition.tsv")
    F.write("fig10_c_psychad_donor_fraction_scz_vs_control.tsv", "c (context)", ["class", "bank", "n_sz", "n_ctl", "median_frac_sz", "median_frac_ctl", "wilcox_p"], rd(pc), [pc],
            ["this repository's own PsychAD whole-tissue composition (psychad_cellxgene/scripts/composition.R, 2026-09-25): median fraction of each donor's nuclei per class, SZ vs control, Wilcoxon p; 8 classes x 2 banks; donors of the HBCC_main and MSSM_main fits"])
    # d. Endresz audit (expression heterogeneity, NOT proportions)
    d = A("endresz_audit_2026-10-04", "results")
    for fn, pan, lab in (("sst_per_dataset.tsv", "audit", "Sst deficit genes (196): median logFC and genes at FDR < 0.1 per dataset"), ("pvalb_per_dataset.tsv", "audit", "Pvalb deficit genes (245) per dataset"),
                         ("pvalb_loo_summary.tsv", "audit", "Pvalb leave-one-dataset-out"), ("sst_cohort_rows.tsv", "audit", "Sst rows per cohort")):
        p = os.path.join(d, fn)
        r = rd(p)
        F.write("fig10_audit_" + fn, pan, list(r[0].keys()), r, [p], [lab + " (endresz_audit_2026-10-04; DE logFC heterogeneity, not proportions)"])
    ps = os.path.join(d, "sst_loo.tsv")
    S = pdr(ps)
    sl = S.groupby("dropped").agg(median_meta_logFC_loo=("meta_logFC_loo", "median"), n_genes=("gene", "nunique"), n_genes_flipped_pos=("sign_flip_pos", lambda s: int((s.astype(str) == "True").sum()))).reset_index()
    F.df("fig10_audit_sst_loo_summary.tsv", "audit", sl, [ps], ["Sst leave-one-dataset-out summary computed from sst_loo.tsv (the folder has the summary only for Pvalb)"])
    pd2 = os.path.join(d, "de_pair_support.tsv")
    DP = pdr(pd2)
    dist = DP.n_datasets_sig.value_counts().sort_index()
    F.write("fig10_audit_de_support_distribution.tsv", "audit", ["datasets_reaching_fdr_0.1", "n_gene_celltype_pairs", "share"], [{"datasets_reaching_fdr_0.1": int(k), "n_gene_celltype_pairs": int(v), "share": v / len(DP)} for k, v in dist.items()], [pd2],
            ["number of datasets that also reach FDR < 0.1 for each meta-analysis FDR < 0.1 (gene, cell type) pair (6,898 pairs); computed from de_pair_support.tsv"])
    # e. Mulvey/Bastos meta-analysis moderators
    m = A("interneuron_meta_moderator_2026-10-04", "results")
    p1, p2, p3 = (os.path.join(m, f) for f in ("moderators.tsv", "pmi_noise_sim.tsv", "datapoints.tsv"))
    F.write("fig10_mulvey_moderators.tsv", "context", list(rd(p1)[0].keys()), rd(p1), [p1], ["OLS with study-clustered SE and leave-one-out: moderators of the Mulvey/Bastos GABA interneuron meta-analysis t' scores (BastosLab/sim-map, commit 62d20dd), 30 models"])
    F.write("fig10_mulvey_pmi_noise_sim.tsv", "context", list(rd(p2)[0].keys()), rd(p2), [p2], ["expected residual PMI confounding under recorded-PMI adjustment, scenario grid (E undocumented hours assumed, not identified)"])
    F.write("fig10_mulvey_datapoints.tsv", "context", list(rd(p3)[0].keys()), rd(p3), [p3], ["330 study x area x cell type x layer datapoints of the meta-analysis (t' = (SZ - HC) / sqrt(SEM_SZ^2 + SEM_HC^2)) with study covariates"])
    F.note("b", "NOT PRODUCED: 'the same comparison after removing the oldest donors' (Figures.md Fig 10b). No table exists that compares proportion results with and without the age cut for Lee or for Endresz. Endresz's Supplemental Table 4 is already restricted to donors < 70 y (HBCC > 20 y), so it cannot show the effect of removing them; "
                "Endresz's released files hold no donor-level cell counts (Supplemental Table 1 lists donor, age, sex, diagnosis, PMI and total cells only; the Zenodo parquet holds DE results). A with/without-age-cut comparison needs donor-by-cell-type counts for every donor: for PsychAD these are in the repo for the HBCC and Mount Sinai donors "
                "(analyses/psychad_pseudobulk/cells_by_donor_celltype.tsv, 31,387 rows; donors.tsv) but not for Lee's full 1,494-donor cohort; for Endresz's other cohorts they are not in the repo. A local R run of a crumblr-type model on the PsychAD donors with and without an age cut would be a new script (estimated cost: low, minutes on the laptop); "
                "Endresz's own with/without comparison needs their donor-level counts (authors or Synapse) and cannot be done from the repository.")
    F.note("a", "Lee-versus-Endresz on one axis does not exist as a table. The nearest plottable tables are fig10_a_lee_psychad_composition_scz.tsv (Lee, PsychAD all cohorts, 92 SCZ rows) and fig10_a_endresz_abundance_by_dataset.tsv (Endresz, per dataset and meta-analysis, 1,103 rows); their cell-type labels differ (Lee subclass/subtype names such as IN_PVALB_CHC; Endresz supertypes such as Sst_25), so a mapping is needed before the two are placed on one panel (celltype_map.tsv in analyses/single_nucleus_study_comparison/ maps Gemma classes, not these).")
    F.note("claims-without-table", "Paragraph 'Claims about cell-type proportions depend on analysis choices' (manuscript.md line 114) relies on: (1) Lee's PVALB-increase claim (Lee Supp Table 3 exported); (2) Mulvey 2025 reports in the other direction (no table: a literature statement; the Mulvey/Bastos moderator tables exported cover PMI and pH moderation of density and mRNA t', not proportions); (3) Endresz: 'if the oldest donors are removed ... consistent only with a decrease in SST' (Endresz Supp Table 4 exported; the with/without-age-cut contrast is in no file); (4) 'Changes in cell proportion can arise from differential capture or identification ... found nothing conclusive' (Fig 10): no analysis file; analyses/endresz_audit_2026-10-04 audits DE heterogeneity of Sst/Pvalb deficit genes across Endresz datasets and the donor-overlap concern, not proportions or capture (the paper logic map, notes/Schizophrenia transcriptome replication - paper logic map.md line 78, lists it as the source and says 'Keep as a limit, not a Results claim'; manuscript.md has the same CHECK comment). The only capture-adjacent table in the repo is the PsychAD donor fraction table (fig10_c).")
    F.note("context", "Mulvey/Bastos tables: moderators.tsv 30 models; pmi_noise_sim.tsv 72 rows; datapoints.tsv 330 rows. These are density and mRNA meta-analysis (counted-cell and expression) results with PMI and pH moderators, not Lee-versus-Endresz proportions.")


def main():
    F8 = Figure("fig08")
    for fn in (fig8_a, fig8_b, fig8_c, fig8_d, fig8_e, fig8_f, fig8_g):
        fn(F8)
    F8.note("STALE / COMPUTE", "No panel needs a rerun. Inputs whose recorded shas differ from today's but whose values are unchanged: combined_controls_pH.tsv (pH Z identical, so ph_pHindex_2026-10-03 and ph_singlecell_2026-10-02 stand). "
            "Population-dependent inputs are current: recurrence_genes.tsv (Dachet direction, 2026-10-08 16:16 after recurrence.tsv 14:31), study_claims.tsv (primary_lists.json e282904dcacb = current). "
            "interneuron_ph_2026-10-01 Task A (enrichment.tsv, not exported) and per_sample_* were built on combined_controls_pH sha 8e029621de28 and the 2026-10-03 population; the exported subject-level scores (marker groups, not population-dependent) are 2026-10-08. "
            "Fig 8f numbers are the remote results of 2026-10-08/09; nothing in the folder marks them as final. Light recomputations done here: 8a gene table (Dachet tar), 8c class rho, medians of 8b pair tables.")
    F8.finish()
    F9 = Figure("fig09")
    for fn in (fig9_a, fig9_b, fig9_c, fig9_d):
        fn(F9)
    F9.note("STALE / COMPUTE", "Nothing stale; nothing needs compute. 9a and 9b depend on sc_split_half_reps.tsv.gz and sc_cross_pairs.tsv (September remote runs, summaries rebuilt 2026-10-07), 9c on 2026-10-03 tables (pH Z unchanged), 9d on effect_sizes built 2026-10-07 from the current primary_lists.json.")
    F9.finish()
    F10 = Figure("fig10")
    fig10(F10)
    F10.note("STALE / COMPUTE", "Not plottable as the legend describes: 10b (with/without the oldest donors) and the differential-capture check have no table. Compute to route if wanted: a PsychAD donor-level proportion model with and without an age cut (local R, minutes; new script); Endresz's own with/without contrast needs their donor-level counts (not in the repository).")
    F10.finish()


if __name__ == "__main__":
    main()
