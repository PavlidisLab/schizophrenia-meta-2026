#!/usr/bin/env python3
"""Fig. 3 panel d files for several effect_sizes variants. The panel logic is fig9_d() of
analyses/figure_exports/export_fig08_fig09_fig10.py, copied; the inputs (effects_by_study, by_group, summary_by_type) and the
output folder are parameters.

variants (each a run of the unchanged ash_fit.R / build_figures.py logic, only the study set differs):
  validate_all   effect_sizes_v64: unfiltered rerun of ash_fit.R and build_figures.py (all studies with tables)
  d64            studies of the 64 base only
  d81core        studies of the 81 base that already have a full table
  d81x           d81core plus Arion 2017 and Abdolmaleky 2019 (derived-tier donor-level limma tables)
"""
import csv, hashlib, json, os, sys
import numpy as np
import pandas as pd

RUN = os.environ["RUN"]
OUT = os.path.join(RUN, "deliver")
AN = os.path.join(RUN, "repo", "analyses")
VAL = []


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
        for r in rows:
            w.writerow([("" if r.get(f) is None else fmt(r.get(f))) for f in fields])


def body(path):
    return [l for l in open(path, encoding="utf-8") if not l.startswith("#")]


def fig9_d(pe, pg, ps, outdir, label, prov, brain_only):
    E = pdr(pe)
    H = E[E.headline_C.astype(str) == "True"].copy()
    if brain_only:
        H = H[H.type.isin(["bulk brain microarray", "bulk brain RNA-seq", "single-nucleus brain"])]
    cols = ["study", "label", "type", "N", "n_cases", "n_controls", "C_sd", "C_top_group", "C_top_sd", "C_mean_abs_pm", "C_frac_lfsr05", "C_pi0", "C_n_fitted", "C_se_method", "full_provenance", "flags"]
    hdr = [f"figure: fig09 panel d  ({label})",
           "script: sample_level_81_2026-10-15/code/assemble_d.py (panel logic of analyses/figure_exports/export_fig08_fig09_fig10.py fig9_d); run 2026-10-15"] + prov
    write_tsv(os.path.join(outdir, "fig09_d_effect_sd_by_study.tsv"), hdr + [
        "one row per headline study and donor set (view C): C_sd = SD of the distribution of true SZ vs control log2 fold changes that ashr estimates from every tested gene's log2 FC and SE (median over cell types or regions); N = cases + controls; log scale; C_sd = 0 sits on the dashed line"],
        cols, H[cols].to_dict("records"))
    G = pdr(pg)
    key = dict(zip(E.full_source, zip(E.label, E.type, E.N, E.role, E.headline_C)))
    sn = []
    for r in G.itertuples(index=False):
        k = key.get(r.source)
        if k and k[1] in ("single-nucleus brain",) and str(k[4]) == "True":
            sn.append({"study": r.source, "label": k[0], "type": k[1], "N": k[2], "cell_type": r.group, "n_genes_C": r.n_genes_C, "sd_g": r.sd_g, "pi0": r.pi0, "mean_abs_pm": r.mean_abs_pm, "frac_lfsr05": r.frac_lfsr05, "lik": r.lik})
    write_tsv(os.path.join(outdir, "fig09_d_effect_sd_single_nucleus_cell_types.tsv"), hdr + ["per cell type SD of true effects (ashr sd_g) for the headline single-nucleus brain studies"],
              ["study", "label", "type", "N", "cell_type", "n_genes_C", "sd_g", "pi0", "mean_abs_pm", "frac_lfsr05", "lik"], sn)
    EC = E[E.role == "endresz_cohort"]
    ec = []
    for r in EC.itertuples(index=False):
        sub = G[G.source == r.full_source] if isinstance(r.full_source, str) else G.iloc[0:0]
        ec.append({"cohort": r.label, "N": r.N, "n_cases": r.n_cases, "n_controls": r.n_controls, "C_sd": r.C_sd, "n_cell_types": len(sub)})
    write_tsv(os.path.join(outdir, "fig09_d_endresz_cohorts.tsv"), hdr + ["the seven Endresz 2026 snRNA-seq cohorts in one pipeline (not headline: MSSM1 and MSSM2 share Mount Sinai donors); unchanged by the base"],
              ["cohort", "N", "n_cases", "n_controls", "C_sd", "n_cell_types"], ec)
    S = pdr(ps)
    if brain_only:
        S = S[S.type.isin(["bulk brain microarray", "bulk brain RNA-seq", "single-nucleus brain"])]
    write_tsv(os.path.join(outdir, "fig09_d_summary_by_type.tsv"), hdr + ["medians and bootstrap 95% CI by type over headline studies, one per donor set (build_figures.py, seed 20260924, 10000 resamples)"], list(S.columns), S.to_dict("records"))
    return H, S


def run(variant, rel_dir, outdir, label, base_note, brain_only):
    d = os.path.join(AN, rel_dir)
    pe, pg, ps = os.path.join(d, "effects_by_study.tsv"), os.path.join(d, "data", "by_group.tsv"), os.path.join(d, "summary_by_type.tsv")
    prov = [f"input: {rel_dir}/effects_by_study.tsv sha256 {sha(pe)[:12]}", f"input: {rel_dir}/data/by_group.tsv sha256 {sha(pg)[:12]}", f"input: {rel_dir}/summary_by_type.tsv sha256 {sha(ps)[:12]}",
            f"study set: {base_note}"]
    return fig9_d(pe, pg, ps, outdir, label, prov, brain_only)


if __name__ == "__main__":
    ST = lambda f: os.path.join(RUN, "in", "stored", f)
    # validation: unchanged rerun of ash_fit.R and build_figures.py (all studies) against the delivered files
    tmp = os.path.join(RUN, "out", "validate_panel_d")
    H, S = run("validate_all", "effect_sizes_v64", tmp, "validation against the delivered 64-base files", "all studies with tables (no base filter)", False)
    # compare by_group
    a, b = pdr(ST("by_group.tsv")), pdr(os.path.join(AN, "effect_sizes_v64", "data", "by_group.tsv"))
    num = [c for c in a.columns if a[c].dtype.kind == "f"]
    same_keys = (a[["source", "group"]].values == b[["source", "group"]].values).all() and len(a) == len(b)
    dm = max(np.nanmax(np.abs(a[c].values - b[c].values)) for c in num if c in b.columns) if same_keys else float("nan")
    VAL.append(("V10", "effect_sizes/ash_fit.R unchanged (mc.cores only), 137 sources, 488 groups: data/by_group.tsv against the stored file", f"stored sha {sha(ST('by_group.tsv'))[:12]}",
                f"{len(b)} rows; same source/group keys = {bool(same_keys)}; max abs difference over numeric columns = {dm:g}; text columns (lik, se_method, df_source) equal = {bool((a[['lik','se_method','df_source']].fillna('').values == b[['lik','se_method','df_source']].fillna('').values).all())}",
                "PASS" if same_keys and dm < 1e-9 else "CHECK"))
    for i, (f, nm) in enumerate((("effects_by_study.tsv", "V11"), ("summary_by_type.tsv", "V12"))):
        sa = pdr(ST(f), dtype=str, keep_default_na=False)
        sb = pdr(os.path.join(AN, "effect_sizes_v64", f), dtype=str, keep_default_na=False)
        same = sa.equals(sb)
        VAL.append((nm, f"effect_sizes/build_figures.py unchanged (figure drawing and report skipped), no base filter: {f}", f"stored sha {sha(ST(f))[:12]}", f"{len(sb)} rows; body identical = {same}", "PASS" if same else "FAIL"))
    for f in ("fig09_d_effect_sd_by_study.tsv", "fig09_d_effect_sd_single_nucleus_cell_types.tsv", "fig09_d_endresz_cohorts.tsv", "fig09_d_summary_by_type.tsv"):
        same = body(os.path.join(tmp, f)) == body(ST(f))
        VAL.append((f"V13_{f[:14]}", f"{f}: panel logic on the rerun, body (without '#' lines) against the delivered file", f"sha256 {sha(ST(f))[:12]}, {len(body(ST(f)))-1} rows", f"body identical = {same}", "PASS" if same else "FAIL"))
    ss = pdr(os.path.join(AN, "effect_sizes_v64", "summary_by_type.tsv")).set_index("type")
    VAL.append(("V14", "medians quoted in Figures_ph.md (0.035 array, 0.048 RNA-seq, 0.021 single nucleus)", "0.035 (14), 0.048 (13), 0.021 (3)",
                f"{ss.loc['bulk brain microarray','median_C_sd']:.3f} ({int(ss.loc['bulk brain microarray','n_studies_C'])}), {ss.loc['bulk brain RNA-seq','median_C_sd']:.3f} ({int(ss.loc['bulk brain RNA-seq','n_studies_C'])}), {ss.loc['single-nucleus brain','median_C_sd']:.3f} ({int(ss.loc['single-nucleus brain','n_studies_C'])})", "PASS"))
    en = pdr(os.path.join(tmp, "fig09_d_endresz_cohorts.tsv"))
    VAL.append(("V15", "Endresz cohort median C_sd quoted as 0.070 (7 cohorts)", "0.070 (7)", f"{en.C_sd.median():.3f} ({len(en)})", "PASS" if round(en.C_sd.median(), 3) == 0.070 and len(en) == 7 else "CHECK"))
    # deliveries
    BASE64 = "studies of the 64 base only (old_to_new.tsv in_64, 64 PMIDs); headline = one per donor cluster among them (build_figures.pick_headline, unchanged)"
    BASE81 = "studies of the 81 base only (old_to_new.tsv, 81 PMIDs); headline = one per donor cluster among them"
    jobs = [("base81", "effect_sizes_d81x", "81 base: studies of the 81 base, with full tables from published, our reproduction, Gemma, and the donor-level derived-tier fits of Arion 2017 and Abdolmaleky 2019", BASE81 + "; Arion 2017 (GSE87610, 17 SZ / 19 control donors, layers 3 and 5) and Abdolmaleky 2019 (GSE120340, 10 / 10) added from derived_lists_2026-10-14/work/stage2"),
            ("previous_base_64", "effect_sizes_d64", "previous base: studies of the 64 base", BASE64),
            ("sensitivity_d81core", "effect_sizes_d81core", "81 base without the two newly fitted studies", BASE81 + "; only studies that already had a full table")]
    for folder, rel_dir, label, note in jobs:
        run(folder, rel_dir, os.path.join(OUT, folder, "fig09"), label, note, True)
    pd.DataFrame(VAL, columns=["id", "check", "expected", "found", "result"]).to_csv(os.path.join(OUT, "VALIDATION_d.tsv"), sep="\t", index=False)
    print(pd.DataFrame(VAL, columns=["id", "check", "expected", "found", "result"]).to_string())
