#!/usr/bin/env python3
"""NEW 2026-10-18. One script regenerates the summary table and the figures of ph_levels_2026-10-18 from the TSVs in this folder.
Reads out/psychad_levels_rho.tsv, out/endresz_levels_rho.tsv, out/nsw_levels_rho.tsv (written by the R scripts); writes
levels_summary_main_rows.tsv, fig_levels_pH.svg/.png, fig_levels_agonal.svg/.png (+ _data.tsv, caption text)."""
import os, sys, hashlib, datetime
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); RUN = os.path.dirname(HERE); OUT = os.path.join(RUN, "out")
ACC = {"all": "#2563eb", "rec383": "#10b981", "rec107": "#f59e0b"}
LAB = {"all": "all shared genes", "rec383": "383 recurrent genes", "rec107": "107 genes in 10+ lists"}
LEVELS = ["total", "total_composition", "class", "subclass"]
XLAB = {"total": "total\n(all nuclei)", "total_composition": "total, composition\nadjusted", "class": "class", "subclass": "subclass"}
GREY, TEXT, SUBTLE, GRID = "#9ca3af", "#1f2937", "#6b7280", "#e5e7eb"
plt.rcParams.update({"svg.fonttype": "none", "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "sans-serif"],
    "font.size": 9, "axes.edgecolor": SUBTLE, "axes.labelcolor": TEXT, "xtick.color": SUBTLE, "ytick.color": SUBTLE,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "axes.grid.axis": "y", "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.titleweight": "normal", "axes.titlesize": 9, "figure.facecolor": "white", "axes.facecolor": "white"})

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()

def load(name):
    return pd.read_csv(os.path.join(OUT, name), sep="\t", comment="#")

def main():
    ps = load("psychad_levels_rho.tsv"); en = load("endresz_levels_rho.tsv"); ns = load("nsw_levels_rho.tsv"); bt = load("batiuk_levels_rho.tsv")
    for d in (ps, en, ns, bt): d["panel"] = d["dataset"] + " | " + d["population"]
    allr = pd.concat([ns, bt, ps, en], ignore_index=True)
    allr.to_csv(os.path.join(RUN, "levels_all_rows.tsv"), sep="\t", index=False)
    main_cell = lambda r: (r.cell_type == "all") or (r.cell_type == "MEDIAN_over_cell_types_weighted_nuclei") or (r.dataset == "Endresz_meta" and r.cell_type == "MEDIAN_over_cell_types_unweighted")
    mk = allr.apply(main_cell, axis=1)
    mainr = allr[mk & allr.level.isin(LEVELS)].copy()
    mainr["row_kind"] = np.where(mainr.cell_type == "all", "single unit", "median over cell types")
    keep = ["dataset", "population", "analysis", "level", "row_kind", "axis", "gene_set", "n_genes", "n_samples", "n_nuclei", "rho", "lo", "hi", "se_fisher", "pct_lo", "pct_hi", "n_boot_ok"]
    hdr = ["NEW 2026-10-18 main rows of the level analysis; built %s by scripts/build_figure.py" % datetime.date.today().isoformat(),
           "level summary rows: single pseudobulk (total, total_composition) or weighted-by-nuclei median over cell types (class, subclass); Endresz: unweighted median over 24 subclasses, no interval",
           "inputs sha256(16): " + "; ".join("%s %s" % (n, sha(os.path.join(OUT, n))[:16]) for n in ("psychad_levels_rho.tsv", "endresz_levels_rho.tsv", "nsw_levels_rho.tsv", "batiuk_levels_rho.tsv")),
           "rho = Spearman across genes; positive pH axis = higher expression at higher pH; interval = Fisher-scale +- 1.96 donor-bootstrap SD (Endresz: none)"]
    with open(os.path.join(RUN, "levels_summary_main_rows.tsv"), "w") as f:
        for h in hdr: f.write("# " + h + "\n")
        mainr[keep].to_csv(f, sep="\t", index=False)
    for axis, tag, panels in (("pH", "pH", None), ("agonal", "agonal", None)):
        figure(allr, axis, tag)

PANELS = [("NSW_OFC", "controls"), ("NSW_OFC", "all donors"), ("Batiuk_2022", "all donors"), ("PsychAD_HBCC", None), ("PsychAD_MSSM", None), ("Endresz_meta", None)]

def figure(allr, axis, tag):
    titles = {"NSW_OFC|controls": "NSW OFC, pH association\ncontrols only", "NSW_OFC|all donors": "NSW OFC, pH association\nall donors, diagnosis adjusted",
              "Batiuk_2022|all donors": "Batiuk 2022 (11 donors)\nall donors, dx adjusted; neurons", "PsychAD_HBCC|": "PsychAD HBCC\nSZ effect", "PsychAD_MSSM|": "PsychAD MSSM\nSZ effect", "Endresz_meta|": "Endresz meta-analysis\nSZ effect, subclass level only"}
    plist = list(PANELS)
    fig, axs = plt.subplots(1, len(plist), figsize=(3.1 * len(plist), 4.3), sharey=True, squeeze=False); axs = axs[0]
    rows_out = []
    for ax, (ds, pop) in zip(axs, plist):
        d = allr[(allr.dataset == ds) & (allr.axis == axis)]
        if pop is not None: d = d[d.population.str.startswith("controls" if pop == "controls" else "all donors")]
        xs = np.arange(len(LEVELS)); off = {"all": -0.18, "rec383": 0.0, "rec107": 0.18}
        if ds == "Endresz_meta": xs_used = [3]
        else: xs_used = list(range(4))
        for gs in ("all", "rec383", "rec107"):
            xx, yy, lo, hi = [], [], [], []
            for i in xs_used:
                lev = LEVELS[i]; q = d[(d.gene_set == gs) & (d.level == lev)]
                if ds == "Endresz_meta": q = q[q.cell_type == "MEDIAN_over_cell_types_unweighted"]
                elif lev in ("class", "subclass"): q = q[q.cell_type == "MEDIAN_over_cell_types_weighted_nuclei"]
                else: q = q[q.cell_type == "all"]
                if len(q) == 0: continue
                r = q.iloc[0]; xx.append(i + off[gs]); yy.append(r.rho); lo.append(r.lo); hi.append(r.hi)
                rows_out.append(dict(panel=ds + ("|" + pop if pop else ""), axis=axis, level=lev, gene_set=gs, x=i + off[gs], rho=r.rho, lo=r.lo, hi=r.hi))
            if not xx: continue
            yy = np.array(yy, float); lo = np.array(lo, float); hi = np.array(hi, float)
            if np.isfinite(lo).any(): ax.errorbar(xx, yy, yerr=[np.nan_to_num(yy - lo), np.nan_to_num(hi - yy)], fmt="none", ecolor=ACC[gs], elinewidth=1.2, capsize=2.5)
            if len(xx) > 1: ax.plot(xx, yy, "-", color=ACC[gs], lw=1.1, alpha=0.7)
            ax.plot(xx, yy, "o", color=ACC[gs], ms=5, label=LAB[gs], zorder=3)
        # per cell type points (all shared genes), class and subclass levels
        for i in (2, 3):
            if ds != "Endresz_meta" and i == 2 or (ds == "Endresz_meta" and i == 3) or (ds != "Endresz_meta" and i == 3):
                q = d[(d.gene_set == "all") & (d.level == LEVELS[i]) & (~d.cell_type.str.startswith("MEDIAN"))]
                if len(q):
                    jit = np.random.RandomState(7).uniform(-0.07, 0.07, len(q))
                    ax.plot(i + off["all"] - 0.32 + jit, q.rho.values, ".", color=GREY, ms=3.5, alpha=0.8, zorder=2)
                    for r in q.itertuples(): rows_out.append(dict(panel=ds + ("|" + pop if pop else ""), axis=axis, level=LEVELS[i], gene_set="all", x=np.nan, rho=r.rho, lo=np.nan, hi=np.nan, cell_type=r.cell_type))
        ax.axhline(0, color=SUBTLE, lw=0.8)
        ax.set_xticks(range(4)); ax.set_xticklabels([XLAB[l] for l in LEVELS], fontsize=7.5); ax.set_xlim(-0.6, 3.6)
        ax.set_title(titles[ds + "|" + (pop or "")], loc="left", color=TEXT)
        for a in list(ax.lines) + list(ax.collections) + list(ax.texts) + list(ax.patches): a.set_clip_on(False)
    axs[0].set_ylabel("Spearman correlation across genes\nwith the bulk %s axis" % ("pH" if axis == "pH" else "agonal (Hardy t)"), color=TEXT)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, frameon=False, fontsize=8, bbox_to_anchor=(0.5, 0.045))
    fig.text(0.01, 0.005, "Source: ph_levels_2026-10-18 (out/*_levels_rho.tsv; scripts/build_figure.py). Bars: rho +/- 1.96 donor-bootstrap SD on the Fisher scale (Endresz: no interval). Grey dots: individual cell types, all-genes set.", fontsize=7, color=SUBTLE, ha="left")
    fig.tight_layout(rect=(0, 0.09, 1, 1))
    for ext in ("svg", "png"): fig.savefig(os.path.join(RUN, "fig_levels_%s.%s" % (tag, ext)), dpi=200)
    pd.DataFrame(rows_out).to_csv(os.path.join(RUN, "fig_levels_%s_data.tsv" % tag), sep="\t", index=False)
    plt.close(fig)

if __name__ == "__main__": main()
