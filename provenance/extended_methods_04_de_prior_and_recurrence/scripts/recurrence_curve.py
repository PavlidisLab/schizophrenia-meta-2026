"""Genes found in at least N studies, observed against chance, for the primary analysis (recurrence.py scenario
'published, brain': one published list per first-look postmortem brain study, each restricted to the genes it
tested; study tier).

Paul, 2026-09-29: "I need the number that goes with our primary analysis. A plot that shows how many genes found in at
least N studies could help." Expected counts use recurrence.study_tier's own probability matrices, so the three
chance models are the ones the recurrence analysis tests against: uniform over tested genes, DE-prior weighted and
expression-matched. Expected genes with at least N lists = sum over genes of P(X_g >= N) (Poisson-binomial tails).

share_chance_* = expected / observed at each N: a set-level estimate of the share of the genes in at least N studies that
the null expects there by chance. It does not say which genes. Above 1, the null expects more genes than are observed.

sd_* (2026-10-01): the total "genes with >= N lists" is a sum over genes of indicator variables 1{X_g >= N}; the model
already treats each gene's own list-membership draws as independent Bernoulli trials (pb_matrix's per-list probability
P[g, k], not a fixed-size draw without replacement), so under that same assumption the indicators across genes are
independent too and Var(sum) = sum over genes of P(X_g >= N) * (1 - P(X_g >= N)); sd_* is its square root. z_* =
(observed - expected) / sd: a Normal-approximation z-score for whether the observed count is distinguishable from
that null, not a per-gene test (recurrence_genes.tsv's q_bh is that). The figure shades expected +/- 1.96 sd.

"5+ studies" (Paul, 2026-09-29: "I see the bulk test is just an estimate of how many there are. That's going to be a
criticism. We should analyze the 388 too"): the genes in 5 or more of the lists. For them: the chance shares at N = 5,
their expression deciles (the expression-matched null's bins, recurrence.expression_weights: LIBD DLPFC mean log CPM,
hippocampus where DLPFC has none, genes in neither in the lowest bin) against every measured gene, and their overlap with
the study-tier genes at q_bh < 0.1 and q_expr_bh < 0.1 (recurrence_genes.tsv).

The DE prior is recurrence.py's DE_PRIOR_SOURCE (our own 2-fold build). Under DE_PRIOR_SOURCE=global
(external files built at |log2FC| >= 2, 4-fold, where 2-fold was meant; Neera, 2026-09-29) the DE-prior columns are
marked provisional.

    python3 analyses/gene_coverage/recurrence_curve.py -> recurrence_curve.tsv, fig_recurrence_curve.svg/.png, recurrence_curve_meta.json,
                                                          five_plus_genes.tsv, five_plus_deciles.tsv, recurrence_curve_report.md
"""
import csv, datetime, hashlib, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
AN = os.path.dirname(HERE)
sys.path.insert(0, AN)
import recurrence as R  # noqa: E402

SCEN = os.environ.get("RECURRENCE_CURVE_SCENARIO", "published, brain")
# output files stay unsuffixed for the default scenario (published, brain, used by report.md); any other scenario
# (e.g. RECURRENCE_CURVE_SCENARIO="published, brain, FDR lists", stringent: input lists restricted to studies whose
# own list applied formal FDR control, 33 of the 60 studies) writes alongside under a slug suffix, so it never
# clobbers the default outputs
SUFFIX = "" if SCEN == "published, brain" else "_" + SCEN.split(", brain, ")[-1].replace(" ", "_")
FIVE = 5
PROVISIONAL = ("" if R.DE_PRIOR_SOURCE == "ours" else
               "provisional: the DE-prior files were built at |log2FC| >= 2 (4-fold) where 2-fold was meant (Neera, 2026-09-29)")
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]


def main():
    units = R.load_units()
    uni = json.load(open(os.path.join(AN, "universes.json")))["lists"]
    lists, M, nk, _, _ = R.load_lists(units, uni)
    reps = R.representatives(units, dict(R.SCENARIOS)[SCEN], lambda k: bool(lists[k]))
    sel = [k for k in sorted(reps, key=lambda k: (0, int(k), "") if k.isdigit() else (1, 0, k)) if k in M]
    st = R.study_tier(sel, M, nk, lists, R.de_prior_weights(), R.expression_levels())
    obs, tail, tailP, tailE = st["obs"], st["tail"], st["tailP"], st["tailE"]
    top = int(obs.max())
    rows = []
    for n in range(1, top + 1):
        o = int((obs >= n).sum())
        row = {"n_studies_at_least": n, "observed_genes": o}
        for key, t in (("uniform", tail), ("de_prior", tailP), ("expression", tailE)):
            p = t[:, n]
            mean = float(p.sum())
            sd = float(np.sqrt((p * (1 - p)).sum()))
            row[f"expected_{key}"] = round(mean, 1)
            row[f"sd_{key}"] = round(sd, 2)
            row[f"share_chance_{key}"] = round(mean / o, 3) if o else ""
            row[f"z_{key}"] = round((o - mean) / sd, 2) if sd > 0 else ""
        rows.append(row)
    genes = st["genes"]
    best = sorted(((int(obs[i]), genes[i]) for i in range(len(genes)) if obs[i] == top), key=lambda x: x[1])
    stamp = datetime.date.today().isoformat()
    ins = {f: sha(os.path.join(AN, f)) for f in ("primary_lists.json", "universes.json")} | {"DE prior": sha(R.DE_PRIOR), "expression": sha(R.EXPR_FILE)}
    with open(os.path.join(HERE, f"recurrence_curve{SUFFIX}.tsv"), "w", newline="") as fh:
        fh.write(f"# script: analyses/gene_coverage/recurrence_curve.py, run {stamp}; scenario '{SCEN}', study tier; {len(sel)} lists; "
                 f"genes measured by at least one list: {len(genes)}; inputs sha256 {ins}\n"
                 f"# share_chance_* = expected / observed (set-level estimate, not per gene); DE prior: {os.path.relpath(R.DE_PRIOR, AN)}"
                 + (f"; expected_de_prior and share_chance_de_prior are {PROVISIONAL}" if PROVISIONAL else "") + "\n")
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    five = five_plus(st, rows, stamp, ins, len(sel))
    json.dump({"generated": stamp, "scenario": SCEN, "n_lists": len(sel), "genes_measured": len(genes), "max_studies": top,
               "genes_at_max": [g for _, g in best], "inputs": ins, "de_prior_file": os.path.relpath(R.DE_PRIOR, AN),
               "de_prior_provisional": PROVISIONAL, "five_plus": five},
              open(os.path.join(HERE, f"recurrence_curve_meta{SUFFIX}.json"), "w"), indent=1)
    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": ["Helvetica", "Arial", "sans-serif"], "svg.fonttype": "none", "font.size": 9})
    x = [r["n_studies_at_least"] for r in rows]
    fig, ax = plt.subplots(figsize=(6.2, 4.2), facecolor="white")
    for key, lab, col, ls in (("observed_genes", "observed", "#1f2937", "-"), ("expected_uniform", "chance: uniform over tested genes", "#2563eb", "--"),
                              ("expected_de_prior", "chance: DE-prior weighted", "#f59e0b", "--"), ("expected_expression", "chance: expression-matched", "#10b981", "--")):
        y = [r[key] if r[key] >= 0.1 else float("nan") for r in rows]  # expectations below 0.1 genes are not drawn
        ax.plot(x, y, ls, color=col, lw=2 if key == "observed_genes" else 1.4, marker="o" if key == "observed_genes" else None, ms=3.5, label=lab)
        if key != "observed_genes":
            null = key.split("_", 1)[1]
            lo = [max(r[key] - 1.96 * r[f"sd_{null}"], 0.05) if r[key] >= 0.1 else float("nan") for r in rows]
            hi = [r[key] + 1.96 * r[f"sd_{null}"] if r[key] >= 0.1 else float("nan") for r in rows]
            ax.fill_between(x, lo, hi, color=col, alpha=0.15, lw=0)
    ax.set_yscale("log"); ax.set_xticks(x)
    ax.set_xlabel("found DE in at least N studies", color="#1f2937"); ax.set_ylabel("genes (log scale)", color="#1f2937")
    ax.tick_params(colors="#6b7280"); ax.grid(axis="y", color="#e5e7eb", lw=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title(f"Genes found in at least N of {len(sel)} postmortem brain studies", loc="left", fontweight="normal")
    ax.annotate(", ".join(g for _, g in best) + f": {top} studies", (top, rows[-1]["observed_genes"]), textcoords="offset points",
                xytext=(-6, -16), ha="right", fontsize=8, color="#1f2937")
    ax.legend(frameon=False, fontsize=8)
    fig.text(0.01, 0.005, f"Source: analyses/gene_coverage/recurrence_curve.tsv (recurrence.py study tier, scenario '{SCEN}')", fontsize=7, color="#6b7280")
    for a in [ax] + list(ax.lines) + list(ax.texts):
        a.set_clip_on(False)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(os.path.join(HERE, f"fig_recurrence_curve{SUFFIX}.svg")); fig.savefig(os.path.join(HERE, f"fig_recurrence_curve{SUFFIX}.png"), dpi=150)
    for r in rows:
        print(r)
    print("most studies:", best)


def five_plus(st, rows, stamp, ins, n_lists):
    """The 5+ studies set: per-gene table, expression deciles against every measured gene, overlap with the 44 and the 20,
    and the report. Returns the summary recorded in recurrence_curve_meta.json."""
    genes, obs, ebins = st["genes"], st["obs"], st["ebins"]
    rg_path = os.path.join(AN, "recurrence_genes.tsv")
    rg = {r["gene"]: r for r in csv.DictReader((l for l in open(rg_path) if not l.startswith("#")), delimiter="\t") if r["scenario"] == SCEN}
    five = [i for i in range(len(genes)) if obs[i] >= FIVE]
    if {genes[i] for i in five} != {g for g, r in rg.items() if int(r["n_lists_observed"] or 0) >= FIVE}:
        raise SystemExit("recurrence_genes.tsv n_lists_observed does not match the study tier recomputed here; rerun recurrence.py")
    s44 = {g for g, r in rg.items() if float(r["q_bh"]) < 0.1}
    s20 = {g for g, r in rg.items() if r.get("q_expr_bh") and float(r["q_expr_bh"]) < 0.1}
    expr = R.expression_levels()
    per = []
    for i in sorted(five, key=lambda i: (-int(obs[i]), genes[i])):
        g = genes[i]
        per.append({"gene": g, "n_studies": int(obs[i]), "expr_decile": int(ebins[i]) + 1,
                    "libd_mean_logcpm": round(expr[g], 3) if g in expr else "", "in_study_tier_q01": g in s44,
                    "in_expression_matched_q01": g in s20, "q_bh": rg[g]["q_bh"], "q_expr_bh": rg[g]["q_expr_bh"],
                    "up": rg[g]["up"], "down": rg[g]["down"]})
    dec = []
    n_all, n5 = len(genes), len(five)
    for b in range(int(ebins.max()) + 1):
        idx = [i for i in range(n_all) if ebins[i] == b]
        f5 = [i for i in five if ebins[i] == b]
        dec.append({"expr_decile": b + 1, "genes_measured": len(idx), "genes_measured_without_libd_value": sum(1 for i in idx if genes[i] not in expr),
                    "genes_5plus": len(f5), "genes_5plus_without_libd_value": sum(1 for i in f5 if genes[i] not in expr),
                    "share_of_measured": round(len(idx) / n_all, 4), "share_of_5plus": round(len(f5) / n5, 4) if n5 else "",
                    "ratio": round((len(f5) / n5) / (len(idx) / n_all), 2) if n5 and idx else ""})
    head = (f"# script: analyses/gene_coverage/recurrence_curve.py, run {stamp}; scenario '{SCEN}', study tier; {n_lists} lists; "
            f"5+ studies = genes in {FIVE} or more lists; recurrence_genes.tsv sha256 {sha(rg_path)}; inputs sha256 {ins}\n")
    for name, rs in ((f"five_plus_genes{SUFFIX}.tsv", per), (f"five_plus_deciles{SUFFIX}.tsv", dec)):
        with open(os.path.join(HERE, name), "w", newline="") as fh:
            fh.write(head)
            w = csv.DictWriter(fh, fieldnames=list(rs[0]), delimiter="\t"); w.writeheader(); w.writerows(rs)
    r5 = next(r for r in rows if r["n_studies_at_least"] == FIVE)
    in5 = {genes[i] for i in five}
    summ = {"n_genes": n5, "observed_at_5": r5["observed_genes"],
            "expected_at_5": {k: r5[f"expected_{k}"] for k in ("uniform", "de_prior", "expression")},
            "share_chance_at_5": {k: r5[f"share_chance_{k}"] for k in ("uniform", "de_prior", "expression")},
            "n_study_tier_q01": len(s44), "n_study_tier_q01_in_5plus": len(s44 & in5),
            "n_expression_matched_q01": len(s20), "n_expression_matched_q01_in_5plus": len(s20 & in5),
            "n_5plus_in_both": len(s44 & s20 & in5), "n_5plus_in_neither": len(in5 - s44 - s20),
            "n_5plus_without_libd_value": sum(1 for g in in5 if g not in expr),
            "recurrence_genes_sha256": sha(rg_path)}
    # report
    L = ["---", "type: detail", f"modified: {stamp}", "---", "", "> [!note] Paul", ">", "",
         "# Genes in at least N studies against chance, and the 5+ studies set", "",
         f"Built by `recurrence_curve.py` ({stamp}); edit the script, not this file. Scenario '{SCEN}' of `recurrence.py`, study tier: "
         f"{n_lists} lists, {n_all} genes measured by at least one. Expected counts are sums of Poisson-binomial tails under the three nulls of "
         "`recurrence.py` (uniform over tested genes, DE-prior weighted, expression-matched).", "",
         (f"**DE-prior columns are {PROVISIONAL}.**" if PROVISIONAL else f"DE prior: `{os.path.relpath(R.DE_PRIOR, AN)}`."), "",
         "## Genes in at least N studies", "",
         "| N | observed | expected: uniform | DE prior | expression-matched | share by chance: uniform | DE prior | expression-matched |",
         "|---:|---:|---:|---:|---:|---:|---:|---:|"]
    L += [f"| {r['n_studies_at_least']} | {r['observed_genes']} | {r['expected_uniform']} | {r['expected_de_prior']} | {r['expected_expression']} | "
          f"{r['share_chance_uniform']} | {r['share_chance_de_prior']} | {r['share_chance_expression']} |" for r in rows]
    sc = summ["share_chance_at_5"]
    ex = summ["expected_at_5"]
    L += ["", "Share by chance = expected / observed: a set-level estimate of how much of the set a null accounts for. It does not say which "
          "genes are the chance ones, and it is not a per-gene error rate. Above 1, the null expects more genes than are observed.", "",
          f"## 5+ studies ({n5} genes in {FIVE} or more of the {n_lists} lists)", "",
          f"- **Chance share (set level).** Genes expected in {FIVE}+ studies: {ex['uniform']} under the uniform null ({sc['uniform']:.0%} of {n5}), "
          f"{ex['expression']} under the expression-matched null ({sc['expression']:.0%}), {ex['de_prior']} under the DE-prior null "
          f"({sc['de_prior']:.0%}" + ("; provisional)." if PROVISIONAL else ")."),
          f"- **Overlap.** {summ['n_study_tier_q01_in_5plus']} of the {summ['n_study_tier_q01']} study-tier genes at q < 0.1 (uniform null) and "
          f"{summ['n_expression_matched_q01_in_5plus']} of the {summ['n_expression_matched_q01']} at q < 0.1 under the expression-matched null are in the set; "
          f"{summ['n_5plus_in_both']} are in both, {summ['n_5plus_in_neither']} in neither.",
          f"- **Expression.** {summ['n_5plus_without_libd_value']} of the {n5} have no LIBD value and sit in decile 1 with the other genes that have none.", "",
          "| expression decile | measured genes (no LIBD value) | 5+ genes (no LIBD value) | share of measured | share of 5+ | ratio |",
          "|---:|---:|---:|---:|---:|---:|"]
    L += [f"| {d['expr_decile']} | {d['genes_measured']} ({d['genes_measured_without_libd_value']}) | {d['genes_5plus']} ({d['genes_5plus_without_libd_value']}) | "
          f"{d['share_of_measured']:.1%} | {d['share_of_5plus']:.1%} | {d['ratio']} |" for d in dec]
    L += ["", "Deciles are the expression-matched null's bins (`recurrence.expression_weights`: ten equal-count bins of LIBD BrainSeq DLPFC mean log CPM, "
          "hippocampus where DLPFC has none, among measured genes with a value; genes with no value join decile 1).", "",
          f"Files: `recurrence_curve{SUFFIX}.tsv`, `five_plus_genes{SUFFIX}.tsv` (per gene), `five_plus_deciles{SUFFIX}.tsv`, `recurrence_curve_meta{SUFFIX}.json`. "
          f"Inputs sha256: {ins}; recurrence_genes.tsv {sha(rg_path)}."]
    open(os.path.join(HERE, f"recurrence_curve_report{SUFFIX}.md"), "w").write("\n".join(L) + "\n")
    return summ


if __name__ == "__main__":
    main()
