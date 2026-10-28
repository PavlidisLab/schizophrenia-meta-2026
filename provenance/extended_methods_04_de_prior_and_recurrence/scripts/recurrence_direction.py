"""Direction agreement of the recurrent genes: when a gene recurs across studies, do the studies agree on up or down?

Paul, 2026-09-29: "Did we talk about agreement of directionality? We need to be clear when we talk about the 44 genes
for example. I think the directions mostly agree." Reads recurrence.py's per-gene output (analyses/recurrence_genes.tsv,
columns up, down, both, other_direction, sign_test_p: counts of the gene's lists by the direction they report, and the
exact two-sided sign test of up vs down) for scenario 'published, brain', study tier, and summarizes it for the gene sets
used in the paper: q_bh < 0.1 under the uniform null (the 44), q_expr_bh < 0.1, genes in 3 or more studies and genes in 5 or
more studies ("5+ studies", Paul 2026-09-29: "We should analyze the 388 too"). Lists
that report the gene without a direction ('other_direction') or in both directions ('both') are counted apart.

    python3 analyses/recurrent_gene_direction/recurrence_direction.py -> direction_genes.tsv, direction_summary.tsv,
                                                                    recurrence_direction_report.md, recurrence_direction_meta.json
"""
import csv, datetime, hashlib, json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
AN = os.path.dirname(HERE)
SRC = os.path.join(AN, "recurrence_genes.tsv")
SCEN = "published, brain"
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
num = lambda x: int(x) if x not in ("", None) else 0


def exp_majority(n):
    """Expected hits in the majority direction for n lists whose directions are independent fair coin flips."""
    return sum(math.comb(n, k) * max(k, n - k) for k in range(n + 1)) / 2 ** n


def main():
    rows = [r for r in csv.DictReader((l for l in open(SRC) if not l.startswith("#")), delimiter="\t") if r["scenario"] == SCEN]
    sets = {"study tier q<0.1 (uniform null)": [r for r in rows if float(r["q_bh"]) < 0.1],
            "expression-matched q<0.1": [r for r in rows if r.get("q_expr_bh") and float(r["q_expr_bh"]) < 0.1],
            "3 or more studies": [r for r in rows if num(r["n_lists_observed"]) >= 3],
            "5+ studies": [r for r in rows if num(r["n_lists_observed"]) >= 5]}
    genes, summ = [], []
    for name, S in sets.items():
        hits = maj = unan = split = sig = 0
        emaj = 0.0
        for r in S:
            u, d = num(r["up"]), num(r["down"])
            m = max(u, d)
            hits += u + d; maj += m; emaj += exp_majority(u + d)
            unan += (u + d > 0 and min(u, d) == 0)
            split += (u + d > 0 and min(u, d) / (u + d) >= 0.4)
            sig += bool(r["sign_test_p"]) and float(r["sign_test_p"]) < 0.05
            if name.startswith("study tier"):
                genes.append({"gene": r["gene"], "n_studies": r["n_lists_observed"], "up": u, "down": d, "both": num(r["both"]),
                              "no_direction": num(r["other_direction"]), "majority": "up" if u > d else "down" if d > u else "tie",
                              "share_majority": round(m / (u + d), 3) if u + d else "", "sign_test_p": r["sign_test_p"]})
        summ.append({"gene_set": name, "genes": len(S), "list_hits_with_direction": hits, "in_majority_direction": maj,
                     "share_majority": round(maj / hits, 3) if hits else "",
                     "share_majority_if_random": round(emaj / hits, 3) if hits else "", "genes_unanimous": unan,
                     "genes_split_40pct_or_more": split, "genes_sign_test_p05": sig})
    genes.sort(key=lambda g: (-int(g["n_studies"]), g["gene"]))
    stamp = datetime.date.today().isoformat()
    hdr = f"# script: analyses/recurrent_gene_direction/recurrence_direction.py, run {stamp}; source: recurrence_genes.tsv sha256 {sha(SRC)}; scenario '{SCEN}', study tier\n"
    for nm, rs in (("direction_genes.tsv", genes), ("direction_summary.tsv", summ)):
        with open(os.path.join(HERE, nm), "w", newline="") as fh:
            fh.write(hdr); w = csv.DictWriter(fh, fieldnames=list(rs[0]), delimiter="\t"); w.writeheader(); w.writerows(rs)
    json.dump({"generated": stamp, "scenario": SCEN, "source": {"recurrence_genes.tsv": sha(SRC)}},
              open(os.path.join(HERE, "recurrence_direction_meta.json"), "w"), indent=1)
    s44 = summ[0]
    s5 = next(r for r in summ if r["gene_set"] == "5+ studies")
    split = [g for g in genes if g["up"] + g["down"] and min(g["up"], g["down"]) / (g["up"] + g["down"]) >= 0.4]
    up = [g["gene"] for g in genes if g["majority"] == "up" and g["share_majority"] == 1]
    dn = [g["gene"] for g in genes if g["majority"] == "down" and g["share_majority"] == 1]
    split_txt = ", ".join("%s (%d up / %d down)" % (g["gene"], g["up"], g["down"]) for g in split)
    L = ["---", "type: detail", f"modified: {stamp}", "---", "", "> [!note] Paul", ">", "",
         "# Do the recurrent genes agree in direction?", "", "## Executive summary", "",
         f"- **Question.** Paul, 2026-09-29: \"I think the directions mostly agree.\" For the genes that recur across the postmortem brain studies (scenario '{SCEN}', study tier), do the studies report the same direction?",
         f"- **Answer. Yes, mostly.** For the {s44['genes']} genes at q < 0.1, {s44['in_majority_direction']} of {s44['list_hits_with_direction']} list hits with a direction ({s44['share_majority']:.0%}) are in the gene's majority direction, against {s44['share_majority_if_random']:.0%} expected if each list's direction were a coin flip; {s44['genes_unanimous']} genes are unanimous and {s44['genes_split_40pct_or_more']} are split 60:40 or closer.",
         f"- **Sign test.** {s44['genes_sign_test_p05']} of {s44['genes']} genes reach p < 0.05; with 7 to 11 lists per gene, 6 of 7 in one direction does not.",
         f"- **Unanimous up:** {', '.join(up)}. **Unanimous down:** {', '.join(dn)}. **Split:** {split_txt}.",
         f"- **5+ studies ({s5['genes']} genes in 5 or more of the lists).** {s5['in_majority_direction']} of {s5['list_hits_with_direction']} list hits with a direction ({s5['share_majority']:.0%}) are in the gene's majority direction, against {s5['share_majority_if_random']:.0%} if directions were coin flips; {s5['genes_unanimous']} genes unanimous, {s5['genes_split_40pct_or_more']} split 60:40 or closer, {s5['genes_sign_test_p05']} at sign test p < 0.05.",
         "- **Caveat.** A shared direction is what a confound with a consistent direction also produces: in controls these genes move with low brain pH the way they move in SZ, and pH is lower in SZ in every dataset tested (`../confounds/confounds_report.md`).", "",
         "| gene set | genes | list hits with a direction | in majority direction | if directions were random | unanimous genes | split ≥ 40% | sign test p < 0.05 |", "|---|---:|---:|---:|---:|---:|---:|---:|"]
    L += [f"| {r['gene_set']} | {r['genes']} | {r['list_hits_with_direction']} | {r['in_majority_direction']} ({r['share_majority']:.0%}) | {r['share_majority_if_random']:.0%} | {r['genes_unanimous']} | {r['genes_split_40pct_or_more']} | {r['genes_sign_test_p05']} |" for r in summ]
    L += ["", "Per gene: `direction_genes.tsv`. Source: `../recurrence_genes.tsv` (sha256 in the TSV headers)."]
    open(os.path.join(HERE, "recurrence_direction_report.md"), "w").write("\n".join(L) + "\n")
    for r in summ:
        print(r)


if __name__ == "__main__":
    main()
