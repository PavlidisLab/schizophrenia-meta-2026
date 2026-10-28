"""Per-gene drug effects for every unit (drug vs vehicle contrasts of the Gemma result sets), on human genes.

    opt/anaconda3/bin/python3 scripts/build_long.py
Rodent genes -> human through HCOP orthologues (../confounds/scripts/build_genesets.map_rodent: Ensembl id first, then
symbol; 1:many resolved by the most supporting sources, ties left unmapped); human genes through HGNC (map_human_symbol).
Per contrast and human gene the probe with the smallest p is kept. z = sign(t) * Phi^-1(1 - p/2).
Writes results/contrasts_long.tsv.gz and results/units_table.tsv (one row per unit x contrast: design, n per group).
"""
import csv, gzip, io, json, math, os, re, sys
from collections import defaultdict
from scipy.stats import norm

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GL = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HERE, "scripts"))
sys.path.insert(0, os.path.join(GL, "confounds", "scripts"))
from units import UNITS, DRUGS  # noqa: E402
import build_genesets as BG  # noqa: E402

DRUG_RE = re.compile(DRUGS, re.I)


def rs_meta(gid, rs):
    j = json.load(open(os.path.join(HERE, "data", "analyses", f"{gid}.json")))
    for a in j["analyses"]["data"]:
        for r in a["resultSets"]:
            if r["id"] == rs:
                sub = (a.get("subsetFactorValue") or {}).get("factorValue")
                vals = {}
                for ef in r["experimentalFactors"]:
                    for v in ef.get("values", []):
                        vals[v["id"]] = v.get("factorValue")
                return vals, (r.get("baselineGroup") or {}).get("id"), (r.get("baselineGroup") or {}).get("factorValue"), sub
    raise SystemExit(f"{gid}: result set {rs} not found")


def sample_counts(gid, fvids, subset):
    d = json.load(open(os.path.join(HERE, "data", "samples", f"{gid}.json")))["data"]
    n, tissues = defaultdict(int), set()
    for s in d:
        if s.get("outlier") or s.get("userFlaggedOutlier"):
            continue
        ch = s["sample"].get("characteristics", [])
        tis = [c["value"] for c in ch if (c.get("category") or "").lower() in ("tissue", "organism part", "organismpart", "region", "biosource")]
        fvs = [str(fv.get("summary") or fv.get("value") or "") for fv in s["sample"].get("fvbasicVOs", [])]
        if subset and not any(subset.lower() in t.lower() or t.lower() in subset.lower() for t in tis + fvs if t):
            continue
        tissues.update(tis)
        ids = {fv["id"] for fv in s["sample"].get("fvbasicVOs", [])}
        for f in fvids:
            if f in ids:
                n[f] += 1
    return n, sorted(tissues)


def main():
    long_rows, table = [], []
    for unit, cls, gid, gse, species, region, rs, duration, dose in UNITS:
        vals, base_id, base_lab, sub = rs_meta(gid, rs)
        with gzip.open(os.path.join(HERE, "data", "resultsets", f"{rs}.tsv.gz"), "rt") as fh:
            rows = list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))
        cons = sorted({int(k.split("_")[1]) for k in rows[0] if k.startswith("contrast_") and k.endswith("_tstat")})
        drug_cons = [c for c in cons if DRUG_RE.search(vals.get(c) or "")]
        n, tissues = sample_counts(gid, drug_cons + [base_id], sub)
        best = {}
        for r in rows:
            sym, ens = r.get("gene_official_symbol") or "", r.get("gene_ensembl_id") or ""
            if not sym or "|" in sym:
                continue
            if species == "human":
                h, route = BG.map_human_symbol(sym)
            else:
                h, route = BG.map_rodent(species, sym, ens.split("|")[0] if ens else None)
            if not h:
                continue
            for c in drug_cons:
                t, p, fc = r.get(f"contrast_{c}_tstat"), r.get(f"contrast_{c}_pvalue"), r.get(f"contrast_{c}_log2fc")
                try:
                    t, p = float(t), float(p)
                except (TypeError, ValueError):
                    continue
                if not (0 < p <= 1) or math.isnan(t):
                    continue
                k = (c, h)
                if k not in best or p < best[k][1]:
                    best[k] = (t, p, fc, route, sym)
        for (c, h), (t, p, fc, route, sym) in best.items():
            z = math.copysign(norm.isf(p / 2), t)
            long_rows.append([unit, cls, gid, gse, rs, c, vals[c], h, sym, route, f"{t:.4g}", f"{p:.4g}", fc, f"{z:.4g}"])
        for c in drug_cons:
            table.append([unit, cls, gid, gse, species, region if not sub else f"{region} (Gemma subset: {sub})", "; ".join(tissues), rs, c,
                          vals[c], base_lab, n.get(c, ""), n.get(base_id, ""), duration, dose, sum(1 for (cc, _) in best if cc == c)])
        print(unit, gse, [vals[c] for c in drug_cons], "genes", len({h for (_, h) in best}))
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with gzip.open(os.path.join(HERE, "results", "contrasts_long.tsv.gz"), "wt", compresslevel=6) as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        fh.write("# medication/scripts/build_long.py: per human gene and drug contrast, the best probe of the Gemma result set\n")
        w.writerow(["unit", "class", "gemma_id", "geo", "result_set", "contrast_fv", "drug", "gene", "source_symbol", "map_route", "t", "p", "log2fc", "z"])
        w.writerows(long_rows)
    with open(os.path.join(HERE, "results", "units_table.tsv"), "w") as fh:
        fh.write("# medication/scripts/build_long.py: one row per drug contrast; n from Gemma sample annotations (outliers excluded; subset units restricted to the subset tissue)\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["unit", "class", "gemma_id", "geo", "species", "region", "sample_tissues", "result_set", "contrast_fv", "drug", "baseline", "n_drug", "n_baseline", "duration", "dose", "n_human_genes"])
        w.writerows(table)


if __name__ == "__main__":
    main()
