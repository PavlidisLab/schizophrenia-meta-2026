#!/usr/bin/env python3
"""Share of the top-200 control-pH genes' proteomic hits that point the low-pH way (2026-10-09).

Same scoring as analyses/brain_genes_in_periphery (sets.tsv ph_top200: top 200 genes by |Z| of control pH in
analyses/ph_and_confound_analyses/results/combined_controls_pH.tsv; predicted SZ direction = the one low pH gives, Z>0 -> down).
Unit = side-studies/brain_proteomics_lists/units.tsv row (union of its lists; a gene's direction = majority over the unit's lists,
ties/no direction dropped). Donor groups: units sharing a bank (units.tsv `bank`; stanley, nsw, pitt, mannheim, upv,
koopmans' banks as one string) count once per gene (majority over the group's units); units with no bank stay alone.
Null: shuffle the direction labels among the unit's own directed genes (5,000x, seed 20261009) and recount concordant
hits, so the unit's up/down balance and the set of proteins it reports (abundance bias) are held fixed.
Outputs: ph_top200_concordance_units.tsv, ph_top200_concordance_groups.tsv, ph_top200_concordance_summary.json
"""
import collections, csv, json, os, random, hashlib
H = os.path.dirname(os.path.abspath(__file__)); R = os.path.normpath(os.path.join(H, "..", ".."))
def rows(p):
    return list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t"))
pred = {}
for l in open(os.path.join(R, "analyses/brain_genes_in_periphery/sets.tsv")):
    if l.startswith("ph_top200\t"):
        f = l.rstrip("\n").split("\t"); pred[f[1]] = f[2]
assert len(pred) == 200, len(pred)
units = rows(os.path.join(H, "units.tsv"))
def unit_dirs(u):
    d = collections.defaultdict(list)
    for lid in u["lists"].split(";"):
        fp = os.path.join(R, "data/lists", lid + ".tsv")
        if not os.path.exists(fp):
            print("missing list file, skipped:", lid); continue
        for r in rows(fp):
            if r["direction"] in ("up", "down") and r["symbol"]:
                d[r["symbol"]].append(r["direction"])
    out = {}
    for g, v in d.items():
        c = collections.Counter(v)
        if c["up"] != c["down"]: out[g] = "up" if c["up"] > c["down"] else "down"
    return out
rnd = random.Random(20261009); N = 5000
res = []; udirs = {}
for u in units:
    d = unit_dirs(u); udirs[u["unit"]] = d
    hits = [g for g in d if g in pred]
    conc = sum(d[g] == pred[g] for g in hits)
    labs = list(d.values()); genes = list(d)
    ge = 0
    exp = 0.0
    if hits:
        for _ in range(N):
            rnd.shuffle(labs)
            m = dict(zip(genes, labs))
            c = sum(m[g] == pred[g] for g in hits)
            exp += c
            ge += c >= conc
    res.append({"unit": u["label"], "bank": u["bank"], "class": u["class"], "directed_genes": len(d), "ph_top200_hits": len(hits),
                "concordant": conc, "share": round(conc / len(hits), 3) if hits else "", "expected_share": round(exp / N / len(hits), 3) if hits else "",
                "p_perm_ge": round((ge + 1) / (N + 1), 4) if hits else ""})
def w(name, rr, head):
    with open(os.path.join(H, name), "w") as f:
        f.write("# built 2026-10-09 by side-studies/brain_proteomics_lists/ph_top200_concordance.py; " + head + "\n")
        wr = csv.DictWriter(f, fieldnames=list(rr[0]), delimiter="\t"); wr.writeheader(); wr.writerows(rr)
w("ph_top200_concordance_units.tsv", res, "per unit; share = concordant / ph_top200 hits with a direction; expected_share and p from 5,000 direction-label shuffles within the unit")
# donor groups
grp = collections.defaultdict(list)
for u in units:
    grp[u["bank"] if u["bank"] else "solo:" + u["unit"]].append(u)
grows = []; tot_h = tot_c = 0; tot_e = 0.0; ge_all = 0
perm_tot = [0] * N
gdat = {}
for b, us in grp.items():
    dd = [udirs[u["unit"]] for u in us]
    genes = set().union(*[set(x) for x in dd])
    call = {}
    for g in genes:
        c = collections.Counter(x[g] for x in dd if g in x)
        if c["up"] != c["down"]: call[g] = "up" if c["up"] > c["down"] else "down"
    hits = [g for g in call if g in pred]; conc = sum(call[g] == pred[g] for g in hits)
    labs = list(call.values()); gl = list(call); exp = 0.0
    for i in range(N):
        rnd.shuffle(labs); m = dict(zip(gl, labs))
        c = sum(m[g] == pred[g] for g in hits); exp += c; perm_tot[i] += c
    gdat[b] = (len(hits), conc)
    tot_h += len(hits); tot_c += conc; tot_e += exp / N
    grows.append({"group": b, "units": ";".join(u["label"] for u in us), "genes_called": len(call), "ph_top200_hits": len(hits), "concordant": conc,
                  "share": round(conc / len(hits), 3) if hits else "", "expected_share": round(exp / N / len(hits), 3) if hits else ""})
pooled_p = (sum(c >= tot_c for c in perm_tot) + 1) / (N + 1)
w("ph_top200_concordance_groups.tsv", grows, "per donor group (units sharing a bank pooled; gene direction = majority over the group's units)")
from math import comb
def binom2(k, n, p=0.5):
    pk = lambda i: comb(n, i) * p**i * (1 - p)**(n - i)
    o = pk(k); return min(1, sum(pk(i) for i in range(n + 1) if pk(i) <= o + 1e-12))
S = {"groups": len(grp), "pooled_hits": tot_h, "pooled_concordant": tot_c, "pooled_share": round(tot_c / tot_h, 3),
     "pooled_expected_share": round(tot_e / tot_h, 3), "p_perm_pooled_ge": round(pooled_p, 4), "binom_p_vs_0.5": binom2(tot_c, tot_h),
     "groups_with_hits": sum(1 for r in grows if r["ph_top200_hits"]), "seed": 20261009, "n_perm": N,
     "input_sha256_sets": hashlib.sha256(open(os.path.join(R, "analyses/brain_genes_in_periphery/sets.tsv"), "rb").read()).hexdigest()[:12]}
json.dump(S, open(os.path.join(H, "ph_top200_concordance_summary.json"), "w"), indent=1)
print(json.dumps(S, indent=1))
for r in grows: print(r["group"][:30].ljust(30), r["ph_top200_hits"], r["concordant"], r["share"], r["expected_share"])
