"""02_membership.py: 5-or-more-lists membership ('published, brain', study tier) recomputed from data/lists rows, independent of recurrence.study_tier.
obs(k, g) = g is a symbol contributed by any row of study k's primary lists under the alternative, and g in M[k] (universe), as in recurrence.study_tier (Lr = lists & M).
Writes out/membership_5plus_summary.tsv, out/membership_5plus_genes.tsv, out/membership_baseline_check.tsv"""
import datetime
from vlib import *
import recurrence as R
sel, M, lists, units = R.scenario_selection("published, brain")
assert sha(CP + "/recurrent_sets.tsv") == FROZEN_SHA
frozen = {}
for r in tsv_rows(CP + "/recurrent_sets.tsv"):
    if r["set"] == "study_brain_5plus": frozen[r["gene"]] = int(r["n_lists_observed"])
rg = {r["gene"]: r for r in tsv_rows(CP + "/analyses/recurrence_genes.tsv") if r["scenario"] == "published, brain"}
obs = {}; stud_genes = {}
for alt in ALTS:
    o = collections.defaultdict(set)
    for k in sel:
        gs = set()
        for n in units["studies"][k]["primary"]:
            for r in list_rows(n): gs |= row_genes(r, alt)
        for g in gs & M[k]: o[g].add(k)
    obs[alt] = o
A = obs["A_current"]
mism = [g for g in set(rg) | set(A) if len(A.get(g, ())) != int(rg[g]["n_lists_observed"] if g in rg else 0)]
f5_A = {g for g, s in A.items() if len(s) >= 5}
T = datetime.date.today().isoformat()
def hdr(fh):
    fh.write(f"# written {T} by scripts/02_membership.py\n")
    fh.write(f"# frozen set: recurrent_sets.tsv sha256 {sha(CP+'/recurrent_sets.tsv')} (copy of analyses/confounds/recurrent_sets.tsv), set study_brain_5plus, {len(frozen)} genes; recurrence_genes.tsv sha256 {sha(CP+'/analyses/recurrence_genes.tsv')}\n")
    fh.write(f"# selection: recurrence.scenario_selection('published, brain') on the copy ({len(sel)} studies; primary_lists.json sha256 {sha(CP+'/analyses/primary_lists.json')}, universes.json sha256 {sha(CP+'/analyses/universes.json')}); lists: work/sztr/data/lists/*.tsv; pipeline HGNC (candidates for previous-dated names) data/reference/hgnc_complete_set.txt sha256 {sha(CP+'/data/reference/hgnc_complete_set.txt')}\n")
with open(OUT + "/membership_baseline_check.tsv", "w") as fh:
    hdr(fh)
    fh.write("# alternative A_current recomputed from rows against the frozen outputs\n")
    w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(["check", "value"])
    w.writerow(["genes in recurrence_genes.tsv (published, brain)", len(rg)]); w.writerow(["genes observed >=1 list, recomputed", len(A)])
    w.writerow(["genes whose n_lists_observed differs from recurrence_genes.tsv", len(mism)])
    w.writerow(["recomputed >=5 set size", len(f5_A)]); w.writerow(["frozen study_brain_5plus size", len(frozen)])
    w.writerow(["symmetric difference recomputed vs frozen", len(f5_A ^ set(frozen))])
    nobs = sum(len(s) for s in A.values()); w.writerow(["observations (study, gene) under A", nobs])
print("baseline mismatches", len(mism), "A>=5", len(f5_A), "frozen", len(frozen), "symdiff", len(f5_A ^ set(frozen)))
rows = []
for alt in ALTS:
    o = obs[alt]; cnt = {g: len(s) for g, s in o.items()}
    S = {g for g, c in cnt.items() if c >= 5}
    lost = sorted(set(frozen) - S); gained = sorted(S - set(frozen))
    # exposure: frozen members with at least one (study, gene) observation under A absent under this alternative
    expo = [g for g in frozen if A[g] - o.get(g, set())]
    lost_obs = sum(len(A[g] - o.get(g, set())) for g in A)
    gain_obs = sum(len(o[g] - A.get(g, set())) for g in o)
    mem_gain = [g for g in frozen if o.get(g, set()) - A[g]]
    # genes gaining at least one observation, not members
    nm_gain = [g for g in o if g not in frozen and o[g] - A.get(g, set())]
    rows.append([alt, len(S), len(lost), len(gained), len(expo), lost_obs, gain_obs, len(mem_gain), len(nm_gain), ";".join(lost[:40]), ";".join(gained[:40])])
    print(alt, len(S), "lost", len(lost), "gained", len(gained), "exposed", len(expo), "obs lost", lost_obs, "obs gained", gain_obs, "members gaining", len(mem_gain), "nonmember gaining", len(nm_gain))
with open(OUT + "/membership_5plus_summary.tsv", "w") as fh:
    hdr(fh)
    fh.write("# alternatives: A_current = the pipeline's symbol column; B1 = rows whose map_route is previous, previous;tag, alias, cleaned->previous, mirna->alias or previous-dated removed (old name left as printed, so not an approved symbol); B2 = B1 plus excel-date, loc->ncbi-gene-info, resolver rows removed; B2p = B2 plus override->approved rows removed (the pilot's 'affected' rule); C = previous-dated rows removed (names resolved from an ambiguous name by the symbol-change date); D0 = A plus every candidate gene of rows the pipeline dropped as ambiguous; D = D0 plus the other candidates (previous or alias of the printed name, pipeline HGNC) of previous-dated rows\n")
    fh.write("# columns: n_ge5 = genes in >=5 study lists; lost = frozen members below 5; gained = non-members at >=5; exposed = frozen members with >=1 observation (study, gene) present under A and absent under the alternative; obs_lost/obs_gained = (study, gene) observations removed/added over all genes\n")
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["alternative", "n_ge5", "lost", "gained", "exposed_members", "obs_lost", "obs_gained", "members_gaining_obs", "nonmembers_gaining_obs", "lost_genes_first40", "gained_genes_first40"]); w.writerows(rows)
with open(OUT + "/membership_5plus_genes.tsv", "w") as fh:
    hdr(fh)
    fh.write("# genes that are frozen members or reach >=5 under any alternative, or differ in count; columns n_<alt> = number of study lists observing the gene\n")
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    gs = sorted({g for g in frozen} | {g for alt in ALTS for g, s in obs[alt].items() if len(s) >= 5})
    w.writerow(["gene", "frozen_member", "frozen_n_lists_observed"] + ["n_" + a for a in ALTS])
    for g in gs:
        if all(len(obs[a].get(g, ())) == len(A.get(g, ())) for a in ALTS) and g in frozen: continue
        w.writerow([g, "yes" if g in frozen else "no", frozen.get(g, "")] + [len(obs[a].get(g, ())) for a in ALTS])
