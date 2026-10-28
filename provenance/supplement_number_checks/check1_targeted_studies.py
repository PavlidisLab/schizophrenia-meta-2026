#!/usr/bin/env python3
"""Check 1: targeted studies, overlap with the 665 profiling corpus, and the total.
Inputs (agent-mail/sz-transcriptome/corpus_additions_2026-10-19): regen/lineage_overlay.json (the 801-node lineage of the 669 run), regen/h4_papers/h4_papers.tsv (delivered h4_papers.tsv 1,506 rows + 47 new rows), and the delivered provenance/extended_methods_01_literature_search/h4_papers.tsv.
Corpus (665) = nodes not category 2 (targeted/one-gene), not the Merikangas landmark ('review'), not category 1m (the 4 panels of EL07.06) and not 'Cited source (not in corpus)'; matched on PMID (4 corpus nodes have no PMID).
Targeted set = PMIDs in h4_papers (1,553) union the lineage category-2 nodes (63). Usage: check1_targeted_studies.py CORPUS_ADDITIONS_DIR OLD_H4_TSV OUT_TSV"""
import sys, json, hashlib, pandas as pd
D, OLD, OUT = sys.argv[1:4]; sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
N = json.load(open(f"{D}/regen/lineage_overlay.json"))["nodes"]; cat = lambda n: str(n.get("category"))
corp = [n for n in N if cat(n) not in ("2", "review", "1m") and not str(n.get("source")).startswith("Cited source")]
cp = {str(n["pmid"]) for n in corp if n.get("pmid")}; c2 = {str(n["pmid"]) for n in N if cat(n) == "2" and n.get("pmid")}; m1 = {str(n["pmid"]) for n in N if cat(n) == "1m"}
H = set(pd.read_csv(f"{D}/regen/h4_papers/h4_papers.tsv", sep="\t", comment="#").pmid.astype(str)); H0 = set(pd.read_csv(OLD, sep="\t", comment="#").pmid.astype(str))
T = H | c2
rows = [("corpus studies (nodes)", len(corp)), ("corpus studies with a PMID", len(cp)), ("h4_papers delivered (provenance file)", len(H0)), ("h4_papers regenerated (+47 additions)", len(H)), ("lineage category-2 nodes (one-gene studies and panels moved out of the corpus before 2026-10-19)", len(c2)),
        ("category-2 nodes not in h4_papers (regenerated)", len(c2 - H)), ("delivered h4_papers papers that are category-2 nodes (funnel E2 note: 35)", len(H0 & c2)), ("targeted set = h4 (1,553) union category 2", len(T)), ("the 4 panels of EL07.06 (category 1m) inside that set", len(m1 & T)),
        ("overlap of the targeted set with the 665 corpus (PMID match)", len(T & cp)), ("total = targeted + 665 - overlap", len(T) + len(corp) - len(T & cp))]
pd.DataFrame(rows, columns=["quantity", "value"]).to_csv(OUT, sep="\t", index=False)
with open(OUT, "r") as fh: body = fh.read()
open(OUT, "w").write(f"# check1_targeted_studies.py (supplement_number_checks); lineage_overlay.json {sha(f'{D}/regen/lineage_overlay.json')}; regen h4_papers.tsv {sha(f'{D}/regen/h4_papers/h4_papers.tsv')}; delivered h4_papers.tsv {sha(OLD)}\n" + body)
print(open(OUT).read())
