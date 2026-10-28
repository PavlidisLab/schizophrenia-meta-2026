"""01_classify.py: classify every printed gene name in the primary lists against current HGNC.
Reads (read only): work/sztr/data/lists (GitHub commit 746d6740)/*.tsv, analyses/primary_lists.json,
corpus/lineage.json, analyses/recurrence.py (imported for scenario_selection).
Writes: out/rows_classified.tsv.gz, out/lists_summary.tsv, out/list_set.tsv"""
import sys, os, csv, json, gzip, collections, datetime
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hlib as lib
H = lib.HGNC()
SZ = lib.SZ; OUT = lib.WD + "/out"
pl = json.load(open(SZ + "/analyses/primary_lists.json"))["studies"]
lin = json.load(open(SZ + "/corpus/lineage.json"))["nodes"]
year = {n["pmid"]: n.get("year") for n in lin if n.get("pmid")}
# recurrence selections
sys.path.insert(0, SZ + "/analyses")
import recurrence as R
sel_lists = {}   # list name -> set of scenario tags
units_all = None
for tag, scen in (("brain", "published, brain"), ("brainFDR", "published, brain, FDR lists")):
    sel, M, lists, units = R.scenario_selection(scen)
    for k in sel:
        for n in units["studies"][k]["primary"]:
            sel_lists.setdefault(n, {"study": k, "label": units["studies"][k]["label"], "scen": set()})["scen"].add(tag)
pm = sorted(set(R.scenario_selection('published, brain')[0]))
set64 = {}
for p in pm:
    for n in pl[p]["primary"]:
        set64.setdefault(n, {"study": p, "label": pl[p]["label"]})
allnames = sorted(set(set64) | set(sel_lists))
with open(OUT + "/list_set.tsv", "w") as fh:
    fh.write("# annotation drift pilot, list inventory; written %s by scripts/01_classify.py\n" % datetime.date.today())
    fh.write("# in64: primary lists of the %d studies of recurrence.scenario_selection('published, brain') on work/sztr (%s)\n" % (len(pm), "GitHub commit 746d6740"))
    fh.write("# in_brain / in_brainFDR: list is a modelled representative in recurrence.scenario_selection('published, brain') / ('published, brain, FDR lists'); primary_lists.json sha256 %s; year from corpus/lineage.json sha256 %s\n" % (lib.sha(SZ + "/analyses/primary_lists.json"), lib.sha(SZ + "/corpus/lineage.json")))
    fh.write("list\tstudy_pmid\tlabel\tyear\tin64\tin_brain\tin_brainFDR\tlist_file_sha256\n")
    for n in allnames:
        s = set64.get(n) or sel_lists[n]
        st = s["study"]
        fh.write("\t".join([n, st, s["label"], str(year.get(st, "")), "yes" if n in set64 else "no",
                            "yes" if "brain" in sel_lists.get(n, {}).get("scen", ()) else "no",
                            "yes" if "brainFDR" in sel_lists.get(n, {}).get("scen", ()) else "no",
                            lib.sha(SZ + "/data/lists/%s.tsv" % n)]) + "\n")
ROUTE_KIND = lambda r: r
rows_out = []
summ = []
with gzip.open(OUT + "/rows_classified.tsv.gz", "wt") as fo:
    fo.write("# annotation drift pilot: one row per row of each list file; written %s by scripts/01_classify.py\n" % datetime.date.today())
    fo.write("# hgnc: hgnc_complete_set.txt sha256 %s; withdrawn.txt sha256 %s; downloaded %s from storage.googleapis.com/public-download-files/hgnc/tsv/tsv/ (HGNC, genenames.org)\n" % (H.sha, H.sha_w, open(lib.HGNC_DIR + "/access_date.txt").read().strip()))
    fo.write("# input lists: work/sztr/data/lists (GitHub commit 746d6740)/<list>.tsv (sha256 per list in list_set.tsv); columns gene_reported, id_reported, symbol, map_route are copied from them\n")
    fo.write("# cls: class of gene_reported against current HGNC (lib.py); hgnc_ids: candidate HGNC ids; gene_flags: readthrough/antisense of the resolved gene; ens_cls: for ENSG ids, found in current HGNC or not; ens_symbol_agree: printed symbol resolves to the same gene as the ENSG id\n")
    w = csv.writer(fo, delimiter="\t", lineterminator="\n")
    w.writerow(["list", "stratum", "gene_reported", "id_reported", "symbol", "map_route", "cls", "cls_flags", "hgnc_ids", "gene_flags", "ens_cls", "ens_symbol_agree"])
    for n in allnames:
        path = SZ + "/data/lists/%s.tsv" % n
        rd = list(csv.DictReader([l for l in open(path, encoding="utf-8") if not l.startswith("#")], delimiter="\t"))
        routes = [r["map_route"] for r in rd]
        st = lib.id_stratum([r["id_reported"] for r in rd], routes)
        for r in rd:
            c, ids, fl = H.lookup(r["gene_reported"])
            gf = ";".join(sorted({H.gene_flags(i) for i in ids if H.gene_flags(i)})) if c in ("current", "previous", "alias") else ""
            ens_c = agree = ""
            idr = (r["id_reported"] or "").strip()
            if idr.startswith("ENSG"):
                base = idr.split(".")[0]
                eh = H.ens.get(base)
                ens_c = "in_hgnc" if eh else "not_in_hgnc"
                if eh and ids:
                    agree = "yes" if eh in ids else "no"
            w.writerow([n, st, r["gene_reported"], idr, r["symbol"], r["map_route"], c, fl, ",".join(ids), gf, ens_c, agree])
print("lists", len(allnames))
