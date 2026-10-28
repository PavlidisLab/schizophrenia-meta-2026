"""Case and control counts for studies whose N the eligibility reading did not give, taken from the Step 2 list
contrasts, which state them ("SZ (8) vs control (9)").

Paul, 2026-09-24: the dashed circles (N unknown) should be fleshed out. The primary list's contrast is used first,
then the study's other first-look lists; the matched text is kept as evidence. Only studies with no N in
lineage.json are filled (build_lineage.py uses this after n_extracted.json and the eligibility text).

    python3 fill_n.py   # -> n_from_lists.json
"""
import datetime
import glob
import hashlib
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEG = HERE.parent / "analyses"
DATA = DEG.resolve().parent / "data"  # repo data/ (lists, manifests, reference, raw)
PAT = re.compile(r"(?:SZ|SCZ|schizophrenia|cases?|patients?)\s*\(?\s*(?:n\s*=\s*)?(\d{1,4})\)?[\s,]*(?:vs\.?|versus)\s*"
                 r"(?:healthy\s+)?(?:controls?|HCs?|CON|CTRL)\s*\(?\s*(?:n\s*=\s*)?(\d{1,4})", re.I)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()[:12]


def main():
    nodes = [n for n in json.loads((HERE / "lineage.json").read_text())["nodes"] if not n.get("external") and not n.get("landmark")]
    specs = json.loads((DATA / "manifest.json").read_text())["lists"]
    for f in sorted(glob.glob(str(DATA / "manifests" / "*.json"))):
        specs += json.loads(Path(f).read_text())["lists"]
    by = {}
    for x in specs:
        by.setdefault(x["pmid"], []).append(x)
    pl = json.loads((DEG / "primary_lists.json").read_text())["studies"]
    out = {}
    for n in nodes:
        if n.get("n") and not str(n.get("n_by", "")).startswith("Step 2 list contrast"):  # stable on rebuilds
            continue
        key = n.get("pmid") or "doi_" + n["id"].split(":", 1)[-1].replace("/", "_")
        prim = set((pl.get(key) or {}).get("primary") or [])
        cands = sorted(by.get(key, []), key=lambda x: (f"{x['pmid']}__{x['list_id']}" not in prim, x.get("role", "first_look") != "first_look"))
        for x in cands:
            m = PAT.search(x.get("contrast", ""))
            if m:
                cases, controls = int(m.group(1)), int(m.group(2))
                out[key] = {"label": n["label"], "cases": cases, "controls": controls, "n": cases + controls,
                            "list": f"{x['pmid']}__{x['list_id']}", "quote": m.group(0)}
                break
    doc = {"generated": datetime.date.today().isoformat(), "script": "fill_n.py",
           "inputs": {"lineage.json": sha(HERE / "lineage.json"), "primary_lists.json": sha(DEG / "primary_lists.json")},
           "note": "N = cases + controls as stated in the list's contrast; the contrast can describe one region or subgroup of the study.",
           "studies": out}
    (HERE / "n_from_lists.json").write_text(json.dumps(doc, indent=1, ensure_ascii=False))
    print(f"{len(out)} studies filled")


if __name__ == "__main__":
    main()
