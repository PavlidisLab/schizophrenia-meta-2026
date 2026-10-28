#!/usr/bin/env python3
"""Keyword pass over the full texts and supplements of the 64 primary brain studies for brain pH and agonal state / manner of death.

    python3 extract_and_scan.py       -> text_cache/<pmid>/*.txt (text extracted from data/raw/<pmid>), scan_files.tsv (what was read), scan_hits.tsv (every matching sentence with source file, cached-text line)

Studies: recurrence.py scenario 'published, brain' (64 representative lists; rerun-full/analyses). Sources under rerun-full/data/raw/<pmid>/: .txt as is; .html/.xml with tags stripped; .pdf through pdftotext -layout; .docx paragraphs and tables (python-docx);
.xlsx/.xls header rows only (a column or cell named pH / brain pH / agonal / manner / cause of death is recorded as a 'table column' hit); .doc, .ppt, .zip, images are not read (listed as unread).
pH hits: sentences with the word pH that also carry a postmortem-tissue cue (brain, tissue, postmortem, PMI, cortex, donor, subject, control, covariate, RIN, matched, autopsy); buffer chemistry sentences are dropped. Agonal hits: agonal, agony, manner / cause / mode / circumstances of death,
terminal illness / coma, hypoxi*, resuscitat*, ventilat*, pyrexia, fever, seizure at death, brain death. Line = line in the cached extraction (text_cache), the source file is named in every row."""
import csv, glob, hashlib, html, os, re, subprocess, sys, datetime, collections
sys.path.insert(0, "rerun-full/analyses")
import recurrence as R
HERE = os.path.dirname(os.path.abspath(__file__)); RAW = "rerun-full/data/raw"; CACHE = os.path.join(HERE, "text_cache"); PDFTXT = "rerun-envs/poppler/bin/pdftotext"
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
sel, M, lists, units = R.scenario_selection("published, brain")
PH = re.compile(r"\bpH\b|\bPH\b(?=[ ,.;:)])")
CUE = re.compile(r"brain|tissue|post-?mortem|PMI|cortex|cortical|donor|subject|control|covariate|RIN|matched|autopsy|cerebell|hippocamp|sample", re.I)
BUF = re.compile(r"buffer|Tris|HEPES|PBS|phosphate|citrate|acetate|saline|medium|solution|incubat|centrifug|lysis|homogeni[sz]|wash", re.I)
AG = re.compile(r"agonal|agony|manner of death|cause of death|causes of death|mode of death|circumstances of death|terminal (?:illness|coma|phase|state)|\bcoma\b|hypoxi|resuscitat|ventilat|pyrexia|\bfever\b|brain death|prolonged dying|suicide|\bsuicid", re.I)
AGSTRICT = re.compile(r"agonal|agony|manner of death|mode of death|circumstances of death|terminal (?:illness|coma|phase|state)|\bcoma\b|hypoxi|resuscitat|ventilat|pyrexia|\bfever\b|brain death|prolonged dying", re.I)
def extract(path):
    e = os.path.splitext(path)[1].lower()
    try:
        if e == ".txt": return open(path, errors="replace").read()
        if e in (".html", ".htm", ".xml"):
            t = open(path, errors="replace").read(); t = re.sub(r"(?is)<(script|style).*?</\1>", " ", t); t = re.sub(r"(?s)<[^>]+>", " ", t); return html.unescape(t)
        if e == ".pdf": return subprocess.run([PDFTXT, "-layout", path, "-"], capture_output=True, text=True, timeout=120).stdout
        if e == ".docx":
            import docx
            d = docx.Document(path); out = [p.text for p in d.paragraphs]
            for t in d.tables:
                for r in t.rows: out.append(" | ".join(c.text.strip() for c in r.cells))
            return "\n".join(out)
    except Exception as ex:
        return None
    return None
files, hits = [], []
for k in sorted(sel, key=lambda x: int(x) if x.isdigit() else 0):
    lab = units["studies"][k]["label"]; d = os.path.join(RAW, k)
    for p in sorted(glob.glob(d + "/**/*", recursive=True)):
        if not os.path.isfile(p): continue
        rel = os.path.relpath(p, RAW); e = os.path.splitext(p)[1].lower(); row = {"pmid": k, "label": lab, "file": rel, "bytes": os.path.getsize(p), "sha256": sha(p), "status": ""}
        if e in (".xlsx", ".xls"):
            heads = []
            try:
                if e == ".xlsx":
                    import openpyxl
                    wb = openpyxl.load_workbook(p, read_only=True, data_only=True)
                    for ws in wb.worksheets:
                        for i, r in enumerate(ws.iter_rows(values_only=True, max_row=6)):
                            heads += [str(c) for c in r if c is not None and re.search(r"\bpH\b|agonal|manner of death|cause of death", str(c), re.I)]
                else:
                    import xlrd
                    wb = xlrd.open_workbook(p)
                    for ws in wb.sheets():
                        for i in range(min(6, ws.nrows)): heads += [str(c) for c in ws.row_values(i) if re.search(r"\bpH\b|agonal|manner of death|cause of death", str(c), re.I)]
            except Exception as ex: row["status"] = "unreadable table"
            row["status"] = row["status"] or ("table header hit" if heads else "table, no header hit")
            for h in dict.fromkeys(heads): hits.append({"pmid": k, "label": lab, "kind": "table_header", "file": rel, "line": "", "sentence": h[:200], "pH": bool(PH.search(h)), "agonal": bool(AG.search(h))})
            files.append(row); continue
        if e not in (".txt", ".html", ".htm", ".xml", ".pdf", ".docx"):
            row["status"] = "not read (" + (e or "no extension") + ")"; files.append(row); continue
        txt = extract(p)
        if txt is None or len(txt.strip()) < 200: row["status"] = "no text extracted"; files.append(row); continue
        os.makedirs(os.path.join(CACHE, k), exist_ok=True); cp = os.path.join(CACHE, k, rel.replace("/", "__") + ".txt"); open(cp, "w").write(txt)
        row["status"] = f"read ({len(txt)} chars)"; files.append(row)
        flat = re.sub(r"[ \t]+", " ", txt)
        for m in list(PH.finditer(flat)) + list(AG.finditer(flat)):
            s0 = max(flat.rfind(". ", 0, m.start()), flat.rfind("\n\n", 0, m.start())); s0 = 0 if s0 < 0 else s0 + 1
            e1 = [x for x in (flat.find(". ", m.end()), flat.find("\n\n", m.end())) if x >= 0]; e1 = min(e1) if e1 else len(flat)
            sent = re.sub(r"\s+", " ", flat[s0:min(e1 + 1, s0 + 700)]).strip()
            isph = bool(PH.search(sent)); isag = bool(AG.search(sent))
            if isph and (BUF.search(sent) and not CUE.search(sent) or not CUE.search(sent)): isph = False
            if not (isph or isag): continue
            hits.append({"pmid": k, "label": lab, "kind": "text", "file": rel, "line": flat.count("\n", 0, m.start()) + 1, "sentence": sent, "pH": isph, "agonal": bool(AGSTRICT.search(sent)) or (isag and not isph and bool(AG.search(sent))), "cached_text": os.path.relpath(cp, HERE)})
seen, uniq = set(), []
for h in hits:
    key = (h["pmid"], h["sentence"][:160])
    if key in seen: continue
    seen.add(key); uniq.append(h)
def wr(name, rows, head):
    import pandas as pd
    with open(os.path.join(HERE, name), "w") as fh:
        fh.write(head + "\n"); pd.DataFrame(rows).to_csv(fh, sep="\t", index=False, line_terminator="\n")
hd = f"# script: extract_and_scan.py, run {datetime.date.today().isoformat()}; 64 studies of recurrence.py scenario 'published, brain'; sources rerun-full/data/raw/<pmid>; line = line in text_cache/<pmid>/<file>.txt"
wr("scan_files.tsv", files, hd); wr("scan_hits.tsv", uniq, hd + f"; {len(uniq)} distinct sentences")
c = collections.Counter(f["status"].split(" (")[0] for f in files); print(len(sel), "studies;", len(files), "files;", dict(c)); print("hits", len(uniq), "pH sentences", sum(1 for h in uniq if h["pH"]), "agonal sentences", sum(1 for h in uniq if h["agonal"]))
