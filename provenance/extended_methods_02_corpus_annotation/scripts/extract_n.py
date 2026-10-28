#!/usr/bin/env python3
"""Read case / control counts from each study's sample description (lineage.json
tissue_n, itself from the biolit eligibility pass) with one Haiku call.
Writes n_extracted.json; build_lineage.py uses it in preference to regex parsing.
Resumable (2026-09-25): only studies not already in n_extracted.json are sent, so a rerun cannot change stored counts."""
import json, re
from datetime import datetime, timezone
from pathlib import Path
import anthropic

HERE = Path(__file__).resolve().parent
MODEL = "claude-haiku-4-5-20251001"
nodes = [n for n in json.loads((HERE / "lineage.json").read_text())["nodes"] if not n.get("external")]
prev = json.loads((HERE / "n_extracted.json").read_text()) if (HERE / "n_extracted.json").exists() else {"counts": {}}
items = [{"id": n["id"], "text": n["tissue_n"]} for n in nodes if n["tissue_n"] and n["id"] not in prev["counts"]]
if not items:
    print("n_extracted.json: nothing new"); raise SystemExit(0)
prompt = ("For each study below, read the sample description and return the number of schizophrenia "
          "cases and the number of controls actually profiled. Count matched pairs as one case and one "
          "control each. Ignore other diagnoses (bipolar, depression) and validation-only cohorts unless "
          "they are the only numbers. Use null when a number is not stated. Return only a JSON list of "
          '{"id": <id exactly as given>, "cases": int|null, "controls": int|null}.\n\n'
          + json.dumps(items, ensure_ascii=False))
client = anthropic.Anthropic()  # repository copy: key from the ANTHROPIC_API_KEY environment variable (credential lookup removed)
resp = client.messages.create(model=MODEL, max_tokens=16000, messages=[{"role": "user", "content": prompt}])
text = "".join(b.text for b in resp.content if b.type == "text")
got = json.loads(re.search(r"\[.*\]", text, re.S).group(0))
sent = {i["id"] for i in items}
valid = {g["id"]: {"cases": g.get("cases"), "controls": g.get("controls")} for g in got if g.get("id") in sent}
bad = [g.get("id") for g in got if g.get("id") not in sent]
(HERE / "n_extracted.json").write_text(json.dumps({
    "model": MODEL, "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    "source": "lineage.json tissue_n (biolit eligibility pass)", "counts": {**prev["counts"], **valid}}, indent=1))
print(f"{len(valid)}/{len(items)} returned with valid ids; rejected ids: {bad}")
