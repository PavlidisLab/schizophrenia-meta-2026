#!/usr/bin/env python3
"""Copy into figures/data/ exactly the data files the figures read, from a delivery folder, and write the manifests.

Usage:
  python3 figures/scripts/sync_figure_data.py [--source DIR]

DIR defaults to $SZ_FIGURE_DATA. The list of files is figures/data/READS.tsv,
written by `rebuild_all.py --audit` with SZ_FIGURE_DATA pointing at DIR; EXTRA below adds files that live elsewhere.
Each copy keeps its relative path, so the scripts read <data>/<folder>/<file> as they read <DIR>/<folder>/<file>.
A file read only by the figure scripts of one supplementary note goes to supplementary_notes/<note>/data/, so a note can
be dropped by deleting its folder; every other file goes to figures/data/. Each of those data folders has its own
MANIFEST.tsv, which records, per file, its sha256, size and the first '# script:' and '# source' header lines,
which name the producing script and inputs in the analysis repo.

Run it after new tables are delivered and before committing, then `rebuild_all.py` against figures/data/.
"""
import argparse
import hashlib
import json
import os
import shutil
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "figures" / "data"
# files the figures read from outside the delivery folder: (path under figures/data, source path)
EXTRA = []  # none: every file the figures read is in the delivery folder and listed in READS.tsv


def header(path, key):
    try:
        with open(path, errors="replace") as fh:
            for i, line in enumerate(fh):
                if i > 40 or not line.startswith("#"):
                    break
                if line[1:].strip().lower().startswith(key):
                    return line[1:].strip()[:300].replace("\t", " ")
    except OSError:
        pass
    return ""


def owner(fig):
    """Supplementary note folder (supplementary_notes/<note>) of a figure's script, or None for a main/extended figure."""
    for m in list((REPO / "figures" / "results").glob(f"{fig}_meta.json")) + list(REPO.glob(f"supplementary_notes/*/results/{fig}_meta.json")):
        s = json.loads(m.read_text()).get("script", "")
        if s.startswith("supplementary_notes/"):
            return "/".join(s.split("/")[:2])
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default=os.environ.get("SZ_FIGURE_DATA"), required="SZ_FIGURE_DATA" not in os.environ)
    src = Path(ap.parse_args().source).expanduser()
    readers = {}
    for l in (DATA / "READS.tsv").read_text().splitlines():
        if l and not l.startswith("#"):
            fig, rel = l.split("\t")[:2]
            readers.setdefault(rel.strip(), set()).add(owner(fig))
    pairs = [(r, src / r, next(iter(o)) if len(o) == 1 else None) for r, o in sorted(readers.items())]
    pairs += [(r, p, None) for r, p in EXTRA]
    rows = {}  # data root -> manifest rows
    for r, p, note in pairs:
        root = REPO / note / "data" if note else DATA
        dst = root / r
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, dst)
        h = hashlib.sha256(dst.read_bytes()).hexdigest()
        rows.setdefault(root, []).append("\t".join([r, h, str(dst.stat().st_size), str(p).replace(str(Path.home()), "~"),
                                                    header(dst, "script"), header(dst, "source") or header(dst, "input")]))
    for root, rs in rows.items():
        (root / "MANIFEST.tsv").write_text("\n".join(["path\tsha256\tbytes\tcopied_from\tscript_header\tsource_header"] + rs) + "\n")
    print(f"copied {len(pairs)} files into {len(rows)} data folders; manifests written")


if __name__ == "__main__":
    main()
