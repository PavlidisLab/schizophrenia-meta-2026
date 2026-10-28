#!/usr/bin/env python3
"""Check every file in the figure data folders against the MANIFEST.tsv of its folder.

Usage: python3 figures/scripts/verify_data.py

The data folders are figures/data/ (main and extended figures) and supplementary_notes/<note>/data/ (one per note);
each has a MANIFEST.tsv. A file that is missing or whose sha256 differs from the manifest is reported and the script
exits 1. Standard library only.
"""
import csv
import hashlib
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "figures" / "data"


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def manifests():
    return sorted([DATA / "MANIFEST.tsv", *REPO.glob("supplementary_notes/*/data/MANIFEST.tsv")])


def verify():
    """Check every file listed in each data folder's MANIFEST.tsv; a note folder that has been deleted is not checked."""
    rows = []
    for m in manifests():
        rows += [{"path": (m.parent / r["path"]).relative_to(REPO).as_posix(), "sha256": r["sha256"]}
                 for r in csv.DictReader(open(m), delimiter="\t")]
    missing = [r["path"] for r in rows if not (REPO / r["path"]).exists()]
    differ = [r["path"] for r in rows if (REPO / r["path"]).exists() and sha256(REPO / r["path"]) != r["sha256"]]
    for p in missing:
        print(f"missing  {p}")
    for p in differ:
        print(f"differs  {p}")
    print(f"{len(rows) - len(missing) - len(differ)} of {len(rows)} files match MANIFEST.tsv")
    return not (missing or differ)



if __name__ == "__main__":
    sys.exit(0 if verify() else 1)
