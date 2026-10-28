#!/usr/bin/env python3
"""Rebuild every figure the manuscript, figure pages and supplement use, from figures/data/.

Usage:
  python3 figures/scripts/rebuild_all.py            rebuild all figures
  python3 figures/scripts/rebuild_all.py --audit    also record every file each script reads (figures/data/READS.tsv)
  python3 figures/scripts/rebuild_all.py --list     list the figures and their scripts, build nothing

The figures and their scripts are listed in figures/scripts/FIGURES.tsv; where no script is named there, it is the one that
makes the fk.save(fig, "<name>", ...) call for that figure (figures/scripts/ and supplementary_notes/<note>/). A figure that
is not listed is not rebuilt. The images are generated and not tracked.
"""
import json
import os
import re
import runpy
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "figures" / "scripts"
FIGURES = SCRIPTS / "FIGURES.tsv"


def figure_scripts():
    """Figure name -> script, found from the fk.save(fig, "<name>", ...) call in each build script."""
    out = {}
    for s in [*SCRIPTS.glob("build_*.py"), *(REPO / "supplementary_notes").glob("*/build_*.py")]:
        for name in re.findall(r'fk\.save\(\s*\w+,\s*"([A-Za-z0-9_]+)"', s.read_text()):
            out[name] = s.relative_to(REPO).as_posix()
    return out


def used_figures():
    """(figure, script) from FIGURES.tsv: the script named there, else the one found from the fk.save() call; empty if neither."""
    detected = figure_scripts(); out = {}
    for l in FIGURES.read_text().splitlines():
        if l.strip() and not l.startswith("#"):
            c = l.split("\t"); out[c[0]] = (c[1] if len(c) > 1 and c[1].strip() else detected.get(c[0], ""))
    return sorted(out.items())


def main():
    if "--list" in sys.argv:
        for fig, script in used_figures():
            print(f"{fig}\t{script or '(no script in this repo)'}")
        return
    audit = "--audit" in sys.argv
    sys.path.insert(0, str(SCRIPTS))
    import figkit as fk
    fd = Path(fk.FD).resolve()
    roots = [fd] if "SZ_FIGURE_DATA" in os.environ else [fd, *(REPO / "supplementary_notes").glob("*/data")]
    reads = {}
    current = [None]
    if audit:
        def hook(event, args):
            if event == "open" and current[0] and isinstance(args[0], (str, Path)):
                p = Path(args[0])
                try:
                    rp = p.resolve()
                except OSError:
                    return
                if (len(args) < 2 or "r" in str(args[1] or "r")):
                    for root in roots:
                        if root.resolve() in rp.parents and rp.exists():
                            reads.setdefault(current[0], set()).add(rp.relative_to(root.resolve()).as_posix())
        sys.addaudithook(hook)
    failed = []
    for fig, script in used_figures():
        if not script:
            print(f"skip {fig}: not built by a script in this repo (the committed image is used)")
            continue
        current[0] = fig
        print(f"build {fig} <- {script}", flush=True)
        try:
            runpy.run_path(str(REPO / script), run_name="__main__")
        except SystemExit:
            pass
        except Exception as e:  # report and continue, so one broken figure does not hide the others
            failed.append((fig, repr(e)))
        current[0] = None
    if audit:
        lines = ["# figure\tfile (relative to figures/data)"]
        for fig in sorted(reads):
            lines += [f"{fig}\t{f}" for f in sorted(reads[fig])]
        (REPO / "figures" / "data" / "READS.tsv").write_text("\n".join(lines) + "\n")
        print(f"recorded {sum(len(v) for v in reads.values())} reads in figures/data/READS.tsv")
    for fig, err in failed:
        print(f"FAILED {fig}: {err}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
