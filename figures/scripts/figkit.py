"""Shared layout and style for the manuscript figures.

Each figure is drawn as ONE matplotlib figure at final print size (180 mm wide,
Nature double column), with every panel an axes in a GridSpec, so text size,
line weight and panel letters are identical across panels and figures.

Palette and font rules follow the global figure contract (CLAUDE.md); lighten()
comes from the plotting skill. Titles and source lines are not drawn on panels:
the legend in Figures.md says what each panel shows, and provenance goes to the
*_meta.json sidecar.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt



def lighten(hex_color: str, frac: float = 0.55) -> str:
    """Blend a hex colour with white (frac = white fraction). Copied from the plotting skill's pavlab_style.lighten so
    that the figures build without the skill installed."""
    h = hex_color.lstrip("#")
    if len(h) != 6:
        raise ValueError(f"lighten expects a 6-digit hex colour, got {hex_color!r}")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f"#{int(r + (255 - r) * frac):02x}{int(g + (255 - g) * frac):02x}{int(b + (255 - b) * frac):02x}"


# Figure data: the tables in figures/data/ (verify_data.py checks them; sync_figure_data.py refreshes them from a
# delivery folder). SZ_FIGURE_DATA points the scripts at another folder.
FD = Path(os.environ.get("SZ_FIGURE_DATA", Path(__file__).resolve().parents[1] / "data"))
OUT_DIR = Path(__file__).resolve().parent.parent / "results"
REPO = Path(__file__).resolve().parents[2]


def note_data(script_file) -> Path:
    """Data folder of a supplementary-note script: supplementary_notes/<note>/data/ beside it (SZ_FIGURE_DATA overrides)."""
    return FD if "SZ_FIGURE_DATA" in os.environ else Path(script_file).resolve().parent / "data"

MM = 1 / 25.4
WIDTH_MM = 180

ACCENT = "#2563eb"    # schizophrenia / primary series
ACCENT_2 = "#10b981"  # age / positive control
ACCENT_3 = "#f59e0b"  # DE prior / pH-related
ACCENT_4 = "#ef4444"  # down / failure, sparingly
NEUTRAL = "#9ca3af"
GRID = "#e5e7eb"
TEXT = "#1f2937"
SUBTLE = "#6b7280"
UP, DOWN = ACCENT, ACCENT_4

SMALL = 6.5  # in-panel annotations, category tick labels
TINY = 5.2   # dense gene labels only

# The six studies that pool donors from several banks without per-bank results, so belong to no cohort cluster
# (Supplementary Methods; recurrence_report.md, Cohort tier). Shared by Extended Fig. 2e and Fig. 2e.
MULTI_BANK = {"Maycox 2009", "Pai 2019", "Mendizabal 2019", "Ruzicka 2024", "Mendez-Victoriano 2026", "Lee 2026 (PsychAD)", "Saetre 2007"}

NULL_COLOURS = {"uniform": ACCENT, "deprior": ACCENT_3, "exprmatched": ACCENT_2}
NULL_LABELS = {"uniform": "random genes", "deprior": "DE-prior-weighted random genes", "exprmatched": "expression-matched random genes"}


def apply_style() -> None:
    plt.rcParams.update({
        "font.family": ["Helvetica", "Arial", "sans-serif"],
        "font.sans-serif": ["Helvetica", "Arial"],
        "svg.fonttype": "none",
        "svg.hashsalt": "figkit",  # fixed clip-path ids, so an unchanged figure rebuilds to the same SVG
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "legend.frameon": False,
        "axes.edgecolor": TEXT,
        "axes.labelcolor": TEXT,
        "axes.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "axes.grid.axis": "y",
        "grid.color": GRID,
        "grid.linewidth": 0.4,
        "xtick.color": SUBTLE,
        "ytick.color": SUBTLE,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "lines.linewidth": 1.0,
        "lines.markersize": 3.5,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    })


def figure(height_mm: float):
    apply_style()
    return plt.figure(figsize=(WIDTH_MM * MM, height_mm * MM))


def letter(ax, s: str, dx: float = -0.02, dy: float = 0.0) -> None:
    """Bold lowercase panel letter at the top-left of the axes' bounding box (incl. tick labels)."""
    fig = ax.figure
    fig.canvas.draw()
    bb = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(max(bb.x0 + dx, 0.004), bb.y1 + dy, s, fontsize=10, fontweight="bold", color="black",
             ha="left", va="bottom")
    # panel box (figure fractions, y up, including tick labels and the letter), written to the *_meta.json by save()
    fig.__dict__.setdefault("_panel_boxes", {}).setdefault(s, []).append(
        [round(float(max(bb.x0 + dx, 0.0)), 4), round(float(bb.y0), 4), round(float(bb.x1), 4), round(float(min(bb.y1 + 0.02, 1.0)), 4)])


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _rel(p: Path) -> str:
    """Path as recorded in *_meta.json: relative to the repository when inside it, else with ~ for the home folder."""
    p = Path(p).resolve()
    repo = REPO
    return p.relative_to(repo).as_posix() if repo in p.parents else str(p).replace(str(Path.home()), "~", 1)


def _plain_log_label(v, _pos=None) -> str:
    """Log-axis tick label: plain numbers below 10,000 (0.1, 1, 10, 1,000), powers of ten from 10^4 up."""
    if v <= 0:
        return ""
    if v < 1e4:
        return f"{v:,.10g}"
    e = np.log10(v)
    return f"$10^{{{int(round(e))}}}$" if abs(e - round(e)) < 1e-9 else ""


def plain_log_ticks(fig) -> None:
    """Replace matplotlib's default 10^n labels on every log axis that still uses the default formatter."""
    from matplotlib.ticker import FuncFormatter, LogFormatter, NullFormatter
    for ax in fig.axes:
        for axis, scale in ((ax.xaxis, ax.get_xscale()), (ax.yaxis, ax.get_yscale())):
            if scale == "log" and isinstance(axis.get_major_formatter(), LogFormatter):
                axis.set_major_formatter(FuncFormatter(_plain_log_label))
                axis.set_minor_formatter(NullFormatter())


def pending_inputs(inputs) -> list[str]:
    """Names of input files whose header says PENDING_375 (delivered unchanged from the 209-gene state)."""
    out = []
    for p in inputs:
        try:
            with open(p, errors="replace") as fh:
                head = "".join(next(fh, "") for _ in range(3))
        except OSError:
            continue
        if "PENDING_375" in head:
            out.append(Path(p).name)
    return out


def save(fig, name: str, script: str, inputs: list[Path], notes: str = "") -> None:
    # a supplementary-note script (supplementary_notes/<note>/...) writes to the results/ folder beside it
    out = REPO / Path(script).parent / "results" if script.startswith("supplementary_notes/") else OUT_DIR
    out.mkdir(exist_ok=True)
    plain_log_ticks(fig)
    panel_boxes = getattr(fig, "_panel_boxes", None)
    titles = []  # axes titles as boxes in figure fractions (y up), so slide crops can leave them out
    fig.canvas.draw()
    _r = fig.canvas.get_renderer()
    for _ax in fig.axes:
        for _t in (_ax._left_title, _ax.title, _ax._right_title):
            if _t.get_text().strip():
                _bb = _t.get_window_extent(_r).transformed(fig.transFigure.inverted())
                titles.append([round(float(_bb.x0), 4), round(float(_bb.y0), 4), round(float(_bb.x1), 4), round(float(_bb.y1), 4)])
    pend = pending_inputs(inputs)
    if pend:  # an input still holds the old (209-gene, 64-study) content: say so on the figure and in the sidecar
        fig.text(0.5, 0.004, "PROVISIONAL: " + ", ".join(pend) + " still hold the previous gene set (209) until the rerun is delivered",
                 ha="center", va="bottom", fontsize=7, color=ACCENT_4)
        notes = (notes + " " if notes else "") + "PROVISIONAL (PENDING_375 inputs): " + ", ".join(pend)
    svg, png = out / f"{name}.svg", out / f"{name}.png"
    fig.savefig(svg, metadata={"Date": None})
    fig.savefig(png, dpi=300)
    plt.close(fig)
    meta = {
        "artifact": [svg.name, png.name],
        "sha256": {svg.name: _sha(svg), png.name: _sha(png)},
        "script": script,
        "inputs": {_rel(p): _sha(p) for p in inputs},
        "input_meta": {_rel(m): json.loads(m.read_text()).get("run_id") if m.exists() else None
                       for m in sorted({p.parent / f"{p.stem}_meta.json" for p in inputs})},
        "built": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "notes": notes,
    }
    if panel_boxes:
        meta["panel_boxes"] = panel_boxes
    if titles:
        meta["title_boxes"] = titles
    (out / f"{name}_meta.json").write_text(json.dumps(meta, indent=2) + "\n")


def n_recurrent() -> int:
    """Size of the recurrent gene set (the 375, study tier q < 0.1), from the fig06 export."""
    d = pd.read_csv(FD / "fig06" / "fig06_b_direction_summary.tsv", sep="\t", comment="#").set_index("gene_set")
    return int(d.loc["study tier q<0.1 (uniform null)", "genes"])


def n_primary() -> int:
    """Number of primary brain studies (one row per study in the fig02 study annotation)."""
    return len(pd.read_csv(FD / "fig02" / "fig02_b_study_annotation.tsv", sep="\t", comment="#"))
