#!/usr/bin/env python3
"""Agreement between donor-separate bulk dataset pairs (SZ t, L0 and L4) after removing a GTEx-derived brain axis built from the Hardy 3 vs 2 or Hardy 4 vs 3 contrast.

    python3 hardy3_axis_removal.py <inputs_root: analyses|analyses_laptop> <axis_dir> <out_tsv>   (run on remote with axes-handoff in place)

Reuses the setup of axes_pritzker_agonal_2026-10-09/decomposition_axes.py by exec (datasets, pairs, pH axis with leave-collection-out, rank-residual function) as single_axis_same_genes.py does,
and replaces the GTEx axis by: GTEx_4v2 (the existing axis, baseline), 3v2 and 4v3 (axis_brain_hardy*.tsv.gz written by gtex_hardy3_contrasts.py; mean t of Cortex, BA9, BA24).
For each target set (ALL 28 datasets; NOPRITZ; NOPRITZ2) and each removal (one axis; pH + axis; 3v2 + 4v3 together; pH + 3v2 + 4v3): median over both-bulk donor-separate pairs of
Spearman rho(t_a, t_b) raw and the correlation of rank residuals after regressing both profiles on the axes. Genes: in both datasets, in the pH axis for the pair, and in every axis named
('own' gene set = genes of the axes removed in that row; 'common' = genes of all four axes pH, 4v2, 3v2, 4v3 for every row)."""
import sys, hashlib, datetime, numpy as np, pandas as pd
from pathlib import Path
import os; HERE = Path(os.environ["AXES_ROOT"]) / "analyses" / "axes_pritzker_agonal_2026-10-09"
src = (HERE / "decomposition_axes.py").read_text().split("# per dataset")[0]
INPUT_ROOT, AXDIR, OUTF = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
sys.argv = ["x", INPUT_ROOT, str(Path.home() / "gtexh3" / "tmp_unused")]
g = {"__name__": "setup", "__file__": str(HERE / "decomposition_axes.py")}
exec(src, g)
from scipy.stats import spearmanr
pairs, SETS, resid_multi, pH_axis, coll_of = (g[k] for k in ("pairs", "SETS", "resid_multi", "pH_axis", "coll_of"))
DE = g["DE"]
def sha(p):
    h = hashlib.sha256(); h.update(Path(p).read_bytes()); return h.hexdigest()
AX = {"4v2": g["GTEX"]}
for nm, fn in (("3v2", "axis_brain_hardy3_vs_hardy2.tsv.gz"), ("4v3", "axis_brain_hardy4_vs_hardy3.tsv.gz")):
    AX[nm] = pd.read_csv(AXDIR / fn, sep="\t", comment="#", index_col=0)["mean_t_3_regions"].to_dict()
rho_axes = {f"{a}_vs_{b}": round(float(spearmanr(*zip(*[(AX[a][x], AX[b][x]) for x in sorted(set(AX[a]) & set(AX[b]))]))[0]), 3) for a, b in (("4v2", "3v2"), ("4v2", "4v3"), ("3v2", "4v3"))}
print("axis-axis Spearman over shared genes", rho_axes)
ONE = [["4v2"], ["3v2"], ["4v3"], ["pH", "4v2"], ["pH", "3v2"], ["pH", "4v3"], ["3v2", "4v3"], ["pH", "3v2", "4v3"], ["pH"]]
rows = []
for sname in ("ALL", "NOPRITZ", "NOPRITZ2"):
    P = [r for r in pairs if r["both_bulk"] == "TRUE" and r["dataset_a"] in SETS[sname] and r["dataset_b"] in SETS[sname]]
    for geneset in ("own", "common"):
        for rem in ONE:
            for lv in ("L0", "L4"):
                raw, par, ng = [], [], []
                for r in P:
                    a, b = r["dataset_a"], r["dataset_b"]; phg = pH_axis(coll_of(a) | coll_of(b))
                    gs = set(x for x in DE[a] if x in DE[b] and lv in DE[a][x] and lv in DE[b][x] and x in phg)
                    need = rem if geneset == "own" else ["4v2", "3v2", "4v3"]
                    for k in need:
                        if k != "pH": gs &= set(AX[k])
                    gl = sorted(gs); ta = np.array([DE[a][x][lv] for x in gl]); tb = np.array([DE[b][x][lv] for x in gl])
                    V = [np.array([phg[x] for x in gl]) if k == "pH" else np.array([AX[k][x] for x in gl]) for k in rem]
                    raw.append(spearmanr(ta, tb)[0]); par.append(np.corrcoef(resid_multi(ta, V), resid_multi(tb, V))[0, 1]); ng.append(len(gl))
                rows.append({"set": sname, "geneset": geneset, "axes_removed": "+".join(rem), "level": lv, "n_pairs": len(P), "median_genes": int(np.median(ng)),
                             "median_raw": round(float(np.median(raw)), 4), "median_after_removal": round(float(np.median(par)), 4), "share_removed": round(float(1 - np.median(par) / np.median(raw)), 3)})
R = pd.DataFrame(rows)
hdr = ["# script: analyses/gtex_death_speed_all_tissues_2026-10-10/hardy3_axis_removal.py", f"# date: {datetime.date.today().isoformat()}", f"# inputs root: {INPUT_ROOT} (remote axes-handoff; b_de and dataset_pairs.tsv, confounds within_controls_*.tsv)",
       f"# source dataset_pairs.tsv sha256 {sha(g['FS'] / 'dataset_pairs.tsv')}"] + [f"# source {fn} sha256 {sha(AXDIR / fn)}" for fn in ("axis_brain_hardy3_vs_hardy2.tsv.gz", "axis_brain_hardy4_vs_hardy3.tsv.gz")] + \
      [f"# source hardy_t_by_region.tsv.gz sha256 {sha('axes-handoff/analyses/gtex_ischemia_2026-10-08/hardy_t_by_region.tsv.gz')}",
       f"# axis-axis Spearman over shared genes: {rho_axes}; genes per axis: " + ", ".join(f"{k} {len(v)}" for k, v in AX.items()),
       "# 4v2 = existing GTEx axis (baseline; ALL/own/4v2 should reproduce 0.1177 -> 0.0717 at L0); median over both-bulk donor-separate pairs; share_removed = 1 - after/raw"]
with open(OUTF, "w") as f:
    f.write("\n".join(hdr) + "\n"); R.to_csv(f, sep="\t", index=False)
print(R[R.geneset == "own"].to_string(index=False))
