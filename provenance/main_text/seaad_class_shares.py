"""SEA-AD MTG class shares (donor-class pseudobulk CPM, EWCE specificity) for the 383 recurrent genes.
Run on the remote environment:  cd ~/second_reference_seaad_2026-10-19 && OMP_NUM_THREADS=1 nice -n 19 python -I seaad_class_shares.py
Stage 1 streams the CSR UMI layer once (10000 rows per chunk) into donor x subclass sums for the gene panel and for all genes (total).
Stage 2 mirrors analyses/interneuron_ph_2026-10-01/scripts/specificity.py (class_means; min 20 nuclei per donor-class; no CPM >= 1 filter).
"""
import os, sys, time, hashlib, resource, datetime
for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[k] = "1"
import numpy as np, pandas as pd, h5py
from scipy.stats import spearmanr, pearsonr

H5 = os.path.expanduser("~/nserc_pilots_2026-10-14/seaad_sst_loss_check/data/SEAAD_MTG_RNAseq_final-nuclei.2026-06-22.h5ad")
H5REL = "nserc_pilots_2026-10-14/seaad_sst_loss_check/data/SEAAD_MTG_RNAseq_final-nuclei.2026-06-22.h5ad"
OUT = os.path.expanduser("~/second_reference_seaad_2026-10-19")
MIN_NUC = 20
FINE_ORDER = ["Astro", "EN", "IN_other", "LAMP5", "Micro", "OPC", "Oligo", "PVALB", "SST", "VIP", "Vascular"]
IN_FINE = ["SST", "PVALB", "VIP", "LAMP5", "IN_other"]
MAP = {"Sst": "SST", "Pvalb": "PVALB", "Chandelier": "PVALB", "Vip": "VIP", "Lamp5": "LAMP5", "Lamp5 Lhx6": "LAMP5",
       "Sncg": "IN_other", "Pax6": "IN_other", "Astrocyte": "Astro", "Oligodendrocyte": "Oligo", "OPC": "OPC",
       "Immune": "Micro", "Endothelial": "Vascular", "VLMC & Perivascular": "Vascular"}
for s in ["L2/3 IT", "L4 IT", "L5 IT", "L5 ET", "L5/6 NP", "L6 IT", "L6 IT Car3", "L6 CT", "L6b"]: MAP[s] = "EN"
MAP_SENS = dict(MAP, **{"Sst Chodl": "SST"})

def dec(a): return np.array([x.decode() if isinstance(x, bytes) else x for x in a], dtype=object)

def cat(o, k):
    v = o[k]; return v["codes"][:], dec(v["categories"][:])

def hgnc_lookup(genes):
    hg = pd.read_csv(os.path.join(OUT, "hgnc_complete_set.txt"), sep="\t", dtype=str, usecols=["symbol", "ensembl_gene_id", "prev_symbol", "alias_symbol"])
    bysym = {r.symbol: r.ensembl_gene_id for r in hg.itertuples() if isinstance(r.ensembl_gene_id, str)}
    prev = {}
    for r in hg.itertuples():
        if isinstance(r.ensembl_gene_id, str) and isinstance(r.prev_symbol, str):
            for p in r.prev_symbol.split("|"): prev.setdefault(p, []).append((r.symbol, r.ensembl_gene_id))
    res = {}
    for g in genes:
        if g in bysym: res[g] = (g, bysym[g], "current_symbol")
        elif g in prev and len(prev[g]) == 1: res[g] = (prev[g][0][0], prev[g][0][1], "previous_symbol->" + prev[g][0][0])
        else: res[g] = (g, None, "unresolved_in_hgnc")
    return res

def stage1():
    t0 = time.time()
    genes = list(pd.read_csv(os.path.join(OUT, "recurrence_genes.tsv"), sep="\t", comment="#").query("scenario == 'published, brain' and q_bh < 0.1").gene)
    ps = pd.read_csv(os.path.join(OUT, "interneuron_specificity_383_v2.tsv"), sep="\t", comment="#")
    assert set(genes) == set(ps.gene), "gene list differs from PsychAD table"
    hl = hgnc_lookup(genes)
    f = h5py.File(H5, "r")
    vid = dec(f["var"]["index"][:]); vens = dec(f["var"]["gene_ids"][:])
    e2i = {}; s2i = {}
    for i, (s, e) in enumerate(zip(vid, vens)):
        e2i.setdefault(e, []).append(i); s2i.setdefault(s, []).append(i)
    rows = []
    for g in genes:
        cur, ens, how = hl[g]
        col, found = None, ""
        if ens and ens in e2i and len(e2i[ens]) == 1: col, found = e2i[ens][0], "ensembl_id" + ("" if how == "current_symbol" else "; " + how)
        elif cur in s2i and len(s2i[cur]) == 1: col, found = s2i[cur][0], "symbol_in_h5ad_var" + ("" if how == "current_symbol" else "; " + how)
        else: found = "not_in_h5ad_var (ensembl %s)" % ens
        rows.append((g, cur, ens or "", col, found))
    gm = pd.DataFrame(rows, columns=["gene", "hgnc_symbol", "ensembl_id", "col", "found_by"])
    o = f["obs"]
    dc, dcat = cat(o, "Donor ID"); sc, scat = cat(o, "Subclass"); cc, ccat = cat(o, "Cognitive Status")
    xc, xcat = cat(o, "Sex"); age = o["Age at Death"][:]
    nd, ns = len(dcat), len(scat)
    don = pd.DataFrame({"donor": dcat[dc], "cog": ccat[cc], "sex": xcat[xc], "age": age}).groupby("donor").agg(
        cog=("cog", "first"), sex=("sex", "first"), age=("age", "first"), ncog=("cog", "nunique"), nsex=("sex", "nunique"))
    assert (don.ncog == 1).all() and (don.nsex == 1).all()
    grp = dc.astype(np.int64) * ns + sc.astype(np.int64)
    ng = nd * ns
    nuc = np.bincount(grp, minlength=ng)
    sel = gm.dropna(subset=["col"]); lut = np.full(len(vid), -1, np.int32); lut[sel.col.astype(int).values] = np.arange(len(sel)); nsel = len(sel)
    mat = f["layers"]["UMIs"]; ind, dat, ptr = mat["indices"], mat["data"], mat["indptr"][:]
    n = len(ptr) - 1; panel = np.zeros(ng * nsel); tot = np.zeros(ng); CH = 10000; nonint = 0
    for r0 in range(0, n, CH):
        r1 = min(n, r0 + CH); lo, hi = ptr[r0], ptr[r1]
        v = dat[lo:hi]; ii = ind[lo:hi]
        nonint += int(np.count_nonzero(v != np.round(v)))
        rl = np.repeat(np.arange(r1 - r0), np.diff(ptr[r0:r1 + 1]))
        g = grp[r0:r1][rl]
        tot += np.bincount(g, weights=v, minlength=ng)
        m = lut[ii]; p = m >= 0
        panel += np.bincount(g[p] * nsel + m[p], weights=v[p], minlength=ng * nsel)
        if (r0 // CH) % 10 == 0: print("rows", r1, n, "%.0fs" % (time.time() - t0), "maxrss_GB %.1f" % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6), flush=True)
    panel = panel.reshape(ng, nsel)
    np.savez_compressed(os.path.join(OUT, "pseudobulk_donor_subclass_panel.npz"), panel=panel, tot=tot, nuc=nuc, donors=dcat, subclasses=scat, panel_genes=sel.gene.values, nonint=nonint)
    gm.to_csv(os.path.join(OUT, "gene_map.tsv"), sep="\t", index=False); don.to_csv(os.path.join(OUT, "donor_info.tsv"), sep="\t")
    print("stage1 done %.0fs nonint=%d" % (time.time() - t0, nonint), flush=True)

def stage2(z, don, gm, donors, mapping):
    """Returns tables with CPM per fine/coarse class for the chosen donors."""
    dcat, scat = list(z["donors"]), list(z["subclasses"])
    ns = len(scat); nd = len(dcat)
    P = z["panel"].reshape(nd, ns, -1); T = z["tot"].reshape(nd, ns); N = z["nuc"].reshape(nd, ns)
    sub2fine = {j: mapping[s] for j, s in enumerate(scat) if s in mapping}
    out = {}
    for level in ("fine", "coarse"):
        lab = {j: (c if level == "fine" else ("IN" if c in IN_FINE else c)) for j, c in sub2fine.items()}
        classes = sorted(set(lab.values())); cm = {}; ndon = {}; nnuc = {}
        for c in classes:
            js = [j for j, l in lab.items() if l == c]; rows = []; nn = 0
            for d in donors:
                di = dcat.index(d)
                if N[di, js].sum() < MIN_NUC: continue
                tt = T[di, js].sum()
                if tt > 0: rows.append(P[di, js].sum(axis=0) / tt * 1e6); nn += int(N[di, js].sum())
            cm[c] = np.mean(rows, axis=0); ndon[c] = len(rows); nnuc[c] = nn
        out[level] = (pd.DataFrame(cm, index=z["panel_genes"]), ndon, nnuc)
    fcm, nd_f, nn_f = out["fine"]; ccm, nd_c, nn_c = out["coarse"]
    fcm = fcm[FINE_ORDER]
    T_ = pd.DataFrame(index=gm.gene)
    spec_f = fcm.div(fcm.sum(axis=1), axis=0); spec_c = ccm.div(ccm.sum(axis=1), axis=0)
    full = pd.concat([fcm.add_prefix("cpm_"), spec_f.add_prefix("spec_"), spec_c[["IN"]].add_prefix("spec_coarse_"), spec_c[["EN"]].add_prefix("spec_coarse_"), ccm[["IN"]].add_prefix("cpm_coarse_")], axis=1)
    full["spec_IN_fine_sum"] = spec_f[IN_FINE].sum(axis=1)
    full["max_class_mean_cpm"] = fcm.max(axis=1)
    full["n_classes_measured"] = (fcm > 0).sum(axis=1)
    for c in FINE_ORDER: full["n_donors_" + c] = nd_f[c]
    full = T_.join(full)
    return full, nd_f, nn_f, nd_c, nn_c

def main():
    if not os.path.exists(os.path.join(OUT, "pseudobulk_donor_subclass_panel.npz")): stage1()
    t0 = time.time()
    z = np.load(os.path.join(OUT, "pseudobulk_donor_subclass_panel.npz"), allow_pickle=True)
    don = pd.read_csv(os.path.join(OUT, "donor_info.tsv"), sep="\t", index_col=0); gm = pd.read_csv(os.path.join(OUT, "gene_map.tsv"), sep="\t")
    pg = list(z["panel_genes"])
    ps = pd.read_csv(os.path.join(OUT, "interneuron_specificity_383_v2.tsv"), sep="\t", comment="#").set_index("gene")
    nodem = list(don.index[don.cog == "No dementia"]); alld = list(don.index)
    runs = {"nodem": stage2(z, don, gm, nodem, MAP), "all": stage2(z, don, gm, alld, MAP), "nodem_sens": stage2(z, don, gm, nodem, MAP_SENS)}
    sha = hashlib.sha256()
    with open(H5, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""): sha.update(b)
    size = os.path.getsize(H5); today = datetime.date.today().isoformat()
    gm = gm.set_index("gene")
    def finish(full, tag):
        full = full.reindex(gm.index)
        full.insert(0, "status", "SEA-AD MTG (%s)" % tag)
        full.insert(1, "ensembl_id", gm.ensembl_id.values); full.insert(2, "found_by", gm.found_by.values)
        full.index.name = "gene"
        return full
    tabs = {k: finish(v[0], "normal-donor nuclei" if k.startswith("nodem") else "all-donor nuclei") for k, v in runs.items()}
    _, nd_f, nn_f, nd_c, nn_c = runs["nodem"]
    dn = don.loc[nodem]
    hdr_common = [
        "# script: second_reference_seaad_2026-10-19/seaad_class_shares.py (class_means logic of analyses/interneuron_ph_2026-10-01/scripts/specificity.py), run %s on the remote environment" % today,
        "# input: %s, %d bytes, sha256 %s (layers/UMIs, raw UMI counts, %d nuclei x %d genes; SEA-AD MTG RNAseq final-nuclei, release file dated 2026-06-22; Gabitto et al. 2024 Nat Neurosci doi 10.1038/s41593-024-01774-5; h5ad uns title: MTG: Seattle Alzheimer's Disease Atlas (SEA-AD); no version string recorded in the file beyond the file date)" % (H5REL, size, sha.hexdigest(), int(z["nuc"].sum()), 36601),
        "# input: hgnc_complete_set.txt sha256 %s (copy of data/reference/hgnc_complete_set.txt); recurrence_genes.tsv (scenario 'published, brain', q_bh < 0.1) sha256 %s; interneuron_specificity_383_v2.tsv sha256 %s" % tuple(hashlib.sha256(open(os.path.join(OUT, f), "rb").read()).hexdigest()[:12] for f in ("hgnc_complete_set.txt", "recurrence_genes.tsv", "interneuron_specificity_383_v2.tsv")),
        "# genes: the 383 genes of recurrence_genes.tsv (the 8 extra genes MT1X SERPINA3 CLDN11 IL33R HILPDA KCNA1 ATP1A3 IFITM1 are already among the 383); matched by Ensembl id via HGNC first, current symbol in the h5ad var index second; see found_by; IL33R is the current HGNC symbol (HGNC:5998; previous symbol IL1RL1); no alias mapping was needed for any gene",
        "# class mapping (SEA-AD Subclass -> fine class): " + "; ".join("%s=%s" % (c, ",".join(s for s, v in MAP.items() if v == c)) for c in FINE_ORDER),
        "# unmapped Subclass labels (nuclei left out): Sst Chodl (%d nuclei; not in the requested mapping; a sensitivity run that adds it to SST is in seaad_agreement_with_psychad.tsv)" % int(z["nuc"].reshape(len(z["donors"]), -1)[:, list(z["subclasses"]).index("Sst Chodl")].sum()),
        "# Immune = SEA-AD label for microglia/PVM and a few other immune cells; used as Micro. 'VLMC & Perivascular' (VLMC + pericyte/SMC) and Endothelial = Vascular. PsychAD left PVM and Adaptive out; SEA-AD has no separate PVM label",
        "# method: per donor and class, summed raw UMIs over all %d genes -> CPM over all genes of that pseudobulk; donor-class pseudobulks with >= %d nuclei; class mean CPM = mean over qualifying donors; spec_<class> = class mean CPM / sum of the 11 fine-class mean CPMs; spec_coarse_IN / EN over 7 coarse classes (IN = SST+PVALB+VIP+LAMP5+IN_other pooled per donor); no CPM >= 1 filter; n_classes_measured = fine classes with mean CPM > 0" % (36601, MIN_NUC),
        "# differences from the PsychAD table: different tissue (middle temporal gyrus, not DLPFC), nuclei sorted/loaded by the SEA-AD protocol (includes NeuN-enriched sorting, so neuronal classes are over-represented relative to tissue; class shares use class means so this does not enter the shares directly), different class labels (see mapping), donors with no dementia (not psychiatric controls), older donors",
    ]
    def don_line(d, ndf, nnf):
        return ["# donors: %d (%s); Age at Death mean %.1f, median %.0f, range %.0f-%.0f (%d with age recorded); sex %d female, %d male" % (len(d), d.cog.iloc[0] if d.cog.nunique() == 1 else "Dementia + No dementia", d.age.mean(), d.age.median(), d.age.min(), d.age.max(), d.age.notna().sum(), (d.sex == "Female").sum(), (d.sex == "Male").sum()),
                "# qualifying donors per fine class: " + ", ".join("%s %d" % (c, ndf[c]) for c in FINE_ORDER),
                "# nuclei in qualifying donor-class pseudobulks per fine class: " + ", ".join("%s %d" % (c, nnf[c]) for c in FINE_ORDER)]
    h1 = hdr_common[:] + don_line(dn, nd_f, nn_f) + ["# donor rule: Cognitive Status == 'No dementia'"]
    ha = hdr_common[:] + don_line(don, runs["all"][1], runs["all"][2]) + ["# donor rule: all donors (Dementia and No dementia)"]
    def write(tab, hdr, name):
        with open(os.path.join(OUT, name), "w") as fh:
            fh.write("\n".join(hdr) + "\n"); tab.round(5).to_csv(fh, sep="\t")
    write(tabs["nodem"], h1, "seaad_class_shares_391.tsv"); write(tabs["all"], ha, "seaad_class_shares_391_alldonors.tsv")
    # agreement
    def agree(a, b, label, ref_filter=None):
        g = [x for x in a.index if x in b.index and a.loc[x, "spec_coarse_IN"] == a.loc[x, "spec_coarse_IN"] and b.loc[x, "spec_coarse_IN"] == b.loc[x, "spec_coarse_IN"]]
        if ref_filter is not None: g = [x for x in g if ref_filter(x)]
        A, B = a.loc[g], b.loc[g]; r = []
        for col in ["spec_coarse_IN"] + ["spec_" + c for c in FINE_ORDER]:
            r.append(dict(comparison=label, quantity=col, n_genes=len(g), spearman=spearmanr(A[col], B[col])[0], pearson=pearsonr(A[col], B[col])[0], frac_same_top_fine_class="", frac_agree_coarse_IN_ge_0_5=""))
        sf = ["spec_" + c for c in FINE_ORDER]
        r.append(dict(comparison=label, quantity="top_fine_class_and_coarse_IN_ge_0.5", n_genes=len(g), spearman="", pearson="",
                      frac_same_top_fine_class=float((A[sf].idxmax(axis=1) == B[sf].idxmax(axis=1)).mean()),
                      frac_agree_coarse_IN_ge_0_5=float(((A.spec_coarse_IN >= .5) == (B.spec_coarse_IN >= .5)).mean())))
        return r
    pf = lambda x: ps.loc[x, "status"] == "in_filtered_table"
    R = agree(tabs["nodem"], ps, "SEA-AD no-dementia vs PsychAD HBCC controls (all genes with values in both)")
    R += agree(tabs["nodem"], ps, "SEA-AD no-dementia vs PsychAD HBCC controls (PsychAD in_filtered_table only)", pf)
    R += agree(tabs["all"], tabs["nodem"], "SEA-AD all donors vs SEA-AD no-dementia")
    R += agree(tabs["nodem_sens"], tabs["nodem"], "SEA-AD no-dementia with Sst Chodl added to SST vs primary")
    ag = pd.DataFrame(R)
    hh = hdr_common[:2] + ["# script: second_reference_seaad_2026-10-19/seaad_class_shares.py, run %s" % today,
        "# comparisons: Spearman and Pearson over genes with values in both tables; frac_same_top_fine_class = fraction of genes whose largest fine-class spec is the same class in both; frac_agree_coarse_IN_ge_0_5 = fraction of genes on the same side of spec_coarse_IN >= 0.5",
        "# PsychAD reference: interneuron_specificity_383_v2.tsv (95 HBCC control donors; 372 genes in the filtered table, 3 computed without the CPM >= 1 filter, 8 absent)", hdr_common[4]]
    with open(os.path.join(OUT, "seaad_agreement_with_psychad.tsv"), "w") as fh:
        fh.write("\n".join(hh) + "\n"); ag.to_csv(fh, sep="\t", index=False, float_format="%.4f")
    print(ag.to_string()); print(tabs["nodem"].loc[["MT1X", "SERPINA3", "CLDN11", "IL33R", "HILPDA", "KCNA1", "ATP1A3", "IFITM1"], ["ensembl_id", "found_by", "spec_coarse_IN", "max_class_mean_cpm"] + ["spec_" + c for c in FINE_ORDER]].to_string())
    print(nd_f, nn_f); print(dn.age.describe(), (dn.sex.value_counts()))
    print("peak_rss_GB %.2f runtime_stage2+ %.0fs" % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6, time.time() - t0))
main()
