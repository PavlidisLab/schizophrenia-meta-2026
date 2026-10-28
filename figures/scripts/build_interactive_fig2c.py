#!/usr/bin/env python3
"""Interactive version of Fig. 2c: one self-contained HTML file (data embedded as JSON, no network needed).

Reads (figures/data/fig02/ and fig07/): fig02_g_gene_by_study_long_80.tsv (calls, 383 genes x 80 studies),
fig02_f_effects_by_study.tsv.gz and fig02_f_effects_study_index.tsv (per-gene effects for 36 studies),
fig02_f_gene_names.tsv, fig02_b_gene_de_prior.tsv, fig02_b_gene_agonal_axis.tsv, fig02_b_gene_interneuron_specificity.tsv,
fig02_b_study_annotation.tsv, fig07/fig07_a2_a4_gene_scatter.tsv, fig01/fig01_b_corpus_studies.tsv.
Writes supplementary/sz-fig2c-interactive.html.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage, optimal_leaf_ordering
from scipy.spatial.distance import pdist

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_fig_agreement import CLUSTERS, cluster_group  # noqa: E402  (cohort cluster names and rule of Fig. 2c)

REPO = Path(__file__).resolve().parents[2]
FD = REPO / "figures" / "data"
OUT = REPO / "supplementary" / "sz-fig2c-interactive.html"
HOLD = set()  # studies whose effects are held back from the tints until checked (Maas 2020 checked 2026-10-19: Group2 is schizophrenia)
CALLS = ["not_tested", "tested_not_hit", "up", "down", "mixed", "hit_no_direction"]


def rd(p, **kw):
    return pd.read_csv(p, sep="\t", comment="#", **kw)


def clean(x, nd=3):
    return None if pd.isna(x) else round(float(x), nd)


def main():
    long = rd(FD / "fig02" / "fig02_g_gene_by_study_long_80.tsv", dtype={"study_id": str})
    genes_rank = long.drop_duplicates("gene").sort_values("rank_within_383")
    genes = genes_rank["gene"].tolist()
    ng = len(genes)
    names = rd(FD / "fig02" / "fig02_f_gene_names.tsv").set_index("gene")
    prior = rd(FD / "fig02" / "fig02_b_gene_de_prior.tsv").set_index("gene")
    ph = rd(FD / "fig07" / "fig07_a2_a4_gene_scatter.tsv").set_index("gene")
    ag = rd(FD / "fig02" / "fig02_b_gene_agonal_axis.tsv").set_index("gene")
    sp = rd(FD / "fig02" / "fig02_b_gene_interneuron_specificity.tsv").set_index("gene")
    sp2 = rd(FD / "fig02" / "fig02_b_gene_specificity_seaad.tsv").set_index("gene")  # genes absent from PsychAD

    # studies: the 72 with a list (as in the static figure) and the 8 with an empty list at FDR 5%
    ann = rd(FD / "fig02" / "fig02_b_study_annotation.tsv", dtype={"study_id": str}).set_index("study_id")
    corp = rd(FD / "fig01" / "fig01_b_corpus_studies.tsv", dtype={"study_id": str}).drop_duplicates("study_id").set_index("study_id")
    sids = list(dict.fromkeys(long["study_id"]))
    empty_nt = rd(FD / "fig02" / "fig02_h_empty_list_studies_n_tested.tsv", dtype={"study_id": str}).set_index("study_id")["n_tested"].to_dict()
    cname = {k: lab for k, lab, _ in CLUSTERS}
    ann["grp_name"] = [cname[g] for g in ann.apply(cluster_group, axis=1)]
    cp = corp.reset_index() if "study_id" not in corp.columns else corp
    cprim = rd(FD / "fig01" / "fig01_b_corpus_studies.tsv", dtype={"study_id": str})
    cprim = cprim[cprim["in_primary_population"] == "yes"]
    by_id = cprim.set_index("study_id")["n_total"]
    by_label = cprim.drop_duplicates("label").set_index("label")["n_total"]
    norm_label = lambda l: l.replace(" (PsychAD)", "").replace("Frohlich", "Fröhlich")

    def n_of(sid, lab):
        v = by_id.get(sid, by_id.get("sc:" + sid, by_label.get(norm_label(lab))))
        if v is None or pd.isna(v):
            v = corp["n_total"].get(sid)
        return None if v is None or pd.isna(v) else float(v)
    mat = np.zeros((len(sids), ng), dtype=np.int8)
    gi = {g: i for i, g in enumerate(genes)}
    si = {s: i for i, s in enumerate(sids)}
    code = {c: i for i, c in enumerate(CALLS)}
    for s, g, c in zip(long["study_id"], long["gene"], long["call"]):
        mat[si[s], gi[g]] = code[c]

    # effects (mean over the groups behind the study's list), matched to study ids by label
    ef = rd(FD / "fig02" / "fig02_f_effects_by_study.tsv.gz", usecols=["study_label", "gene", "effect", "p", "is_fig2c_group"])
    ix = rd(FD / "fig02" / "fig02_f_effects_study_index.tsv", dtype={"pmid": str})
    ef = ef[ef["is_fig2c_group"].astype(str).str.upper() == "TRUE"]
    grp = ef.groupby(["study_label", "gene"])
    eff, pv = grp["effect"].mean(), grp["p"].apply(lambda s: s.iloc[0] if len(s) == 1 else np.nan)
    ng_groups = grp.size()
    label_to_id = {}
    for sid in sids:
        lab = ann["label"].get(sid, corp["label"].get(sid))
        label_to_id[lab] = sid
    meta = ix.set_index("study_label")
    studies = []
    unmatched = []
    for sid in sids:
        lab = ann["label"].get(sid, corp["label"].get(sid))
        e = [None] * ng
        p = [None] * ng
        elab = lab if lab in meta.index else lab + " (own list)"
        if elab in HOLD:
            elab = None  # effect held back until the source is checked  # the effects index labels Pai 2019 "(own list)"
        has = elab is not None and elab in meta.index
        if has:
            for g in genes:
                if (elab, g) in eff.index:
                    e[gi[g]] = clean(eff[(elab, g)], 3)
                    pp = pv[(elab, g)]
                    p[gi[g]] = None if pd.isna(pp) else float(f"{pp:.2g}")
        c = corp.loc[sid] if sid in corp.index else None
        row = {
            "id": sid, "label": lab,
            "year": int(ann["year"].get(sid, c["year"] if c is not None else 0) or 0),
            "cluster": str(ann["grp_name"].get(sid, "no list")),
            "tissue": None if c is None or pd.isna(c["tissue_short"]) else c["tissue_short"],
            "bank": None if c is None or pd.isna(c["bank_label"]) else c["bank_label"],
            "tech": str(ann["technology"].get(sid, c["technology"] if c is not None else "")),
            "n": n_of(sid, lab),
            "ncase": clean(c["n_cases"], 0) if c is not None else None,
            "nctl": clean(c["n_controls"], 0) if c is not None else None,
            "nlist": clean(ann["n_list_genes"].get(sid, 0), 0),
            "ntested": clean(ann["n_tested"].get(sid, empty_nt.get(sid)), 0),
            "haslist": sid in ann.index,
        }
        if has:
            m = meta.loc[elab]
            row["eff"] = e
            row["p"] = p
            row["method"] = str(m["method"])
            row["source"] = str(m["source_kind"])
            row["groups"] = int(m["n_groups"])
            row["access"] = None if pd.isna(m["access_note"]) else str(m["access_note"])
        else:
            unmatched.append(lab)
        studies.append(row)
    print("studies with effects:", sum("eff" in s for s in studies), "of", len(studies))

    # genes
    G = []
    for g in genes:
        r = prior.loc[g]
        spr = sp.loc[g] if g in sp.index else None
        shares = None
        note = None if spr is None else str(spr["status"])
        if (spr is None or pd.isna(spr["spec_IN_fine_sum"])) and g in sp2.index and not pd.isna(sp2.loc[g, "spec_IN_fine_sum"]) and sp2.loc[g, "spec_IN_fine_sum"] is not None:
            r2 = sp2.loc[g]
            gl2 = sum(r2[k] for k in ("spec_Astro", "spec_Micro", "spec_OPC", "spec_Oligo"))
            shares = [clean(r2["spec_IN_fine_sum"], 3), clean(r2["spec_EN"], 3), clean(gl2, 3), clean(r2["spec_Vascular"], 3)]
            note = "second_reference"
        elif spr is not None and not pd.isna(spr["spec_IN_fine_sum"]):
            gl = sum(spr[k] for k in ("spec_Astro", "spec_Micro", "spec_OPC", "spec_Oligo"))
            shares = [clean(spr["spec_IN_fine_sum"], 3), clean(spr["spec_EN"], 3), clean(gl, 3), clean(spr["spec_Vascular"], 3)]
        G.append({
            "g": g, "name": str(names["approved_name"].get(g, "")),
            "nup": int((mat[:, gi[g]] == 2).sum()), "ndown": int((mat[:, gi[g]] == 3).sum()),
            "ph": clean(ph["z_ph_controls"].get(g), 2), "phq": clean(ph["q_ph_controls"].get(g), 4),
            "phc": None if pd.isna(ph["ph_confound_consistent"].get(g)) else str(ph["ph_confound_consistent"].get(g)).upper() == "TRUE",
            "ag": clean(ag["agonal_t_axis"].get(g), 2) if g in ag.index else None,
            "agq": clean(ag["q_agonal"].get(g), 4) if g in ag.index else None,
            "agc": None if g not in ag.index or pd.isna(ag["agonal_consistent"].get(g)) else str(ag["agonal_consistent"].get(g)).upper() == "TRUE",
            "top": str(r["top10pct"]) == "yes", "nde": int(r["n_de_total"]),
            "spec": shares,
            "specNote": note,
        })

    # clustering of rows (studies) and columns (genes): the direction of the calls (up +1, down -1) plus, for
    # genes a study measured but did not list, half the effect (log2 fold change / 0.5, limited to +-1) where held;
    # correlation distance, average linkage, leaves ordered to put similar neighbours together
    v = np.where(mat == 2, 1.0, np.where(mat == 3, -1.0, 0.0))
    for k, st in enumerate(studies):
        if "eff" in st:
            e = np.array([0.0 if x is None else x for x in st["eff"]])
            v[k] += 0.5 * np.clip(e / 0.5, -1, 1) * (mat[k] == 1)

    def order(x):
        x = x + 1e-6 * np.random.RandomState(0).randn(*x.shape)  # all-zero rows have no correlation
        d = pdist(x, metric="correlation")
        return [int(i) for i in leaves_list(optimal_leaf_ordering(linkage(d, "average"), d))]
    row_cluster = order(v)
    col_cluster = order(v.T)

    DATA = {"genes": G, "studies": studies, "calls": ["".join(str(x) for x in r) for r in mat.tolist()],
            "rowCluster": row_cluster, "colCluster": col_cluster}
    html = TEMPLATE.replace("__DATA__", json.dumps(DATA, separators=(",", ":"), allow_nan=False))
    OUT.write_text(html, encoding="utf-8")
    print(f"{OUT}: {OUT.stat().st_size / 1e6:.2f} MB; {len(G)} genes, {len(studies)} studies; unmatched effect labels: {unmatched[:5]}")


TEMPLATE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Recurrent genes by study</title>
<style>
:root{--bg:#fff;--fg:#1f2937;--sub:#6b7280;--grid:#e5e7eb;--card:#fff}
@media (prefers-color-scheme:dark){:root{--bg:#111827;--fg:#e5e7eb;--sub:#9ca3af;--grid:#374151;--card:#1f2937}}
body{margin:0;background:var(--bg);color:var(--fg);font:13px Arial,Helvetica,sans-serif}
header{padding:12px 16px;border-bottom:1px solid var(--grid)}
h1{font-size:15px;font-weight:600;margin:0 0 6px}
.controls{display:flex;flex-wrap:wrap;gap:14px;align-items:center;margin:6px 0}
.controls label{display:flex;gap:6px;align-items:center;color:var(--sub)}
select,input,button{font:inherit;color:var(--fg);background:var(--bg);border:1px solid var(--grid);border-radius:4px;padding:2px 6px}
.note{color:var(--sub);margin:4px 0 0}
#legend{display:flex;flex-wrap:wrap;gap:4px 14px;margin-top:8px}
#legend span{display:flex;align-items:center;gap:5px;color:var(--sub)}
#legend i{width:12px;height:12px;border:1px solid var(--grid);display:inline-block}
#app{display:flex;padding:8px 16px}
#left{flex:0 0 150px}
#pane{flex:1;overflow-x:auto;overflow-y:hidden;border-left:1px solid var(--grid);border-right:1px solid var(--grid)}
#right{flex:0 0 110px}
.lab{white-space:nowrap;overflow:hidden;text-overflow:ellipsis;color:var(--fg);font-size:11px;line-height:12px;height:12px}
.slab{cursor:pointer;color:var(--sub);text-align:right;padding-right:6px}
.slab:hover{color:var(--fg);text-decoration:underline}
#tip{position:fixed;display:none;max-width:340px;background:var(--card);border:1px solid var(--sub);border-radius:6px;padding:8px 10px;pointer-events:none;box-shadow:0 2px 8px rgba(0,0,0,.2);z-index:5;line-height:1.35}
#tip b{font-weight:600}#tip .s{color:var(--sub)}#tip hr{border:0;border-top:1px solid var(--grid);margin:5px 0}
</style></head><body>
<header>
<h1>The recurrent genes in each study</h1>
<details id="ui" open><summary style="cursor:pointer;color:var(--sub)">Controls and legend</summary>
<div class="controls">
<label>Studies <select id="rowOrder"><option value="waterfall">waterfall (studies listing the first gene first)</option><option value="cluster">clustered by similarity of results</option><option value="cohort">by cohort cluster</option><option value="n">by sample size</option><option value="year">by year</option><option value="gene">by the clicked gene</option><option value="dens">by how many genes are listed or have a fold change</option></select></label>
<label>Genes <select id="colOrder"><option value="rank">by recurrence</option><option value="cluster">clustered by similarity across studies</option><option value="study">by the clicked study</option></select></label>
<label>Minimum log<sub>2</sub> fold change to colour <input id="thr" type="range" min="0" max="0.5" step="0.01" value="0.05"><b id="thrv">0.05</b></label>
<label>Unlisted cells <select id="style"><option value="subtle">pale tint</option><option value="solid">full colour</option></select></label>
<label><input type="checkbox" id="showEmpty" checked> show the 8 studies with no listed genes</label>
<label>Hide studies with less than <input id="minD" type="range" min="0" max="100" step="1" value="0"><b id="minDv">0</b>% of the genes listed or with a fold change <span id="nshown" style="color:var(--sub)"></span></label>
<label>Find gene <input id="find" size="9" placeholder="symbol"></label>
<button id="reset">Reset order</button>
</div>
<div class="note">Each column is one of the 383 recurrent genes and each row one study. Clustering puts studies with the most similar results next to each other (the genes each lists up or down, and the fold changes where held), and likewise for genes across studies. Click a gene (its column or its name) to sort the studies by it, or click a study name (or shift-click a cell) to sort the genes by it; click again to reverse. Click a label on the left of the strips to sort the genes by that strip. Hover over a cell, strip or gene for details. Scroll sideways for more genes.</div>
<div id="legend"></div>
</details>
</header>
<div id="app"><div id="left"><canvas id="cl"></canvas></div><div id="pane"><canvas id="cv"></canvas></div><div id="right"><canvas id="cr"></canvas></div></div>
<div id="tip"></div>
<script id="data" type="application/json">__DATA__</script>
<script>
const D=JSON.parse(document.getElementById('data').textContent);
const G=D.genes,S=D.studies,NG=G.length,NS=S.length;
const CW=8,TOPH=56,STH=10,NSTRIP=4,LABH=46,PADT=4,LW=150,RW=110;
let RH=10;
const strips=[['studies up','nup'],['studies down','ndown'],['pH','ph'],['slow death','ag'],['interneuron','spec0'],['DE prior','nde']];
const UP='#2563eb',DOWN='#ef4444',GREY='#f3f4f6';
const CALLCOL=['#ffffff',GREY,UP,DOWN,'#7c3aed','#4b5563'];
const CALLTXT=['The study did not measure this gene','The study measured this gene and did not list it','Listed as higher in schizophrenia','Listed as lower in schizophrenia','Listed, with the direction differing between lists','Listed, with no direction stated'];
const cl={}; // cluster label colours
['CMC, Mount Sinai, Pitt, LIBD','Stanley','NSW / Victoria','Harvard / McLean','Pritzker','one-study cluster','several banks (no cluster)','no list'].forEach((k,i)=>cl[k]=['#10b981','#0891b2','#8b5cf6','#ec4899','#84cc16','#9ca3af','#1f2937','#d1d5db'][i]);
let showEmpty=true,minD=0;
const dens=new Array(NS).fill(0);
function recomputeDens(){const t=parseFloat(document.getElementById('thr').value);S.forEach((st,i)=>{let n=0;for(let g=0;g<NG;g++){const c=calls[i][g];if(c>=2)n++;else if(c===1&&st.eff&&st.eff[g]!=null&&Math.abs(st.eff[g])>=t)n++}dens[i]=n/NG*100})}
const visList=()=>[...Array(NS).keys()].filter(i=>(showEmpty||S[i].haslist)&&dens[i]>=minD);
let colOrd=[...Array(NG).keys()],rowOrd=visList();
const calls=D.calls.map(s=>Array.from(s,Number));
const rank=[...Array(NG).keys()];
let sortKey=null,sortDir=1,rowRev=false,rowGene=null,rowGeneDir=1,colStudy=null,colStudyDir=1;
function cellScore(si,gi){const c=calls[si][gi];if(c===2)return 2;if(c===3)return -2;if(c===4||c===5)return 0.5;if(c===0)return null;const st=S[si],e=st.eff?st.eff[gi]:null;if(e==null)return 0;return Math.abs(e)>=thr()?Math.max(-1,Math.min(1,e/0.5)):0}
const cv=document.getElementById('cv'),ctx=cv.getContext('2d'),cr=document.getElementById('cr'),rctx=cr.getContext('2d');
const tip=document.getElementById('tip');
let dpr=window.devicePixelRatio||1;
const cl2=document.getElementById('cl'),lctx=cl2.getContext('2d');
const MATT=TOPH+NSTRIP*STH+PADT*2;
let H=MATT+visList().length*RH+LABH;
function size(c,w,h){c.width=w*dpr;c.height=h*dpr;c.style.width=w+'px';c.style.height=h+'px';c.getContext('2d').setTransform(dpr,0,0,dpr,0,0)}
function sizeAll(){dpr=window.devicePixelRatio||1;const hh=document.querySelector('header').offsetHeight;RH=Math.max(5,Math.min(12,Math.floor((innerHeight-hh-MATT-LABH-30)/visList().length)));H=MATT+rowOrd.length*RH+LABH;size(cv,NG*CW,H);size(cr,RW,H);size(cl2,LW,H)}
sizeAll();addEventListener('resize',()=>{sizeAll();draw()});
function thr(){return parseFloat(document.getElementById('thr').value)}
function effColor(e,solid){const a=solid?1:0.08+0.14*Math.min(Math.abs(e)/0.5,1);const t=e>0?[37,99,235]:[239,68,68];const b=[243,244,246];const m=b.map((x,i)=>Math.round((1-a)*x+a*t[i]));return `rgb(${m})`}
function strength(k){return strips[k]}
function phCls(g){if(g.ph==null)return 0;const sig=g.phq!=null&&g.phq<0.05;return g.phc?(sig?3:2):1}
function agCls(g){if(g.ag==null||g.agc==null)return 0;return g.agc?(g.agq<0.05?3:2):1}
const PHC=['#ffffff','#9ca3af','#fcd9a0','#f59e0b'],AGC=['#ffffff','#9ca3af','#e2a9e6','#a21caf'];
function draw(){
 document.getElementById('nshown').textContent='('+rowOrd.length+' of '+NS+' studies shown)';
 ctx.clearRect(0,0,NG*CW,H);
 const t=thr(),solid=document.getElementById('style').value==='solid';
 // bars: studies up above, down below
 const mx=Math.max(...G.map(g=>Math.max(g.nup,g.ndown)));
 const mid=TOPH*0.6;
 ctx.fillStyle=getComputedStyle(document.body).color;ctx.fillRect(0,mid,NG*CW,0.5);
 colOrd.forEach((gi,x)=>{const g=G[gi];ctx.fillStyle=UP;const hu=g.nup/mx*mid*0.95;ctx.fillRect(x*CW+1,mid-hu,CW-2,hu);ctx.fillStyle=DOWN;const hd=g.ndown/mx*(TOPH-mid)*0.95;ctx.fillRect(x*CW+1,mid,CW-2,hd)});
 // strips
 const y0=TOPH+PADT;
 colOrd.forEach((gi,x)=>{const g=G[gi];
  ctx.fillStyle=PHC[phCls(g)];ctx.fillRect(x*CW,y0,CW,STH-1);
  ctx.fillStyle=AGC[agCls(g)];ctx.fillRect(x*CW,y0+STH,CW,STH-1);
  const s=g.spec?g.spec[0]:null;ctx.fillStyle=s==null?'#fff':`rgba(15,118,110,${Math.min(s/0.6,1)})`;ctx.fillRect(x*CW,y0+2*STH,CW,STH-1);
  if(g.top){ctx.fillStyle=getComputedStyle(document.body).color;ctx.beginPath();ctx.arc(x*CW+CW/2,y0+3*STH+STH/2-1,2.2,0,6.3);ctx.fill()}});
 // matrix
 rowOrd.forEach((si,y)=>{const st=S[si];colOrd.forEach((gi,x)=>{let c=calls[si][gi];let col=CALLCOL[c];
   if(c===1&&st.eff){const e=st.eff[gi];if(e!=null&&Math.abs(e)>=t&&t>=0)col=effColor(e,solid)}
   ctx.fillStyle=col;ctx.fillRect(x*CW,MATT+y*RH,CW,RH)})});
 // gene labels
 ctx.fillStyle=getComputedStyle(document.body).color;ctx.font='8px Arial';ctx.textBaseline='middle';
 colOrd.forEach((gi,x)=>{ctx.save();ctx.translate(x*CW+CW/2,MATT+rowOrd.length*RH+4);ctx.rotate(Math.PI/2);ctx.fillText(G[gi].g,0,0);ctx.restore()});
 // right bars: genes tested and listed per study
 rctx.clearRect(0,0,RW,H);
 const bw=RW-30;const lg=v=>Math.log10(Math.max(v,1))/4.3*bw;
 rowOrd.forEach((si,y)=>{const st=S[si],yy=MATT+y*RH;
  if(st.ntested){rctx.fillStyle='#9ca3af';rctx.fillRect(4,yy+1,lg(st.ntested),RH-2);rctx.fillStyle='#1f2937';rctx.fillRect(4,yy+1,lg(st.nlist),RH-2)}
  else{rctx.fillStyle='#9ca3af';rctx.font=Math.min(9,RH)+'px Arial';rctx.textBaseline='middle';rctx.fillText('*',4,yy+RH/2)}});
 // log axis under the bars, with white lines across them
 const yb=MATT+rowOrd.length*RH;rctx.strokeStyle='#ffffff';rctx.lineWidth=1;rctx.fillStyle='#6b7280';rctx.font='8px Arial';rctx.textBaseline='top';rctx.textAlign='center';
 [1,10,100,1000,10000].forEach((v,k)=>{const x=4+lg(v);if(k>0){rctx.beginPath();rctx.moveTo(x,MATT);rctx.lineTo(x,yb);rctx.stroke()}rctx.fillText(v>=1000?(v/1000)+'k':String(v),x,yb+3)});
 rctx.textAlign='left';rctx.textBaseline='alphabetic';rctx.fillText('genes tested (grey)',4,MATT-14);rctx.fillText('and listed (dark)',4,MATT-5);
 buildLeft();
 if(typeof drawHL==='function')drawHL();
}
function buildLeft(){
 const fg=getComputedStyle(document.body).color,sub='#6b7280';
 lctx.clearRect(0,0,LW,H);lctx.textBaseline='middle';
 const y0=TOPH+PADT;
 lctx.font='9px Arial';lctx.textAlign='right';
 lctx.fillStyle=sub;lctx.fillText('studies up / down'+(sortKey==='nup'?(sortDir>0?' \u25BE':' \u25B4'):''),LW-4,TOPH/2);
 [['pH association','ph'],['slow death','ag'],['interneuron specificity','spec0'],['DE prior (dot)','nde']].forEach(([n,k],i)=>{lctx.fillStyle=sub;lctx.fillText(n+(sortKey===k?(sortDir>0?' \u25BE':' \u25B4'):''),LW-4,y0+i*STH+STH/2-0.5)});
 lctx.textAlign='left';lctx.font=Math.min(9,RH-1)+'px Arial';
 rowOrd.forEach((si,y)=>{const st=S[si],yy=MATT+y*RH;lctx.fillStyle=cl[st.cluster]||'#ccc';lctx.fillRect(2,yy,8,RH);
  const gv=st.n?Math.round(255*(1-Math.min(Math.max((Math.log10(st.n)-1)/1.5,0),1))):255;lctx.fillStyle=st.n?`rgb(${gv},${gv},${gv})`:'#fff';lctx.fillRect(11,yy,8,RH);lctx.fillStyle=fg;
  let t=st.label;while(t.length>3&&lctx.measureText(t).width>LW-26)t=t.slice(0,-1);lctx.fillText(t,23,yy+RH/2+0.5)});
 lctx.font='8px Arial';lctx.fillStyle=sub;lctx.textAlign='left';lctx.textBaseline='middle';lctx.save();lctx.translate(6,MATT-2);lctx.rotate(-Math.PI/2);lctx.fillText('cohort',0,0);lctx.restore();lctx.save();lctx.translate(15,MATT-2);lctx.rotate(-Math.PI/2);lctx.fillText('N',0,0);lctx.restore();
}
function key(g,k){if(k==='spec0')return g.spec?g.spec[0]:null;return g[k]}
function sortBy(k){if(sortKey===k)sortDir=-sortDir;else{sortKey=k;sortDir=1}
 document.getElementById('colOrder').value='rank';
 colOrd=[...Array(NG).keys()].sort((a,b)=>{const x=key(G[a],k),y=key(G[b],k);if(x==null&&y==null)return a-b;if(x==null)return 1;if(y==null)return -1;return sortDir*(y-x)||a-b});
 applyRows();draw()}
function applyCols(){const v=document.getElementById('colOrder').value;sortKey=null;
 if(v==='cluster')colOrd=D.colCluster.slice();
 else if(v==='study'&&colStudy!=null){colOrd=[...Array(NG).keys()].sort((a,b)=>{const x=cellScore(colStudy,a),y=cellScore(colStudy,b);if(x==null&&y==null)return a-b;if(x==null)return 1;if(y==null)return -1;return colStudyDir*(y-x)||a-b})}
 else colOrd=[...Array(NG).keys()]}
function applyRows(){const v=document.getElementById('rowOrder').value;
 if(v==='cluster'){const ok=new Set(visList());rowOrd=D.rowCluster.filter(i=>ok.has(i))}
 else if(v==='cohort')rowOrd=visList().sort((a,b)=>(S[a].cluster>S[b].cluster?1:S[a].cluster<S[b].cluster?-1:0)||a-b);
 else if(v==='n')rowOrd=visList().sort((a,b)=>(S[b].n||0)-(S[a].n||0));
 else if(v==='gene'&&rowGene!=null){rowOrd=visList().sort((a,b)=>{const x=cellScore(a,rowGene),y=cellScore(b,rowGene);if(x==null&&y==null)return a-b;if(x==null)return 1;if(y==null)return -1;return rowGeneDir*(y-x)||a-b})}
 else if(v==='dens')rowOrd=visList().sort((a,b)=>dens[b]-dens[a]||a-b);
 else if(v==='year')rowOrd=visList().sort((a,b)=>S[a].year-S[b].year||a-b);
 else{const r=c=>c===2?3:c===3?2:(c===4||c===5)?1:0;
  rowOrd=visList().sort((a,b)=>{for(const gi of colOrd){const d=r(calls[b][gi])-r(calls[a][gi]);if(d)return d}return (S[b].nlist||0)-(S[a].nlist||0)})}
 if(rowRev)rowOrd.reverse()}
function sortRows(v){const sel=document.getElementById('rowOrder');rowRev=(sel.value===v)?!rowRev:false;sel.value=v;applyRows();draw()}
document.getElementById('rowOrder').onchange=()=>{rowRev=false;applyRows();draw()};
document.getElementById('colOrder').onchange=()=>{applyCols();applyRows();draw()};
function sortRowsByGene(gi){const sel=document.getElementById('rowOrder');rowGeneDir=(sel.value==='gene'&&rowGene===gi)?-rowGeneDir:1;rowGene=gi;rowRev=false;sel.value='gene';applyRows();hl=gi;draw();drawHL()}
function sortGenesBy(si){const sel=document.getElementById('colOrder');colStudyDir=(sel.value==='study'&&colStudy===si)?-colStudyDir:1;colStudy=si;sel.value='study';applyCols();applyRows();draw()}
document.getElementById('thr').oninput=e=>{document.getElementById('thrv').textContent=(+e.target.value).toFixed(2);recomputeDens();applyRows();sizeAll();draw()};
document.getElementById('minD').oninput=e=>{minD=+e.target.value;document.getElementById('minDv').textContent=minD;applyRows();sizeAll();draw()};
document.getElementById('style').onchange=draw;
document.getElementById('showEmpty').onchange=e=>{showEmpty=e.target.checked;applyRows();sizeAll();draw()};
document.getElementById('reset').onclick=()=>{document.getElementById('rowOrder').value='waterfall';document.getElementById('colOrder').value='rank';sortKey=null;applyCols();applyRows();draw()};
document.getElementById('find').onkeydown=e=>{if(e.key!=='Enter')return;const q=e.target.value.trim().toUpperCase();const gi=G.findIndex(g=>g.g===q);if(gi<0)return;const x=colOrd.indexOf(gi);const p=document.getElementById('pane');p.scrollLeft=Math.max(0,x*CW-p.clientWidth/2);hl=gi;draw();drawHL()};
let hl=null;
function drawHL(){if(hl==null)return;const x=colOrd.indexOf(hl);ctx.strokeStyle='#f59e0b';ctx.lineWidth=1.5;ctx.strokeRect(x*CW,TOPH,CW,H-TOPH-LABH+LABH)}
const dim=(n,d)=>n==null?'not recorded':n.toLocaleString();
function geneCard(g){
 const q=v=>v==null?'':(v<0.001?'q < 0.001':'q = '+v.toFixed(3));
 let h=`<b>${g.g}</b> <span class="s">${g.name}</span><hr>`;
 h+=`Listed as higher by ${g.nup} studies and lower by ${g.ndown}.<br>`;
 h+=g.ph==null?'<span class="s">pH association: not measured</span><br>':`pH association in controls: z = ${g.ph} (${q(g.phq)}); ${g.phc?'goes the way low pH goes in schizophrenia':'goes the other way'}<br>`;
 h+=g.ag==null?'<span class="s">Slow-death association: not measured</span><br>':`Slow against fast death (GTEx): t = ${g.ag} (${q(g.agq)}); ${g.agc?'goes the way slow death goes in schizophrenia':'goes the other way'}<br>`;
 h+=g.top?'In the top 10% of genes most often differentially expressed across Gemma experiments.<br>':'';
 if(g.spec){const nm=['interneurons','excitatory neurons','glia','vascular'];const col=['#0f766e','#2563eb','#a16207','#9ca3af'];
  h+='<div style="margin-top:4px;display:flex;gap:8px;align-items:center">'+g.spec.map((s,i)=>{if(s==null)return `<span style="display:flex;flex-direction:column;align-items:center;font-size:10px;color:#9ca3af"><svg width="34" height="34"></svg>${nm[i].split(' ')[0]} n/a</span>`;const r=3+Math.sqrt(s)*14;return `<span title="${nm[i]}" style="display:flex;flex-direction:column;align-items:center;font-size:10px"><svg width="34" height="34"><circle cx="17" cy="17" r="${r}" fill="${col[i]}" opacity=".8"/></svg>${nm[i].split(' ')[0]} ${(s*100).toFixed(0)}%</span>`}).join('')+'</div><span class="s">Share of the gene’s mean expression across cell classes ('+(g.specNote==='second_reference'?'SEA-AD middle temporal gyrus, 42 donors without dementia; the gene is not in PsychAD':'PsychAD controls'+(g.specNote&&g.specNote.indexOf('below_filter')>=0?'; expression is very low, so the shares are noisy':''))+')</span>'}
 else h+='<span class="s">Cell-type share: '+(g.specNote==='not_in_psychad_gene_list'?'this gene is not in the PsychAD gene list and is not in the second reference':'not available')+'</span>';
 return h}
function studyCard(st){
 let h=`<b>${st.label}</b> <span class="s">${st.year||''}</span><br>`;
 const bits=[st.tissue,st.bank,st.tech].filter(x=>x&&x!=='nan');if(bits.length)h+=`<span class="s">${bits.join(' · ')}</span><br>`;
 if(st.n)h+=`${st.ncase!=null?st.ncase+' cases and '+st.nctl+' controls':st.n+' subjects'}<br>`;
 h+=st.haslist?`${st.nlist} genes listed of ${st.ntested!=null?st.ntested.toLocaleString():'?'} tested`:'No genes listed at FDR 5%';
 return h}
function cellCard(si,gi){const st=S[si],g=G[gi],c=calls[si][gi];
 let h=studyCard(st)+'<hr>'+`<b>${g.g}</b> <span class="s">${g.name}</span><br>${CALLTXT[c]}`;
 if(st.eff&&c!==0){const e=st.eff[gi];if(e!=null){const p=st.p[gi];h+=`<br>log<sub>2</sub> fold change ${e>0?'+':''}${e.toFixed(3)}`+(p!=null?` (P = ${p})`:'')+`<br><span class="s">${st.source==='stage'||st.source==='raw_supplement'?'our reanalysis':st.source==='gemma'?'Gemma':'published table'}${st.groups>1?', mean of '+st.groups+' groups':''}; ${st.method}</span>`}}
 else if(!st.eff&&c!==0)h+='<br><span class="s">No table of all genes is held for this study.</span>';
 return h}
function show(html,ev){tip.innerHTML=html;tip.style.display='block';const w=tip.offsetWidth,h=tip.offsetHeight;let x=ev.clientX+14,y=ev.clientY+14;if(x+w>innerWidth-8)x=ev.clientX-w-14;if(y+h>innerHeight-8)y=Math.max(8,innerHeight-h-8);tip.style.left=x+'px';tip.style.top=y+'px'}
cv.onmousemove=ev=>{const r=cv.getBoundingClientRect();const x=ev.clientX-r.left,y=ev.clientY-r.top;const col=Math.floor(x/CW);if(col<0||col>=NG)return tip.style.display='none';const gi=colOrd[col];
 if(y>=MATT&&y<MATT+rowOrd.length*RH){show(cellCard(rowOrd[Math.floor((y-MATT)/RH)],gi),ev)}
 else if(y<MATT||y>=MATT+rowOrd.length*RH){show(geneCard(G[gi]),ev)}else tip.style.display='none'};
cv.onmouseleave=()=>tip.style.display='none';
cv.onclick=ev=>{const r=cv.getBoundingClientRect();const x=ev.clientX-r.left,y=ev.clientY-r.top;const col=Math.floor(x/CW);if(col<0||col>=NG)return;const gi=colOrd[col];
 if(y>=MATT&&y<MATT+rowOrd.length*RH){if(ev.shiftKey||ev.altKey)sortGenesBy(rowOrd[Math.floor((y-MATT)/RH)]);else sortRowsByGene(gi)}else sortRowsByGene(gi)};
cr.onmousemove=ev=>{const r=cr.getBoundingClientRect();const y=ev.clientY-r.top;if(y>=MATT&&y<MATT+rowOrd.length*RH)show(studyCard(S[rowOrd[Math.floor((y-MATT)/RH)]]),ev);else tip.style.display='none'};
cr.onmouseleave=()=>tip.style.display='none';
cl2.onmousemove=ev=>{const r=cl2.getBoundingClientRect();const y=ev.clientY-r.top;const i=Math.floor((y-MATT)/RH);if(i>=0&&i<rowOrd.length)show(studyCard(S[rowOrd[i]]),ev);else tip.style.display='none';cl2.style.cursor=(y<MATT&&y>0)?'pointer':'default'};
cl2.onmouseleave=()=>tip.style.display='none';
cl2.onclick=ev=>{const r=cl2.getBoundingClientRect();const x=ev.clientX-r.left,y=ev.clientY-r.top;if(y>=MATT-34&&y<MATT&&x<20){sortRows(x<10?'cohort':'n');return}
 if(y>=MATT&&y<MATT+rowOrd.length*RH){sortGenesBy(rowOrd[Math.floor((y-MATT)/RH)]);return}
 if(y<TOPH)sortBy('nup');else{const i=Math.floor((y-TOPH-PADT)/STH);const k=['ph','ag','spec0','nde'][i];if(k)sortBy(k)}};
const leg=[['#fff','not measured'],[GREY,'measured, not listed'],[UP,'up in SZ (hit list)'],[DOWN,'down in SZ (hit list)'],['#7c3aed','mixed'],['#4b5563','no direction stated'],['rgb(214,224,248)','not listed, higher (tint)'],['rgb(248,222,222)','not listed, lower (tint)']];
document.getElementById('legend').innerHTML=leg.map(([c,t])=>`<span><i style="background:${c}"></i>${t}</span>`).join('')+'<span><i style="background:#f59e0b"></i>pH: low-pH direction, q &lt; 0.05 (pale: not significant; grey: opposite)</span><span><i style="background:#a21caf"></i>slow death: same direction, q &lt; 0.05</span><span><i style="background:#0f766e"></i>interneuron share of expression</span>';
const cll=Object.keys(cl).map(k=>`<span><i style="background:${cl[k]}"></i>${k}</span>`).join('');
document.getElementById('legend').innerHTML+='<span style="flex-basis:100%"></span><span style="color:var(--fg)">Left strips: cohort cluster</span>'+cll+'<span style="color:var(--fg)">sample size (cases plus controls, log)</span><span><i style="background:linear-gradient(90deg,#fff,#000);width:60px"></i>10 to 300</span>';
document.getElementById('nshown').textContent='';
document.getElementById('ui').addEventListener('toggle',()=>{sizeAll();draw()});
recomputeDens();applyCols();applyRows();draw();
</script></body></html>
"""

if __name__ == "__main__":
    main()
