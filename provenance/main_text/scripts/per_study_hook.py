# NEW 2026-10-17 (RA). Executed inside the namespace of other_disorders_176.py (line "exec(per_study_hook)" inserted before the Hagihara section; the script exits afterwards).
# One row per DE table or Gemma list, for the current set (OD_SETNAME): panels a, b, c. Same functions as the pooled rows (enrich_row, dir_rows, axis_rho), same matched sets (seed SEED):
#   panel a, b: first NPERM_LIST (2000) of the matched sets, as the stored per-list rows; panel c: first NPERM_AXIS (5000), as the stored axis rows.
import csv as _csv, re as _re
_SET = os.environ.get("OD_SETNAME", "study_brain")
_OUT = os.environ["PS_OUT"]; os.makedirs(_OUT, exist_ok=True)
# pH axis v1 (stored combined_controls_pH.tsv) alongside the v2 axis the run was started with (PH_Z_FILE)
_v1 = {}
for _l in open(os.path.join(AN, "confounds/results/combined_controls_pH.tsv")):
    if _l.startswith("#") or _l.startswith("gene\t"):
        continue
    _f = _l.rstrip("\n").split("\t")
    if _f[1] == "pH" and _f[0] in GI:
        try:
            _v1[_f[0]] = float(_f[3])
        except ValueError:
            pass
_AX = {"pH_v2": axes["control pH Z (low pH = negative)"], "pH_v1": (_v1, -1), "agonal": axes["GTEx agonal t (slow vs fast death)"]}
assert os.environ.get("PH_Z_FILE", "").endswith("axis_v2_controls.tsv.gz"), "run with the v2 axis"

# ---- pooled membership
_ad_cap = list(_keep)           # AD lists under the 18.3% DE cap (from ad_tables.py)
_member = collections.defaultdict(list)
for _dis in ("BD", "MDD", "ASD"):
    for _drop in (False, True):
        _wl = (("restricted", "restricted"), ("dropped", "all series"), ("all", "all series")) if (MODE209 and _dis in ("BD", "ASD")) else (("all", ""),)
        for _which, _ss in _wl:
            if _dis == "ASD" and _which != "restricted" and MODE209:
                continue
            _tp, _kp = pooled(_dis, _drop, _which)
            for _k in _kp:
                _member[_k.name].append(_tp.name)
for _k in ad_tabs:
    _member[_k.name].append(_all.name)
    if _k in _ad_cap:
        _member[_k.name].append(_cap.name)

# ---- metadata
_asd = {}
for _r in tsv(ASD_SEL):
    _asd[f"ASD__{_r['short_name']}__" + ("all" if _r["subset"] == "all" else _r["subset"].replace(" ", "_"))] = _r
_adm = {}
for _r in tsv(os.environ["AD_SEL"]):
    _adm[f"{_r['geo_series']}|{_r['result_set']}"] = _r
def _region_from_name(n):
    m = _re.findall(r"\(([^()]*)\)", n or "")
    return "; ".join(m) if m else ""
meta = []   # (table, dict)
for _t in tables:       # Gandal
    meta.append((_t, dict(disorder=_t.disorder, label=_t.name, region="meta-analysis (cortical regions pooled by the authors)", platform="microarray" if "array" in _t.name else "RNA-seq", n_samples="", pooled_into="none (single table)", independent="yes", shared_donor_rule="", same_series_as="")))
_seen = {}
for _L, _t in gem_rows:
    _dis = _L["disorder"]; _lid = _L["list_id"]
    if _dis == "ASD":
        _a = _asd[_lid]; _geo = _L["short_name"]; _n = _a["n_samples"]; _pl = _a["platforms"]; _reg = ""; _desc = _a["name"]
        _plat = _pl
    else:
        _geo = _L["geo_series"] if "geo_series" in _L else _L["short_name"]; _n = (int(_L["n_case"]) + int(_L["n_control"])) if _L["n_case"] not in ("", None) else ""
        _plat = (_L.get("platforms", "") + (" (microarray)" if str(_L.get("is_microarray")) == "True" else " (RNA-seq)")); _desc = _L.get("name", ""); _reg = ""
    _reg = _region_from_name(_desc)
    _key = (_dis, _geo)
    _same = _seen.get(_key, "")
    _seen.setdefault(_key, _lid)
    meta.append((_t, dict(disorder=_dis, label=f"{_geo} [{_lid}]", region=_reg if _reg else "see description", platform=_plat, n_samples=_n, pooled_into=" | ".join(_member.get(_t.name, [])),
                          independent="yes" if not _same else f"no (same GEO series as {_same})", shared_donor_rule=(("retain" if RETAIN.get(_lid) else "drop") if _lid in RETAIN else "not classified"), same_series_as=_same, description=_desc)))
for _t in ad_tabs:
    _r = [x for x in tsv(os.environ["AD_SEL"]) if x["geo_series"] == _t.geo and x["independent"] == "yes" and f"Gemma result set {x['result_set']}" in _t.source][0]
    meta.append((_t, dict(disorder="AD", label=_t.geo, region=_r["region_label"], platform=_r["platform_class"] + " (" + _r["platforms"] + ")", n_samples=_r["n_samples"], pooled_into=" | ".join(_member.get(_t.name, [])), independent=_r["independent"],
                          shared_donor_rule="", same_series_as="" if _r["same_geo_series_as_another_row"] == "no" else _r["same_geo_series_as_another_row"], description=_r["note"])))
_cols_meta = ["disorder", "study_or_list_label", "table", "region", "platform", "n_samples", "set", "pooled_into", "independent", "shared_donor_rule", "same_series_as", "description"]
def _base(t, m):
    return dict(disorder=m["disorder"], study_or_list_label=m["label"], table=t.name, region=m["region"], platform=m["platform"], n_samples=m["n_samples"], set=len(S176),
                pooled_into=m["pooled_into"], independent=m["independent"], shared_donor_rule=m["shared_donor_rule"], same_series_as=m["same_series_as"], description=m.get("description", ""))
A, B, C = [], [], []
for t, m in meta:
    for thr in (0.05, 0.10):
        r = enrich_row(t, thr, DRl, use_logit=False)
        d = _base(t, m); d.update(fdr=thr, n_set_measured=r["n_176_measured"], observed=r["obs_176_DE"], matched_mean=r["null_mean"], matched_sd=r["null_sd"],
                                  effect=round(r["obs_176_DE"] - r["null_mean"], 2), n_matched_sets=len(DRl), table_DE_fraction_q05=r["pool_DE_frac"] if thr == 0.05 else "")
        A.append(d)
    for r in dir_rows(t, 0.05, DRl):
        if r["set"] == "nominal p<0.05":
            continue
        d = _base(t, m); d.update(subset=r["set"], n_set_with_SZ_sign=r["n_176_with_SZ_sign_in_set"], concordant=r["concordant"], n_with_direction=r["n_176_with_SZ_sign_in_set"],
                                  share=r["frac_concordant"], matched_mean=r["null_frac_mean"], matched_sd=r["null_frac_sd"],
                                  effect=(round(r["frac_concordant"] - r["null_frac_mean"], 4) if r["frac_concordant"] != "" else ""), n_matched_sets=len(DRl))
        B.append(d)
    d = _base(t, m); R = {}
    for k, (ax, orient) in _AX.items():
        R[k] = axis_rho(t, ax, orient, DR, NPERM_AXIS)
    d.update(n_set_measured=R["pH_v2"][0] if R["pH_v2"] else "", n_matched_sets=NPERM_AXIS,
             rho_pH_v2=round(R["pH_v2"][1], 4), matched_mean_rho_pH_v2=round(R["pH_v2"][2], 4), matched_sd_rho_pH_v2=round(R["pH_v2"][3], 4), effect_pH_v2=round(R["pH_v2"][1] - R["pH_v2"][2], 4),
             rho_agonal=round(R["agonal"][1], 4), matched_mean_rho_agonal=round(R["agonal"][2], 4), matched_sd_rho_agonal=round(R["agonal"][3], 4), effect_agonal=round(R["agonal"][1] - R["agonal"][2], 4),
             rho_pH_v1=round(R["pH_v1"][1], 4), matched_mean_rho_pH_v1=round(R["pH_v1"][2], 4), matched_sd_rho_pH_v1=round(R["pH_v1"][3], 4), effect_pH_v1=round(R["pH_v1"][1] - R["pH_v1"][2], 4))
    C.append(d)
for nm, rows_ in (("a", A), ("b", B), ("c", C)):
    cols = list(rows_[0])
    with open(os.path.join(_OUT, f"ps_{nm}_{len(S176)}.tsv"), "w", newline="") as fh:
        w = _csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows_)
print("per-study rows: tables", len(meta), "a", len(A), "b", len(B), "c", len(C))
sys.exit(0)
