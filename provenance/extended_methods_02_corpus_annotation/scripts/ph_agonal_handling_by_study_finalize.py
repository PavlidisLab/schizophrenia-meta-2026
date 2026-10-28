#!/usr/bin/env python3
"""Final pH and agonal-state handling table for the 64 primary brain studies: the keyword classification (prelim_classification.tsv) with manual calls where the keyword pass and h2_studies.tsv (cov_ph) disagreed.

    python3 finalize.py -> ph_agonal_handling_64.tsv, ph_calls_summary.txt

basis: 'keyword + h2 agree' (rule-based class equals the h2 class by the mapping in classify_prelim.py; quote from the scan), 'reviewed' (I read the sentences of the study; the call and the quote are mine, h2 may differ), 'unread' (no extractable text on remote; h2 label or RA's cov_corrections note quoted).
Categories (RA brief): matched, adjusted covariate, sensitivity only, discussed as biology, discussed as confound without adjustment, silent; plus 'reported only' (a value or table column, no handling) and 'considered' (screened as candidate covariate).
Quotes carry file (under rerun-full/data/raw/) and line in the cached extraction (text_cache/<pmid>/)."""
import re, pandas as pd
import os, json
D = os.path.dirname(os.path.abspath(__file__)) + "/"
EXCL = {k for k in json.load(open(D + "../../corpus/corpus_exclusions.json")) if not k.startswith("_")}  # left the corpus (2026-10-10)
P = pd.read_csv(D + "prelim_classification.tsv", sep="\t", dtype=str).fillna(""); P = P[~P.pmid.isin(EXCL)]; H = pd.read_csv(D + "scan_hits.tsv", sep="\t", comment="#", dtype=str).fillna("")
O = {  # pmid: (call, note, regex selecting the quote among the study's pH hits)
 "16038679": ("matched", "samples balanced for sample pH; genes correlated with brain pH examined (pH contributed to variance of some genes, mostly independent of disease)", r"balanced across patient demographics|Brain pH contributed"),
 "18778695": ("sensitivity only", "global ANCOVA of pH and PMI on gene-expression variance; pH confound found in one cohort and gone after removing 4 subjects; correlations of pH with qPCR values", r"confounding effect of pH"),
 "18985160": ("adjusted covariate", "brain pH one of the fixed effects in the ANOVA with diagnosis; cases and controls also approximately matched on pH", r"fixed effects \(brain pH|major non-diagnosis"),
 "21091092": ("matched", "CSF pH (not brain tissue pH) did not differ", r"did not differ from controls"),
 "22558445": ("sensitivity only", "matched for age, PMI and pH; ANCOVA showed no pH effect; genes correlated with pH covaried", r"matched for age, PMI and pH|ANCOVA"),
 "22675524": ("reported only", "pH given in the demographics table; no handling sentence found in the extracted text (h2 says sensitivity: not confirmed here)", r"Brain pH 6\.\d"),
 "24169640": ("matched", "groups matched for brain pH; Pearson correlation of pH with expression examined (pH correlated with APOL1)", r"matched for descriptive variables|brain pH was additionally correlated"),
 "24594777": ("discussed as confound without adjustment", "Introduction and Discussion name tissue pH as a factor affecting results; pH only in the demographics", r"can affect results|lack of robustness"),
 "27668389": ("adjusted covariate", "pH in the covariate list of the differential-expression and eQTL models (h2 missed it)", r"covariates were used in the dif"),
 "29112193": ("matched", "groups did not differ in pH; pH tested as a potential mediator and not included as a covariate", r"did not differ in mean age|Tissue pH was not included"),
 "29463818": ("matched", "no significant difference in CNS pH between diagnostic groups", r"no significant differences in age"),
 "31174959": ("discussed as confound without adjustment", "RNA degradation as a confounder; adjusting for observed quality metrics (RIN, pH, mito rate) shown insufficient; quality-surrogate approach used instead", r"not fully addressed by adjusting|adequately control"),
 "31288836": ("reported only", "supplementary donor table has a pH column; no handling in the extracted text", r"^pH$"),
 "33184259": ("silent", "the only pH hit is primer text", r"zzz"),
 "33558674": ("considered as candidate covariate", "tissue pH among the technical and demographic variables screened for covariate selection (supplement)", r"tissue pH"),
 "34059796": ("adjusted covariate", "supplement: brain pH included as a covariate when correlated with the target gene; excluded as a covariate for macrophage markers (acidosis accompanies inflammation); reports lower brain pH in SZ and BD", r"included as covariates if a significant correlation|Brain pH was excluded"),
 "36223459": ("matched", "supplementary figure compares age, PMI, pH and neuroleptic treatment between SZ and controls", r"Comparison of age, PMI"),
 "41567988": ("not covariate in primary; sensitivity with pH", "supplement: brain pH not included as a covariate because pH is lower in SZ and acidosis is typical of inflamed tissue; separate comparisons including brain pH as a covariate (Table S7)", r"not included as a covariate|including brain pH as a covariate"),
 "42461238": ("adjusted covariate", "CNS pH and RIN as covariates in the primary analyses; groups also did not differ in CNS pH", r"covariates in this, our primary|no significant differences in sex ratio"),
}
UNREAD = {"21538462": ("unread", "no extractable text on remote (h2: reported_only)"), "39526318": ("matched", "text not extracted on remote; h2 cov_corrections.tsv (checked 2026-10-06 against the h2 text): 'pH matched between groups; covariates entered were the top two by variance (illness duration, age at onset); pH not among them for SZ'")}
MAP = {"adjusted": "adjusted covariate", "sensitivity": "sensitivity only", "matched": "matched", "confound_unadjusted": "discussed as confound without adjustment", "reported": "reported only", "silent": "silent"}
rows = []
for _, r in P.iterrows():
    k = r.pmid; hh = H[(H.pmid == k) & (H.pH == "True")]
    if k in UNREAD: call, note = UNREAD[k]; basis = "unread"; q = ""
    elif k in O:
        call, note, rx = O[k]; basis = "reviewed"; m = hh[hh.sentence.str.contains(rx, regex=True, case=False)]
        q = f"{m.iloc[0].file}:{m.iloc[0].line} | {m.iloc[0].sentence[:300]}" if len(m) else ""
    else:
        call, note, basis, q = MAP[r.scan_ph_class], "", "keyword + h2 agree", r.ph_quote
    rows.append({"pmid": k, "label": r.label, "year": r.year, "ph_call": call, "ph_basis": basis, "ph_note": note, "ph_quote_file_line": q, "h2_cov_ph": r.h2_cov_ph, "agonal_keyword_class": r.scan_agonal_class, "agonal_sentences": r.n_agonal_sentences,
                 "agonal_quote_file_line": r.agonal_quote, "h2_cov_agonal": r.h2_cov_agonal, "files_read": r.files_read, "files_total": r.files_total, "unread_files": r.unread_files})
T = pd.DataFrame(rows)
with open(D + "ph_agonal_handling_64.tsv", "w") as fh:
    fh.write("# finalize.py: pH and agonal-state handling in the 64 primary brain studies (recurrence.py 'published, brain'); pH calls: keyword scan (extract_and_scan.py) checked against analyses/h2/h2_studies.tsv cov_ph, manual review where they disagreed; agonal column = keyword pass only (agonal, manner / cause / mode of death, hypoxia, coma ... in text), not reviewed\n")
    T.to_csv(fh, sep="\t", index=False, line_terminator="\n")
c = T.ph_call.value_counts()
agstrict = T[T.agonal_quote_file_line != ""]
txt = [f"pH calls over {len(T)} studies: " + "; ".join(f"{k} {v}" for k, v in c.items()), f"basis: " + "; ".join(f"{k} {v}" for k, v in T.ph_basis.value_counts().items()),
       f"h2 cov_ph vs final: agree on {(T.apply(lambda r: {'adjusted':'adjusted covariate','matched_or_compared':'matched','sensitivity':'sensitivity only','reported_only':'reported only','not_mentioned':'silent'}.get(r.h2_cov_ph, '?') == r.ph_call, axis=1)).sum()} of {len(T)}",
       f"studies with an agonal / manner-of-death keyword sentence in their text: {int((T.agonal_sentences.astype(int) > 0).sum())}; with a keyword class other than silent: {int((T.agonal_keyword_class != 'silent').sum())}"]
open(D + "ph_calls_summary.txt", "w").write("\n".join(txt) + "\n"); print("\n".join(txt))
