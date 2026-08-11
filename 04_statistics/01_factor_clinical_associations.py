# -*- coding: utf-8 -*-
"""Manuscript step 2: model postprocess -> supplementary tables + association stats."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import config as cfg
import os, json, re
import h5py
import numpy as np
import pandas as pd
from scipy import stats

DIR = str(cfg.PROJECT_DIR)
MS = os.path.join(DIR, "manuscript")
TB = os.path.join(MS, "tables")

f = h5py.File(os.path.join(MS, "mofa_ms_K8.hdf5"), "r")
samples = [s.decode() for s in f["samples/all"][:]]
views = [v.decode() for v in f["views/views"][:]]
Z = pd.DataFrame(f["expectations/Z/all"][:].T, index=samples,
                 columns=[f"Factor{i+1}" for i in range(f["expectations/Z/all"].shape[0])])
W = {}
for v in views:
    feats = [x.decode() for x in f[f"features/{v}"][:]]
    W[v] = pd.DataFrame(f[f"expectations/W/{v}"][:].T, index=feats, columns=Z.columns)
r2 = pd.DataFrame(f["variance_explained/r2_per_factor/all"][:], index=views, columns=Z.columns)
if r2.values.max() <= 1.0:
    r2 *= 100
f.close()

meta = pd.read_csv(os.path.join(DIR, "sample_metadata.tsv"), sep="\t").set_index("patient").loc[samples]
cmpd = pd.read_csv(os.path.join(DIR, "metabolome_compound_names.tsv"), sep="\t", index_col=0)["compound_name"].astype(str).to_dict()

# ---- clinical: CRP / calprotectin from SAAデータ.xlsx ----
import openpyxl
canon_re = re.compile(r"(?:[gn]om[-_ ]?)?(sap|saga|asic|okh|okn)[-_ ]?0*(\d+)", re.I)
def canon(s):
    m = canon_re.search(str(s).lower())
    if not m:
        return None
    fam, num = m.group(1), int(m.group(2))
    if fam == "okh":
        fam = "okn"
    return {"sap": f"NOMSap_{num:03d}", "saga": f"Saga_{num:03d}",
            "asic": f"ASIC{num:03d}", "okn": f"OKN{num:03d}"}[fam]

wb = openpyxl.load_workbook(str(cfg.CLINICAL_XLSX), read_only=True, data_only=True)
ws = wb.worksheets[0]
rows = list(ws.iter_rows(values_only=True))
hdr = [str(x) if x is not None else "" for x in rows[0]]
def col(name):
    return hdr.index(name)
clin = {}
for r in rows[2:]:
    pid = canon(r[0]) if r[0] else None
    if pid is None:
        continue
    def num(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return np.nan
    clin[pid] = {"CRP": num(r[col("CRP")]), "calprotectin": num(r[col("便中カルプロテクチン")])}
clin = pd.DataFrame(clin).T
meta = meta.join(clin, how="left")
print("CRP available:", int(meta["CRP"].notna().sum()), "| calprotectin:", int(meta["calprotectin"].notna().sum()))

# ---- factor scores table (S7) ----
Z.join(meta[["diagnosis"]]).to_csv(os.path.join(TB, "TableS7_factor_scores.csv"))
r2.round(2).to_csv(os.path.join(TB, "TableS5_variance_explained_K8.csv"))

# ---- factor associations ----
rows_a = []
for fac in Z.columns:
    z = Z[fac]
    uc = z[meta.diagnosis == "UC"]; cd = z[meta.diagnosis == "CD"]
    p_dx = stats.mannwhitneyu(uc, cd).pvalue
    r_crp, p_crp = stats.spearmanr(z, meta["CRP"], nan_policy="omit")
    r_cal, p_cal = stats.spearmanr(z, meta["calprotectin"], nan_policy="omit")
    mb = meta["microbiome_batch"].astype(str)
    gb_ = [z[mb == g] for g in mb.unique() if (mb == g).sum() > 2]
    p_mb = stats.kruskal(*gb_).pvalue if len(gb_) > 1 else np.nan
    mtb = meta["metabolome_batch"].astype(str)
    gb2 = [z[mtb == g] for g in mtb.unique() if (mtb == g).sum() > 2]
    p_mtb = stats.kruskal(*gb2).pvalue if len(gb2) > 1 else np.nan
    rows_a.append({"factor": fac, "p_UCvsCD": p_dx, "rho_CRP": r_crp, "p_CRP": p_crp,
                   "rho_calprotectin": r_cal, "p_calprotectin": p_cal,
                   "p_micro_batch": p_mb, "p_metab_batch": p_mtb})
assoc = pd.DataFrame(rows_a)
for c in ["p_UCvsCD", "p_CRP", "p_calprotectin"]:
    ranks = assoc[c].rank()
    assoc[c + "_fdr"] = (assoc[c] * len(assoc) / ranks).clip(upper=1).round(4)
assoc.round(4).to_csv(os.path.join(TB, "TableS4_factor_clinical_associations.csv"), index=False)
print(assoc.round(3).to_string(index=False))

# ---- top-15 loadings tables per dominant view factor ----
dom = {v: r2.loc[v].idxmax() for v in views}
print("dominant factor per view:", dom)
for v, name in [("microbiome", "TableS1_microbiome_loadings"),
                ("mutation", "TableS2_mutation_loadings"),
                ("metabolome", "TableS3_metabolome_loadings"),
                ("transcriptome", "TableS8_transcriptome_loadings")]:
    fac = dom[v]
    w = W[v][fac]
    t = w.reindex(w.abs().sort_values(ascending=False).head(15).index)
    lab = [cmpd.get(str(i), str(i)) if v == "metabolome" else str(i) for i in t.index]
    pd.DataFrame({"feature": t.index, "label": lab, "factor": fac, "weight": t.values}).to_csv(
        os.path.join(TB, f"{name}.csv"), index=False)

# ---- factor correlation matrix ----
Z.corr().round(3).to_csv(os.path.join(TB, "factor_correlations.csv"))
cmax = Z.corr().abs().where(~np.eye(len(Z.columns), dtype=bool)).max().max()
print("max |r| between factors:", round(cmax, 3))

# ---- complete-case counts ----
n3 = int((meta.has_mutation & meta.has_microbiome & meta.has_metabolome).sum())
n4 = int(meta[["has_mutation", "has_microbiome", "has_metabolome", "has_transcriptome"]].all(axis=1).sum())
summary = {
    "n_views_data": {v: int(pd.read_csv(os.path.join(DIR, fn), sep="\t", index_col=0).notna().any().sum())
                     for v, fn in [("mutation", "view_mutation_gene_rare.tsv"), ("microbiome", "view_microbiome_clr.tsv"),
                                   ("metabolome", "view_metabolome_log2.tsv"), ("transcriptome", "view_transcriptome_logcpm.tsv")]},
    "complete_3omics": n3, "complete_4omics": n4,
    "diagnosis": meta.diagnosis.value_counts().to_dict(),
    "dominant_factor_per_view": dom,
    "max_factor_correlation": float(round(cmax, 3)),
    "r2_per_factor": {v: r2.loc[v].round(2).to_dict() for v in views},
}
meta.to_csv(os.path.join(MS, "meta_with_clinical.tsv"), sep="\t")
with open(os.path.join(TB, "stats_step2.json"), "w", encoding="utf-8") as fh:
    json.dump(summary, fh, indent=1, ensure_ascii=False)
print(json.dumps(summary, indent=1))
