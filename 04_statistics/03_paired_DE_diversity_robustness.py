# -*- coding: utf-8 -*-
"""Manuscript step 5: additional analyses to substantiate an expanded Results section.
1) paired inflamed vs non-inflamed differential expression (n=21 pairs)
2) alpha diversity by diagnosis and MEFV E148Q
3) characterisation of the transcriptome-dominant Factor 1 (inflammation score)
4) K=6 vs K=8 factor-matching robustness
5) complete-case (four-omics, n=18) sensitivity model
6) detailed loadings for narrative
"""
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
res = {}

meta = pd.read_csv(os.path.join(MS, "meta_with_clinical.tsv"), sep="\t", index_col=0)
with h5py.File(os.path.join(MS, "mofa_ms_K8.hdf5"), "r") as f:
    samples = [s.decode() for s in f["samples/all"][:]]
    views = [v.decode() for v in f["views/views"][:]]
    Z = pd.DataFrame(f["expectations/Z/all"][:].T, index=samples,
                     columns=[f"Factor{i+1}" for i in range(f["expectations/Z/all"].shape[0])])
    W = {v: pd.DataFrame(f[f"expectations/W/{v}"][:].T,
                         index=[x.decode() for x in f[f"features/{v}"][:]], columns=Z.columns) for v in views}
meta = meta.loc[samples]

# ---------- 1) paired inflamed vs non-inflamed DE ----------
lc = pd.read_csv(os.path.join(DIR, "rna_logcpm_all_samples.tsv"), sep="\t", index_col=0)
inf_c = [c for c in lc.columns if c.endswith("__inf")]
nor_c = [c for c in lc.columns if c.endswith("__nor")]
pairs = sorted(set(c[:-5] for c in inf_c) & set(c[:-5] for c in nor_c))
pairs_in = [p for p in pairs if p in meta.index]
A = lc[[p + "__inf" for p in pairs_in]]
B = lc[[p + "__nor" for p in pairs_in]]
expressed = ((A.values >= 1).sum(1) >= 5) | ((B.values >= 1).sum(1) >= 5)
A, B = A[expressed], B[expressed]
lfc = (A.values - B.values).mean(1)
pv = np.array([stats.wilcoxon(A.values[i], B.values[i])[1] if np.any(A.values[i] != B.values[i]) else 1.0
               for i in range(A.shape[0])])
de = pd.DataFrame({"gene": A.index, "log2FC_inf_vs_nor": lfc, "p": pv})
de["q"] = (de.p * len(de) / de.p.rank()).clip(upper=1)
anno = pd.read_csv(cfg.GENE_ANNOTATION, sep="\t", index_col=0)
sym_map = {str(k): str(v) for k, v in anno["Gene_Symbol"].items()}
de["symbol"] = [sym_map.get(str(g), str(g)) for g in de.gene]
de = de.sort_values("p")
de.to_csv(os.path.join(TB, "TableS9_paired_DE_inf_vs_nor.csv"), index=False)
up = de[(de.q < 0.05) & (de.log2FC_inf_vs_nor > 0.5)]
dn = de[(de.q < 0.05) & (de.log2FC_inf_vs_nor < -0.5)]
res["paired_DE"] = {
    "n_pairs": len(pairs_in), "n_genes_tested": int(len(de)),
    "n_up_q05_lfc05": int(len(up)), "n_down_q05_lfc05": int(len(dn)),
    "top_up": up.head(12)[["symbol", "log2FC_inf_vs_nor", "q"]].round(3).to_dict("records"),
    "top_down": dn.head(12)[["symbol", "log2FC_inf_vs_nor", "q"]].round(3).to_dict("records"),
}
print(f"paired DE: {len(pairs_in)} pairs, {len(de)} genes, up={len(up)}, down={len(dn)}")
print("  top up:", [r["symbol"] for r in res["paired_DE"]["top_up"]])
print("  top down:", [r["symbol"] for r in res["paired_DE"]["top_down"]])

# inflammation score = mean logCPM of top-100 up genes, inflamed sample
top100 = up.head(100).gene.tolist() if len(up) >= 20 else de.head(100).gene.tolist()
inflam = lc.loc[top100, inf_c].mean()
inflam.index = [c[:-5] for c in inflam.index]
inflam = inflam.reindex(meta.index)

# ---------- 2) alpha diversity ----------
rel = pd.read_csv(os.path.join(MS, "genus_relabundance_35pts.tsv"), sep="\t", index_col=0)
p_ = rel.values
shannon = pd.Series(-(np.where(p_ > 0, p_ * np.log(p_), 0)).sum(0), index=rel.columns)
observed = pd.Series((rel > 0).sum(0), index=rel.columns)
alpha = pd.DataFrame({"shannon": shannon, "observed_genera": observed})
alpha = alpha.join(meta[["diagnosis", "MEFV_E148Q_dosage", "CRP", "microbiome_batch"]])
alpha.to_csv(os.path.join(TB, "TableS10_alpha_diversity.csv"))
uc = alpha[alpha.diagnosis == "UC"].shannon
cd = alpha[alpha.diagnosis == "CD"].shannon
p_dx = stats.mannwhitneyu(uc, cd).pvalue
car = alpha[alpha.MEFV_E148Q_dosage > 0].shannon
non = alpha[alpha.MEFV_E148Q_dosage == 0].shannon
p_mefv = stats.mannwhitneyu(car, non).pvalue
r_crp, p_crp = stats.spearmanr(alpha.shannon, alpha.CRP, nan_policy="omit")
res["alpha"] = {
    "n": len(alpha),
    "shannon_UC_median": round(float(uc.median()), 3), "shannon_CD_median": round(float(cd.median()), 3),
    "p_UCvsCD": round(float(p_dx), 4),
    "shannon_E148Q_carrier": round(float(car.median()), 3), "shannon_noncarrier": round(float(non.median()), 3),
    "p_MEFV": round(float(p_mefv), 4),
    "rho_CRP": round(float(r_crp), 3), "p_CRP": round(float(p_crp), 4),
    "observed_median": int(observed.median()),
}
print("alpha diversity:", res["alpha"])

# ---------- 3) Factor 1 characterisation ----------
f1 = {}
for fac in Z.columns:
    r, p = stats.spearmanr(Z[fac], inflam, nan_policy="omit")
    f1[fac] = (round(float(r), 3), round(float(p), 4))
res["factor_vs_inflammation_score"] = f1
print("factor vs mucosal inflammation score:", {k: v for k, v in f1.items() if v[1] < 0.1})

# also: factor vs alpha diversity
fa = {}
for fac in Z.columns:
    r, p = stats.spearmanr(Z.loc[alpha.index, fac], alpha.shannon)
    if p < 0.1:
        fa[fac] = (round(float(r), 3), round(float(p), 4))
res["factor_vs_shannon"] = fa
print("factor vs Shannon:", fa)

# ---------- 4) K=6 vs K=8 robustness ----------
with h5py.File(os.path.join(MS, "mofa_ms_K6.hdf5"), "r") as f:
    Z6 = pd.DataFrame(f["expectations/Z/all"][:].T, index=[s.decode() for s in f["samples/all"][:]],
                      columns=[f"K6_Factor{i+1}" for i in range(f["expectations/Z/all"].shape[0])])
C = pd.DataFrame(np.corrcoef(Z.T.values, Z6.loc[Z.index].T.values)[:8, 8:],
                 index=Z.columns, columns=Z6.columns).abs()
best = C.max(axis=1).round(3)
res["K6_vs_K8_matching"] = {"max_abs_r_per_K8_factor": best.to_dict(),
                            "n_factors_matched_r_gt_0.8": int((best > 0.8).sum())}
C.round(3).to_csv(os.path.join(TB, "TableS11_K6_vs_K8_factor_matching.csv"))
print("K6-K8 matching (max |r| per K8 factor):", best.to_dict())

# ---------- 5) complete-case sensitivity (four-omics, n=18) ----------
cc = meta[meta[["has_mutation", "has_microbiome", "has_metabolome", "has_transcriptome"]].all(axis=1)].index.tolist()
res["complete_case_n"] = len(cc)
views_files = {"mutation": "view_mutation_gene_rare.tsv", "microbiome": "view_microbiome_clr.tsv",
               "metabolome": "view_metabolome_log2.tsv", "transcriptome": "view_transcriptome_logcpm.tsv"}
mats, fnames = [], []
for name, fn in views_files.items():
    df = pd.read_csv(os.path.join(DIR, fn), sep="\t", index_col=0)[cc]
    mats.append([df.T.to_numpy(dtype=float)])
    fnames.append([str(i) for i in df.index])
from mofapy2.run.entry_point import entry_point
ent = entry_point()
ent.set_data_options(scale_views=True)
ent.set_data_matrix(mats, likelihoods=["gaussian"] * 4, views_names=list(views_files),
                    groups_names=["all"], samples_names=[cc], features_names=fnames)
ent.set_model_options(factors=8, spikeslab_weights=True, ard_weights=True)
ent.set_train_options(iter=2000, convergence_mode="slow", seed=42, verbose=False, freqELBO=5, quiet=True)
ent.build(); ent.run()
out_cc = os.path.join(MS, "mofa_ms_completecase.hdf5")
ent.save(out_cc)
with h5py.File(out_cc, "r") as f:
    r2cc = pd.DataFrame(f["variance_explained/r2_per_factor/all"][:],
                        index=[v.decode() for v in f["views/views"][:]])
    Zcc = pd.DataFrame(f["expectations/Z/all"][:].T, index=cc)
if r2cc.values.max() <= 1.0:
    r2cc *= 100
Ccc = pd.DataFrame(np.corrcoef(Z.loc[cc].T.values, Zcc.T.values)[:8, 8:], index=Z.columns).abs()
res["complete_case"] = {"r2_total_per_view": r2cc.sum(axis=1).round(2).to_dict(),
                        "max_abs_r_vs_main": Ccc.max(axis=1).round(3).to_dict()}
r2cc.round(2).to_csv(os.path.join(TB, "TableS12_completecase_variance.csv"))
print("complete-case n=18 r2:", res["complete_case"]["r2_total_per_view"])
print("complete-case factor recovery:", res["complete_case"]["max_abs_r_vs_main"])

# ---------- 6) loadings detail ----------
cmpd = pd.read_csv(os.path.join(DIR, "metabolome_compound_names.tsv"), sep="\t", index_col=0)["compound_name"].astype(str).to_dict()
det = {}
for v, fac in [("transcriptome", "Factor1"), ("metabolome", "Factor2"), ("microbiome", "Factor3"),
               ("transcriptome", "Factor3"), ("mutation", "Factor7"), ("mutation", "Factor4"),
               ("transcriptome", "Factor2")]:
    w = W[v][fac]
    t = w.reindex(w.abs().sort_values(ascending=False).head(10).index)
    lab = [(cmpd.get(str(i), str(i)) if v == "metabolome" else str(i).split("|")[0].replace("g__", ""))
           for i in t.index]
    det[f"{v}|{fac}"] = [{"f": l, "w": round(float(x), 2)} for l, x in zip(lab, t.values)]
res["loadings"] = det
for k, vv in det.items():
    print(k, "->", [(d["f"], d["w"]) for d in vv[:6]])

def _d(o):
    if isinstance(o, (np.integer,)): return int(o)
    if isinstance(o, (np.floating,)): return float(o)
    return str(o)
with open(os.path.join(TB, "stats_step5.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, indent=1, ensure_ascii=False, default=_d)
print("\nsaved stats_step5.json")
