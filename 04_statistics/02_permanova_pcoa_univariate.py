# -*- coding: utf-8 -*-
"""Manuscript step 3: PERMANOVA / PCoA / univariate genotype associations / subtype tests."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import config as cfg
import os, json, re
import numpy as np
import pandas as pd
from scipy import stats
from collections import defaultdict

BASE = str(cfg.DATA_ROOT)
DIR = str(cfg.PROJECT_DIR)
MS = os.path.join(DIR, "manuscript")
TB = os.path.join(MS, "tables")
rng = np.random.default_rng(42)

canon_re = re.compile(r"(?:[gn]om[-_ ]?)?(sap|saga|asic|okh|okn)[-_ ]?0*(\d+)", re.I)
def canon(s):
    m = canon_re.search(str(s).lower())
    fam, num = m.group(1), int(m.group(2))
    if fam == "okh":
        fam = "okn"
    return {"sap": f"NOMSap_{num:03d}", "saga": f"Saga_{num:03d}",
            "asic": f"ASIC{num:03d}", "okn": f"OKN{num:03d}"}[fam]

meta = pd.read_csv(os.path.join(MS, "meta_with_clinical.tsv"), sep="\t", index_col=0)

# ---------- genus relative abundance (35 patients) ----------
import json as _json
MUT = cfg.MUTATION_DIR
m1 = _json.load(open(os.path.join(MUT, "microbiome1.json"), encoding="utf-8"))
m2 = _json.load(open(os.path.join(MUT, "microbiome2.json"), encoding="utf-8"))
def to_genus(mjson):
    rows = {}
    for taxon, vals in mjson.items():
        g = re.search(r"(g__[^;]+)", taxon)
        if not g:
            continue
        for s, v in vals.items():
            if v is None:
                continue
            rows.setdefault(g.group(1), {}).setdefault(s, 0.0)
            rows[g.group(1)][s] += float(v)
    return pd.DataFrame(rows).T
g1, g2 = to_genus(m1), to_genus(m2)
common = sorted(set(g1.index) & set(g2.index))
g1, g2 = g1.loc[common].fillna(0), g2.loc[common].fillna(0)
sample_of = {}
for s in g1.columns:
    sample_of.setdefault(canon(s), (s, 1))
for s in g2.columns:
    sample_of.setdefault(canon(s), (s, 2))
patients_m = [p for p in meta.index if p in sample_of]
counts = pd.DataFrame({p: (g1 if sample_of[p][1] == 1 else g2)[sample_of[p][0]] for p in patients_m})
rel = counts / counts.sum(axis=0)
rel.to_csv(os.path.join(MS, "genus_relabundance_35pts.tsv"), sep="\t")

# Bray-Curtis
X = rel.T.to_numpy()
n = X.shape[0]
BC = np.zeros((n, n))
for i in range(n):
    for j in range(i + 1, n):
        BC[i, j] = BC[j, i] = np.abs(X[i] - X[j]).sum() / (X[i] + X[j]).sum()
BCdf = pd.DataFrame(BC, index=patients_m, columns=patients_m)

def permanova(D, labels, n_perm=999):
    labels = np.asarray(labels)
    uniq = np.unique(labels)
    N = len(labels)
    def pseudo_f(lab):
        SST = (D ** 2).sum() / (2 * N)
        SSW = 0.0
        for u in uniq:
            idx = np.where(lab == u)[0]
            if len(idx) > 1:
                sub = D[np.ix_(idx, idx)]
                SSW += (sub ** 2).sum() / (2 * len(idx))
        SSA = SST - SSW
        a = len(uniq)
        return (SSA / (a - 1)) / (SSW / (N - a))
    F0 = pseudo_f(labels)
    cnt = 0
    for _ in range(n_perm):
        if pseudo_f(rng.permutation(labels)) >= F0:
            cnt += 1
    return F0, (cnt + 1) / (n_perm + 1)

# ---------- gene carrier matrix (any missense, canonical tx) ----------
d = _json.load(open(os.path.join(MUT, r"Client_aa_mut\all_genes_all_patients.json"), encoding="utf-8"))
tx_by_gene = defaultdict(set)
for gt in d:
    g, tx = gt.split("-", 1)
    tx_by_gene[g].add(tx)
canon_tx = {g: (sorted([t for t in txs if t.startswith("NP_")]) or sorted(txs))[0] for g, txs in tx_by_gene.items()}
carrier = pd.DataFrame(False, index=sorted(tx_by_gene), columns=meta.index)
for gt, patmap in d.items():
    g, tx = gt.split("-", 1)
    if tx != canon_tx[g]:
        continue
    for pat_raw, rec in patmap.items():
        pat = canon(pat_raw)
        if any(m["ref_aa"] != m["mut_aa"] for m in rec.get("mutations", [])):
            carrier.loc[g, pat] = True

# PERMANOVA per gene (>=5 in each group among the 35)
perma_rows = []
sub_car = carrier[patients_m]
for g in carrier.index:
    lab = sub_car.loc[g].to_numpy()
    if lab.sum() < 5 or (~lab).sum() < 5:
        continue
    F0, p = permanova(BC, lab.astype(int))
    perma_rows.append({"gene": g, "n_carrier": int(lab.sum()), "n_noncarrier": int((~lab).sum()),
                       "pseudo_F": round(float(F0), 3), "p": p})
# covariates
for name, lab in [("diagnosis_UCvsCD", (meta.loc[patients_m, "diagnosis"] == "CD").to_numpy().astype(int)),
                  ("microbiome_batch", meta.loc[patients_m, "microbiome_batch"].to_numpy())]:
    F0, p = permanova(BC, lab)
    perma_rows.append({"gene": f"[{name}]", "n_carrier": "", "n_noncarrier": "",
                       "pseudo_F": round(float(F0), 3), "p": p})
perma = pd.DataFrame(perma_rows).sort_values("p")
perma.to_csv(os.path.join(TB, "permanova_gene_microbiome.csv"), index=False)
print("=== PERMANOVA (top 12) ===")
print(perma.head(12).to_string(index=False))

# MEFV E148Q carrier PERMANOVA + PCoA
e148q = (pd.read_csv(os.path.join(DIR, "view_mutation_rareburden_mefv.tsv"), sep="\t", index_col=0)
         .loc["MEFV_E148Q"] > meta["MEFV_E148Q_dosage"].min())
e148q = meta["MEFV_E148Q_dosage"] > 0
lab = e148q.loc[patients_m].to_numpy().astype(int)
F0, p = permanova(BC, lab)
print(f"MEFV E148Q PERMANOVA: F={F0:.3f} p={p:.3f} (carriers {lab.sum()}/{len(lab)})")
mefv_perma = {"F": round(float(F0), 3), "p": p, "n_carrier": int(lab.sum()), "n": len(lab)}

# PCoA (classical MDS)
J = np.eye(n) - np.ones((n, n)) / n
B = -0.5 * J @ (BC ** 2) @ J
evals, evecs = np.linalg.eigh(B)
order = np.argsort(evals)[::-1]
coords = evecs[:, order[:2]] * np.sqrt(np.abs(evals[order[:2]]))
pcoa = pd.DataFrame(coords, index=patients_m, columns=["PCo1", "PCo2"])
pcoa["MEFV_E148Q"] = ["carrier" if x else "non-carrier" for x in lab]
pcoa["var1"] = evals[order[0]] / evals[evals > 0].sum()
pcoa["var2"] = evals[order[1]] / evals[evals > 0].sum()
pcoa.to_csv(os.path.join(MS, "pcoa_mefv.tsv"), sep="\t")

# ---------- Spearman gene rare-burden x genus CLR ----------
gb = pd.read_csv(os.path.join(DIR, "view_mutation_gene_rare.tsv"), sep="\t", index_col=0)
clr = pd.read_csv(os.path.join(DIR, "view_microbiome_clr.tsv"), sep="\t", index_col=0)[patients_m]
res = []
for g in gb.index:
    x = gb.loc[g, patients_m]
    if x.std() == 0:
        continue
    for tax in clr.index:
        r, pv = stats.spearmanr(x, clr.loc[tax])
        res.append((g, tax, r, pv))
sp = pd.DataFrame(res, columns=["gene", "genus", "rho", "p"])
sp["q"] = sp["p"] * len(sp) / sp["p"].rank()
sp["q"] = sp["q"].clip(upper=1)
sp.sort_values("p").head(500).to_csv(os.path.join(TB, "spearman_gene_genus.csv"), index=False)
print(f"gene x genus pairs: {len(sp)}, min q = {sp['q'].min():.3f}, nominal p<0.05: {(sp.p < 0.05).sum()}")

# Mantel test (euclid on gene burden vs BC), 999 perms
gbm = gb[patients_m].T.to_numpy()
Dg = np.sqrt(((gbm[:, None, :] - gbm[None, :, :]) ** 2).sum(-1))
iu = np.triu_indices(n, 1)
r0 = stats.spearmanr(Dg[iu], BC[iu]).statistic
cnt = 0
for _ in range(999):
    perm = rng.permutation(n)
    if abs(stats.spearmanr(Dg[np.ix_(perm, perm)][iu], BC[iu]).statistic) >= abs(r0):
        cnt += 1
mantel = {"r": round(float(r0), 3), "p": (cnt + 1) / 1000}
print("Mantel:", mantel)

# ---------- univariate genotype associations (MEFV E148Q, CASP9) ----------
metab = pd.read_csv(os.path.join(DIR, "view_metabolome_log2.tsv"), sep="\t", index_col=0)
cmpd = pd.read_csv(os.path.join(DIR, "metabolome_compound_names.tsv"), sep="\t", index_col=0)["compound_name"].astype(str).to_dict()

def mw_table(feat_df, group_mask, label_map=None):
    out = []
    cols = feat_df.columns
    gm = group_mask.loc[cols]
    for feat in feat_df.index:
        x = feat_df.loc[feat]
        a, b = x[gm], x[~gm]
        if a.notna().sum() < 3 or b.notna().sum() < 3:
            continue
        pv = stats.mannwhitneyu(a.dropna(), b.dropna()).pvalue
        out.append({"feature": feat, "label": (label_map or {}).get(str(feat), str(feat)),
                    "median_carrier": round(float(a.median()), 3),
                    "median_noncarrier": round(float(b.median()), 3), "p": pv})
    t = pd.DataFrame(out).sort_values("p")
    t["q"] = (t["p"] * len(t) / t["p"].rank()).clip(upper=1)
    return t

geno_stats = {}
for gname, mask in [("MEFV_E148Q", meta["MEFV_E148Q_dosage"] > 0),
                    ("CASP9", carrier.loc["CASP9"])]:
    mm = clr.columns
    t_mic = mw_table(clr, mask)
    mcols = metab.columns[metab.notna().any()]
    t_met = mw_table(metab[mcols].dropna(axis=1, how="all"), mask, cmpd)
    t_mic.to_csv(os.path.join(TB, f"univariate_{gname}_microbiome.csv"), index=False)
    t_met.to_csv(os.path.join(TB, f"univariate_{gname}_metabolome.csv"), index=False)
    nmask_mic = int(mask.loc[clr.columns].sum())
    metab_pat = [c for c in metab.columns if metab[c].notna().any()]
    nmask_met = int(mask.loc[metab_pat].sum())
    geno_stats[gname] = {
        "micro": {"n_carrier": nmask_mic, "n_total": len(clr.columns),
                  "top": t_mic.head(5)[["label", "p", "q"]].to_dict("records")},
        "metab": {"n_carrier": nmask_met, "n_total": len(metab_pat),
                  "top": t_met.head(5)[["label", "p", "q"]].to_dict("records")},
        "min_q_micro": float(t_mic["q"].min()), "min_q_metab": float(t_met["q"].min()),
    }
    print(f"--- {gname}: micro carriers {nmask_mic}/{len(clr.columns)}, metab {nmask_met}/{len(metab_pat)}")
    print(t_mic.head(5).to_string(index=False))
    print(t_met.head(5).to_string(index=False))

# ---------- CASP7 R6P (isoform NP_001253986.1) ----------
c7 = d.get("CASP7-NP_001253986.1", {})
r6p = {}
for pat_raw, rec in c7.items():
    pat = canon(pat_raw)
    r6p[pat] = any(m["ref_aa"] == "R" and m["position"] == 6 and m["mut_aa"] == "P" for m in rec.get("mutations", []))
r6p = pd.Series(r6p).reindex(meta.index).fillna(False)
tab = pd.crosstab(r6p, meta["diagnosis"])
odds, p_c7 = stats.fisher_exact(tab.values)
c7_stats = {"carriers_total": int(r6p.sum()),
            "UC": f"{int((r6p & (meta.diagnosis=='UC')).sum())}/{int((meta.diagnosis=='UC').sum())}",
            "CD": f"{int((r6p & (meta.diagnosis=='CD')).sum())}/{int((meta.diagnosis=='CD').sum())}",
            "fisher_p": round(float(p_c7), 4), "odds": round(float(odds), 3)}
print("CASP7 R6P:", c7_stats)

# ---------- cytidine ----------
cyt = metab.loc["6175"] if "6175" in metab.index else metab.loc[6175]
cyt = cyt.dropna()
uc = cyt[meta.loc[cyt.index, "diagnosis"] == "UC"]
cd = cyt[meta.loc[cyt.index, "diagnosis"] == "CD"]
p_cyt = stats.mannwhitneyu(uc, cd).pvalue
crp = meta.loc[cyt.index, "CRP"]
ok = crp.notna()
r_crp, p_crp = stats.spearmanr(cyt[ok], crp[ok])
cyt_stats = {"n": len(cyt), "n_UC": len(uc), "n_CD": len(cd), "p_UCvsCD": round(float(p_cyt), 4),
             "median_UC": round(float(uc.median()), 2), "median_CD": round(float(cd.median()), 2),
             "rho_CRP": round(float(r_crp), 3), "p_CRP": round(float(p_crp), 4), "n_CRP": int(ok.sum())}
print("cytidine:", cyt_stats)

out = {"permanova_top": perma.head(12).to_dict("records"), "mefv_permanova": mefv_perma,
       "mantel": mantel, "spearman_min_q": float(sp["q"].min()),
       "spearman_n_pairs": len(sp), "genotype": geno_stats, "CASP7_R6P": c7_stats, "cytidine": cyt_stats}
def _default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return str(o)
with open(os.path.join(TB, "stats_step3.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=1, ensure_ascii=False, default=_default)
print("saved stats_step3.json")
