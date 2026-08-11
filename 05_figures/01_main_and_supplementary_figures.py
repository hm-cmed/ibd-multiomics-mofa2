# -*- coding: utf-8 -*-
"""Manuscript step 4: generate all main and supplementary figures."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import config as cfg
import os, json
import h5py
import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, leaves_list
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DIR = str(cfg.PROJECT_DIR)
MS = os.path.join(DIR, "manuscript")
TB = os.path.join(MS, "tables")
FG = os.path.join(MS, "figures")

plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": "#94A3B8", "axes.linewidth": 0.8,
                     "axes.labelcolor": "#334155", "text.color": "#334155",
                     "xtick.color": "#64748B", "ytick.color": "#64748B", "figure.dpi": 300})
C_UC, C_CD = "#2563EB", "#F59E0B"
C_POS, C_NEG = "#B91C1C", "#1D4ED8"
GRAY = "#94A3B8"

# ---- load model ----
with h5py.File(os.path.join(MS, "mofa_ms_K8.hdf5"), "r") as f:
    samples = [s.decode() for s in f["samples/all"][:]]
    views = [v.decode() for v in f["views/views"][:]]
    Z = pd.DataFrame(f["expectations/Z/all"][:].T, index=samples,
                     columns=[f"Factor{i+1}" for i in range(f["expectations/Z/all"].shape[0])])
    W = {v: pd.DataFrame(f[f"expectations/W/{v}"][:].T,
                         index=[x.decode() for x in f[f"features/{v}"][:]], columns=Z.columns) for v in views}
    r2 = pd.DataFrame(f["variance_explained/r2_per_factor/all"][:], index=views, columns=Z.columns)
if r2.values.max() <= 1.0:
    r2 *= 100

meta = pd.read_csv(os.path.join(MS, "meta_with_clinical.tsv"), sep="\t", index_col=0).loc[samples]
cmpd = pd.read_csv(os.path.join(DIR, "metabolome_compound_names.tsv"), sep="\t", index_col=0)["compound_name"].astype(str).to_dict()
VIEW_LABEL = {"mutation": "Mutation", "microbiome": "Microbiome", "metabolome": "Metabolome",
              "transcriptome": "Transcriptome"}

def label_feat(v, feat):
    if v == "metabolome":
        return cmpd.get(str(feat), str(feat))[:26]
    if v == "microbiome":
        return str(feat).replace("g__", "")[:26]
    return str(feat).split("|")[0][:26]

# ================= Figure 1: cohort & availability =================
avail = pd.DataFrame({
    "Mutation": meta.has_mutation, "Microbiome": meta.has_microbiome,
    "Metabolome": meta.has_metabolome, "Transcriptome": meta.has_transcriptome}).astype(int)
order = avail.sum(axis=1).sort_values(ascending=False).index
fig, axes = plt.subplots(1, 2, figsize=(9, 3.1), gridspec_kw={"width_ratios": [1, 2.2]})
ax = axes[0]
ns = avail.sum()
D = {"Mutation": 49, "Microbiome": 69, "Metabolome": 62, "Transcriptome": 2000}
bars = ax.bar(range(4), ns.values, color=["#7C3AED", "#059669", "#D97706", "#DC2626"], width=0.6)
ax.axhline(40, ls="--", color=GRAY, lw=0.9)
ax.text(3.45, 40.6, "N = 40", fontsize=7.5, color="#64748B", ha="right")
for i, (v, n_) in enumerate(ns.items()):
    ax.text(i, n_ + 0.6, f"{n_}", ha="center", fontsize=8)
    ax.text(i, 2.0, f"D={D[v]}", ha="center", fontsize=6.8, color="white", rotation=90)
ax.set_xticks(range(4), ns.index, rotation=20, ha="right")
ax.set_ylabel("Patients with data")
ax.set_ylim(0, 45)
ax.set_title("a", loc="left", fontsize=12, fontweight="bold")
ax.spines[["top", "right"]].set_visible(False)

ax = axes[1]
M = avail.loc[order].T.values
im = ax.imshow(M, cmap=matplotlib.colors.ListedColormap(["#E2E8F0", "#2563EB"]), aspect="auto")
ax.set_yticks(range(4), avail.columns)
ax.set_xticks([])
for i, p in enumerate(order):
    c = C_UC if meta.loc[p, "diagnosis"] == "UC" else C_CD
    ax.add_patch(plt.Rectangle((i - 0.5, -0.95), 1, 0.35, color=c, clip_on=False))
ax.text(-0.6, -0.78, "Dx", ha="right", fontsize=7)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=C_UC, label="UC"), Patch(color=C_CD, label="CD"),
                   Patch(color="#2563EB", label="available"), Patch(color="#E2E8F0", label="missing")],
          loc="lower right", frameon=False, fontsize=6.8, ncols=4, bbox_to_anchor=(1, -0.28))
ax.set_title("b", loc="left", fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(FG, "Fig1_cohort_availability.png"), bbox_inches="tight")
plt.close(fig)

# ================= Figure 2: MOFA overview =================
fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.0), gridspec_kw={"width_ratios": [1.35, 1]})
ax = axes[0]
im = ax.imshow(r2.values, cmap="Blues", aspect="auto", vmin=0, vmax=r2.values.max())
ax.set_xticks(range(r2.shape[1]), r2.columns, rotation=30, ha="right")
ax.set_yticks(range(4), [VIEW_LABEL[v] for v in r2.index])
for i in range(r2.shape[0]):
    for j in range(r2.shape[1]):
        v_ = r2.values[i, j]
        if v_ >= 0.5:
            ax.text(j, i, f"{v_:.1f}", ha="center", va="center", fontsize=7,
                    color="white" if v_ > r2.values.max() * 0.55 else "#334155")
fig.colorbar(im, ax=ax, shrink=0.85, label="Variance explained (%)")
ax.set_title("a", loc="left", fontsize=12, fontweight="bold")

ax = axes[1]
for g, c in [("UC", C_UC), ("CD", C_CD)]:
    sub = Z[meta.diagnosis == g]
    ax.scatter(sub["Factor1"], sub["Factor2"], s=34, color=c, edgecolor="white", lw=0.7, label=g)
ax.set_xlabel("Factor 1 (transcriptome-dominant)")
ax.set_ylabel("Factor 2 (metabolome / CRP-associated)")
ax.legend(frameon=False)
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("b", loc="left", fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(FG, "Fig2_mofa_overview.png"), bbox_inches="tight")
plt.close(fig)

# ================= Figure 3: PERMANOVA / PCoA / Spearman =================
perma = pd.read_csv(os.path.join(TB, "permanova_gene_microbiome.csv"))
perma_g = perma[~perma.gene.str.startswith("[")].sort_values("pseudo_F", ascending=False).head(20)
pcoa = pd.read_csv(os.path.join(MS, "pcoa_mefv.tsv"), sep="\t", index_col=0)
sp = pd.read_csv(os.path.join(TB, "spearman_gene_genus.csv"))

fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.4), gridspec_kw={"width_ratios": [1.15, 0.9, 1.25]})
ax = axes[0]
hl = {"CASP7", "RIGI", "TICAM1", "MEFV", "CASP9"}
colors = ["#D97706" if g in hl else GRAY for g in perma_g.gene]
ax.bar(range(len(perma_g)), perma_g.pseudo_F, color=colors, width=0.65)
for i, (_, row) in enumerate(perma_g.iterrows()):
    if row.p < 0.06:
        ax.text(i, row.pseudo_F + 0.04, f"p={row.p:.3f}", rotation=90, fontsize=6, ha="center")
ax.set_xticks(range(len(perma_g)), perma_g.gene, rotation=75, fontsize=6.5)
ax.set_ylabel("PERMANOVA pseudo-F")
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("a", loc="left", fontsize=12, fontweight="bold")

ax = axes[1]
for g, c, mk in [("carrier", "#D97706", "o"), ("non-carrier", "#64748B", "s")]:
    sub = pcoa[pcoa.MEFV_E148Q == g]
    ax.scatter(sub.PCo1, sub.PCo2, s=30, color=c, marker=mk, edgecolor="white", lw=0.6, label=f"E148Q {g}")
ax.set_xlabel(f"PCo1 ({pcoa['var1'].iloc[0]*100:.0f}%)")
ax.set_ylabel(f"PCo2 ({pcoa['var2'].iloc[0]*100:.0f}%)")
ax.legend(frameon=False, fontsize=7)
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("b", loc="left", fontsize=12, fontweight="bold")

ax = axes[2]
top_genes = sp.groupby("gene")["p"].min().sort_values().head(10).index
top_gen = sp[sp.gene.isin(top_genes)].groupby("genus")["p"].min().sort_values().head(12).index
mat = sp[sp.gene.isin(top_genes) & sp.genus.isin(top_gen)].pivot(index="genus", columns="gene", values="rho")
im = ax.imshow(mat.values, cmap="RdBu_r", vmin=-0.7, vmax=0.7, aspect="auto")
ax.set_xticks(range(mat.shape[1]), mat.columns, rotation=75, fontsize=6.5)
ax.set_yticks(range(mat.shape[0]), [g.replace("g__", "") for g in mat.index], fontsize=6.5)
fig.colorbar(im, ax=ax, shrink=0.8, label="Spearman \u03c1")
ax.set_title("c", loc="left", fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(FG, "Fig3_gene_microbiome.png"), bbox_inches="tight")
plt.close(fig)

# ================= Figure 4: genotype univariate boxplots =================
clr = pd.read_csv(os.path.join(DIR, "view_microbiome_clr.tsv"), sep="\t", index_col=0)
metab = pd.read_csv(os.path.join(DIR, "view_metabolome_log2.tsv"), sep="\t", index_col=0)
metab.index = metab.index.astype(str)
carrier_masks = {"MEFV E148Q": meta["MEFV_E148Q_dosage"] > 0}
# CASP9 carrier from univariate tables' groups — rebuild via burden file
c9tab = pd.read_csv(os.path.join(TB, "univariate_CASP9_microbiome.csv"))

def panel_box(ax, feat_df, mask, feats, labels, ylab):
    data, xt = [], []
    for ft in feats:
        x = feat_df.loc[ft]
        data.append((x[mask.reindex(x.index).fillna(False)].dropna(),
                     x[~mask.reindex(x.index).fillna(True)].dropna()))
    pos = 0
    rng_ = np.random.default_rng(0)
    for (a, b), lb in zip(data, labels):
        for k, (vals, c) in enumerate([(a, "#D97706"), (b, "#64748B")]):
            ax.boxplot([vals], positions=[pos + k * 0.7], widths=0.5, showfliers=False,
                       medianprops=dict(color="#334155"), boxprops=dict(color=GRAY),
                       whiskerprops=dict(color=GRAY), capprops=dict(color=GRAY))
            ax.scatter(rng_.normal(pos + k * 0.7, 0.07, len(vals)), vals, s=12, color=c,
                       edgecolor="white", lw=0.4, zorder=3)
        xt.append((pos + 0.35, lb))
        pos += 2.2
    ax.set_xticks([p for p, _ in xt], [l for _, l in xt], rotation=25, ha="right", fontsize=7)
    ax.set_ylabel(ylab)
    ax.spines[["top", "right"]].set_visible(False)

uni = {}
for gname in ["MEFV_E148Q", "CASP9"]:
    uni[gname] = {
        "mic": pd.read_csv(os.path.join(TB, f"univariate_{gname}_microbiome.csv")),
        "met": pd.read_csv(os.path.join(TB, f"univariate_{gname}_metabolome.csv"))}

with open(os.path.join(TB, "stats_step3.json"), encoding="utf-8") as fh:
    s3 = json.load(fh)
casp9_carrier = pd.Series({p: False for p in meta.index})
# reconstruct CASP9 carrier: use file from step3? easiest: recompute via mutation gene rare? CASP9 carrier = any missense
# stored implicitly; rebuild from all_genes json quickly
import re as _re
from collections import defaultdict as _dd
d_ = json.load(open(os.path.join(str(cfg.MUTATION_DIR), "Client_aa_mut", "all_genes_all_patients.json"), encoding="utf-8"))
canon_re = _re.compile(r"(?:[gn]om[-_ ]?)?(sap|saga|asic|okh|okn)[-_ ]?0*(\d+)", _re.I)
def canon(s):
    m = canon_re.search(str(s).lower())
    fam, num = m.group(1), int(m.group(2))
    fam = "okn" if fam == "okh" else fam
    return {"sap": f"NOMSap_{num:03d}", "saga": f"Saga_{num:03d}", "asic": f"ASIC{num:03d}", "okn": f"OKN{num:03d}"}[fam]
c9 = d_["CASP9-" + sorted([t.split("-", 1)[1] for t in d_ if t.startswith("CASP9-") and t.split("-", 1)[1].startswith("NP_")])[0]]
for pat_raw, rec in c9.items():
    casp9_carrier[canon(pat_raw)] = any(m["ref_aa"] != m["mut_aa"] for m in rec.get("mutations", []))
carrier_masks["CASP9"] = casp9_carrier.astype(bool)

fig, axes = plt.subplots(2, 2, figsize=(10, 6.2))
for row, gname, key in [(0, "MEFV E148Q", "MEFV_E148Q"), (1, "CASP9", "CASP9")]:
    tm = uni[key]["mic"].head(4)
    panel_box(axes[row][0], clr, carrier_masks[gname], tm.feature.tolist(),
              [f"{f.replace('g__','')}\np={p:.3f}" for f, p in zip(tm.feature, tm.p)], "CLR abundance")
    axes[row][0].set_title(f"{'ac'[row]}", loc="left", fontsize=12, fontweight="bold")
    tt = uni[key]["met"].head(4)
    panel_box(axes[row][1], metab, carrier_masks[gname], [str(f) for f in tt.feature],
              [f"{l}\np={p:.3f}" for l, p in zip(tt.label, tt.p)], "log2 intensity")
    axes[row][1].set_title(f"{'bd'[row]}", loc="left", fontsize=12, fontweight="bold")
from matplotlib.patches import Patch
fig.legend(handles=[Patch(color="#D97706", label="carrier"), Patch(color="#64748B", label="non-carrier")],
           loc="upper right", frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FG, "Fig4_genotype_univariate.png"), bbox_inches="tight")
plt.close(fig)

# ================= Figure 5: loadings =================
dom = [("mutation", "Factor7"), ("microbiome", "Factor3"), ("metabolome", "Factor2"), ("transcriptome", "Factor1")]
fig, axes = plt.subplots(1, 4, figsize=(13.5, 4.0))
for ax, (v, fac) in zip(axes, dom):
    w = W[v][fac]
    t = w.reindex(w.abs().sort_values(ascending=False).head(12).index)[::-1]
    ax.barh(range(len(t)), t.values, color=[C_POS if x > 0 else C_NEG for x in t.values], height=0.6)
    ax.set_yticks(range(len(t)), [label_feat(v, i) for i in t.index], fontsize=6.6)
    ax.axvline(0, color=GRAY, lw=0.7)
    ax.set_title(f"{'abcd'[dom.index((v, fac))]}", loc="left", fontsize=12, fontweight="bold")
    ax.text(0.5, 1.02, f"{VIEW_LABEL[v]} \u2013 {fac}", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=8.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlabel("MOFA2 weight")
fig.tight_layout()
fig.savefig(os.path.join(FG, "Fig5_loadings.png"), bbox_inches="tight")
plt.close(fig)

# ================= Figure 6: subtype & CRP =================
cyt = metab.loc["6175"].dropna()
crp = meta["CRP"]
fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.2))
ax = axes[0]
counts = np.array([[21, 27 - 21], [6, 13 - 6]], dtype=float)
for i, (g, c) in enumerate([("UC", C_UC), ("CD", C_CD)]):
    frac = counts[i, 0] / counts[i].sum() * 100
    ax.bar(i, frac, color=c, width=0.55)
    ax.text(i, frac + 2, f"{int(counts[i,0])}/{int(counts[i].sum())}", ha="center", fontsize=8)
ax.set_xticks([0, 1], ["UC", "CD"])
ax.set_ylabel("CASP7 R6P carriers (%)")
ax.set_ylim(0, 100)
ax.text(0.5, 92, "Fisher p=0.072", ha="center", fontsize=7.5)
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("a", loc="left", fontsize=12, fontweight="bold")

ax = axes[1]
rng_ = np.random.default_rng(1)
for i, (g, c) in enumerate([("UC", C_UC), ("CD", C_CD)]):
    vals = cyt[meta.loc[cyt.index, "diagnosis"] == g]
    ax.boxplot([vals], positions=[i], widths=0.5, showfliers=False,
               medianprops=dict(color="#334155"), boxprops=dict(color=GRAY),
               whiskerprops=dict(color=GRAY), capprops=dict(color=GRAY))
    ax.scatter(rng_.normal(i, 0.06, len(vals)), vals, s=22, color=c, edgecolor="white", lw=0.6, zorder=3)
ax.set_xticks([0, 1], [f"UC (n={(meta.loc[cyt.index,'diagnosis']=='UC').sum()})",
                       f"CD (n={(meta.loc[cyt.index,'diagnosis']=='CD').sum()})"])
ax.set_ylabel("Cytidine, log2 intensity")
ax.set_title("b", loc="left", fontsize=12, fontweight="bold")
ax.spines[["top", "right"]].set_visible(False)

ax = axes[2]
ok = cyt.index.intersection(crp.dropna().index)
for g, c in [("UC", C_UC), ("CD", C_CD)]:
    idx = [p for p in ok if meta.loc[p, "diagnosis"] == g]
    ax.scatter(cyt[idx], crp[idx], s=26, color=c, edgecolor="white", lw=0.6, label=g)
ax.set_xlabel("Cytidine, log2 intensity")
ax.set_ylabel("Serum CRP (mg/dL)")
ax.legend(frameon=False, fontsize=7.5)
ax.text(0.97, 0.92, "\u03c1=\u22120.56, p=0.007", transform=ax.transAxes, ha="right", fontsize=7.5)
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("c", loc="left", fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(FG, "Fig6_subtype_crp.png"), bbox_inches="tight")
plt.close(fig)

# ================= Supp figures =================
# S1: grid K
grid = {}
for K in [5, 6, 8, 10]:
    p = os.path.join(MS, "tables", f"variance_explained_K{K}.csv")
    if os.path.exists(p):
        t = pd.read_csv(p, index_col=0)
        grid[K] = t.sum(axis=1)
gd = pd.DataFrame(grid)
fig, ax = plt.subplots(figsize=(5.5, 3), dpi=300)
x = np.arange(len(gd.columns))
w_ = 0.19
for i, v in enumerate(gd.index):
    ax.bar(x + i * w_, gd.loc[v], width=w_, label=VIEW_LABEL[v],
           color=["#7C3AED", "#059669", "#D97706", "#DC2626"][i])
ax.set_xticks(x + 1.5 * w_, [f"K={k}" for k in gd.columns])
ax.set_ylabel("Total variance explained (%)")
ax.legend(frameon=False, fontsize=7)
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("", loc="left", fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(FG, "FigS1_grid_search.png"), bbox_inches="tight")
plt.close(fig)

# S2: factor correlations
C = Z.corr()
fig, ax = plt.subplots(figsize=(4.6, 3.8), dpi=300)
im = ax.imshow(C.values, cmap="RdBu_r", vmin=-1, vmax=1)
ax.set_xticks(range(len(C)), C.columns, rotation=45, ha="right", fontsize=7)
ax.set_yticks(range(len(C)), C.columns, fontsize=7)
for i in range(len(C)):
    for j in range(len(C)):
        ax.text(j, i, f"{C.values[i,j]:.2f}", ha="center", va="center", fontsize=6,
                color="white" if abs(C.values[i, j]) > 0.6 else "#334155")
fig.colorbar(im, ax=ax, shrink=0.8)
ax.set_title("", loc="left", fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(FG, "FigS2_factor_correlations.png"), bbox_inches="tight")
plt.close(fig)

# S3: patient x factor heatmap
Zz = (Z - Z.mean()) / Z.std()
lk = linkage(Zz.values, method="ward")
order3 = leaves_list(lk)
fig, ax = plt.subplots(figsize=(9, 3.2), dpi=300)
im = ax.imshow(Zz.iloc[order3].T.values, cmap="RdBu_r", vmin=-2.5, vmax=2.5, aspect="auto")
ax.set_yticks(range(Z.shape[1]), Z.columns, fontsize=7)
ax.set_xticks([])
for i, p in enumerate(Zz.index[order3]):
    c = C_UC if meta.loc[p, "diagnosis"] == "UC" else C_CD
    ax.add_patch(plt.Rectangle((i - 0.5, -1.15), 1, 0.45, color=c, clip_on=False))
fig.colorbar(im, ax=ax, shrink=0.8, label="z-scored factor score")
ax.set_title("", loc="left", fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(FG, "FigS3_patient_heatmap.png"), bbox_inches="tight")
plt.close(fig)

# TableS6 from grid csvs
gd.T.round(2).assign(total=gd.sum().round(2)).to_csv(os.path.join(TB, "TableS6_grid_search.csv"))
print("figures done:", os.listdir(FG))
