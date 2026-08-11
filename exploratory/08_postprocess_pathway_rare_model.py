# -*- coding: utf-8 -*-
"""Post-process the trained MOFA model: factors, variance explained, associations, plots."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import config as cfg
import os, h5py
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DIR = str(cfg.PROJECT_DIR)
RES = os.path.join(DIR, "results_rareburden")
os.makedirs(RES, exist_ok=True)

f = h5py.File(os.path.join(DIR, "mofa_model_rareburden.hdf5"), "r")
samples = [s.decode() for s in f["samples/all"][:]]
view_names = [v.decode() for v in f["views/views"][:]]
Z = f["expectations/Z/all"][:]                      # factors x samples
factors = [f"Factor{i+1}" for i in range(Z.shape[0])]
Zdf = pd.DataFrame(Z.T, index=samples, columns=factors)

W = {}
for v in view_names:
    feats = [x.decode() for x in f[f"features/{v}"][:]]
    W[v] = pd.DataFrame(f[f"expectations/W/{v}"][:].T, index=feats, columns=factors)

r2 = f["variance_explained/r2_per_factor/all"][:]   # views x factors (fraction or %)
r2df = pd.DataFrame(r2, index=view_names, columns=factors)
if r2df.values.max() <= 1.0:
    r2df *= 100

meta = pd.read_csv(os.path.join(DIR, "sample_metadata.tsv"), sep="\t").set_index("patient").loc[samples]

# ---------- associations ----------
rows = []
for fac in factors:
    z = Zdf[fac]
    uc = z[meta["diagnosis"] == "UC"].dropna()
    cd = z[meta["diagnosis"] == "CD"].dropna()
    p_dx = stats.mannwhitneyu(uc, cd).pvalue if len(uc) > 2 and len(cd) > 2 else np.nan
    r_mb = np.nan
    mb = meta["microbiome_batch"].astype(str)
    if mb.nunique() > 1:
        groups = [z[mb == g].dropna() for g in mb.unique() if (mb == g).sum() > 2]
        if len(groups) > 1:
            r_mb = stats.kruskal(*groups).pvalue
    metb = meta["metabolome_batch"].astype(str)
    r_mtb = np.nan
    groups = [z[metb == g].dropna() for g in metb.unique() if (metb == g).sum() > 2]
    if len(groups) > 1:
        r_mtb = stats.kruskal(*groups).pvalue
    rows.append({"factor": fac, "p_UCvsCD": p_dx, "p_microbiome_batch": r_mb,
                 "p_metabolome_batch": r_mtb})
assoc = pd.DataFrame(rows)

# ---------- write tables ----------
Zdf.join(meta[["diagnosis"]]).to_csv(os.path.join(RES, "factors.tsv"), sep="\t")
r2df.round(2).to_csv(os.path.join(RES, "variance_explained.tsv"), sep="\t")
assoc.to_csv(os.path.join(RES, "factor_associations.tsv"), sep="\t", index=False)
cmpd = pd.read_csv(os.path.join(DIR, "metabolome_compound_names.tsv"), sep="\t",
                   index_col=0)["compound_name"].astype(str).to_dict()
with pd.ExcelWriter(os.path.join(RES, "top_weights.xlsx")) as xw:
    for v in view_names:
        top = []
        for fac in factors:
            w = W[v][fac]
            t = w.reindex(w.abs().sort_values(ascending=False).head(15).index)
            for feat, val in t.items():
                name = cmpd.get(str(feat), "") if v == "metabolome" else ""
                top.append({"factor": fac, "feature": feat, "name": name, "weight": val})
        pd.DataFrame(top).to_excel(xw, sheet_name=v, index=False)

print("=== variance explained (%) ===")
print(r2df.round(1).to_string())
print("\n=== factor associations (p-values) ===")
print(assoc.round(4).to_string(index=False))

# ---------- plots ----------
plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#94A3B8",
                     "axes.labelcolor": "#334155", "text.color": "#334155",
                     "xtick.color": "#64748B", "ytick.color": "#64748B"})
C_UC, C_CD = "#2563EB", "#F59E0B"   # blue / amber (CVD-safe pair)

# (a) variance explained heatmap (sequential blues)
fig, ax = plt.subplots(figsize=(7.5, 2.6), dpi=150)
im = ax.imshow(r2df.values, cmap="Blues", aspect="auto",
               vmin=0, vmax=max(1e-9, r2df.values.max()))
ax.set_xticks(range(len(factors)), factors)
ax.set_yticks(range(len(view_names)), view_names)
for i in range(r2df.shape[0]):
    for j in range(r2df.shape[1]):
        val = r2df.values[i, j]
        ax.text(j, i, f"{val:.1f}", ha="center", va="center",
                color="white" if val > r2df.values.max() * 0.6 else "#334155", fontsize=8)
ax.set_title("Variance explained per factor (%)", loc="left")
fig.colorbar(im, ax=ax, shrink=0.8, label="%")
fig.tight_layout()
fig.savefig(os.path.join(RES, "variance_explained.png"), bbox_inches="tight")
plt.close(fig)

# (b) diagnosis-associated factor: strip + box
best = assoc.sort_values("p_UCvsCD").iloc[0]
fac = best["factor"]
fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), dpi=150)
ax = axes[0]
for k, (g, c) in enumerate([("UC", C_UC), ("CD", C_CD)]):
    v = Zdf.loc[meta["diagnosis"] == g, fac].dropna()
    ax.boxplot([v], positions=[k], widths=0.5, showfliers=False,
               medianprops=dict(color="#334155"), boxprops=dict(color="#94A3B8"),
               whiskerprops=dict(color="#94A3B8"), capprops=dict(color="#94A3B8"))
    ax.scatter(np.random.default_rng(1).normal(k, 0.06, len(v)), v, s=28, color=c,
               edgecolor="white", linewidth=0.8, zorder=3, label=g)
ax.set_xticks([0, 1], ["UC", "CD"])
ax.set_ylabel(f"{fac} score")
ax.set_title(f"{fac} vs diagnosis (Mann-Whitney p={best['p_UCvsCD']:.3g})", loc="left")

ax = axes[1]
f1, f2 = factors[0], factors[1]
for g, c in [("UC", C_UC), ("CD", C_CD)]:
    sub = Zdf[meta["diagnosis"] == g]
    ax.scatter(sub[f1], sub[f2], s=30, color=c, edgecolor="white", linewidth=0.8, label=g)
ax.set_xlabel(f1); ax.set_ylabel(f2)
ax.set_title("Factor 1 vs Factor 2", loc="left")
ax.legend(frameon=False)
for a in axes:
    a.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(RES, "factors_by_diagnosis.png"), bbox_inches="tight")
plt.close(fig)

# (c) top weights of the diagnosis factor, per view (diverging by sign)
fig, axes = plt.subplots(1, len(view_names), figsize=(4.2 * len(view_names), 4.4), dpi=150)
for ax, v in zip(axes, view_names):
    w = W[v][fac]
    t = w.reindex(w.abs().sort_values(ascending=False).head(10).index)[::-1]
    labels = []
    for feat in t.index:
        lab = cmpd.get(str(feat), str(feat)) if v == "metabolome" else str(feat)
        lab = lab.split("|")[0]
        lab = lab.replace("g__", "") if v == "microbiome" else lab
        labels.append(lab[:28])
    colors = ["#B91C1C" if x > 0 else "#1D4ED8" for x in t.values]
    ax.barh(range(len(t)), t.values, color=colors, height=0.62)
    ax.set_yticks(range(len(t)), labels, fontsize=7.5)
    ax.axvline(0, color="#94A3B8", lw=0.8)
    ax.set_title(v, loc="left", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
fig.suptitle(f"Top weights on {fac} (red = higher in {fac}+, blue = lower)", x=0.01, ha="left", fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(os.path.join(RES, "top_weights_diagnosis_factor.png"), bbox_inches="tight")
plt.close(fig)

print("\nplots & tables saved to", RES)
