# -*- coding: utf-8 -*-
"""Figure 7: paired mucosal DE, Factor 1 as inflammation axis, alpha diversity by MEFV E148Q."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import config as cfg
import os, json
import h5py
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DIR = str(cfg.PROJECT_DIR)
MS = os.path.join(DIR, "manuscript")
TB, FG = os.path.join(MS, "tables"), os.path.join(MS, "figures")

plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": "#94A3B8", "axes.linewidth": 0.8,
                     "axes.labelcolor": "#334155", "text.color": "#334155",
                     "xtick.color": "#64748B", "ytick.color": "#64748B", "figure.dpi": 300})
C_UC, C_CD, GRAY = "#2563EB", "#F59E0B", "#94A3B8"

de = pd.read_csv(os.path.join(TB, "TableS9_paired_DE_inf_vs_nor.csv"))
meta = pd.read_csv(os.path.join(MS, "meta_with_clinical.tsv"), sep="\t", index_col=0)
alpha = pd.read_csv(os.path.join(TB, "TableS10_alpha_diversity.csv"), index_col=0)
with h5py.File(os.path.join(MS, "mofa_ms_K8.hdf5"), "r") as f:
    Z = pd.DataFrame(f["expectations/Z/all"][:].T, index=[s.decode() for s in f["samples/all"][:]],
                     columns=[f"Factor{i+1}" for i in range(f["expectations/Z/all"].shape[0])])

lc = pd.read_csv(os.path.join(DIR, "rna_logcpm_all_samples.tsv"), sep="\t", index_col=0)
lc.index = lc.index.astype(str)
up = de[(de.q < 0.05) & (de.log2FC_inf_vs_nor > 0.5)]
top100 = [str(g) for g in up.head(100).gene]
inf_c = [c for c in lc.columns if c.endswith("__inf")]
inflam = lc.loc[top100, inf_c].mean()
inflam.index = [c[:-5] for c in inflam.index]

fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.4))

# a: volcano
ax = axes[0]
x = de.log2FC_inf_vs_nor.values
y = -np.log10(np.clip(de.p.values, 1e-12, None))
sig_up = (de.q < 0.05) & (de.log2FC_inf_vs_nor > 0.5)
sig_dn = (de.q < 0.05) & (de.log2FC_inf_vs_nor < -0.5)
ax.scatter(x[~(sig_up | sig_dn)], y[~(sig_up | sig_dn)], s=3, color="#CBD5E1", rasterized=True)
ax.scatter(x[sig_up], y[sig_up], s=4, color="#B91C1C", rasterized=True)
ax.scatter(x[sig_dn], y[sig_dn], s=4, color="#1D4ED8", rasterized=True)
lab_pos = {"MUC5AC": (18, 6), "OAS2": (10, 14), "S100A8": (24, -2), "OSM": (30, -12),
           "NDUFS3": (-26, 10), "CPT2": (-30, -2), "HADH": (-28, -14)}
for sym, off in lab_pos.items():
    row = de[de.symbol == sym]
    if len(row):
        r = row.iloc[0]
        ax.annotate(sym, (r.log2FC_inf_vs_nor, -np.log10(max(r.p, 1e-12))),
                    textcoords="offset points", xytext=off, fontsize=6.5, style="italic",
                    ha="center", arrowprops=dict(arrowstyle="-", lw=0.4, color="#64748B"))
ax.axvline(0, color=GRAY, lw=0.6)
ax.set_xlabel("log$_2$ fold change (inflamed vs non-inflamed)")
ax.set_ylabel("$-$log$_{10}$ $P$")
ax.set_title("a", loc="left", fontsize=12, fontweight="bold")
ax.text(0.98, 0.05, f"{int(sig_up.sum())} up / {int(sig_dn.sum())} down\n(q < 0.05, n = 20 pairs)",
        transform=ax.transAxes, ha="right", fontsize=6.8)
ax.spines[["top", "right"]].set_visible(False)

# b: Factor1 vs inflammation score
ax = axes[1]
idx = [p for p in inflam.index if p in Z.index]
for g, c in [("UC", C_UC), ("CD", C_CD)]:
    sel = [p for p in idx if meta.loc[p, "diagnosis"] == g]
    ax.scatter(Z.loc[sel, "Factor1"], inflam[sel], s=30, color=c, edgecolor="white", lw=0.6, label=g)
r, p = stats.spearmanr(Z.loc[idx, "Factor1"], inflam[idx])
ax.set_xlabel("Factor 1 score")
ax.set_ylabel("Mucosal inflammation score\n(mean log$_2$CPM, top-100 up genes)")
ax.text(0.97, 0.93, f"$\\rho$ = {r:.2f}, $P$ = {p:.1e}", transform=ax.transAxes, ha="right", fontsize=7.5)
ax.legend(frameon=False, fontsize=7.5, loc="lower left")
ax.set_title("b", loc="left", fontsize=12, fontweight="bold")
ax.spines[["top", "right"]].set_visible(False)

# c: alpha diversity by MEFV E148Q
ax = axes[2]
rng = np.random.default_rng(3)
groups = [("carrier", alpha[alpha.MEFV_E148Q_dosage > 0].shannon, "#D97706"),
          ("non-carrier", alpha[alpha.MEFV_E148Q_dosage == 0].shannon, "#64748B")]
for i, (lbl, vals, c) in enumerate(groups):
    ax.boxplot([vals], positions=[i], widths=0.5, showfliers=False,
               medianprops=dict(color="#334155"), boxprops=dict(color=GRAY),
               whiskerprops=dict(color=GRAY), capprops=dict(color=GRAY))
    ax.scatter(rng.normal(i, 0.07, len(vals)), vals, s=26, color=c, edgecolor="white", lw=0.6, zorder=3)
pm = stats.mannwhitneyu(groups[0][1], groups[1][1]).pvalue
ax.set_xticks([0, 1], [f"E148Q carrier\n(n = {len(groups[0][1])})", f"non-carrier\n(n = {len(groups[1][1])})"])
ax.set_ylabel("Shannon diversity (genus level)")
ax.text(0.5, 0.95, f"$P$ = {pm:.3f}", transform=ax.transAxes, ha="center", fontsize=7.5)
ax.set_title("c", loc="left", fontsize=12, fontweight="bold")
ax.spines[["top", "right"]].set_visible(False)

fig.tight_layout()
fig.savefig(os.path.join(FG, "Fig7_mucosal_DE_diversity.png"), bbox_inches="tight")
plt.close(fig)
print("saved Fig7; volcano up/down:", int(sig_up.sum()), int(sig_dn.sum()),
      "| Factor1 rho:", round(r, 3), "| alpha p:", round(pm, 4))
