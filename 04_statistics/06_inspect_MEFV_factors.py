# -*- coding: utf-8 -*-
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import config as cfg
import pandas as pd, numpy as np
from scipy import stats
DIR = str(cfg.PROJECT_DIR)
RES = DIR + r"\results_rareburden_mefv"

Z = pd.read_csv(RES + r"\factors.tsv", sep="\t", index_col=0).drop(columns=["diagnosis"])
meta = pd.read_csv(DIR + r"\sample_metadata.tsv", sep="\t").set_index("patient").loc[Z.index]
view = pd.read_csv(DIR + r"\view_mutation_rareburden_mefv.tsv", sep="\t", index_col=0)

print("=== mutation weights on Factor4 (all 14 features) ===")
tw = pd.read_excel(RES + r"\top_weights.xlsx", sheet_name="mutation")
t4 = tw[tw.factor == "Factor4"]
print(t4[["feature", "weight"]].to_string(index=False))

print("\n=== Factor4 co-loading: metabolome / transcriptome top weights ===")
for v in ["metabolome", "transcriptome"]:
    tv = pd.read_excel(RES + r"\top_weights.xlsx", sheet_name=v)
    t = tv[tv.factor == "Factor4"].head(6)
    for _, r in t.iterrows():
        nm = r["name"] if isinstance(r["name"], str) and r["name"] not in ("", "nan") else str(r["feature"])
        print(f"  [{v}] {nm[:40]:42s} w={r['weight']:+.2f}")

print("\n=== MEFV features vs all factors (spearman, p<0.1) ===")
for feat in ["MEFV_E148Q", "MEFV_L110P", "MEFV_P369S_R408Q", "MEFV_rare_missense"]:
    x = view.loc[feat].reindex(Z.index)
    for fac in Z.columns:
        r, p = stats.spearmanr(Z[fac], x)
        if p < 0.1:
            print(f"  {feat:20s} ~ {fac}: rho={r:+.2f} p={p:.4f}")

print("\n=== MEFV variants vs diagnosis (UC/CD) ===")
dx = meta["diagnosis"]
for feat in ["MEFV_E148Q", "MEFV_L110P", "MEFV_P369S_R408Q"]:
    x = view.loc[feat].reindex(Z.index)
    tab = pd.crosstab(x > 0, dx)
    if tab.shape == (2, 2):
        odds, p = stats.fisher_exact(tab.values)
        print(f"  {feat}: carriers UC={int(tab.loc[True].get('UC',0))}/{int((dx=='UC').sum())}, "
              f"CD={int(tab.loc[True].get('CD',0))}/{int((dx=='CD').sum())}, Fisher p={p:.3f}")
