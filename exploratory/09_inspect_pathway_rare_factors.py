# -*- coding: utf-8 -*-
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import config as cfg
import pandas as pd, numpy as np
from scipy import stats
DIR = str(cfg.PROJECT_DIR)
RES = DIR + r"\results_rareburden"

Z = pd.read_csv(RES + r"\factors.tsv", sep="\t", index_col=0).drop(columns=["diagnosis"])
meta = pd.read_csv(DIR + r"\sample_metadata.tsv", sep="\t").set_index("patient").loc[Z.index]

for fac in ["Factor2", "Factor4"]:
    print(f"===== {fac} =====")
    for view in ["mutation", "microbiome", "transcriptome"]:
        tw = pd.read_excel(RES + r"\top_weights.xlsx", sheet_name=view)
        t = tw[tw.factor == fac].head(6)
        print(f"-- {view}")
        for _, r in t.iterrows():
            print(f"   {str(r['feature'])[:40]:42s} w={r['weight']:+.2f}")

print("\n=== factor vs MEFV covariates (spearman, p) ===")
for cov in ["MEFV_E148Q_dosage", "MEFV_L110P_dosage"]:
    for fac in Z.columns:
        r, p = stats.spearmanr(Z[fac], meta[cov])
        if p < 0.1:
            print(f"  {cov} ~ {fac}: rho={r:+.2f} p={p:.4f}")
print("(p<0.1のみ表示)")

# rare total burden vs factors
burden = pd.read_csv(DIR + r"\mutation_rare_missense_pathway_raw.tsv", sep="\t", index_col=0)
tot = burden.loc["Total_rare"].reindex(Z.index)
for fac in Z.columns:
    r, p = stats.spearmanr(Z[fac], tot)
    if p < 0.05:
        print(f"Total_rare_burden ~ {fac}: rho={r:+.2f} p={p:.4f}")
