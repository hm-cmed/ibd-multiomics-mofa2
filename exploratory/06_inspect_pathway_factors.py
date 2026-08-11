# -*- coding: utf-8 -*-
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import config as cfg
import pandas as pd, numpy as np
from scipy import stats
DIR = str(cfg.PROJECT_DIR)
RES = DIR + r"\results_pathway"

Z = pd.read_csv(RES + r"\factors.tsv", sep="\t", index_col=0)
tw = pd.read_excel(RES + r"\top_weights.xlsx", sheet_name="mutation")
print("=== mutation weights on Factor4 ===")
print(tw[tw.factor == "Factor4"][["feature", "weight"]].to_string(index=False))

fam = Z.index.str.extract(r"^(NOMSap|Saga|ASIC|OKN)")[0]
print("\n=== Factor4 score by patient family ===")
print(Z.groupby(fam.values)["Factor4"].agg(["mean", "count"]).round(2).to_string())
kw = stats.kruskal(*[Z.loc[fam.values == g, "Factor4"] for g in ["NOMSap", "Saga"]])
print("NOMSap vs Saga kruskal p =", round(kw.pvalue, 5))

burden = pd.read_csv(DIR + r"\mutation_pathway_burden_raw.tsv", sep="\t", index_col=0)
tot = burden.loc["Total_panel"]
r = Z["Factor4"].corr(tot.reindex(Z.index), method="spearman")
print("\nFactor4 vs Total_panel burden: spearman =", round(r, 3))
