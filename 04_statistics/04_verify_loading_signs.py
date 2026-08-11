# -*- coding: utf-8 -*-
"""Verify loading signs used in the manuscript narrative."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import config as cfg
import os, h5py
import pandas as pd
from scipy import stats

DIR = str(cfg.PROJECT_DIR)
MS = os.path.join(DIR, "manuscript")
with h5py.File(os.path.join(MS, "mofa_ms_K8.hdf5"), "r") as f:
    views = [v.decode() for v in f["views/views"][:]]
    samples = [s.decode() for s in f["samples/all"][:]]
    Z = pd.DataFrame(f["expectations/Z/all"][:].T, index=samples,
                     columns=[f"Factor{i+1}" for i in range(f["expectations/Z/all"].shape[0])])
    W = {v: pd.DataFrame(f[f"expectations/W/{v}"][:].T,
                         index=[x.decode() for x in f[f"features/{v}"][:]], columns=Z.columns) for v in views}
cmpd = pd.read_csv(os.path.join(DIR, "metabolome_compound_names.tsv"), sep="\t", index_col=0)["compound_name"].astype(str).to_dict()
meta = pd.read_csv(os.path.join(MS, "meta_with_clinical.tsv"), sep="\t", index_col=0).loc[samples]

wm = W["metabolome"]["Factor2"]
print("cytidine (6175) weight on Factor2:", round(float(wm.get("6175", float('nan'))), 3))
print("Factor2 metabolome extremes:")
t = wm.reindex(wm.abs().sort_values(ascending=False).head(10).index)
for i, v in t.items():
    print(f"   {cmpd.get(str(i), i)[:28]:30s} {v:+.2f}")

print("\nFactor2 vs CRP:", [round(x, 3) for x in stats.spearmanr(Z["Factor2"], meta["CRP"], nan_policy="omit")])
print("Factor1 vs CRP:", [round(x, 3) for x in stats.spearmanr(Z["Factor1"], meta["CRP"], nan_policy="omit")])

# transcriptome Factor1 direction check on a few genes
wt = W["transcriptome"]["Factor1"]
for g in wt.index[:0]:
    pass
sel = [i for i in wt.index if i.split("|")[0] in ("IL11", "OSM", "S100A8", "CLCA1", "SELE", "CSF3")]
print("\nFactor1 transcriptome:")
for i in sel:
    print(f"   {i.split('|')[0]:10s} {wt[i]:+.2f}")

wt3 = W["transcriptome"]["Factor3"]
sel3 = [i for i in wt3.index if i.split("|")[0] in ("IGHG3", "IGHG1", "IGKV3-7")]
print("Factor3 transcriptome (Ig):", {i.split("|")[0]: round(float(wt3[i]), 2) for i in sel3})
wmic3 = W["microbiome"]["Factor3"]
print("Factor3 microbiome:", {i.replace("g__", ""): round(float(wmic3[i]), 2)
                              for i in ["g__Faecalibacterium", "g__Escherichia-Shigella", "g__Fusicatenibacter"]
                              if i in wmic3.index})
