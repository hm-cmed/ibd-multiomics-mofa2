# -*- coding: utf-8 -*-
"""Train one MOFA model for the manuscript. Usage: ms_run_one_k.py <K> [save]"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import config as cfg
import os, sys, json
import pandas as pd

K = int(sys.argv[1])
SAVE = len(sys.argv) > 2 and sys.argv[2] == "save"
DIR = str(cfg.PROJECT_DIR)
MS = os.path.join(DIR, "manuscript")

views_files = {"mutation": "view_mutation_gene_rare.tsv", "microbiome": "view_microbiome_clr.tsv",
               "metabolome": "view_metabolome_log2.tsv", "transcriptome": "view_transcriptome_logcpm.tsv"}
mats, feat_names, samples = [], [], None
for name, f in views_files.items():
    df = pd.read_csv(os.path.join(DIR, f), sep="\t", index_col=0)
    if samples is None:
        samples = list(df.columns)
    mats.append([df[samples].T.to_numpy(dtype=float)])
    feat_names.append([str(i) for i in df.index])

from mofapy2.run.entry_point import entry_point
ent = entry_point()
ent.set_data_options(scale_views=True)
ent.set_data_matrix(mats, likelihoods=["gaussian"] * 4, views_names=list(views_files),
                    groups_names=["all"], samples_names=[samples], features_names=feat_names)
ent.set_model_options(factors=K, spikeslab_weights=True, ard_weights=True)
ent.set_train_options(iter=2000, convergence_mode="slow", seed=42, verbose=False, freqELBO=5, quiet=True)
ent.build()
ent.run()
out = os.path.join(MS, f"mofa_ms_K{K}.hdf5")
ent.save(out)

import h5py
with h5py.File(out, "r") as f:
    r2 = f["variance_explained/r2_per_factor/all"][:]
    vnames = [v.decode() for v in f["views/views"][:]]
r2df = pd.DataFrame(r2, index=vnames, columns=[f"Factor{i+1}" for i in range(r2.shape[1])])
if r2df.values.max() <= 1.0:
    r2df *= 100
r2df.round(3).to_csv(os.path.join(MS, "tables", f"variance_explained_K{K}.csv"))
tot = r2df.sum(axis=1)
res = {"K": K, "per_view_total": {k: round(v, 2) for k, v in tot.items()},
       "total": round(float(tot.sum()), 2),
       "n_active_factors": int((r2df.sum(axis=0) > 0.5).sum())}
print(json.dumps(res))
if not SAVE:
    os.remove(out)
