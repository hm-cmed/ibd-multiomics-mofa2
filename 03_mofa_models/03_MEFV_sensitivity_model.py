# -*- coding: utf-8 -*-
"""Run MOFA2 (mofapy2) on the 4-view dataset."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import config as cfg
import os
import numpy as np
import pandas as pd

DIR = str(cfg.PROJECT_DIR)

views = {
    "mutation":      "view_mutation_rareburden_mefv.tsv",
    "microbiome":    "view_microbiome_clr.tsv",
    "metabolome":    "view_metabolome_log2.tsv",
    "transcriptome": "view_transcriptome_logcpm.tsv",
}
likelihoods = ["gaussian", "gaussian", "gaussian", "gaussian"]

mats, feat_names = [], []
samples = None
for name, f in views.items():
    df = pd.read_csv(os.path.join(DIR, f), sep="\t", index_col=0)
    if samples is None:
        samples = list(df.columns)
    df = df[samples]
    mats.append([df.T.to_numpy(dtype=float)])  # [group][samples x features]
    feat_names.append([f"{name}_{i}" if not str(i).strip() else str(i) for i in df.index])
    print(f"{name}: {df.shape[0]} features x {df.shape[1]} samples, "
          f"missing samples: {int(df.isna().all().sum())}")

from mofapy2.run.entry_point import entry_point
ent = entry_point()
ent.set_data_options(scale_views=True)
ent.set_data_matrix(
    mats,
    likelihoods=likelihoods,
    views_names=list(views.keys()),
    groups_names=["all"],
    samples_names=[samples],
    features_names=feat_names,
)
ent.set_model_options(factors=8, spikeslab_weights=True, ard_weights=True)
ent.set_train_options(iter=2000, convergence_mode="slow", seed=42, verbose=False, freqELBO=5)
ent.build()
ent.run()
out = os.path.join(DIR, "mofa_model_rareburden_mefv.hdf5")
ent.save(out)
print("saved:", out)
