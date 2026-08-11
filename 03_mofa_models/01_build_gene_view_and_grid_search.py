# -*- coding: utf-8 -*-
"""Manuscript pipeline step 1: gene-level rare-missense mutation view + K grid search + final model."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import config as cfg
import os, json
import numpy as np
import pandas as pd

DIR = str(cfg.PROJECT_DIR)
MS = os.path.join(DIR, "manuscript")
os.makedirs(MS, exist_ok=True)
os.makedirs(os.path.join(MS, "figures"), exist_ok=True)
os.makedirs(os.path.join(MS, "tables"), exist_ok=True)

# --- mutation view: gene-level rare missense dosage, genes with >=2 carriers, z-scored ---
gb = pd.read_csv(os.path.join(DIR, "mutation_rare_missense_gene_burden.tsv"), sep="\t", index_col=0)
carriers = (gb > 0).sum(axis=1)
keep = carriers[carriers >= 2].index
mut = gb.loc[keep]
z = mut.sub(mut.mean(axis=1), axis=0).div(mut.std(axis=1), axis=0)
z.to_csv(os.path.join(DIR, "view_mutation_gene_rare.tsv"), sep="\t")
print(f"mutation gene view: {z.shape[0]} genes (>=2 carriers of rare missense) x {z.shape[1]} patients")

# --- load all views ---
views_files = {
    "mutation":      "view_mutation_gene_rare.tsv",
    "microbiome":    "view_microbiome_clr.tsv",
    "metabolome":    "view_metabolome_log2.tsv",
    "transcriptome": "view_transcriptome_logcpm.tsv",
}
likelihoods = ["gaussian"] * 4
mats, feat_names = [], []
samples = None
for name, f in views_files.items():
    df = pd.read_csv(os.path.join(DIR, f), sep="\t", index_col=0)
    if samples is None:
        samples = list(df.columns)
    df = df[samples]
    mats.append([df.T.to_numpy(dtype=float)])
    feat_names.append([str(i) for i in df.index])
    print(f"{name}: D={df.shape[0]}, n_with_data={int(df.notna().any().sum())}")

from mofapy2.run.entry_point import entry_point

def run_model(K, outfile=None):
    ent = entry_point()
    ent.set_data_options(scale_views=True)
    ent.set_data_matrix(mats, likelihoods=likelihoods,
                        views_names=list(views_files.keys()), groups_names=["all"],
                        samples_names=[samples], features_names=feat_names)
    ent.set_model_options(factors=K, spikeslab_weights=True, ard_weights=True)
    ent.set_train_options(iter=2000, convergence_mode="slow", seed=42, verbose=False, freqELBO=5, quiet=True)
    ent.build()
    ent.run()
    if outfile:
        ent.save(outfile)
    # variance explained
    import h5py, tempfile
    tmp = outfile or os.path.join(MS, f"_tmp_K{K}.hdf5")
    if not outfile:
        ent.save(tmp)
    with h5py.File(tmp, "r") as f:
        r2 = f["variance_explained/r2_per_factor/all"][:]
        vnames = [v.decode() for v in f["views/views"][:]]
    if not outfile:
        os.remove(tmp)
    r2df = pd.DataFrame(r2, index=vnames)
    if r2df.values.max() <= 1.0:
        r2df *= 100
    return r2df

grid = {}
for K in [5, 6, 8, 10]:
    r2 = run_model(K)
    tot = r2.sum(axis=1)
    grid[K] = {"per_view_total": tot.to_dict(), "total": float(tot.sum()),
               "n_active_factors": int((r2.sum(axis=0) > 0.5).sum())}
    print(f"K={K}: total R2={tot.sum():.1f} per-view={tot.round(1).to_dict()} active={grid[K]['n_active_factors']}")

with open(os.path.join(MS, "tables", "grid_search_K.json"), "w") as f:
    json.dump(grid, f, indent=1)

# final model: K chosen below by inspection; train K=6 and K=8 saved
for K in [6, 8]:
    out = os.path.join(MS, f"mofa_ms_K{K}.hdf5")
    r2 = run_model(K, outfile=out)
    r2.round(2).to_csv(os.path.join(MS, "tables", f"variance_explained_K{K}.csv"))
    print(f"saved model K={K}")
print("DONE")
