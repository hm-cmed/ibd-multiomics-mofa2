# -*- coding: utf-8 -*-
"""Merge per-sample counts -> logCPM matrices; QC vs StringTie TPM; write patient-level RNA view."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import config as cfg
import os, glob, re
import numpy as np
import pandas as pd

BASE = str(cfg.DATA_ROOT)
QDIR = str(cfg.QUANT_DIR)
OUT = str(cfg.PROJECT_DIR)

# ---- merge counts ----
files = sorted(glob.glob(os.path.join(QDIR, "counts_*.tsv")))
cols = {}
for f in files:
    name = os.path.basename(f)[len("counts_"):-len(".tsv")]
    s = pd.read_csv(f, sep="\t", index_col=0).iloc[:, 0]
    cols[name] = s
counts = pd.DataFrame(cols)
print("counts matrix:", counts.shape)

# ---- sample name normalization ----
def norm_sample(s):
    if s == "FMF_SS":
        return None, None
    m = re.match(r"(?i)(nomsap|okh|asic|saga)_?0*(\d+)(inf|nor)$", s)
    if not m:
        return None, None
    fam, num, frac = m.group(1).lower(), int(m.group(2)), m.group(3)
    fam_map = {"nomsap": "NOMSap_{:03d}", "okh": "OKN{:03d}", "asic": "ASIC{:03d}", "saga": "Saga_{:03d}"}
    return fam_map[fam].format(num), frac

renames, patients_frac = {}, {}
for c in counts.columns:
    pid, frac = norm_sample(c)
    if pid is None:
        continue
    renames[c] = f"{pid}__{frac}"
counts = counts[list(renames)].rename(columns=renames)
# de-duplicate (old vs new BAM of same sample: NOMSap_001inf from both deliveries? check)
dup = counts.columns[counts.columns.duplicated()].tolist()
if dup:
    print("duplicated sample columns (keep first):", dup)
    counts = counts.loc[:, ~counts.columns.duplicated()]
print("normalized samples:", counts.shape[1])

# ---- logCPM ----
cpm = counts / counts.sum(axis=0) * 1e6
logcpm = np.log2(cpm + 1)

# ---- QC vs StringTie TPM (16-sample overlap) ----
try:
    tpm_old = pd.read_csv(cfg.VENDOR_TPM, sep="\t", index_col=0)
    tpm_old.index = tpm_old.index.astype(str)
    common_g = logcpm.index.astype(str).intersection(tpm_old.index)
    lc = logcpm.copy(); lc.index = lc.index.astype(str)
    print("QC: genes in common with StringTie profile:", len(common_g))
    for old_c in tpm_old.columns[:16]:
        pid = re.sub(r"(inf|nor)$", "", old_c)
        frac = old_c[-3:]
        m = re.match(r"(?i)nomsap_?0*(\d+)", pid)
        key = f"NOMSap_{int(m.group(1)):03d}__{frac}"
        if key in lc.columns:
            a = np.log2(tpm_old.loc[common_g, old_c] + 1)
            b = lc.loc[common_g, key]
            r = a.corr(b, method="spearman")
            print(f"  {old_c} vs {key}: spearman={r:.3f}")
except Exception as e:
    print("QC failed:", e)

# ---- write full matrices ----
counts.to_csv(os.path.join(OUT, "rna_counts_all_samples.tsv"), sep="\t")
logcpm.to_csv(os.path.join(OUT, "rna_logcpm_all_samples.tsv"), sep="\t")

# ---- patient-level view (inf) ----
inf_cols = [c for c in logcpm.columns if c.endswith("__inf")]
sub_cpm = cpm[inf_cols]
expressed = (sub_cpm >= 1).sum(axis=1) >= int(np.ceil(0.25 * len(inf_cols)))
sub = logcpm.loc[expressed, inf_cols]
hvg = sub.var(axis=1).sort_values(ascending=False).head(2000).index
view = sub.loc[hvg]
view.columns = [c.replace("__inf", "") for c in view.columns]

# gene symbols from expression profile annotation
try:
    anno = pd.read_csv(cfg.GENE_ANNOTATION, sep="\t", index_col=0)
    sym = {str(k): v for k, v in anno["Gene_Symbol"].items()}
    view.index = [f"{sym.get(str(i), str(i))}|{i}" for i in view.index]
except Exception as e:
    print("symbol map failed:", e)

view.to_csv(os.path.join(OUT, "view_transcriptome_logcpm_v2.tsv"), sep="\t")
print("\nRNA v2 view:", view.shape, "patients:", sorted(view.columns.tolist()))
