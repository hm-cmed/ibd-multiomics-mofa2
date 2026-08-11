# -*- coding: utf-8 -*-
"""Assemble 4-view MOFA2 input (mutation / microbiome / metabolome / transcriptome)
at patient level for the 40-patient cohort."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import config as cfg
import json, os, re
import numpy as np
import pandas as pd

BASE = str(cfg.DATA_ROOT)
MUT = cfg.MUTATION_DIR
OUT = str(cfg.PROJECT_DIR)
os.makedirs(OUT, exist_ok=True)

# ---------- canonical patient IDs ----------
PAT = re.compile(r"(?:[gn]om[-_ ]?)?(sap|saga|asic|okh|okn)[-_ ]?0*(\d+)", re.I)

def canon(idstr):
    m = PAT.search(str(idstr).lower())
    if not m:
        return None
    fam, num = m.group(1), int(m.group(2))
    if fam == "okh":
        fam = "okn"
    return {"sap": f"NOMSap_{num:03d}", "saga": f"Saga_{num:03d}",
            "asic": f"ASIC{num:03d}", "okn": f"OKN{num:03d}"}[fam]

# ---------- cohort & diagnosis ----------
import openpyxl
wb = openpyxl.load_workbook(cfg.COHORT_XLSX, data_only=True)
ws = wb.worksheets[0]
diag = {}
for row in ws.iter_rows(min_row=2, values_only=True):
    if row[0]:
        diag[canon(row[0])] = row[1]
mut_raw = json.load(open(os.path.join(MUT, "client_gene_mut.json"), encoding="utf-8"))
patients = sorted({canon(k) for k in mut_raw})
print("cohort:", len(patients), "patients")

# ---------- view 1: mutation ----------
genes = sorted({g for v in mut_raw.values() for g in v})
mut_count = pd.DataFrame(index=genes, columns=patients, dtype=float)
for pid_raw, gd in mut_raw.items():
    pid = canon(pid_raw)
    for g, d in gd.items():
        v = d.get("variants", "N/A")
        mut_count.loc[g, pid] = np.nan if str(v) == "N/A" else float(v)
mut_bin = (mut_count > 0).astype(float).where(mut_count.notna())
print(f"mutation: {mut_count.shape}, NA rate {mut_count.isna().mean().mean():.3f}")

# ---------- view 2: microbiome (genus level, CLR) ----------
m1 = json.load(open(os.path.join(MUT, "microbiome1.json"), encoding="utf-8"))
m2 = json.load(open(os.path.join(MUT, "microbiome2.json"), encoding="utf-8"))

def to_genus_table(mjson):
    """taxonomy-string -> {sample: count}  ==> genus x sample count matrix"""
    rows = {}
    for taxon, vals in mjson.items():
        m = re.search(r"(g__[^;]+)", taxon)
        if not m:
            continue
        genus = m.group(1)
        for s, v in vals.items():
            if v is None:
                continue
            rows.setdefault(genus, {}).setdefault(s, 0.0)
            rows[genus][s] += float(v)
    return pd.DataFrame(rows).T  # genus x sample

g1 = to_genus_table(m1)
g2 = to_genus_table(m2)
common = sorted(set(g1.index) & set(g2.index))
g1, g2 = g1.loc[common].fillna(0), g2.loc[common].fillna(0)
print(f"microbiome: batch1 {g1.shape}, batch2 {g2.shape}, common genera {len(common)}")

micro_sample = {}   # patient -> (sample_key, batch)
for s in g1.columns:
    micro_sample.setdefault(canon(s), (s, 1))          # batch1 preferred
for s in g2.columns:
    micro_sample.setdefault(canon(s), (s, 2))

def clr(counts):
    rel = (counts + 0.5) / (counts + 0.5).sum(axis=0)
    logrel = np.log(rel)
    return logrel - logrel.mean(axis=0)

c1, c2 = clr(g1), clr(g2)
micro = pd.DataFrame(index=common, columns=patients, dtype=float)
for pid, (s, b) in micro_sample.items():
    if pid in micro.columns:
        micro[pid] = (c1 if b == 1 else c2)[s]
# prevalence filter: genus with nonzero count in >=25% of used samples
used = [micro_sample[p][0] for p in patients if p in micro_sample]
raw = pd.concat([g1, g2], axis=1)[[s for s, b in micro_sample.values()]]
prev = (raw > 0).mean(axis=1)
keep = prev[prev >= 0.25].index
micro = micro.loc[sorted(keep)]
n_micro = micro.notna().any(axis=0).sum()
print(f"microbiome view: {micro.shape}, patients with data {n_micro}")

# ---------- view 3: metabolome ----------
mm = json.load(open(os.path.join(MUT, "m3metabol.json"), encoding="utf-8"))

# batch3: SPMUM-HMT-020 (2026-08 発見分) 表3 定量結果, 17 new samples
HMT020 = cfg.METABOLOME_BATCH3_XLSX
wb3 = openpyxl.load_workbook(HMT020, read_only=True, data_only=True)
ws3 = wb3["表3 定量結果"]
rows3 = list(ws3.iter_rows(values_only=True))
samp3 = [str(x) for x in rows3[3][4:21] if x]  # 17 sample columns
b3 = {}  # cid -> {sample: value}
for r in rows3[4:]:
    if not r[0] or r[2] is None:
        continue
    cid = str(r[2]).split(".")[0]
    vals = {}
    for j, s in enumerate(samp3):
        v = r[4 + j]
        if v is None or str(v).strip() in ("N.D.", "N.A.", "-", ""):
            vals[s] = 0.0
        else:
            try:
                vals[s] = float(v)
            except ValueError:
                vals[s] = np.nan
    b3[cid] = vals
common_cids = sorted(set(mm.keys()) & set(b3.keys()))
print(f"metabolome: batch1+2 {len(mm)} CIDs, batch3(HMT-020) {len(b3)} CIDs x {len(samp3)} samples, common {len(common_cids)}")

metab = pd.DataFrame(index=common_cids, columns=patients, dtype=float)
metab_batch = {}
for cid in common_cids:
    for s, v in mm[cid].items():
        pid = canon(s)
        if pid in metab.columns and v is not None:
            metab.loc[cid, pid] = np.log2(float(v) + 1)
            metab_batch.setdefault(pid, "HMT_1_2")
    for s, v in b3[cid].items():
        pid = canon(s)
        if pid in metab.columns and not np.isnan(v):
            if metab_batch.get(pid) == "HMT_1_2":
                continue  # keep earlier batch value if both exist
            metab.loc[cid, pid] = np.log2(v + 1)
            metab_batch.setdefault(pid, "HMT_020")
# compound names
cid2name = {}
for f in [os.path.join(BASE, "metobolome", "metabolome_content1.csv"),
          os.path.join(BASE, "metobolome", "metabolome_content2.csv")]:
    try:
        t = pd.read_csv(f)
        namecol = "Compound name"
        cidcol = "PubChem CID"
        for _, r in t.iterrows():
            cid = str(r[cidcol]).split(".")[0]
            if cid and cid != "nan":
                cid2name.setdefault(cid, str(r[namecol]))
    except Exception as e:
        print("name map fail", f, e)
n_metab = metab.notna().any(axis=0).sum()
print(f"metabolome view: {metab.shape}, patients with data {n_metab}")

# ---------- view 4: transcriptome (v2: uniform re-quantification of 49 BAMs, inf, logCPM) ----------
rna_v2 = pd.read_csv(os.path.join(OUT, "view_transcriptome_logcpm_v2.tsv"),
                     sep="\t", index_col=0)
rna_view = pd.DataFrame(index=rna_v2.index, columns=patients, dtype=float)
rna_sample = {}
for c in rna_v2.columns:
    if c in rna_view.columns:
        rna_view[c] = rna_v2[c]
        rna_sample[c] = c + "__inf"
n_rna = rna_view.notna().any(axis=0).sum()
print(f"transcriptome view (v2): {rna_view.shape}, patients with data {n_rna}"
      f" (out-of-cohort dropped: {sorted(set(rna_v2.columns) - set(patients))})")

# ---------- metadata & availability ----------
meta = pd.DataFrame({"patient": patients})
meta["diagnosis"] = meta["patient"].map(diag)
meta["has_mutation"] = meta["patient"].isin(mut_count.columns[mut_count.notna().any()])
meta["has_microbiome"] = meta["patient"].map(lambda p: micro[p].notna().any())
meta["microbiome_sample"] = meta["patient"].map(lambda p: micro_sample.get(p, ("", ""))[0])
meta["microbiome_batch"] = meta["patient"].map(lambda p: micro_sample.get(p, ("", ""))[1])
meta["has_metabolome"] = meta["patient"].map(lambda p: metab[p].notna().any())
meta["metabolome_batch"] = meta["patient"].map(lambda p: metab_batch.get(p, ""))
meta["has_transcriptome"] = meta["patient"].map(lambda p: rna_view[p].notna().any())
meta["rna_sample"] = meta["patient"].map(lambda p: rna_sample.get(p, ""))
meta["n_views"] = meta[["has_mutation", "has_microbiome", "has_metabolome", "has_transcriptome"]].sum(axis=1)

# ---------- write ----------
mut_count.to_csv(os.path.join(OUT, "view_mutation_count.tsv"), sep="\t")
mut_bin.to_csv(os.path.join(OUT, "view_mutation_binary.tsv"), sep="\t")
micro.to_csv(os.path.join(OUT, "view_microbiome_clr.tsv"), sep="\t")
metab.to_csv(os.path.join(OUT, "view_metabolome_log2.tsv"), sep="\t")
rna_view.to_csv(os.path.join(OUT, "view_transcriptome_logcpm.tsv"), sep="\t")
meta.to_csv(os.path.join(OUT, "sample_metadata.tsv"), sep="\t", index=False)
pd.Series(cid2name, name="compound_name").rename_axis("PubChem_CID").to_csv(
    os.path.join(OUT, "metabolome_compound_names.tsv"), sep="\t")

print("\n=== availability summary ===")
print(meta.groupby("n_views").size())
print("\n4 views complete:", meta.loc[meta.n_views == 4, "patient"].tolist())
print("3+ views:", meta.loc[meta.n_views >= 3, "patient"].tolist())
print("\npairwise overlap:")
for a in ["has_mutation", "has_microbiome", "has_metabolome", "has_transcriptome"]:
    row = []
    for b in ["has_mutation", "has_microbiome", "has_metabolome", "has_transcriptome"]:
        row.append(int((meta[a] & meta[b]).sum()))
    print(a.replace("has_", "").ljust(14), row)
print("\nsaved to", OUT)
