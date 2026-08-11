# -*- coding: utf-8 -*-
"""Rare functional (missense/LoF) allele-dosage burden from Client_aa_mut/all_genes_all_patients.json.
- canonical transcript per gene (smallest NP_, else smallest XP_)
- dosage: hetero=1, homo=2
- rare = cohort carrier count <= 4/40 (10%)
Also writes MEFV known-variant covariates (E148Q, L110P dosage) into sample metadata.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import config as cfg
import json, os, re
import numpy as np
import pandas as pd
from collections import defaultdict, Counter

BASE = str(cfg.DATA_ROOT)
DIR = str(cfg.PROJECT_DIR)
SRC = cfg.AA_MUTATION_JSON
d = json.load(open(SRC, encoding="utf-8"))

canon_re = re.compile(r"^(?:[gn]om[-_ ]?)?(sap|saga|asic|okh|okn)[-_ ]?0*(\d+)", re.I)
def canon(s):
    m = canon_re.search(str(s).lower())
    fam, num = m.group(1), int(m.group(2))
    if fam == "okh":
        fam = "okn"
    return {"sap": f"NOMSap_{num:03d}", "saga": f"Saga_{num:03d}",
            "asic": f"ASIC{num:03d}", "okn": f"OKN{num:03d}"}[fam]

# canonical transcript per gene
tx_by_gene = defaultdict(set)
for gt in d:
    g, tx = gt.split("-", 1)
    tx_by_gene[g].add(tx)
canon_tx = {g: (sorted([t for t in txs if t.startswith("NP_")]) or sorted(txs))[0]
            for g, txs in tx_by_gene.items()}

carrier = defaultdict(dict)   # (gene, aa) -> {patient: dosage}
stats_class = Counter()
trunc = []
for gt, patmap in d.items():
    g, tx = gt.split("-", 1)
    if tx != canon_tx[g]:
        continue
    for pat_raw, rec in patmap.items():
        pat = canon(pat_raw)
        if rec.get("ref_aa_length") != rec.get("patient_aa_length"):
            trunc.append((g, pat, rec.get("ref_aa_length"), rec.get("patient_aa_length")))
        for m in rec.get("mutations", []):
            ref, mut = m["ref_aa"], m["mut_aa"]
            if ref == mut:
                stats_class["same"] += 1
                continue
            cls = "stopgain" if mut in ("*", "X") else ("stoploss" if ref in ("*", "X") else "missense")
            stats_class[cls] += 1
            key = (g, f"{ref}{m['position']}{mut}")
            dosage = 2 if str(m.get("zygosity", "")).startswith("hom") else 1
            carrier[key][pat] = dosage

patients = sorted({canon(p) for patmap in d.values() for p in patmap})
print("patients:", len(patients), "| variant classes:", dict(stats_class),
      "| truncated-length records:", len(trunc))

freq = {k: len(v) for k, v in carrier.items()}
rare = {k for k, n in freq.items() if n <= 4}
common = sorted(((k, n) for k, n in freq.items() if n > 4), key=lambda x: -x[1])
print(f"variants (canonical tx): {len(carrier)} | rare(<=4 carriers): {len(rare)}")
print("most common (filtered out):", [(f"{g}:{aa}", n) for (g, aa), n in common[:10]])

# MEFV detail
mefv_all = sorted((aa, len(v)) for (g, aa), v in carrier.items() if g == "MEFV")
print("MEFV variants (carriers):", mefv_all)

# gene-level rare dosage burden
genes_all = sorted(tx_by_gene)
gene_burden = pd.DataFrame(0.0, index=genes_all, columns=patients)
for (g, aa), pats in carrier.items():
    if (g, aa) in rare:
        for pat, dos in pats.items():
            gene_burden.loc[g, pat] += dos
tot = gene_burden.sum()
print("per-patient rare dosage burden: min/median/max:", tot.min(), tot.median(), tot.max())

PATHWAYS = {
    "Inflammasome_sensor": ["NLRP1", "NLRP3", "NLRP6", "NLRP7", "NLRP9", "NLRP12", "NLRC4",
                            "AIM2", "IFI16", "MEFV", "NLRX1", "NLRC5", "CARD8", "NOD1", "NOD2",
                            "TXNIP", "NEK7"],
    "Pyroptosis_IL1": ["PYCARD", "CASP1", "CASP4", "CASP5", "GSDMD", "IL1B", "IL18", "IL33",
                       "MLKL", "ZBP1"],
    "Apoptosis": ["CASP2", "CASP3", "CASP7", "CASP8", "CASP9", "BCL2", "BCL2L1", "BIRC2",
                  "BIRC3", "XIAP", "FADD", "TRAF2", "TRAF3"],
    "NFkB_signaling": ["IKBKB", "IKBKG", "CHUK", "RELA", "NFKB1", "TNFAIP3", "RIPK2",
                       "MAP3K7", "TAB1", "TAB2", "TAB3", "TRAF6", "MYD88", "IRAK1",
                       "TLR4", "TICAM1", "TNF", "IL6"],
    "NucleicAcid_sensing_IFN": ["STING1", "TBK1", "IKBKE", "IRF3", "IRF7", "MAVS", "RIGI",
                                "OAS1", "RNASEL", "TREX1", "RNASEH2A", "RNASEH2B", "RNASEH2C",
                                "POLR3A", "POLR3B"],
    "IFN_JAK_STAT": ["IFNA1", "IFNA2", "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8", "IFNA10",
                     "IFNA13", "IFNA14", "IFNA16", "IFNA17", "IFNA21", "IFNB1",
                     "JAK1", "TYK2", "STAT1", "STAT2"],
    "MAPK_PI3K": ["MAPK1", "MAPK3", "MAPK8", "MAPK14", "AKT1", "PIK3CA", "PIK3R1"],
    "Cytoskeleton_invasion": ["ACTB", "RHOA", "RAC1", "CDC42", "WASL", "CRK", "NCK1",
                              "BCAR1", "PTK2", "PTK2B", "DNM1L", "CYBB"],
    "Autophagy": ["IRGM", "ATG16L1", "ATG16L2", "NOD2"],
}
rows = {}
for pw, gs in PATHWAYS.items():
    gs = [g for g in gs if g in gene_burden.index]
    rows[pw] = gene_burden.loc[gs].sum(axis=0)
rows["Total_rare"] = gene_burden.sum(axis=0)
pw_burden = pd.DataFrame(rows).T
print("\npathway rare-dosage burden (first 8 patients):")
print(pw_burden.iloc[:, :8].to_string())

z = pw_burden.sub(pw_burden.mean(axis=1), axis=0).div(pw_burden.std(axis=1).replace(0, np.nan), axis=0)
z = z.dropna(how="all")
gene_burden.to_csv(os.path.join(DIR, "mutation_rare_missense_gene_burden.tsv"), sep="\t")
pw_burden.to_csv(os.path.join(DIR, "mutation_rare_missense_pathway_raw.tsv"), sep="\t")
z.to_csv(os.path.join(DIR, "view_mutation_rareburden.tsv"), sep="\t")
print("view:", z.shape)

# MEFV known-variant covariates into metadata
meta = pd.read_csv(os.path.join(DIR, "sample_metadata.tsv"), sep="\t")
for var in ["E148Q", "L110P"]:
    col = {}
    pats = carrier.get(("MEFV", var), {})
    for p in meta["patient"]:
        col[p] = pats.get(p, 0)
    meta[f"MEFV_{var}_dosage"] = meta["patient"].map(col)
meta.to_csv(os.path.join(DIR, "sample_metadata.tsv"), sep="\t", index=False)
print("MEFV covariates added: E148Q carriers =", (meta.MEFV_E148Q_dosage > 0).sum(),
      ", L110P carriers =", (meta.MEFV_L110P_dosage > 0).sum())
