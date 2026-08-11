# -*- coding: utf-8 -*-
"""Mutation view = pathway rare burden + explicit MEFV features (no frequency filter on MEFV)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import config as cfg
import json, os, re
import numpy as np
import pandas as pd

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

# MEFV canonical transcript NP_000234.1
mefv = d["MEFV-NP_000234.1"]
patients = sorted({canon(p) for p in mefv})
rows = {"MEFV_E148Q": {}, "MEFV_L110P": {}, "MEFV_P369S_R408Q": {}, "MEFV_rare_missense": {}}
detail = {}
RARE_SET = {"E84K", "G304R", "G320R", "S503C", "R202Q"}
for pat_raw, rec in mefv.items():
    pat = canon(pat_raw)
    e148q = l110p = p369 = rare = 0
    muts = []
    for m in rec.get("mutations", []):
        aa = f"{m['ref_aa']}{m['position']}{m['mut_aa']}"
        dos = 2 if str(m.get("zygosity", "")).startswith("hom") else 1
        muts.append(f"{aa}({m.get('zygosity')})")
        if aa == "E148Q":
            e148q = dos
        elif aa == "L110P":
            l110p = dos
        elif aa in ("P369S", "R408Q"):
            p369 += dos
        elif aa in RARE_SET:
            rare += dos
    rows["MEFV_E148Q"][pat] = e148q
    rows["MEFV_L110P"][pat] = l110p
    rows["MEFV_P369S_R408Q"][pat] = p369
    rows["MEFV_rare_missense"][pat] = rare
    detail[pat] = "; ".join(muts) if muts else "-"

mefv_df = pd.DataFrame(rows).T[patients]
print("MEFV feature matrix:")
print(mefv_df.T.assign(genotype=[detail[p] for p in patients]).to_string())

# combine with pathway rare burden (raw), z-score everything together
pw = pd.read_csv(os.path.join(DIR, "mutation_rare_missense_pathway_raw.tsv"), sep="\t", index_col=0)
combined = pd.concat([pw, mefv_df.reindex(columns=pw.columns)])
z = combined.sub(combined.mean(axis=1), axis=0).div(combined.std(axis=1).replace(0, np.nan), axis=0)
z = z.dropna(how="all")
z.to_csv(os.path.join(DIR, "view_mutation_rareburden_mefv.tsv"), sep="\t")
pd.Series(detail, name="MEFV_genotype").rename_axis("patient").to_csv(
    os.path.join(DIR, "MEFV_genotypes.tsv"), sep="\t")
print("\nview:", z.shape, "features:", list(z.index))
