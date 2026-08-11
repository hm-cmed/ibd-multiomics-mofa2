# -*- coding: utf-8 -*-
"""Aggregate the 114-gene mutation view into pathway-level burden (mean variants/gene, z-scored)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import config as cfg
import os
import numpy as np
import pandas as pd

DIR = str(cfg.PROJECT_DIR)
cnt = pd.read_csv(os.path.join(DIR, "view_mutation_count.tsv"), sep="\t", index_col=0)

PATHWAYS = {
    # inflammasome sensors / pyrin (FMF関連含む)
    "Inflammasome_sensor": ["NLRP1", "NLRP3", "NLRP6", "NLRP7", "NLRP9", "NLRP12", "NLRC4",
                            "AIM2", "IFI16", "MEFV", "NLRX1", "NLRC5", "CARD8", "NOD1", "NOD2",
                            "TXNIP", "NEK7"],
    # inflammasome effectors / pyroptosis
    "Pyroptosis_IL1": ["PYCARD", "CASP1", "CASP4", "CASP5", "GSDMD", "IL1B", "IL18", "IL33",
                       "MLKL", "ZBP1"],
    # apoptosis
    "Apoptosis": ["CASP2", "CASP3", "CASP7", "CASP8", "CASP9", "BCL2", "BCL2L1", "BIRC2",
                  "BIRC3", "XIAP", "FADD", "TRAF2", "TRAF3"],
    # NF-kB axis (TLR/NOD downstream)
    "NFkB_signaling": ["IKBKB", "IKBKG", "CHUK", "RELA", "NFKB1", "TNFAIP3", "RIPK2",
                       "MAP3K7", "TAB1", "TAB2", "TAB3", "TRAF6", "MYD88", "IRAK1",
                       "TLR4", "TICAM1", "TNF", "IL6"],
    # cGAS-STING / RLR / nucleic-acid sensing & type-I IFN induction
    "NucleicAcid_sensing_IFN": ["STING1", "TBK1", "IKBKE", "IRF3", "IRF7", "MAVS", "RIGI",
                                "OAS1", "RNASEL", "TREX1", "RNASEH2A", "RNASEH2B", "RNASEH2C",
                                "POLR3A", "POLR3B"],
    # IFN ligands + JAK/STAT
    "IFN_JAK_STAT": ["IFNA1", "IFNA2", "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8", "IFNA10",
                     "IFNA13", "IFNA14", "IFNA16", "IFNA17", "IFNA21", "IFNB1",
                     "JAK1", "TYK2", "STAT1", "STAT2"],
    # MAPK / PI3K-AKT
    "MAPK_PI3K": ["MAPK1", "MAPK3", "MAPK8", "MAPK14", "AKT1", "PIK3CA", "PIK3R1"],
    # cytoskeleton / bacterial invasion (Yersinia adhesion module)
    "Cytoskeleton_invasion": ["ACTB", "RHOA", "RAC1", "CDC42", "WASL", "CRK", "NCK1",
                              "BCAR1", "PTK2", "PTK2B", "DNM1L", "CYBB"],
    # autophagy
    "Autophagy": ["IRGM", "ATG16L2", "NOD2"],
}

used = sorted({g for gs in PATHWAYS.values() for g in gs})
missing = [g for g in used if g not in cnt.index]
uncovered = [g for g in cnt.index if g not in used]
print("panel genes:", len(cnt.index), "| in pathways:", len(used),
      "| not in count matrix:", missing, "| panel genes unassigned:", uncovered)

rows = {}
for pw, gs in PATHWAYS.items():
    gs = [g for g in gs if g in cnt.index]
    sub = cnt.loc[gs]
    n_avail = sub.notna().sum(axis=0)
    burden = sub.mean(axis=0, skipna=True)          # mean variants per gene
    burden[n_avail < 0.5 * len(gs)] = np.nan        # too few genes measured
    rows[pw] = burden
rows["Total_panel"] = cnt.mean(axis=0, skipna=True)
burden = pd.DataFrame(rows).T

# z-score each pathway across patients
z = burden.sub(burden.mean(axis=1), axis=0).div(burden.std(axis=1), axis=0)
z.to_csv(os.path.join(DIR, "view_mutation_pathway.tsv"), sep="\t")
burden.round(3).to_csv(os.path.join(DIR, "mutation_pathway_burden_raw.tsv"), sep="\t")
print("\nburden matrix:", z.shape)
print(burden.round(2).iloc[:, :6].to_string())
