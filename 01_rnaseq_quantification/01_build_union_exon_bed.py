# -*- coding: utf-8 -*-
"""Union-exon BED (chr, start0, end, gene_id) from the NCBI GTF, for samtools bedcov."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import config as cfg
import re, os
from collections import defaultdict

GTF = str(cfg.GTF)
OUT = str(cfg.EXON_BED)
os.makedirs(os.path.dirname(OUT), exist_ok=True)

gid_re = re.compile(r'gene_id "([^"]+)"')
exons = defaultdict(list)  # (chr, gid) -> [(s,e)]
with open(GTF, encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.split("\t", 9)
        if len(p) < 9 or p[2] != "exon":
            continue
        m = gid_re.search(p[8])
        if not m:
            continue
        exons[(p[0], m.group(1))].append((int(p[3]) - 1, int(p[4])))

n_genes = 0
with open(OUT, "w", newline="\n") as w:
    for (chrom, gid), ivs in sorted(exons.items()):
        ivs.sort()
        merged = []
        for s, e in ivs:
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        for s, e in merged:
            w.write(f"{chrom}\t{s}\t{e}\t{gid}\n")
        n_genes += 1
print("genes:", n_genes, "->", OUT)
