# -*- coding: utf-8 -*-
"""Fast gene-level read counting: samtools view stream -> pandas chunks -> searchsorted.
Assigns each primary aligned read (MAPQ>=10) to a gene by its start position within
union-exon intervals. Ambiguous (interval of 2 genes overlapping) resolved by 'first'.
Usage: python fast_count.py <bam_windows_path> <out_counts_tsv>
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import config as cfg
import sys, os, subprocess
import numpy as np
import pandas as pd
from collections import defaultdict

BED = str(cfg.EXON_BED)
bam_win, out_tsv = sys.argv[1], sys.argv[2]

def to_wsl(p):
    p = p.replace("\\", "/")
    return "/mnt/" + p[0].lower() + p[2:]

# --- load intervals per chromosome ---
bed = pd.read_csv(BED, sep="\t", header=None, names=["chr", "s", "e", "gid"])
chrom_iv = {}
gene_ids = []
gid_index = {}
for g in bed["gid"].unique():
    gid_index[g] = len(gene_ids)
    gene_ids.append(g)
for c, sub in bed.groupby("chr", sort=False):
    sub = sub.sort_values("s")
    starts = sub["s"].to_numpy(np.int64)
    ends = sub["e"].to_numpy(np.int64)
    gids = sub["gid"].map(gid_index).to_numpy(np.int64)
    # drop intervals fully contained in previous (keep max end monotone)
    keep_e = np.maximum.accumulate(ends)
    chrom_iv[c] = (starts, keep_e, ends, gids)

counts = np.zeros(len(gene_ids), dtype=np.int64)

cmd = ["wsl.exe", "-d", "Ubuntu-22.04", "--", "samtools", "view",
       "-q", "10", "-F", "0x904", to_wsl(bam_win)]
proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=1 << 22)

reader = pd.read_csv(proc.stdout, sep="\t", header=None, usecols=[2, 3],
                     names=["chr", "pos"], chunksize=2_000_000,
                     dtype={"chr": "category", "pos": np.int64},
                     engine="c", quoting=3, on_bad_lines="skip")
total = 0
for chunk in reader:
    total += len(chunk)
    for c, sub in chunk.groupby("chr", observed=True):
        iv = chrom_iv.get(str(c))
        if iv is None:
            continue
        starts, keep_e, ends, gids = iv
        pos0 = sub["pos"].to_numpy(np.int64) - 1  # BAM POS is 1-based
        idx = np.searchsorted(starts, pos0, side="right") - 1
        ok = idx >= 0
        idxc = np.clip(idx, 0, None)
        inside = ok & (pos0 < ends[idxc])
        hit_gids = gids[idxc[inside]]
        np.add.at(counts, hit_gids, 1)
proc.wait()
if proc.returncode not in (0, None):
    print("WARN samtools exit", proc.returncode, file=sys.stderr)

s = pd.Series(counts, index=gene_ids, name=os.path.basename(bam_win).replace("_sorted.bam", ""))
s.rename_axis("Gene_ID").to_csv(out_tsv, sep="\t")
print(f"{os.path.basename(bam_win)}: reads={total:,} assigned={counts.sum():,} ({counts.sum()/max(total,1):.1%})")
