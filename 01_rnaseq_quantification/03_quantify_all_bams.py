# -*- coding: utf-8 -*-
"""Quantify all RNA-seq BAMs (34 new + 16 re-extracted) with fast_count.py, 3 in parallel."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import config as cfg
import os, glob, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
FAST = os.path.join(HERE, "fast_count.py")
QDIR = str(cfg.QUANT_DIR)
NEW = str(cfg.BAM_DIR_PRIMARY)
OLD = os.path.join(QDIR, "HN00212591_BAM", "BAM_files")

bams = sorted(glob.glob(os.path.join(NEW, "*", "*_sorted.bam"))) + \
       sorted(glob.glob(os.path.join(OLD, "*", "*_sorted.bam")))
print(f"total BAMs: {len(bams)}", flush=True)

def run(bam):
    name = os.path.basename(bam).replace("_sorted.bam", "")
    out = os.path.join(QDIR, f"counts_{name}.tsv")
    if os.path.exists(out) and os.path.getsize(out) > 100_000:
        return f"{name}: SKIP (done)"
    t0 = time.time()
    r = subprocess.run([sys.executable, FAST, bam, out], capture_output=True, text=True)
    msg = r.stdout.strip() or r.stderr.strip()[:200]
    return f"{name}: {msg} [{(time.time()-t0)/60:.1f} min]"

with ThreadPoolExecutor(max_workers=3) as ex:
    for res in ex.map(run, bams):
        print(res, flush=True)
print("ALL_DONE", flush=True)
