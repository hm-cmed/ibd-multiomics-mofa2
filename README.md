# Four-layer multi-omics integration of an inflammatory bowel disease cohort

Analysis code accompanying *"Multi-omic factors link systemic inflammation and host genotype in
inflammatory bowel disease"*.

No new algorithm was developed in this study. This repository contains the custom code that
prepares the data, configures and runs published software, performs the statistical analyses,
and produces every figure in the paper.

## Data protection

This repository contains **no patient data and no patient identifiers**. Raw whole-genome,
RNA-sequencing, 16S and metabolomics data are under controlled access; see the Data availability
statement of the paper. Seven scripts from the internal analysis tree were withheld from this
release because they contain patient identifiers or perform sample re-identification; none of
them contributes to any reported result. See [`EXCLUDED.md`](EXCLUDED.md) for the full list and
justification.

All filesystem locations are supplied through `config.yaml`, which is git-ignored, so no local
or institutional path is committed.

## Requirements

| Component | Version |
|---|---|
| Python | 3.12.7 |
| numpy | 1.26.4 |
| scipy | 1.13.1 |
| pandas | 2.2.2 |
| matplotlib | 3.9.2 |
| h5py | 3.11.0 |
| mofapy2 | 0.7.4 |
| openpyxl | 3.1.5 |
| PyYAML | any (optional; only if using `config.yaml` rather than environment variables) |
| SAMtools | 1.13 (read counting only) |

```bash
pip install -r requirements.txt
```

### Which MOFA2 component is used

MOFA2 is distributed by the bioFAM group as two components: **`mofapy2`**
(https://github.com/bioFAM/mofapy2), the reference Python implementation containing the
variational inference engine, and **`MOFA2`** (Bioconductor), an R interface that calls
`mofapy2` internally through reticulate/basilisk. All models reported in the paper were trained
by calling `mofapy2.run.entry_point` directly from Python — the same inference engine the R
package uses. `mofapy2` is installed from PyPI and is not redistributed here.

## Configuration

```bash
cp config.example.yaml config.yaml   # then edit the paths
```

Alternatively set the environment variables listed in `config.py`.

**To reproduce the published results you only need `project_dir`.** Point it at a directory
containing the processed matrices distributed as Supplementary Data 1 of the paper, and every
step from `03_mofa_models` onward will run. The remaining paths are required only to regenerate
those matrices from controlled-access raw data.

### Bundled reference data

`reference/gene_annotation.tsv` (46,428 genes; Gene_ID, symbol, description, biotype, HGNC and
Ensembl identifiers) is derived from the public NCBI RefSeq GRCh38 release 109 annotation and is
included in this repository because several scripts require it to map gene identifiers to
symbols. It contains no patient data. To use a different annotation build, place your own
`gene_annotation.tsv` in `project_dir`; it takes precedence over the bundled copy.

## Reproducing the analysis

From the processed matrices (Supplementary Data 1):

```bash
python 03_mofa_models/01_build_gene_view_and_grid_search.py
python 03_mofa_models/02_train_single_K.py 8 save
python 04_statistics/01_factor_clinical_associations.py
python 04_statistics/02_permanova_pcoa_univariate.py
python 04_statistics/03_paired_DE_diversity_robustness.py
python 05_figures/01_main_and_supplementary_figures.py
python 05_figures/02_figure7_mucosal_DE_diversity.py
```

From controlled-access raw data, prepend:

```bash
python 01_rnaseq_quantification/01_build_union_exon_bed.py
python 01_rnaseq_quantification/03_quantify_all_bams.py
python 01_rnaseq_quantification/04_merge_counts_and_qc.py
python 02_view_construction/01_build_four_views.py
python 02_view_construction/03_mutation_rare_missense_burden.py
```

MOFA2 models are trained with a fixed random seed (42), and permutation-based tests use a fixed
NumPy generator seed (42), so all reported values are reproducible.

## Directory structure

### `01_rnaseq_quantification/`
| Script | Purpose |
|---|---|
| `01_build_union_exon_bed.py` | Builds union-exon models per gene from the RefSeq GTF |
| `02_count_reads_per_gene.py` | Counts primary alignments (MAPQ >= 10, secondary/supplementary excluded) per gene for one BAM, streaming from `samtools view` |
| `03_quantify_all_bams.py` | Runs the counter over all BAMs, three in parallel |
| `04_merge_counts_and_qc.py` | Merges per-sample counts, normalizes to log2(CPM+1), verifies concordance with the vendor quantification, writes the transcriptome view |

### `02_view_construction/`
| Script | Purpose |
|---|---|
| `01_build_four_views.py` | Assembles the microbiome (genus CLR), metabolome (log2) and mutation views at patient level, harmonizes identifiers across omics layers, writes sample metadata |
| `03_mutation_rare_missense_burden.py` | **Mutation view used in the primary model**: rare (cohort carrier frequency <= 10%) missense allele-dosage burden per gene |
| `04_mutation_view_with_MEFV.py` | Mutation view with explicit *MEFV* dosage features, for the reported sensitivity model |

### `03_mofa_models/`
| Script | Purpose |
|---|---|
| `01_build_gene_view_and_grid_search.py` | Builds the gene-level rare-burden view and runs the K grid search |
| `02_train_single_K.py` | Trains and saves one MOFA2 model for a given K (used for K = 5, 6, 8, 10) |
| `03_MEFV_sensitivity_model.py` | Sensitivity model with explicit *MEFV* dosage features |

### `04_statistics/`
| Script | Purpose |
|---|---|
| `01_factor_clinical_associations.py` | Factor associations with diagnosis, CRP and measurement batch; loading tables |
| `02_permanova_pcoa_univariate.py` | PERMANOVA, PCoA, exhaustive gene-by-genus Spearman testing, Mantel test, genotype-stratified comparisons |
| `03_paired_DE_diversity_robustness.py` | Paired inflamed versus non-inflamed differential expression, mucosal inflammation score, alpha diversity, K = 6 versus K = 8 factor matching, complete-case sensitivity model |
| `04_verify_loading_signs.py` | Verifies the sign conventions quoted in the manuscript |
| `05_postprocess_MEFV_model.py`, `06_inspect_MEFV_factors.py` | Post-processing and inspection of the *MEFV* sensitivity model |

### `05_figures/`
| Script | Purpose |
|---|---|
| `01_main_and_supplementary_figures.py` | Figures 1–6 and Supplementary Figures 1–3 |
| `02_figure7_mucosal_DE_diversity.py` | Figure 7 |

### `exploratory/`
Alternative encodings of the mutation view evaluated during model selection that produced **no
result reported in the paper**. Retained for transparency; nothing in the main analysis depends
on them. See [`exploratory/README.md`](exploratory/README.md).

### `reference/`
`gene_annotation.tsv` — see "Bundled reference data" above.

### Root
`config.py`, `config.example.yaml`, `requirements.txt`, and
`report_software_versions.py` (prints the versions of all dependencies actually installed).

## Licence

[MIT](LICENSE) — see `LICENSE`. `mofapy2` is a separate work distributed by its authors under
LGPL-3.0 and is not included in this repository.

## Citation

If you use this code, please cite the accompanying paper and the MOFA2 methodology papers:

- Argelaguet, R. et al. Multi-Omics Factor Analysis — a framework for unsupervised integration of multi-omics data sets. *Mol. Syst. Biol.* **14**, e8124 (2018).
- Argelaguet, R. et al. MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data. *Genome Biol.* **21**, 111 (2020).
