# Scripts excluded from this public release

Seven scripts from the internal analysis tree are deliberately **not** included here.
All of them fall into one of two categories: they contain patient identifiers or encode the
internal structure of the institutional data archive, and **none of them contributes to any
result reported in the paper**. Excluding them therefore removes every patient identifier from
this repository without affecting reproducibility.

| Excluded script | Reason |
|---|---|
| `01_data_curation/01_locate_omics_files.py` | Contains patient identifiers and encodes the internal layout of the institutional data archive. Used only to locate files on disk; contributes to no reported result. |
| `01_data_curation/02_build_data_inventory.py` | Contains patient identifiers. Produces an internal data-availability inventory; contributes to no reported result. |
| `01_data_curation/03_genotype_concordance_ST_vs_cohort.py` | Performs genotype-based re-identification of anonymously coded whole-genome samples and contains patient identifiers. The samples it resolved were **not** included in the analysed cohort, so it contributes to no reported result. Withheld on data-protection grounds. |
| `01_data_curation/04_check_metabolome_sample_ids.py` | Audits sample identifiers within raw vendor deliveries; contributes to no reported result. |
| `01_rnaseq_quantification/05_legacy_prep_ST_cohort.py` | Prepares a vendor-quantified delivery from a separate sample set that was not analysed in this study. |
| `01_rnaseq_quantification/06_legacy_prep_nomsap_cohort.py` | Contains patient identifiers. Superseded by the uniform re-quantification pipeline used for all reported transcriptome results. |
| `02_view_construction/99_superseded_rare_burden_v1.py` | An early implementation of the rare-variant burden built from an incomplete annotation file that omitted *MEFV*. Detected before any model was fitted and never used for any reported result. |

## What the exclusions do and do not affect

The reproducibility chain reported in the paper is:

```
01_rnaseq_quantification/01 → 03 → 04
02_view_construction/01 → 03
03_mofa_models/01 → 02
04_statistics/01 → 02 → 03
05_figures/01 → 02
```

None of the excluded scripts appears in this chain. Every figure, table and statistic in the
paper can be regenerated from the included scripts together with the processed matrices
provided as Supplementary Data 1.

## Residual disclosure in the included code

The included scripts contain no patient identifiers. They do contain the **cohort sample-naming
convention** (regular expressions matching institutional prefixes such as `sap`, `saga`, `asic`
and `okn`), which is required to harmonize identifiers across omics layers. These are naming
patterns rather than individual identifiers, and no individual code appears. This residual
disclosure should nevertheless be confirmed as acceptable by the institutional review board
before the repository is made public.
