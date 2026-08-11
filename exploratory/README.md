# Exploratory model-selection code (not reported in the paper)

The scripts in this directory implement **alternative encodings of the mutation view** that were
evaluated during model selection but that produced **no result reported in the paper**. They are
included for transparency about how the final mutation view was chosen, and because reviewers
and readers may reasonably ask why a rare-missense allele-dosage burden was preferred over
simpler alternatives.

Nothing in the main analysis depends on this directory. Deleting it does not affect
reproduction of any figure, table or statistic in the paper.

## What was tried, and why it was not used

| Encoding | Scripts | Outcome |
|---|---|---|
| **Binary gene-level carrier status** (variant present / absent per gene, Bernoulli likelihood) | `02`, `03` | The mutation view explained essentially no variance in any factor. Too sparse and too coarse at this sample size. |
| **Pathway-level burden over all missense variants** (mean variants per gene, aggregated into nine innate-immunity modules) | `01`, `04`, `05`, `06` | The mutation view did load onto a factor, but that factor tracked *total* variant burden rather than any specific pathway, was unrelated to diagnosis, and shared almost no variance with the other omics layers. Dominated by common polymorphisms. |
| **Pathway-level burden restricted to rare missense variants** | `07`, `08`, `09` | Improved on the previous encoding and motivated the final approach, but aggregating to nine modules discarded gene-level resolution. |

These observations led to the encoding used in the paper: a **gene-level rare-missense
allele-dosage burden**, restricted to variants with a cohort carrier frequency of 10% or less,
implemented in `../02_view_construction/03_mutation_rare_missense_burden.py`. Because that
frequency filter necessarily excludes common *MEFV* variants of interest, an additional
sensitivity model with explicit *MEFV* dosage features was fitted; that model **is** reported in
the paper and lives in the main tree
(`../03_mofa_models/03_MEFV_sensitivity_model.py`).

## Running these scripts

They use the same `config.py` as the main tree and write to separate output directories
(`results/`, `results_pathway/`, `results_rareburden/` inside `project_dir`) so that they cannot
overwrite the reported results.
