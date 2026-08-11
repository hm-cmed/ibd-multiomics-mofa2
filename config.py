# -*- coding: utf-8 -*-
"""Central path configuration.

All filesystem locations used by the analysis scripts are defined here so that no
absolute path is hard-coded in the analysis code itself.

Set the locations in one of two ways:

1. Copy ``config.example.yaml`` to ``config.yaml`` and edit it (requires PyYAML), or
2. Export the corresponding environment variables (see ENV_KEYS below).

Raw sequencing data are under controlled access and are not distributed with this
repository; see the Data availability statement of the accompanying paper. The
processed matrices needed to rerun everything from ``03_mofa_models`` onward are
distributed as Supplementary Data 1 — point ``project_dir`` at the directory holding
them and the downstream scripts will run without any raw data.
"""
from __future__ import annotations

import os
from pathlib import Path

_HERE = Path(__file__).resolve().parent

ENV_KEYS = {
    "project_dir": "MOFA_IBD_PROJECT_DIR",
    "data_root": "MOFA_IBD_DATA_ROOT",
    "mutation_dir": "MOFA_IBD_MUTATION_DIR",
    "bam_dir_primary": "MOFA_IBD_BAM_DIR",
    "gtf": "MOFA_IBD_GTF",
    "clinical_xlsx": "MOFA_IBD_CLINICAL_XLSX",
    "cohort_xlsx": "MOFA_IBD_COHORT_XLSX",
    "metabolome_batch3_xlsx": "MOFA_IBD_METABOLOME_XLSX",
}

_DEFAULTS = {
    "project_dir": str(_HERE.parent / "analysis_output"),
    "data_root": "",
    "mutation_dir": "",
    "bam_dir_primary": "",
    "gtf": "",
    "clinical_xlsx": "",
    "cohort_xlsx": "",
    "metabolome_batch3_xlsx": "",
}


def _load() -> dict:
    cfg = dict(_DEFAULTS)
    yaml_path = _HERE / "config.yaml"
    if yaml_path.exists():
        try:
            import yaml  # type: ignore
            loaded = yaml.safe_load(yaml_path.read_text(encoding="utf-8")) or {}
            cfg.update({k: v for k, v in loaded.items() if v})
        except ImportError:
            raise SystemExit(
                "config.yaml is present but PyYAML is not installed. "
                "Run 'pip install pyyaml', or set the environment variables listed in config.py."
            )
    for key, env in ENV_KEYS.items():
        if os.environ.get(env):
            cfg[key] = os.environ[env]
    return cfg


_CFG = _load()

PROJECT_DIR = Path(_CFG["project_dir"])
DATA_ROOT = Path(_CFG["data_root"]) if _CFG["data_root"] else Path("")
MUTATION_DIR = Path(_CFG["mutation_dir"]) if _CFG["mutation_dir"] else DATA_ROOT / "mutation_analysis"
BAM_DIR_PRIMARY = Path(_CFG["bam_dir_primary"]) if _CFG["bam_dir_primary"] else DATA_ROOT / "bam"
GTF = Path(_CFG["gtf"]) if _CFG["gtf"] else DATA_ROOT / "GRCh38_NCBI_109.gtf"
CLINICAL_XLSX = Path(_CFG["clinical_xlsx"]) if _CFG["clinical_xlsx"] else DATA_ROOT / "clinical.xlsx"
COHORT_XLSX = Path(_CFG["cohort_xlsx"]) if _CFG["cohort_xlsx"] else DATA_ROOT / "cohort_ids.xlsx"
METABOLOME_BATCH3_XLSX = (
    Path(_CFG["metabolome_batch3_xlsx"]) if _CFG["metabolome_batch3_xlsx"]
    else DATA_ROOT / "metabolome_batch3.xlsx"
)

# Derived locations inside the project directory
MANUSCRIPT_DIR = PROJECT_DIR / "manuscript"
TABLES_DIR = MANUSCRIPT_DIR / "tables"
FIGURES_DIR = MANUSCRIPT_DIR / "figures"
QUANT_DIR = PROJECT_DIR / "_quant"
EXON_BED = QUANT_DIR / "gene_exons.bed"
AA_MUTATION_JSON = MUTATION_DIR / "Client_aa_mut" / "all_genes_all_patients.json"

# Gene annotation table (Gene_ID -> Gene_Symbol, description, biotype, HGNC, Ensembl),
# derived from the public NCBI RefSeq GRCh38 release 109 annotation. It contains no patient
# data and is therefore distributed inside this repository under reference/. A copy placed in
# PROJECT_DIR takes precedence, so users may substitute their own annotation build.
GENE_ANNOTATION = (
    PROJECT_DIR / "gene_annotation.tsv"
    if (PROJECT_DIR / "gene_annotation.tsv").exists()
    else _HERE / "reference" / "gene_annotation.tsv"
)

# Optional: vendor (StringTie) TPM matrix for the 16 samples that were quantified twice,
# used only for the concordance check in 01_rnaseq_quantification/04. If absent, that
# check is skipped and the pipeline proceeds.
VENDOR_TPM = PROJECT_DIR / "vendor_stringtie_tpm.tsv"

for _d in (PROJECT_DIR, MANUSCRIPT_DIR, TABLES_DIR, FIGURES_DIR, QUANT_DIR):
    _d.mkdir(parents=True, exist_ok=True)


def require(*names: str) -> None:
    """Raise a clear error if a script needs a raw-data path that has not been configured."""
    missing = [n for n in names if not str(globals().get(n, "")).strip()
               or not Path(str(globals().get(n))).exists()]
    if missing:
        raise SystemExit(
            "The following paths are not configured or do not exist: "
            + ", ".join(missing)
            + "\nEdit config.yaml (see config.example.yaml) or set the corresponding "
              "environment variables. Raw sequencing data are under controlled access; "
              "see the Data availability statement of the paper."
        )
