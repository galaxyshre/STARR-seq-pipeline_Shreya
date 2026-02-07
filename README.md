# STARR-seq analysis pipeline

This repository describes a step-by-step pipeline for STARR-seq data analysis,
developed in the Taipale group and intended for automation by HGI.

## Scope
- STARR-seq activity quantification (RNA vs input)
- Genomic and random STARR-seq libraries
- Transcription start site (TSS) identification from template-switch libraries
- Motif discovery
- Downstream machine learning modelling

## Pipeline overview
1. Data intake and QC (FastQC / MultiQC)
2. Read processing (adapter trimming, UMI handling, deduplication)
3. Insert identity
   - Genomic libraries: alignment to hg38
   - Random libraries: Starcode clustering
4. Activity quantification (signal vs input)
5. TSS identification (Sahu et al. pipeline)
6. Motif discovery
7. Machine learning modelling (HGI)

## Repository structure
- `docs/`     : step-by-step documentation for each pipeline stage
- `scripts/`  : runnable scripts (bash)
- `configs/`  : configuration files and sample sheets

## External dependencies
- **Sahu et al. code** (Zenodo): https://zenodo.org/records/5159644  
  (Downloaded separately; not re-hosted in this repository.  
  The local path should be specified via `SAHU_CODE_DIR` in `configs/config.yaml`.)

## Status
This pipeline is under active development.

