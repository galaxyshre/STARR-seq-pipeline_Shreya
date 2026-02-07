# 05 — TSS identification (template-switch libraries; Sahu et al. framework)
## Reference implementation (Sahu et al.)

The TSS identification pipeline follows the framework described in
Sahu et al. (2022). The original reference code is available at:

https://zenodo.org/records/5159644

The code is not re-hosted in this repository; users should download it
separately and provide the local path when running the pipeline.

## Purpose
Identify transcription start sites (TSSs) captured by template-switch (TS) libraries and associate TSS positions with inserts.

## Inputs
- TS library FASTQs (R1)
- Reference genome (hg38) or input file for genomic library , input R1 file for  Random library 
- Sahu/TS pipeline code location

## Outputs
- TSS calls (BED / TSV)
- Summary stats: number of TSS reads, adapter-containing fraction, unique TSS counts
- Insert-associated TSS annotations (where applicable)

## Pipeline location
TSS pipeline code:
- Sahu's code present in Zenodo

## Notes / TODO
- Specify read structure assumptions (TS adapter, primers, UMI if any)
- Define final output formats expected by downstream modelling

