# 04 — Activity quantification (signal vs background)

## Purpose
Quantify regulatory activity by comparing RNA (signal) abundance to DNA input (background) abundance per insert.

## Inputs
- Insert identity (Step 03): genomic coordinates or random clusters
- Counts per insert from:
  - STARR RNA libraries (signal)
  - Input DNA libraries (background)
- Replicate metadata (cell line, batch, replicate)

## Outputs
- Activity table per insert:
  - insert_id
  - RNA_count (or UMI count)
  - DNA_count (or UMI count)
  - normalized RNA and DNA
  - activity score (log2 fold change vs input)
- Replicate summaries and correlation metrics

## Method (framework)
Activity is computed as enrichment of RNA over DNA input after library-size normalization.
(TBD: finalize exact pseudocount and filtering thresholds.)

## QC checks / expected results
- Activity distributions sensible (centered ~0 with active tail)
- Replicate agreement (correlation / IDR-like reproducibility where relevant)
- Low-count inserts filtered to avoid unstable ratios

## Notes / TODO decisions
- Decide filtering thresholds: minimum DNA count, minimum RNA count
- Decide whether to use UMI-based counts vs read counts
- Define how replicates are merged (mean, median, model-based)

