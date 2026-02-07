# 02 — Read processing (adapter trimming, UMI handling, deduplication)

## Purpose
Prepare raw FASTQs for insert identification and quantification by removing adapters/primers, handling UMIs (if present), and removing technical duplicates where appropriate.

## Inputs
- Raw FASTQ files (R1/R2) from Step 01
- (Optional) UMI definition: read (R1 or R2) and positions (e.g., first 5 bp)

## Outputs
- Trimmed FASTQs (R1/R2)
- UMI-extracted read headers or sidecar UMI files (if used)
- Deduplicated read sets (method depends on UMI availability)
- Summary stats (reads kept, trimmed %, duplicate %)

## Tools
- cutadapt or fastp (adapter/primer trimming)
- (Optional) umi_tools (UMI extraction/dedup)
- seqkit (sanity checks)
- samtools (after mapping)

## Step-by-step commands (representative)
### A) Adapter/primer trimming (example with cutadapt)
(TBD: fill exact adapter/primer sequences used in our STARR-seq libraries)

### B) UMI extraction (if applicable)
(TBD: specify UMI length and location; e.g. R1 positions 1–5)

### C) Deduplication approach
- If UMIs are present: prefer UMI-aware dedup (umi_tools) after mapping
- If no UMIs: coordinate-based dedup (Picard) after mapping (genomic libraries)

## QC checks / expected results
- Adapter content reduced in FastQC
- Read length distribution matches expected insert/read design
- Reasonable retention after trimming
- Duplication interpreted relative to depth and library complexity

## Notes / TODO decisions
- Finalize adapter/primer sequences for trimming.
- Confirm UMI length and position for each library type (STARR RNA, input DNA, TSS libraries).
- Decide whether dedup is done pre- or post-mapping for each library.

