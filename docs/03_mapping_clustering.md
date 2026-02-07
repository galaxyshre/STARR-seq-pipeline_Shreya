# 03 — Insert identity (hg38 mapping for genomic libraries; Starcode clustering for random libraries)

## Purpose
Define insert identity for downstream quantification:
- Genomic libraries: map inserts to hg38 and define genomic coordinates
- Random libraries: cluster sequences to define unique inserts

## Inputs
- Processed FASTQs from Step 02 (or raw FASTQs if processing is not applied yet)
- Reference genome: hg38 FASTA + bowtie2 index


## Outputs
### Genomic libraries
- BAM files (sorted + indexed)
- Mapping stats: % mapped, properly paired, MAPQ distribution, duplicate rate
### Random libraries
- Cluster table: cluster_id → representative sequence → counts
- Summary: number of clusters, cluster size distribution

## Tools
- bowtie2 + samtools (+ picard if dedup needed)
- starcode (random library clustering)
- seqkit (optional pre-checks)
- Preseq tools for the library complexity analysis
## Step-by-step commands (representative)
### A) Genomic mapping (hg38; paired-end)
(TBD: paste final bowtie2 command and filtering thresholds used)

### B) Random library clustering (Starcode)
(TBD: paste starcode parameters; e.g. distance threshold, input format)

## QC checks / expected results
- Mapping rate and unique mapping within expected range for hg38
- MAPQ distribution consistent with correct mapping
- Random library: cluster count consistent with library complexity expectations
- No unexpected inflation of very large clusters (possible PCR artefacts)

## Notes / TODO decisions
- Confirm whether mapping is to hg38 primary assembly and which bowtie2 settings are final.
- Confirm Starcode distance threshold for clustering (Sahu et al. style).

