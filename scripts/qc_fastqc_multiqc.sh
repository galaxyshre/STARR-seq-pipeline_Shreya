#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./qc_fastqc_multiqc.sh /path/to/fastqs out_qc 8
FASTQ_DIR="${1:-.}"
OUTDIR="${2:-qc}"
THREADS="${3:-8}"

mkdir -p "${OUTDIR}/fastqc"

# Run FastQC on all fastq.gz in FASTQ_DIR
fastqc "${FASTQ_DIR}"/*.fastq.gz -t "${THREADS}" -o "${OUTDIR}/fastqc"

# Aggregate with MultiQC
multiqc "${OUTDIR}/fastqc" -o "${OUTDIR}"
