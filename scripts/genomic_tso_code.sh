#!/bin/bash
#BSUB -n 16
#BSUB -M 200GB
#BSUB -o STARR-TSS.%J.out.log
#BSUB -e STARR-TSS.%J.err.log
#BSUB -J STARR-TSS
#BSUB -q long
#BSUB -G team373
#BSUB -R "select[mem>200GB] rusage[mem=200GB] span[hosts=1]"

# =============================================================================
# STARR-seq TSS Calling Pipeline (config-driven)
# =============================================================================
# Identifies transcription start sites (TSSs) from a TSO/template-switch
# library, restricted to input-DNA fragments whose 20bp barcode also appears
# in the reporter library.
#
# USAGE
# -----
#   Direct:        ./starr_tss.sh --conf my_sample.conf
#   LSF (command): bsub -n 16 -M 200GB -q long -G team373 \
#                       -R "select[mem>200GB] rusage[mem=200GB] span[hosts=1]" \
#                       -J tss_<sample> -o tss_<sample>.%J.out.log \
#                       -e tss_<sample>.%J.err.log \
#                       "bash starr_tss.sh --conf my_sample.conf"
#   LSF (stdin):   bsub -env "all, CONF=my_sample.conf" < starr_tss.sh
#
usage() {
  cat <<EOF
Usage: $(basename "$0") --conf <config.conf>
   or: CONF=<config.conf> bsub < $(basename "$0")

Runs the STARR-seq TSS calling pipeline using settings from the given
config file. Edit the .conf file per sample; never edit this script.
EOF
}

CONF="${CONF:-}"   # env-var fallback, useful with:  bsub < starr_tss.sh
while [[ $# -gt 0 ]]; do
  case "$1" in
    --conf=*) CONF="${1#*=}"; shift ;;
    --conf)
      if [[ $# -lt 2 ]]; then
        echo "ERROR: --conf requires a file argument" >&2; usage; exit 1
      fi
      CONF="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "${CONF}" ]]; then
  echo "ERROR: no config given. Use --conf <file> or set \$CONF." >&2
  usage; exit 1
fi
if [[ ! -f "${CONF}" ]]; then
  echo "ERROR: config file not found: ${CONF}" >&2; exit 1
fi

set -euo pipefail

# -------------------------------
THREADS=16
UMI_LEN=5
UPSTREAM=60
DOWNSTREAM=60

REF_DIR="/nfs/team373/dz10/ref/hs1"

INPUT_DNA_DIR_RAW="/lustre/scratch125/gengen/projects_v2/starr-seq/GP_input_aviti"

CONDA_ENV="/software/conda/users/dz10/encode-atac-seq"

source "${CONF}"

BT2_INDEX="${BT2_INDEX:-${REF_DIR}/bt2/hs1/hs1}"
GENOME="${GENOME:-${REF_DIR}/hs1.fa}"
INPUT_DNA_R1_RAW="${INPUT_DNA_R1_RAW:-${INPUT_DNA_DIR_RAW}/Adept_SS_L1_R1.fastq.gz}"
INPUT_DNA_R2_RAW="${INPUT_DNA_R2_RAW:-${INPUT_DNA_DIR_RAW}/Adept_SS_L1_R2.fastq.gz}"


missing_var=0
for v in SAMPLE REPORTER_R1 REPORTER_R2 TSO_R1_RAW TSO_R2_RAW OUT_DIR; do
  if [[ -z "${!v:-}" ]]; then
    echo "ERROR: required config variable '${v}' is not set in ${CONF}" >&2
    missing_var=1
  fi
done
[[ "${missing_var}" -eq 0 ]] || exit 1

missing_file=0
for f in "${REPORTER_R1}" "${REPORTER_R2}" \
         "${TSO_R1_RAW}" "${TSO_R2_RAW}" \
         "${INPUT_DNA_R1_RAW}" "${INPUT_DNA_R2_RAW}" \
         "${GENOME}"; do
  if [[ ! -r "${f}" ]]; then
    echo "ERROR: input file not found / not readable: ${f}" >&2
    missing_file=1
  fi
done
[[ "${missing_file}" -eq 0 ]] || exit 1

if [[ ! -e "${BT2_INDEX}.1.bt2" && ! -e "${BT2_INDEX}.1.bt2l" ]]; then
  echo "ERROR: bowtie2 index not found at prefix '${BT2_INDEX}'" >&2
  echo "       (expected ${BT2_INDEX}.1.bt2 or ${BT2_INDEX}.1.bt2l)" >&2
  exit 1
fi

# Warn if the threads requested exceed the slots LSF actually allocated
if [[ -n "${LSB_DJOB_NUMPROC:-}" && "${THREADS}" -gt "${LSB_DJOB_NUMPROC}" ]]; then
  echo "WARNING: THREADS=${THREADS} but LSF allocated ${LSB_DJOB_NUMPROC} slots" \
       "(#BSUB -n). Keep them aligned to avoid oversubscription." >&2
fi

module load ISG/conda
conda activate "${CONDA_ENV}"


WORK_DIR="${OUT_DIR}/input_dna"
BC_DIR="${OUT_DIR}/barcodes"
mkdir -p "${WORK_DIR}" "${BC_DIR}"

UMI_PATTERN=$(printf 'N%.0s' $(seq 1 ${UMI_LEN}))


strip_fq_ext() { local n; n=$(basename "$1"); n="${n%.gz}"; n="${n%.fastq}"; n="${n%.fq}"; printf '%s' "${n}"; }
INPUT_DNA_R1_NAME=$(strip_fq_ext "${INPUT_DNA_R1_RAW}")
INPUT_DNA_R2_NAME=$(strip_fq_ext "${INPUT_DNA_R2_RAW}")
TSO_R1_NAME=$(strip_fq_ext "${TSO_R1_RAW}")
TSO_R2_NAME=$(strip_fq_ext "${TSO_R2_RAW}")

echo "=== Pipeline config ==="
echo "  CONF        = ${CONF}"
echo "  SAMPLE      = ${SAMPLE}"
echo "  UMI_LEN     = ${UMI_LEN}  (pattern: ${UMI_PATTERN})"
echo "  THREADS     = ${THREADS}"
echo "  Flanks      = ${UPSTREAM}bp upstream + ${DOWNSTREAM}bp downstream"
echo "  REF_DIR     = ${REF_DIR}"
echo "  REPORTER_R1 = ${REPORTER_R1}"
echo "  REPORTER_R2 = ${REPORTER_R2}"
echo "  TSO_R1_RAW  = ${TSO_R1_RAW}"
echo "  TSO_R2_RAW  = ${TSO_R2_RAW}"
echo "  OUT_DIR     = ${OUT_DIR}"
echo "  WORK_DIR    = ${WORK_DIR}"
echo


cutadapt \
    -j ${THREADS} \
    -g "^CCTCCTGTGAGTTTGGTTGGTGTACAGTAGCTTCCACCA...AGATGGGAAGAGGACACCTCTGAACTCC" \
    -G "^NGAGTTCAGAGGTGTCCTCTTCCCATCT...TGGTGGAAGCTACTGTACACCAA" \
    -e 0.15 --pair-filter=any --discard-untrimmed \
    --minimum-length 18 --maximum-length 22 \
    -o "${BC_DIR}/R1_barcodes.fastq.gz" \
    -p "${BC_DIR}/R2_barcodes.fastq.gz" \
    "${REPORTER_R1}" "${REPORTER_R2}"

seqkit fx2tab "${BC_DIR}/R1_barcodes.fastq.gz" \
  | awk -F'\t' '{split($1, a, " "); print a[1] "\t" $2}' \
  > "${WORK_DIR}/bc_from_R1.txt"

seqkit fx2tab "${BC_DIR}/R2_barcodes.fastq.gz" \
  | awk -F'\t' '{split($1, a, " "); print a[1] "\t" $2}' \
  > "${WORK_DIR}/bc_from_R2.txt"

# ---- A2. Pre-filter input DNA: drop reads with vector-only ends -------------
INPUT_DNA_R1_FILT="${INPUT_DNA_DIR_RAW}/${INPUT_DNA_R1_NAME}.filtered.fastq.gz"
INPUT_DNA_R2_FILT="${INPUT_DNA_DIR_RAW}/${INPUT_DNA_R2_NAME}.filtered.fastq.gz"

cutadapt -j ${THREADS} -g "^ACCGGT" \
    --discard-trimmed --no-indels -e 0 \
    -o "${INPUT_DNA_R1_FILT}" "${INPUT_DNA_R1_RAW}"

cutadapt -j ${THREADS} -g "^GGAGTTCAGAGGTGTCCTCTTCCCATCTGTCGAC" \
    --discard-trimmed --no-indels -e 0 \
    -o "${INPUT_DNA_R2_FILT}" "${INPUT_DNA_R2_RAW}"

seqkit pair -j ${THREADS} -1 "${INPUT_DNA_R1_FILT}" -2 "${INPUT_DNA_R2_FILT}"

INPUT_DNA_R1_PAIRED="${INPUT_DNA_DIR_RAW}/${INPUT_DNA_R1_NAME}.filtered.paired.fastq.gz"
INPUT_DNA_R2_PAIRED="${INPUT_DNA_DIR_RAW}/${INPUT_DNA_R2_NAME}.filtered.paired.fastq.gz"

INPUT_DNA_R1_UMI="${INPUT_DNA_DIR_RAW}/${INPUT_DNA_R1_NAME}.umi.fastq.gz"
INPUT_DNA_R2_UMI="${INPUT_DNA_DIR_RAW}/${INPUT_DNA_R2_NAME}.umi.fastq.gz"

umi_tools extract --bc-pattern="${UMI_PATTERN}" \
    -I "${INPUT_DNA_R1_PAIRED}" -S "${INPUT_DNA_R1_UMI}" \
    --read2-in="${INPUT_DNA_R2_PAIRED}" --read2-out="${INPUT_DNA_R2_UMI}"

cutadapt -j ${THREADS} \
    -g "^GGAGTTCAGAGGTGTCCTCTTCCCATCT...TGGTGGAAGCTACTGTACACCAA" \
    -e 0.15 --no-indels \
    --minimum-length 18 --maximum-length 22 \
    --discard-untrimmed \
    -o "${WORK_DIR}/R2_barcodes.fastq.gz" \
    "${INPUT_DNA_R2_UMI}" \
    > "${WORK_DIR}/cutadapt_R2_bc.log" 2>&1

seqkit fx2tab "${WORK_DIR}/R2_barcodes.fastq.gz" \
  | awk -F'\t' '{split($1, a, " "); print a[1] "\t" $2}' \
  > "${WORK_DIR}/barcode_table.txt"

echo "[A5] Input-DNA barcode entries:    $(wc -l < ${WORK_DIR}/barcode_table.txt)"
echo "[A5] Input-DNA unique barcodes:    $(cut -f2 ${WORK_DIR}/barcode_table.txt | sort -u | wc -l)"
echo "[A5] Reporter unique barcodes R2:  $(cut -f2 ${WORK_DIR}/bc_from_R2.txt | sort -u | wc -l)"

awk 'NR==FNR { a[$2]; next } $2 in a' \
    "${WORK_DIR}/bc_from_R2.txt" \
    "${WORK_DIR}/barcode_table.txt" \
  > "${WORK_DIR}/filtered_barcode_table.txt"

awk '{print $1}' "${WORK_DIR}/filtered_barcode_table.txt" \
  > "${WORK_DIR}/keep_names.txt"

echo "[A6] Reads matched between input DNA and reporter: $(wc -l < ${WORK_DIR}/keep_names.txt)"

MATCHED_R1="${WORK_DIR}/matched_R1.fastq.gz"
MATCHED_R2="${WORK_DIR}/matched_R2.fastq.gz"

seqkit grep -j ${THREADS} -f "${WORK_DIR}/keep_names.txt" \
    "${INPUT_DNA_R1_UMI}" -o "${MATCHED_R1}"

seqkit grep -j ${THREADS} -f "${WORK_DIR}/keep_names.txt" \
    "${INPUT_DNA_R2_UMI}" -o "${MATCHED_R2}"

seqkit stats "${MATCHED_R1}"

INPUT_BAM="${WORK_DIR}/matched_R1.bam"

bowtie2 -x "${BT2_INDEX}" \
    -U "${MATCHED_R1}" \
    --very-sensitive --no-unal \
    --threads ${THREADS} \
  | samtools view -bS -q 20 -F 2304 - \
  | samtools sort -@ ${THREADS} -o "${INPUT_BAM}" -

samtools index "${INPUT_BAM}"

echo "[A8] Input fragments aligned MAPQ >= 20, primary only: $(samtools view -c ${INPUT_BAM})"

TSO_R1_FILT="${WORK_DIR}/${TSO_R1_NAME}.filtered.fastq.gz"

cutadapt -j ${THREADS} -g "^ACCGGT" \
    --discard-trimmed --no-indels -e 0 \
    -o "${TSO_R1_FILT}" "${TSO_R1_RAW}"

seqkit pair -j ${THREADS} -1 "${TSO_R1_FILT}" -2 "${TSO_R2_RAW}" -o "${WORK_DIR}/"

TSO_R1_PAIRED="${WORK_DIR}/${TSO_R1_NAME}.filtered.paired.fastq.gz"
TSO_R2_PAIRED="${WORK_DIR}/${TSO_R2_NAME}.paired.fastq.gz"

TSO_R2_FILT="${WORK_DIR}/${TSO_R2_NAME}.paired.filtered.fastq.gz"

cutadapt -j ${THREADS} -g "^NTCGAC" \
    --discard-trimmed --no-indels -e 0 \
    -o "${TSO_R2_FILT}" "${TSO_R2_PAIRED}"

seqkit pair -j ${THREADS} -1 "${TSO_R2_FILT}" -2 "${TSO_R1_PAIRED}" -o "${WORK_DIR}/"

TSO_R1_FINAL="${WORK_DIR}/${TSO_R1_NAME}.filtered.paired.paired.fastq.gz"

PARSED_FQ="${WORK_DIR}/${SAMPLE}.parsed.fq.gz"

zcat "${TSO_R1_FINAL}" \
  | awk -v U="${UMI_LEN}" 'BEGIN{OFS="\n"}
    {
      m=NR%4;
      if (m==1) { sub(/^@/,""); split($0,a," "); rid=a[1]; }
      else if (m==2) { seq=$0; }
      else if (m==3) { plus=$0; }
      else if (m==0) {
        qual=$0;
        umi    = substr(seq,  1, U);
        rest   = substr(seq,  U+1);
        rest_q = substr(qual, U+1);
        if      (rest ~ /^GGGGG/) numg=5;
        else if (rest ~ /^GGGG/)  numg=4;
        else if (rest ~ /^GGG/)   numg=3;
        else                      next;
        out_seq  = substr(rest,   4);
        out_qual = substr(rest_q, 4);
        if (length(out_seq) < 20) next;
        printf "@%s_UMI:%s_NG:%d\n%s\n+\n%s\n", rid, umi, numg, out_seq, out_qual;
      }
    }' \
  | gzip > "${PARSED_FQ}"

echo "[B2] Parsed TSO reads: $(zcat ${PARSED_FQ} | awk 'NR%4==1' | wc -l)"

echo "[B2] Leading-G distribution after UMI strip:"
zcat "${TSO_R1_FINAL}" | awk 'NR%4==2' | awk -v U="${UMI_LEN}" '{
    rest=substr($0, U+1);
    if      (rest ~ /^GGGGG/) g="5G";
    else if (rest ~ /^GGGG/)  g="4G";
    else if (rest ~ /^GGG/)   g="3G";
    else                      g="noG";
    counts[g]++; total++;
  }
  END {
    for (k in counts) printf "       %-4s %10d  %.2f%%\n", k, counts[k], 100*counts[k]/total;
    printf "       %-4s %10d\n", "TOTAL", total;
  }'

TSO_BAM="${WORK_DIR}/${SAMPLE}.aligned.bam"

bowtie2 -x "${BT2_INDEX}" \
    -U "${PARSED_FQ}" \
    --local --very-sensitive-local --no-unal \
    -p ${THREADS} \
  | samtools view -b -F 2304 - \
  | samtools sort -@ ${THREADS} -o "${TSO_BAM}" -

samtools index "${TSO_BAM}"
samtools flagstat "${TSO_BAM}"

echo "[B4] MAPQ distribution top 10:"
samtools view -F 4 "${TSO_BAM}" | awk '{print $5}' \
  | sort | uniq -c | sort -rn | head

echo "[B4] numG vs leading soft-clip:"
samtools view -F 4 "${TSO_BAM}" | awk '{
    if (match($1, /NG:[0-9]+/)) ng = substr($1, RSTART+3, RLENGTH-3) + 0;
    flag    = $2 + 0;
    is_rev  = (and(flag, 16) > 0);
    cigar   = $6;
    leading = 0;
    if (!is_rev && match(cigar, /^[0-9]+S/)) leading = substr(cigar, 1, RLENGTH-1) + 0;
    if ( is_rev && match(cigar, /[0-9]+S$/)) leading = substr(cigar, RSTART, RLENGTH-1) + 0;
    counts[ng "\t" leading]++;
  }
  END {
    print "numg\tleading_softclip\tcount";
    for (k in counts) print k "\t" counts[k];
  }' | sort -k1,1n -k2,2n

TSS_BED="${WORK_DIR}/${SAMPLE}.tss.bed"

samtools view -F 4 -q 20 "${TSO_BAM}" \
  | awk 'BEGIN{OFS="\t"}
    {
      flag   = $2 + 0;
      is_rev = (and(flag, 16) > 0);
      chrom  = $3;
      pos    = $4 + 0;
      cigar  = $6;
      name   = $1;

      umi="."; ng=0;
      if (match(name, /UMI:[ACGTN]+/)) umi = substr(name, RSTART+4, RLENGTH-4);
      if (match(name, /NG:[0-9]+/))    ng  = substr(name, RSTART+3, RLENGTH-3) + 0;

      leading = 0;
      if (!is_rev && match(cigar, /^[0-9]+S/)) leading = substr(cigar, 1, RLENGTH-1) + 0;
      if ( is_rev && match(cigar, /[0-9]+S$/)) leading = substr(cigar, RSTART, RLENGTH-1) + 0;

      if (leading > ng - 3) next;

      n=0; ref_len=0;
      for (i=1; i<=length(cigar); i++) {
        c = substr(cigar, i, 1);
        if (c ~ /[0-9]/) { n = n*10 + (c+0); }
        else {
          if (c=="M" || c=="D" || c=="N" || c=="=" || c=="X") ref_len += n;
          n = 0;
        }
      }

      if (is_rev) { tss = pos - 1 + ref_len - 1; strand = "-"; }
      else        { tss = pos - 1;               strand = "+"; }

      print chrom, tss, tss+1, name, $5, strand, umi, ng, leading;
    }' \
  | sort -k1,1 -k2,2n > "${TSS_BED}"

echo "[B5] TSS BED records: $(wc -l < ${TSS_BED})"
awk '{strands[$6]++} END {for (s in strands) printf "       strand %s: %d\n", s, strands[s]}' "${TSS_BED}"


FRAG_BED="${WORK_DIR}/input_fragments.bed"

samtools view -F 4 -q 20 -b "${INPUT_BAM}" \
  | bedtools bamtobed -i - \
  | sort -k1,1 -k2,2n > "${FRAG_BED}"

echo "[C1] Input fragments: $(wc -l < ${FRAG_BED})"

TSS_IN_FRAG="${WORK_DIR}/${SAMPLE}.tss_in_fragments.tsv"

bedtools intersect -a "${TSS_BED}" -b "${FRAG_BED}" -wa -wb -sorted \
  > "${TSS_IN_FRAG}"

echo "[C2] TSS x fragment overlaps: $(wc -l < ${TSS_IN_FRAG})"
echo "[C2] Distinct TSS reads with fragment hit: $(awk '{print $4}' ${TSS_IN_FRAG} | sort -u | wc -l)"

FRAG_TSS_TSV="${WORK_DIR}/${SAMPLE}.fragment_tss.tsv"

awk 'BEGIN{OFS="\t"; FS="\t"}
{
  tss_chrom=$1; tss_pos=$2+0; tss_strand=$6; umi=$7;
  frag_start=$11+0; frag_end=$12+0; frag_strand=$15;
  frag_id = $1 ":" $11 "-" $12;

  if (tss_strand == "+") off = tss_pos - frag_start;
  else                   off = (frag_end - 1) - tss_pos;

  key = frag_id "\t" frag_start "\t" frag_end "\t" frag_strand "\t" tss_pos "\t" tss_strand;
  reads[key]++;
  if (!((key, umi) in seen)) { seen[key,umi]=1; umis[key]++; }
  off_of[key]   = off;
  chrom_of[key] = tss_chrom;
}
END {
  print "frag_id\tfrag_start\tfrag_end\tfrag_strand\ttss_chrom\ttss_pos\ttss_strand\toffset_from_frag5p\tn_unique_umi\tn_reads";
  for (k in reads) {
    split(k, a, "\t");
    print a[1], a[2], a[3], a[4], chrom_of[k], a[5], a[6], off_of[k], umis[k], reads[k];
  }
}' "${TSS_IN_FRAG}" \
  | (read header; printf "%s\n" "$header"; sort -k1,1 -k6,6n) \
  > "${FRAG_TSS_TSV}"

echo "[C3] Fragment-TSS records: $(tail -n +2 ${FRAG_TSS_TSV} | wc -l)"

FLANK_BED="${WORK_DIR}/${SAMPLE}.tss_flanks.bed"

tail -n +2 "${FRAG_TSS_TSV}" \
  | awk -v U=${UPSTREAM} -v D=${DOWNSTREAM} 'BEGIN{OFS="\t"}
    {
      tss_chrom=$5; tss_pos=$6+0; tss_strand=$7;
      frag_id=$1; n_umi=$9; n_reads=$10;
      if (tss_strand == "+") { start = tss_pos - U; end = tss_pos + 1 + D; }
      else                   { start = tss_pos - D; end = tss_pos + 1 + U; }
      if (start < 0) next;
      name = frag_id "_TSS:" tss_chrom ":" tss_pos ":" tss_strand "_nUMI:" n_umi "_nReads:" n_reads;
      print tss_chrom, start, end, name, ".", tss_strand;
    }' \
  | sort -k1,1 -k2,2n -u > "${FLANK_BED}"

FLANK_FA="${WORK_DIR}/${SAMPLE}.tss_flanks.fasta"

bedtools getfasta -fi "${GENOME}" -bed "${FLANK_BED}" -name -s > "${FLANK_FA}"

awk '/^>/{print; next} {print toupper($0)}' "${FLANK_FA}" > "${FLANK_FA}.tmp" \
  && mv "${FLANK_FA}.tmp" "${FLANK_FA}"
