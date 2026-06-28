#!/usr/bin/env bash
# =============================================================================
# All user-tunable parameters live in K562_pipeline.conf.
# DO NOT EDIT this script for normal runs.
#
# Usage:
#   ./K562_pipeline.sh                                    # uses ./K562_pipeline.conf

# =============================================================================

set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONF_FILE="${SCRIPT_DIR}/K562_pipeline.conf"

# -----------------------------------------------------------------------------
# CLI parsing
# -----------------------------------------------------------------------------
usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --conf FILE         Configuration file (default: ./K562_pipeline.conf)
  --resume            Skip steps whose outputs already exist non-empty.
  --threads N         Override THREADS from config.
  --work-dir DIR      Override WORK_DIR from config.
  --treat-name NAME   Override TREAT_NAME from config.
  --ctrl-name NAME    Override CTRL_NAME from config.
  -h, --help          Show this message.
EOF
}

CLI_RESUME=""
CLI_THREADS=""
CLI_WORK_DIR=""
CLI_TREAT_NAME=""
CLI_CTRL_NAME=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --conf)       CONF_FILE="$2";        shift 2 ;;
        --resume)     CLI_RESUME=1;          shift   ;;
        --threads)    CLI_THREADS="$2";      shift 2 ;;
        --work-dir)   CLI_WORK_DIR="$2";     shift 2 ;;
        --treat-name) CLI_TREAT_NAME="$2";   shift 2 ;;
        --ctrl-name)  CLI_CTRL_NAME="$2";    shift 2 ;;
        -h|--help)    usage; exit 0                  ;;
        *)            echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

# -----------------------------------------------------------------------------
# Source config
# -----------------------------------------------------------------------------
if [[ ! -f "${CONF_FILE}" ]]; then
    echo "ERROR: config file not found: ${CONF_FILE}" >&2
    echo "       Pass --conf /path/to/file or create ${CONF_FILE}" >&2
    exit 1
fi
source "${CONF_FILE}"

# CLI overrides win
[[ -n "${CLI_RESUME}"     ]] && RESUME="${CLI_RESUME}"
[[ -n "${CLI_THREADS}"    ]] && THREADS="${CLI_THREADS}"
[[ -n "${CLI_WORK_DIR}"   ]] && WORK_DIR="${CLI_WORK_DIR}"
[[ -n "${CLI_TREAT_NAME}" ]] && TREAT_NAME="${CLI_TREAT_NAME}"
[[ -n "${CLI_CTRL_NAME}"  ]] && CTRL_NAME="${CLI_CTRL_NAME}"

THREADS="${THREADS:-16}"

# Defaults for anything the conf may omit
RESUME="${RESUME:-0}"
TRIM_LEN="${TRIM_LEN:-37}"
BOWTIE2_MAXINS="${BOWTIE2_MAXINS:-1000}"
DEDUP_MIN_MAPQ="${DEDUP_MIN_MAPQ:-20}"
BIGWIG_BIN_SIZE="${BIGWIG_BIN_SIZE:-25}"
MACS_GENOME="${MACS_GENOME:-hs}"
PSR_SEED="${PSR_SEED:-42}"
JAVA_HEAP_PICARD="${JAVA_HEAP_PICARD:-20g}"
JAVA_HEAP_INSERTSIZE="${JAVA_HEAP_INSERTSIZE:-8g}"
IDR_SCORE_10="${IDR_SCORE_10:-540}"
IDR_SCORE_05="${IDR_SCORE_05:-625}"

SORT_THREADS=$(( THREADS - 1 > 0 ? THREADS - 1 : 1 ))

# Shortcuts
TREAT="${TREAT_NAME}"
CTRL="${CTRL_NAME}"
SAMPLES=("${TREAT}" "${CTRL}")

declare -A SRC_R1=( ["${TREAT}"]="${TREAT_R1}" ["${CTRL}"]="${CTRL_R1}" )
declare -A SRC_R2=( ["${TREAT}"]="${TREAT_R2}" ["${CTRL}"]="${CTRL_R2}" )
declare -A PRE_PAIRED=(
    ["${TREAT}"]="${TREAT_PRE_PAIRED:-0}"
    ["${CTRL}"]="${CTRL_PRE_PAIRED:-0}"
)
declare -A REUSE_BAM=(
    ["${TREAT}"]="${TREAT_REUSE_BAM:-}"
    ["${CTRL}"]="${CTRL_REUSE_BAM:-}"
)
declare -A SKIP_TRIM=(
    ["${TREAT}"]="${TREAT_SKIP_TRIM:-0}"
    ["${CTRL}"]="${CTRL_SKIP_TRIM:-0}"
)

# -----------------------------------------------------------------------------
# Logging + helpers
# -----------------------------------------------------------------------------
_ts()        { date '+%Y-%m-%d %H:%M:%S'; }
log()        { printf '[%s] %s\n'             "$(_ts)" "$*"; }
log_step()   { printf '\n[%s] ===== %s =====\n' "$(_ts)" "$*"; }
log_warn()   { printf '[%s] WARNING: %s\n'    "$(_ts)" "$*" >&2; }
log_error()  { printf '[%s] ERROR: %s\n'      "$(_ts)" "$*" >&2; }

require_tool() {
    if ! command -v "$1" >/dev/null 2>&1; then
        log_error "required tool not on PATH: $1"
        exit 1
    fi
}
require_file() {
    if [[ ! -e "$1" ]]; then
        log_error "required file not found: $1"
        exit 1
    fi
}

# True if RESUME is on AND every listed output exists non-empty.
have_outputs() {
    [[ "${RESUME}" == "1" ]] || return 1
    local f
    for f in "$@"; do
        [[ -s "$f" ]] || return 1
    done
    return 0
}

# -----------------------------------------------------------------------------
# Environment
# -----------------------------------------------------------------------------
if [[ -n "${LMOD_CONDA:-}" ]]; then
    module load "${LMOD_CONDA}"
fi
if [[ -n "${CONDA_ENV:-}" ]]; then
    # shellcheck disable=SC1091
    conda activate "${CONDA_ENV}"
fi
if [[ -n "${LMOD_IDR:-}" ]]; then
    module load "${LMOD_IDR}" || log_warn "could not load ${LMOD_IDR}; IDR step may fail"
fi

for tool in seqkit cutadapt bowtie2 samtools picard bedtools \
            bamCoverage bamCompare macs3 idr awk; do
    require_tool "$tool"
done

require_file "${TREAT_R1}"
require_file "${TREAT_R2}"
require_file "${CTRL_R1}"
require_file "${CTRL_R2}"
require_file "${BT2_INDEX}.1.bt2"

mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}"
mkdir -p fastq trim37 bam bw peaks logs

# -----------------------------------------------------------------------------
# Stage reused BAMs (if any) and build PROCESS_SAMPLES
# -----------------------------------------------------------------------------
# PROCESS_SAMPLES = samples that need raw->BAM processing (Steps 1-4).
# SAMPLES         = both samples; downstream steps (5+) iterate this so QC
#                   and analysis still happen for reused BAMs.
PROCESS_SAMPLES=()
for s in "${SAMPLES[@]}"; do
    reuse="${REUSE_BAM[$s]}"
    if [[ -z "${reuse}" ]]; then
        PROCESS_SAMPLES+=("${s}")
        continue
    fi

    require_file "${reuse}"
    target_bam="bam/${s}.final.bam"
    target_bai="bam/${s}.final.bam.bai"

    # Re-link every run so the symlink target always reflects the conf.
    if [[ -L "${target_bam}" || -e "${target_bam}" ]]; then
        rm -f "${target_bam}" "${target_bai}"
    fi
    ln -s "${reuse}" "${target_bam}"

    if [[ -e "${reuse}.bai" ]]; then
        ln -s "${reuse}.bai" "${target_bai}"
    elif [[ -e "${reuse%.bam}.bai" ]]; then
        ln -s "${reuse%.bam}.bai" "${target_bai}"
    else
        log "Indexing reused BAM for ${s} (no .bai found alongside source)"
        samtools index "${target_bam}"
    fi
    log "Reusing pre-processed BAM for ${s}: ${reuse}"
done

log "=== Pipeline config (from ${CONF_FILE}) ==="
log "  WORK_DIR    = ${WORK_DIR}"
log "  TREAT       = ${TREAT}  (R1=${TREAT_R1})"
log "  CTRL        = ${CTRL}   (R1=${CTRL_R1})"
log "  REF_DIR     = ${REF_DIR}"
log "  BT2_INDEX   = ${BT2_INDEX}"
log "  THREADS     = ${THREADS}  (sort threads = ${SORT_THREADS})"
log "  TRIM_LEN    = ${TRIM_LEN}"
log "  MAX_INS     = ${BOWTIE2_MAXINS}"
log "  RESUME      = ${RESUME}"


log_step "Step 1: repair FASTQs"
if [[ ${#PROCESS_SAMPLES[@]} -eq 0 ]]; then
    log "  All samples have REUSE_BAM set — skipping Steps 1-4 entirely."
fi
for s in "${PROCESS_SAMPLES[@]}"; do
    out_r1="fastq/${s}.R1.fastq.gz"
    out_r2="fastq/${s}.R2.fastq.gz"

    if have_outputs "${out_r1}" "${out_r2}"; then
        log "  -> ${s} already in fastq/, skipping"
        continue
    fi

    # If a previous run left broken symlinks/files behind, clear them.
    rm -f "${out_r1}" "${out_r2}"

    if [[ "${PRE_PAIRED[$s]}" == "1" ]]; then
        log "  -> ${s}: marked PRE_PAIRED — staging via symlink (no re-pairing)"
        ln -s "${SRC_R1[$s]}" "${out_r1}"
        ln -s "${SRC_R2[$s]}" "${out_r2}"
    else
        log "  -> ${s}: running seqkit pair"
        seqkit pair \
            -1 "${SRC_R1[$s]}" \
            -2 "${SRC_R2[$s]}" \
            -O fastq/ \
            --id-regexp '^(\S+)' \
            -j "${THREADS}" \
            2> "logs/${s}.seqkit.log"
        mv "fastq/$(basename "${SRC_R1[$s]}")" "${out_r1}"
        mv "fastq/$(basename "${SRC_R2[$s]}")" "${out_r2}"
    fi

    n1=$(zcat "${out_r1}" | awk 'END{print NR/4}')
    n2=$(zcat "${out_r2}" | awk 'END{print NR/4}')
    log "     R1=${n1}  R2=${n2}"
    if [[ "${n1}" != "${n2}" ]]; then
        log_error "${s}: R1 and R2 read counts differ (${n1} vs ${n2})"
        if [[ "${PRE_PAIRED[$s]}" == "1" ]]; then
            log_error "  -> PRE_PAIRED=1 but inputs are NOT actually paired."
            log_error "     Set ${s^^}_PRE_PAIRED=0 in the conf to run seqkit pair."
        fi
        exit 1
    fi
done

log_step "Step 2: hard-trim to ${TRIM_LEN} bp"
for s in "${PROCESS_SAMPLES[@]}"; do
    out_r1="trim37/${s}.R1.fq.gz"
    out_r2="trim37/${s}.R2.fq.gz"

    if have_outputs "${out_r1}" "${out_r2}"; then
        log "  -> ${s} already trimmed, skipping"
        continue
    fi


    if [[ "${SKIP_TRIM[$s]}" == "1" ]]; then
        log "  -> ${s}: marked SKIP_TRIM — passing fastq through untrimmed"
        rm -f "${out_r1}" "${out_r2}"
        ln -s "../fastq/${s}.R1.fastq.gz" "${out_r1}"
        ln -s "../fastq/${s}.R2.fastq.gz" "${out_r2}"
        continue
    fi

    log "  -> ${s}"
    cutadapt \
        -l "${TRIM_LEN}" -j "${THREADS}" \
        -o "${out_r1}" -p "${out_r2}" \
        "fastq/${s}.R1.fastq.gz" \
        "fastq/${s}.R2.fastq.gz" \
        > "logs/${s}.cutadapt.log"
done

log_step "Step 3: align"
for s in "${PROCESS_SAMPLES[@]}"; do
    if have_outputs "bam/${s}.sorted.bam" "bam/${s}.sorted.bam.bai"; then
        log "  -> ${s} already aligned, skipping"
        continue
    fi
    log "  -> ${s}"
    bowtie2 -x "${BT2_INDEX}" \
        -1 "trim37/${s}.R1.fq.gz" \
        -2 "trim37/${s}.R2.fq.gz" \
        --maxins "${BOWTIE2_MAXINS}" \
        -p "${THREADS}" \
        --rg-id "${s}" \
        --rg "SM:${s}" --rg "LB:${s}" --rg "PL:ILLUMINA" --rg "PU:unit1" \
        2> "logs/${s}.bt2.log" \
      | samtools sort -@ "${SORT_THREADS}" -m 2G -o "bam/${s}.sorted.bam" -
    samtools index "bam/${s}.sorted.bam"
    grep "overall alignment rate" "logs/${s}.bt2.log" || true
done


log_step "Step 4: dedup + filter"
for s in "${PROCESS_SAMPLES[@]}"; do
    if have_outputs "bam/${s}.final.bam" "bam/${s}.final.bam.bai"; then
        log "  -> ${s} already filtered, skipping"
        continue
    fi
    log "  -> Picard MarkDuplicates for ${s}"

    _JAVA_OPTIONS="-Xmx${JAVA_HEAP_PICARD}" picard MarkDuplicates \
        I="bam/${s}.sorted.bam" \
        O="bam/${s}.dedup.bam" \
        M="logs/${s}.dup.metrics" \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=logs \
        > "logs/${s}.picard.stdout.log" \
        2> "logs/${s}.picard.stderr.log"

    samtools index "bam/${s}.dedup.bam"

    samtools view -@ "${SORT_THREADS}" -b -f 2 -q "${DEDUP_MIN_MAPQ}" \
        "bam/${s}.dedup.bam" \
        -o "bam/${s}.final.bam"
    samtools index "bam/${s}.final.bam"
    samtools flagstat "bam/${s}.final.bam" > "logs/${s}.flagstat.txt"
done

log_step "Step 5: insert-size metrics"
for s in "${SAMPLES[@]}"; do
    if have_outputs "logs/${s}.insertsize.metrics" "logs/${s}.insertsize.pdf"; then
        log "  -> ${s} already done, skipping"
        continue
    fi
    log "  -> ${s}"
    _JAVA_OPTIONS="-Xmx${JAVA_HEAP_INSERTSIZE}" picard CollectInsertSizeMetrics \
        I="bam/${s}.final.bam" \
        O="logs/${s}.insertsize.metrics" \
        H="logs/${s}.insertsize.pdf" \
        VALIDATION_STRINGENCY=LENIENT \
        2> "logs/${s}.insertsize.log"
done

# =============================================================================
# 6. bigWig coverage tracks
# =============================================================================
log_step "Step 6: bigWig tracks"
for s in "${SAMPLES[@]}"; do
    if have_outputs "bw/${s}.cpm.bw"; then
        log "  -> ${s} CPM bw exists, skipping"
        continue
    fi
    log "  -> ${s} CPM bigWig"
    bamCoverage \
        -b "bam/${s}.final.bam" \
        -o "bw/${s}.cpm.bw" \
        --binSize "${BIGWIG_BIN_SIZE}" \
        --normalizeUsing CPM \
        --extendReads \
        -p "${THREADS}" \
        2> "logs/${s}.bamcoverage.log"
done

LOG2_BW="bw/${TREAT}.log2ratio.bw"
if have_outputs "${LOG2_BW}"; then
    log "  -> log2 ratio bigWig exists, skipping"
else
    log "  -> log2(treat/input) bigWig"
    bamCompare \
        -b1 "bam/${TREAT}.final.bam" \
        -b2 "bam/${CTRL}.final.bam" \
        -o "${LOG2_BW}" \
        --operation log2 \
        --binSize "${BIGWIG_BIN_SIZE}" \
        --scaleFactorsMethod readCount \
        --extendReads \
        -p "${THREADS}" \
        2> "logs/${TREAT}.bamcompare.log"
fi


callpeaks () {
    local tag=$1 t_bam=$2 c_bam=$3 outdir=$4
    if have_outputs "${outdir}/${tag}.peaks.filt.bed"; then
        log "  -> ${tag} peaks exist, skipping"
        return
    fi
    log "==== MACS3: ${tag} ===="
    mkdir -p "${outdir}"
    macs3 callpeak \
        -t "${t_bam}" \
        -c "${c_bam}" \
        -f BAMPE -g "${MACS_GENOME}" \
        -n "${tag}" --outdir "${outdir}" \
        --keep-dup all \
        2> "logs/${tag}.macs2.log"

    if [[ -f "${BLACKLIST}" ]]; then
        bedtools intersect -v \
            -a "${outdir}/${tag}_peaks.narrowPeak" \
            -b "${BLACKLIST}" \
            > "${outdir}/${tag}.peaks.filt.bed"
    else
        log_warn "blacklist not found at ${BLACKLIST} — using unfiltered peaks"
        cp "${outdir}/${tag}_peaks.narrowPeak" "${outdir}/${tag}.peaks.filt.bed"
    fi
    log "  $(wc -l < "${outdir}/${tag}.peaks.filt.bed") peaks after blacklist"
}

log_step "Step 7: peak calling"
callpeaks "${TREAT}" "bam/${TREAT}.final.bam" "bam/${CTRL}.final.bam" "peaks/${TREAT}"


log_step "Step 8: pseudo-replicate IDR"

PSR1="bam/${TREAT}.psr1.bam"
PSR2="bam/${TREAT}.psr2.bam"

if have_outputs "${PSR1}" "${PSR2}" "${PSR1}.bai" "${PSR2}.bai"; then
    log "  -> pseudo-replicate BAMs already present"
else
    nreads=$(samtools view -c "bam/${TREAT}.final.bam")
    half=$(( nreads / 2 ))
    log "  total reads=${nreads}  half=${half}"

    samtools sort -n -@ "${SORT_THREADS}" -m 2G \
        -o "bam/${TREAT}.nsort.bam" "bam/${TREAT}.final.bam"

    # Pair-preserving random split (seed for reproducibility)
    samtools view -h "bam/${TREAT}.nsort.bam" \
      | awk -v seed="${PSR_SEED}" 'BEGIN{srand(seed); OFS="\t"}
             /^@/ {print; next}
             {
               if ($1 != last) { r = rand(); last = $1 }
               print r"\t"$0
             }' \
      | sort -k1,1 \
      | cut -f2- \
      | awk -v half="${half}" -v out1=psr1.sam -v out2=psr2.sam '
             BEGIN{ n=0 }
             /^@/ { print > out1; print > out2; next }
             {
               if (n < half) print > out1; else print > out2;
               if ($1 != last) n++;
               last = $1;
             }'

    samtools sort -@ "${SORT_THREADS}" -m 2G -o "${PSR1}" psr1.sam && samtools index "${PSR1}"
    samtools sort -@ "${SORT_THREADS}" -m 2G -o "${PSR2}" psr2.sam && samtools index "${PSR2}"
    rm -f psr1.sam psr2.sam "bam/${TREAT}.nsort.bam"
fi

callpeaks "${TREAT}_psr1" "${PSR1}" "bam/${CTRL}.final.bam" "peaks/${TREAT}_psr1"
callpeaks "${TREAT}_psr2" "${PSR2}" "bam/${CTRL}.final.bam" "peaks/${TREAT}_psr2"

log "Building master peak list (union of pseudo-rep peaks)"
cat "peaks/${TREAT}_psr1/${TREAT}_psr1.peaks.filt.bed" \
    "peaks/${TREAT}_psr2/${TREAT}_psr2.peaks.filt.bed" \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - -c 7 -o max \
  | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"peak_"NR,0,".",$4,-1,-1,-1}' \
  > "peaks/${TREAT}.master.narrowPeak"

log "  master peaks: $(wc -l < "peaks/${TREAT}.master.narrowPeak")"

if have_outputs "peaks/${TREAT}.idr.txt"; then
    log "  -> IDR output exists, skipping"
else
    idr --samples "peaks/${TREAT}_psr1/${TREAT}_psr1.peaks.filt.bed" \
                  "peaks/${TREAT}_psr2/${TREAT}_psr2.peaks.filt.bed" \
        --peak-list "peaks/${TREAT}.master.narrowPeak" \
        --input-file-type narrowPeak \
        --output-file "peaks/${TREAT}.idr.txt" \
        --rank signal.value \
        --soft-idr-threshold 0.1 \
        --plot \
        2> "logs/${TREAT}.idr.log"
fi

# score = -125 * log2(IDR)
awk -v t="${IDR_SCORE_10}" 'BEGIN{OFS="\t"} $5 >= t' \
    "peaks/${TREAT}.idr.txt" > "peaks/${TREAT}.idr0.1.bed"
awk -v t="${IDR_SCORE_05}" 'BEGIN{OFS="\t"} $5 >= t' \
    "peaks/${TREAT}.idr.txt" > "peaks/${TREAT}.idr0.05.bed"

# =============================================================================
# 9. Summary tables
# =============================================================================
log_step "Read accounting"
{
  printf "%-20s %12s %12s %12s %12s %8s\n" \
    "sample" "raw_pairs" "aligned" "dedup" "final_MAPQ${DEDUP_MIN_MAPQ}" "%kept"
  for s in "${SAMPLES[@]}"; do
    final=$(samtools view -c "bam/${s}.final.bam")

    # Reused BAM: intermediate files don't exist; show "(reused)" for them.
    if [[ -n "${REUSE_BAM[$s]}" ]]; then
        printf "%-20s %12s %12s %12s %12d %8s\n" \
            "${s}" "(reused)" "(reused)" "(reused)" "${final}" "-"
        continue
    fi

    raw=$(grep -E "^[0-9]+ reads" "logs/${s}.bt2.log" 2>/dev/null | awk '{print $1}' || echo 0)
    [[ -z "${raw}" ]] && raw=0
    aligned=$(samtools view -c -F 4 "bam/${s}.sorted.bam")
    dedup=$(samtools view -c -F 4 "bam/${s}.dedup.bam")
    pct=$(awk -v a="${final}" -v b="${raw}" \
        'BEGIN{ if(b>0) printf "%.1f", 100*a/(2*b); else print 0}')
    printf "%-20s %12d %12d %12d %12d %8s\n" "${s}" "${raw}" "${aligned}" "${dedup}" "${final}" "${pct}"
  done
} | tee logs/map_summary.txt

log_step "Peak counts"
{
  printf "%-35s %12s\n" "set" "n_peaks"
  for f in "peaks/${TREAT}/${TREAT}.peaks.filt.bed" \
           "peaks/${TREAT}_psr1/${TREAT}_psr1.peaks.filt.bed" \
           "peaks/${TREAT}_psr2/${TREAT}_psr2.peaks.filt.bed" \
           "peaks/${TREAT}.master.narrowPeak" \
           "peaks/${TREAT}.idr0.1.bed" \
           "peaks/${TREAT}.idr0.05.bed"; do
      [[ -f "${f}" ]] && printf "%-35s %12d\n" "$(basename "${f}")" "$(wc -l < "${f}")"
  done
} | tee logs/peak_summary.txt


