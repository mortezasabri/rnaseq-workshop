# Shared settings and helpers for the upstream step scripts. Sourced, not run.
# Loads config/config.sh, reads and checks the sample sheet, and defines the folder layout.

set -euo pipefail

REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$REPO"
source config/config.sh

# --- Folder layout ----------------------------------------------------------------------
# data/     large files, never committed: reads, reference, alignments
# results/  small outputs that are worth keeping and showing
REF_DIR=data/reference
FASTQ_DIR=data/fastq
TRIM_DIR=data/trimmed
BAM_DIR=data/bam
UP=results/upstream

GTF="$REF_DIR/$(basename "$GTF_URL")"
CDNA="$REF_DIR/$(basename "$CDNA_URL")"
HISAT2_INDEX="$REF_DIR/grch38/genome"
SALMON_INDEX="$REF_DIR/salmon_index"
TX2GENE="$REF_DIR/tx2gene.tsv"

log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
die() { printf '[%s] ERROR: %s\n' "$(date +%H:%M:%S)" "$*" >&2; exit 1; }

# need <tool>... : stop with a clear message if a program is missing
need() {
  local t
  for t in "$@"; do
    command -v "$t" >/dev/null 2>&1 ||
      die "'$t' not found. Activate the environment first:  micromamba activate rnaseq-workshop  (see README.md)"
  done
}

# --- Sample sheet -----------------------------------------------------------------------
# Fills the arrays SAMPLE_IDS, CONDITIONS, BATCHES, FQ1, FQ2 and sets PAIRED=yes|no.
read_samples() {
  [[ -f $SAMPLES ]] || die "sample sheet not found: $SAMPLES"
  SAMPLE_IDS=(); CONDITIONS=(); BATCHES=(); FQ1=(); FQ2=()
  local header_seen=no line n_paired=0
  while IFS= read -r line || [[ -n $line ]]; do
    line=${line%$'\r'}
    [[ -z ${line//[[:space:]]/} || $line == \#* ]] && continue
    local -a f
    IFS=$'\t' read -r -a f <<< "$line"
    if [[ $header_seen == no ]]; then
      [[ "${f[*]}" == "sample condition batch fastq_1 fastq_2" ]] ||
        die "$SAMPLES: the header must be the five TAB-separated columns: sample condition batch fastq_1 fastq_2 (found: ${f[*]})"
      header_seen=yes; continue
    fi
    (( ${#f[@]} == 5 )) || die "$SAMPLES: this row does not have 5 TAB-separated columns: $line"
    [[ ${f[0]} =~ ^[A-Za-z0-9_-]+$ ]] || die "$SAMPLES: sample name '${f[0]}' may only use letters, digits, _ and -"
    SAMPLE_IDS+=("${f[0]}"); CONDITIONS+=("${f[1]}"); BATCHES+=("${f[2]}"); FQ1+=("${f[3]}"); FQ2+=("${f[4]}")
    [[ ${f[4]} != "-" ]] && n_paired=$((n_paired + 1))
  done < "$SAMPLES"

  local n=${#SAMPLE_IDS[@]}
  (( n >= 2 )) || die "$SAMPLES: needs at least 2 samples, found $n"
  local dup
  dup=$(printf '%s\n' "${SAMPLE_IDS[@]}" | sort | uniq -d | head -1)
  [[ -z $dup ]] || die "$SAMPLES: sample name '$dup' appears more than once"
  local levels
  levels=$(printf '%s\n' "${CONDITIONS[@]}" | sort -u | tr '\n' ' ')
  (( $(printf '%s\n' "${CONDITIONS[@]}" | sort -u | wc -l) == 2 )) ||
    die "$SAMPLES: the condition column needs exactly 2 different values, found: $levels"
  [[ " $levels" == *" $REFERENCE_LEVEL "* ]] ||
    die "REFERENCE_LEVEL is '$REFERENCE_LEVEL' but the conditions are: $levels. Set REFERENCE_LEVEL in config/config.sh"
  if (( n_paired == n )); then PAIRED=yes
  elif (( n_paired == 0 )); then PAIRED=no
  else die "$SAMPLES: mix of paired-end and single-end rows; use one layout for all samples"
  fi
}

# fastq paths after each step, by sample name
raw_r1()  { echo "$FASTQ_DIR/$1_R1.fastq.gz"; }
raw_r2()  { echo "$FASTQ_DIR/$1_R2.fastq.gz"; }
trim_r1() { echo "$TRIM_DIR/$1_trimmed_R1.fastq.gz"; }   # "trimmed" mid-name: MultiQC strips it from name ends
trim_r2() { echo "$TRIM_DIR/$1_trimmed_R2.fastq.gz"; }
bam_of()  { echo "$BAM_DIR/$1.bam"; }

# Salmon 2.x writes info.json into its index, Salmon 1.x versionInfo.json
salmon_index_ready() { [[ -s $SALMON_INDEX/info.json || -s $SALMON_INDEX/versionInfo.json ]]; }
