#!/usr/bin/env bash
# Step 3: trim adapters and low-quality ends with Trimmomatic, then re-check with FastQC.
#   ILLUMINACLIP   remove Illumina adapter sequence
#   SLIDINGWINDOW  cut the read once the mean quality in a window drops too low (TRIM_WINDOW)
#   MINLEN         drop reads that became too short (TRIM_MINLEN)
# Paired-end: a pair is kept only if both reads survive; the orphans go to *.unpaired.fastq.gz.
source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
need trimmomatic fastqc
read_samples
mkdir -p "$TRIM_DIR" "$UP/trimmomatic" "$UP/fastqc_trimmed"

# Adapter files ship with Trimmomatic; ADAPTER_DIR overrides the location if needed.
ADAPTER_DIR=${ADAPTER_DIR:-$(dirname "$(readlink -f "$(command -v trimmomatic)")")/adapters}
adapters="$ADAPTER_DIR/$([[ $PAIRED == yes ]] && echo TruSeq3-PE-2.fa || echo TruSeq3-SE.fa)"
[[ -s $adapters ]] || die "adapter file not found: $adapters. Set ADAPTER_DIR to Trimmomatic's adapters folder."
steps=("ILLUMINACLIP:$adapters:2:30:10:2:True" "SLIDINGWINDOW:$TRIM_WINDOW" "MINLEN:$TRIM_MINLEN")
# Trimmomatic runs out of Java memory with many threads (its default limit is 1 GB);
# 4 threads with 4 GB is plenty for a few million reads.
trim_threads=$(( THREADS < 4 ? THREADS : 4 ))
java_mem=-Xmx4g

for s in "${SAMPLE_IDS[@]}"; do
  if [[ -s $(trim_r1 "$s") ]]; then log "  $s already trimmed"; continue; fi
  log "  trimming $s"
  if [[ $PAIRED == yes ]]; then
    trimmomatic "$java_mem" PE -threads "$trim_threads" -phred33 \
      "$(raw_r1 "$s")" "$(raw_r2 "$s")" \
      "$TRIM_DIR/$s.r1.part.fastq.gz" "$TRIM_DIR/${s}_R1.unpaired.fastq.gz" \
      "$TRIM_DIR/$s.r2.part.fastq.gz" "$TRIM_DIR/${s}_R2.unpaired.fastq.gz" \
      "${steps[@]}" 2> "$UP/trimmomatic/$s.trimmomatic.log" ||
      die "Trimmomatic failed for $s; see $UP/trimmomatic/$s.trimmomatic.log"
    mv "$TRIM_DIR/$s.r2.part.fastq.gz" "$(trim_r2 "$s")"
  else
    trimmomatic "$java_mem" SE -threads "$trim_threads" -phred33 \
      "$(raw_r1 "$s")" "$TRIM_DIR/$s.r1.part.fastq.gz" \
      "${steps[@]}" 2> "$UP/trimmomatic/$s.trimmomatic.log" ||
      die "Trimmomatic failed for $s; see $UP/trimmomatic/$s.trimmomatic.log"
  fi
  mv "$TRIM_DIR/$s.r1.part.fastq.gz" "$(trim_r1 "$s")"
  grep -E "^Input" "$UP/trimmomatic/$s.trimmomatic.log" | sed 's/^/      /'
done

todo=()
for s in "${SAMPLE_IDS[@]}"; do
  for fq in "$(trim_r1 "$s")" $([[ $PAIRED == yes ]] && trim_r2 "$s"); do
    [[ -s $UP/fastqc_trimmed/$(basename "${fq%.fastq.gz}")_fastqc.zip ]] || todo+=("$fq")
  done
done
if (( ${#todo[@]} > 0 )); then
  log "FastQC on ${#todo[@]} trimmed FASTQ files"
  fastqc --quiet --threads "$THREADS" --outdir "$UP/fastqc_trimmed" "${todo[@]}"
fi
log "Trimmed reads in $TRIM_DIR, reports in $UP/trimmomatic and $UP/fastqc_trimmed"
