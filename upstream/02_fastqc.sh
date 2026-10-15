#!/usr/bin/env bash
# Step 2: quality control of the raw reads with FastQC.
# One HTML report per FASTQ file in results/upstream/fastqc_raw/. MultiQC (step 7) combines them.
source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
need fastqc
read_samples
out="$UP/fastqc_raw"
mkdir -p "$out"

todo=()
for s in "${SAMPLE_IDS[@]}"; do
  for fq in "$(raw_r1 "$s")" $([[ $PAIRED == yes ]] && raw_r2 "$s"); do
    [[ -s $fq ]] || die "missing $fq; run step 1 (01_download.sh) first"
    [[ -s $out/$(basename "${fq%.fastq.gz}")_fastqc.zip ]] || todo+=("$fq")
  done
done

if (( ${#todo[@]} == 0 )); then log "FastQC reports already present"; exit 0; fi
log "FastQC on ${#todo[@]} raw FASTQ files"
fastqc --quiet --threads "$THREADS" --outdir "$out" "${todo[@]}"
log "Reports in $out"
