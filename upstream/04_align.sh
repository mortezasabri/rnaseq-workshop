#!/usr/bin/env bash
# Step 4: align the trimmed reads to the genome with HISAT2 (splice-aware: a read may span an
# intron), then sort and index the alignments with samtools.
# Output: data/bam/<sample>.bam (+ .bai) - open these in IGV to look at the reads.
# Alignment rates go to results/upstream/hisat2/, samtools statistics to results/upstream/samtools/.
source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
need hisat2 samtools
read_samples
[[ -s $HISAT2_INDEX.1.ht2 ]] || die "genome index missing; run step 0 (00_reference.sh) first"
mkdir -p "$BAM_DIR" "$UP/hisat2" "$UP/samtools"

for s in "${SAMPLE_IDS[@]}"; do
  bam=$(bam_of "$s")
  if [[ -s $bam.bai ]]; then log "  $s already aligned"; continue; fi
  [[ -s $(trim_r1 "$s") ]] || die "missing trimmed reads for $s; run step 3 (03_trim.sh) first"
  log "  aligning $s"
  if [[ $PAIRED == yes ]]; then reads=(-1 "$(trim_r1 "$s")" -2 "$(trim_r2 "$s")"); else reads=(-U "$(trim_r1 "$s")"); fi
  hisat2 -p "$THREADS" -x "$HISAT2_INDEX" "${reads[@]}" \
      --new-summary --summary-file "$UP/hisat2/$s.hisat2.txt" 2> "$BAM_DIR/$s.hisat2.stderr" |
    samtools sort -@ 4 -m 1G -o "$bam.part" - ||
    die "alignment failed for $s; see $BAM_DIR/$s.hisat2.stderr"
  mv "$bam.part" "$bam"
  samtools index "$bam"
  samtools flagstat "$bam" > "$UP/samtools/$s.flagstat.txt"
  echo "      $(tail -1 "$UP/hisat2/$s.hisat2.txt")"
done
log "Alignments in $BAM_DIR"
