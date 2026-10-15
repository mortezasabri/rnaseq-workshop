#!/usr/bin/env bash
# Step 6 (optional, RUN_SALMON=yes): the alignment-free route. Salmon maps reads directly to
# transcript sequences and estimates how many reads came from each transcript. It is much faster
# than alignment + counting and is what the DESeq2 authors recommend; the R step 0_upstream
# compares its gene totals with featureCounts. It also reports the library strandedness it detected.
# Output: results/upstream/salmon/<sample>/quant.sf (one row per transcript)
source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
if [[ $RUN_SALMON != yes ]]; then log "RUN_SALMON=$RUN_SALMON: skipping Salmon"; exit 0; fi
need salmon
read_samples
salmon_index_ready || die "Salmon index missing; run step 0 (00_reference.sh) first"
mkdir -p "$UP/salmon"

for s in "${SAMPLE_IDS[@]}"; do
  out="$UP/salmon/$s"
  if [[ -s $out/quant.sf ]]; then log "  $s already quantified"; continue; fi
  log "  quantifying $s"
  if [[ $PAIRED == yes ]]; then reads=(-1 "$(trim_r1 "$s")" -2 "$(trim_r2 "$s")"); else reads=(-r "$(trim_r1 "$s")"); fi
  salmon quant -i "$SALMON_INDEX" -l A "${reads[@]}" -p "$THREADS" -o "$out" \
    > "$UP/salmon/$s.salmon.log" 2>&1 || die "Salmon failed for $s; see $UP/salmon/$s.salmon.log"
  awk -F'"' '/"percent_mapped"|"expected_format"/ { print }' "$out/aux_info/meta_info.json" "$out/lib_format_counts.json" |
    tr -s ' ' | sed 's/^/     /'
done
log "Salmon results in $UP/salmon"
