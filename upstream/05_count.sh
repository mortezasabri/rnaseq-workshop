#!/usr/bin/env bash
# Step 5: count reads per gene with featureCounts (Subread package).
# A read (pair) is counted for a gene when it overlaps the gene's exons; reads hitting several
# genes or no gene are reported as unassigned (see counts.txt.summary).
# Output: results/upstream/featureCounts/counts.txt         (featureCounts' own format)
#         results/upstream/featureCounts/gene_counts.tsv    (clean matrix: gene_id + one column per sample)
source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
need featureCounts
read_samples
[[ -s $GTF ]] || die "annotation missing; run step 0 (00_reference.sh) first"
out="$UP/featureCounts"
mkdir -p "$out"

bams=()
for s in "${SAMPLE_IDS[@]}"; do
  [[ -s $(bam_of "$s").bai ]] || die "missing alignment for $s; run step 4 (04_align.sh) first"
  bams+=("$(bam_of "$s")")
done
pe=(); [[ $PAIRED == yes ]] && pe=(-p --countReadPairs)

log "featureCounts on ${#bams[@]} BAM files (strandedness -s $STRANDED)"
featureCounts -T "$THREADS" "${pe[@]}" -s "$STRANDED" -a "$GTF" -o "$out/counts.txt" "${bams[@]}" \
  > "$out/featureCounts.log" 2>&1 || die "featureCounts failed; see $out/featureCounts.log"

# Clean matrix: drop the comment line and the coordinate columns, name columns by sample.
{
  printf 'gene_id'; printf '\t%s' "${SAMPLE_IDS[@]}"; printf '\n'
  awk -F'\t' 'NR > 2 { printf "%s", $1; for (i = 7; i <= NF; i++) printf "\t%s", $i; printf "\n" }' "$out/counts.txt"
} > "$out/gene_counts.tsv"

awk -F'\t' 'NR == 1 { next } { for (i = 2; i <= NF; i++) tot[i] += $i } $1 == "Assigned" { for (i = 2; i <= NF; i++) a[i] = $i }
  END { for (i = 2; i <= NF; i++) printf "      %5.1f%% assigned to genes\n", 100 * a[i] / tot[i] }' "$out/counts.txt.summary" |
  paste <(printf '      %s\n' "${SAMPLE_IDS[@]}") -
log "Count matrix: $out/gene_counts.tsv ($(($(wc -l < "$out/gene_counts.tsv") - 1)) genes x ${#SAMPLE_IDS[@]} samples)"
