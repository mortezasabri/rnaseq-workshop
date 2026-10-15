#!/usr/bin/env bash
# Step 7: one report for everything. MultiQC collects the FastQC, Trimmomatic, HISAT2, samtools,
# featureCounts and Salmon outputs into results/upstream/multiqc/multiqc_report.html and also
# exports the plots as images (multiqc_plots/) for the slides.
source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
need multiqc
out="$UP/multiqc"
log "MultiQC over $UP"
rm -rf "$out"    # generated folder; start clean so no plot from an earlier run lingers
multiqc "$UP" --outdir "$out" --force --export --no-ai --title "RNA-seq workshop: upstream QC" \
  --filename multiqc_report.html > "$UP/multiqc.log" 2>&1 || die "MultiQC failed; see $UP/multiqc.log"

# The report gets shared: replace the local folder path MultiQC writes into it with a relative one.
report="$out/multiqc_report.html"
sed "s|$REPO/||g" "$report" > "$report.tmp" && mv "$report.tmp" "$report"
log "Report: $out/multiqc_report.html"
