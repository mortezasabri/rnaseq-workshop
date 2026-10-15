#!/usr/bin/env bash
# Run the whole workshop analysis hands-off:
#   FASTQ files -> quality control -> trimming -> alignment -> count matrix        (upstream)
#   -> exploration -> differential expression -> GO / pathway / disease analysis   (downstream)
#   -> presentation slides                                                         (slides)
#
# Usage:  ./run_all.sh [all|upstream|downstream|slides] [--dry-run]
#   all         everything, in order (default)
#   upstream    upstream/00_reference.sh ... 07_multiqc.sh
#   downstream  downstream/00_upstream_summary.R ... 12_ReactomePA.R
#   slides      slides/rnaseq_intro.qmd -> docs/rnaseq_intro.html
#   --dry-run   check tools, R packages, settings and the sample sheet; run nothing
#
# Settings: config/config.sh. Override any of them for one run:  THREADS=4 ./run_all.sh
# Upstream steps skip work that is already done; delete a step's output to redo it.
#
# Exit codes: 0 = finished, 1 = setup problem (stopped before running anything), 2 = a step failed
# Logs: logs/run_all.log (one line per step), logs/<part>_<step>.log (all output of that step)

set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source config/config.sh

part=all; dry_run=no
for arg in "$@"; do
  case $arg in
    all|upstream|downstream|slides) part=$arg ;;
    --dry-run) dry_run=yes ;;
    -h|--help) sed -n '2,19p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "Unknown argument '$arg'. Usage: ./run_all.sh [all|upstream|downstream|slides] [--dry-run]" >&2; exit 1 ;;
  esac
done
case $part in
  all) parts=(upstream downstream slides) ;;
  *)   parts=("$part") ;;
esac
wants() { [[ " ${parts[*]} " == *" $1 "* ]]; }

upstream_steps=(upstream/0[0-9]_*.sh)
downstream_steps=(downstream/[0-9][0-9]_*.R)
# Steps 01 and 02 produce what every later downstream step needs; if one fails, stop there.
core_steps=(downstream/01_eda.R downstream/02_dea.R)

# --- Checks: everything a run needs, before anything runs ------------------------------
problems=()
problem() { problems+=("$1"); }

if wants upstream; then
  tools=(curl gzip tar fastqc trimmomatic hisat2 hisat2-inspect samtools featureCounts multiqc)
  [[ $RUN_SALMON == yes ]] && tools+=(salmon)
  for t in "${tools[@]}"; do
    command -v "$t" >/dev/null 2>&1 || problem "program '$t' not found"
  done
fi
if wants downstream || wants slides; then
  if ! command -v Rscript >/dev/null 2>&1; then
    problem "Rscript not found"
  elif ! Rscript install.R --check >/dev/null 2>&1; then
    problem "$(Rscript install.R --check 2>&1 | head -1)"
  fi
fi
wants slides && ! command -v quarto >/dev/null 2>&1 && problem "program 'quarto' not found (needed for the slides)"

# The sample sheet, checked by the same code the steps use
sheet_report=$(bash -c 'source upstream/_common.sh && read_samples &&
  printf "%s\n" "${#SAMPLE_IDS[@]}" "$PAIRED" "$(printf "%s\n" "${BATCHES[@]}" | sort -u | wc -l)" &&
  printf "%s\n" "${FQ1[@]}" "${FQ2[@]}"' 2>&1)
if [[ $? -ne 0 ]]; then
  problem "sample sheet: ${sheet_report##*ERROR: }"
  n_samples=?; paired=?; n_batches=?
else
  mapfile -t lines <<< "$sheet_report"
  n_samples=${lines[0]}; paired=${lines[1]}; n_batches=${lines[2]}
  if wants upstream; then
    missing=0
    for f in "${lines[@]:3}"; do
      [[ $f == - || $f =~ ^(https?|ftp):// ]] && continue
      [[ -f $f ]] || { missing=$((missing + 1)); (( missing <= 3 )) && problem "FASTQ file not found: $f"; }
    done
    (( missing > 3 )) && problem "... and $((missing - 3)) more FASTQ files not found"
  fi
fi

if wants downstream && ! wants upstream && [[ $COUNTS_SOURCE == pipeline ]] &&
   [[ ! -s results/upstream/featureCounts/gene_counts.tsv ]]; then
  problem "no count matrix yet: run './run_all.sh upstream' first, or use the published counts: COUNTS_SOURCE=airway ./run_all.sh downstream"
fi
[[ $COUNTS_SOURCE == pipeline || $COUNTS_SOURCE == airway ]] || problem "COUNTS_SOURCE must be 'pipeline' or 'airway', not '$COUNTS_SOURCE'"

echo "RNA-seq workshop pipeline"
echo "  parts:         ${parts[*]}"
echo "  sample sheet:  $SAMPLES ($n_samples samples, paired-end: $paired, batches: $n_batches)"
echo "  comparison:    condition vs. reference level '$REFERENCE_LEVEL'"
wants upstream && echo "  reads:         $([[ $N_READS == 0 ]] && echo "all" || echo "first $N_READS") per FASTQ file, $THREADS threads"
wants downstream && echo "  counts for R:  $COUNTS_SOURCE;  significant = padj < $PADJ_CUTOFF and |log2FC| >= $LFC_CUTOFF"

if (( ${#problems[@]} > 0 )); then
  echo
  echo "Cannot run; fix these first:"
  printf '  - %s\n' "${problems[@]}"
  if printf '%s\n' "${problems[@]}" | grep -qE "not found \(needed|program '|Rscript|R packages"; then
    echo "Missing programs or packages: activate the environment first:  micromamba activate rnaseq-workshop"
  fi
  exit 1
fi
if [[ $dry_run == yes ]]; then
  echo
  echo "Dry run: all checks passed. Steps that would run:"
  wants upstream && printf '  %s\n' "${upstream_steps[@]}"
  wants downstream && printf '  %s\n' "${downstream_steps[@]}"
  wants slides && echo "  slides/render.sh"
  exit 0
fi

# --- Run ----------------------------------------------------------------------------------
mkdir -p logs
failed=()
run_step() {   # run_step <part> <script>; returns the script's exit code
  local part=$1 script=$2 name rc start=$SECONDS
  name=$(basename "${script%.*}")
  local log="logs/${part}_${name}.log"
  printf '[%s] %-11s %-24s ' "$(date +%H:%M:%S)" "$part" "$name"
  case $script in
    *.sh) bash "$script" > "$log" 2>&1 ;;
    *.R)  Rscript "$script" > "$log" 2>&1 ;;
  esac
  rc=$?
  local mins
  mins=$(awk -v s=$((SECONDS - start)) 'BEGIN { printf "%.1f", s / 60 }')
  if (( rc == 0 )); then
    echo "OK      (${mins} min)"
    printf '%s  OK      %s/%s  (%s min)\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$part" "$name" "$mins" >> logs/run_all.log
  else
    echo "FAILED  (${mins} min)  full log: $log"
    tail -n 12 "$log" | sed 's/^/      /'
    printf '%s  FAILED  %s/%s  (%s min)  see %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$part" "$name" "$mins" "$log" >> logs/run_all.log
    failed+=("$part/$name")
  fi
  return $rc
}

echo
if wants upstream; then
  for s in "${upstream_steps[@]}"; do
    # each upstream step needs the previous one, so stop at the first failure
    run_step upstream "$s" || { echo "Stopped: the upstream part needs every step to succeed."; exit 2; }
  done
fi
if wants downstream; then
  for s in "${downstream_steps[@]}"; do
    if ! run_step downstream "$s" && [[ " ${core_steps[*]} " == *" $s "* ]]; then
      echo "Stopped: later steps need the output of $s."; exit 2
    fi
  done
fi
if wants slides; then
  run_step slides slides/render.sh
fi

echo
if (( ${#failed[@]} > 0 )); then
  echo "Finished with ${#failed[@]} failed step(s): ${failed[*]}"
  echo "The other steps ran normally. Fix the cause (shown above and in logs/) and rerun that part."
  exit 2
fi
echo "All done."
wants downstream && echo "  Results:  results/downstream/ (one folder per step), upstream QC: results/upstream/multiqc/multiqc_report.html"
wants slides && echo "  Slides:   docs/rnaseq_intro.html (open in a web browser; press F for full screen)"
exit 0
