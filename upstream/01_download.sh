#!/usr/bin/env bash
# Step 1: get the reads.
# Takes the first N_READS reads of each FASTQ in the sample sheet (a URL or a local file) and
# writes them to data/fastq/<sample>_R1.fastq.gz (and _R2 for paired-end). With a URL only that
# part of the file is downloaded, so a 2-million-read subset of a 20-million-read run is quick.
# N_READS=0 keeps every read.
source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
need curl gzip
read_samples
mkdir -p "$FASTQ_DIR"
JOBS=${DOWNLOAD_JOBS:-4}   # samples fetched at the same time

# fetch <source> <destination.fastq.gz>
fetch() {
  local src=$1 dest=$2 tmp="$2.part" lines=$(( N_READS * 4 ))
  if [[ -s $dest ]]; then log "  $(basename "$dest") already present"; return; fi
  local -a reader
  if [[ $src =~ ^(https?|ftp):// ]]; then
    reader=(curl -fsL --retry 3 "$src")   # errors are reported below by exit code
  else
    [[ -f $src ]] || die "FASTQ file not found: $src"
    reader=(cat "$src")
  fi

  local -a st
  if (( N_READS > 0 )); then
    # head stops after N_READS reads and closes the pipe, which deliberately cuts the download
    # short. The "broken pipe" this gives curl/cat/gzip is expected, so it is not an error here.
    set +o pipefail
    "${reader[@]}" | gzip -dcf 2>/dev/null | head -n "$lines" | gzip -c > "$tmp"
    st=("${PIPESTATUS[@]}")
    set -o pipefail
  else
    set +o pipefail
    "${reader[@]}" | gzip -dcf | gzip -c > "$tmp"
    st=("${PIPESTATUS[@]}")
    set -o pipefail
  fi
  # curl: 0 = finished, 23 = stopped by head on purpose; cat/gzip: 141 = stopped by head on purpose
  [[ ${st[0]} == 0 || ${st[0]} == 23 || ${st[0]} == 141 ]] ||
    die "reading $src failed (exit ${st[0]}): network problem or wrong URL/path. Rerun to retry."

  local n_lines n
  n_lines=$(gzip -dc "$tmp" | wc -l)
  n=$(( n_lines / 4 ))
  (( n > 0 && n_lines % 4 == 0 )) || die "$src gave no complete FASTQ records (is it a FASTQ file?)"
  (( N_READS == 0 || n == N_READS )) || log "  note: $(basename "$src") has only $n reads (asked for $N_READS)"
  mv "$tmp" "$dest"
  log "  $(basename "$dest"): $n reads"
}

fetch_sample() {
  local i=$1 s=${SAMPLE_IDS[$1]}
  fetch "${FQ1[$i]}" "$(raw_r1 "$s")"
  if [[ $PAIRED == yes ]]; then
    fetch "${FQ2[$i]}" "$(raw_r2 "$s")"
    local n1 n2
    n1=$(gzip -dc "$(raw_r1 "$s")" | wc -l); n2=$(gzip -dc "$(raw_r2 "$s")" | wc -l)
    if (( n1 != n2 )); then
      rm -f "$(raw_r1 "$s")" "$(raw_r2 "$s")"
      die "$s: read 1 and read 2 files have different numbers of reads ($((n1/4)) vs $((n2/4)))"
    fi
  fi
}

log "Fetching ${#SAMPLE_IDS[@]} samples (paired-end: $PAIRED, reads per file: $([[ $N_READS == 0 ]] && echo all || echo "$N_READS"))"
for i in "${!SAMPLE_IDS[@]}"; do
  fetch_sample "$i" &
  while (( $(jobs -rp | wc -l) >= JOBS )); do wait -n 2>/dev/null || true; done
done
wait || true

# A sample counts as done only when all of its files are in place.
missing=()
for s in "${SAMPLE_IDS[@]}"; do
  [[ -s $(raw_r1 "$s") && ( $PAIRED == no || -s $(raw_r2 "$s") ) ]] || missing+=("$s")
done
(( ${#missing[@]} == 0 )) ||
  die "could not fetch: ${missing[*]}. See the messages above; rerun to retry (finished files are kept)."
log "Reads ready in $FASTQ_DIR"
