#!/usr/bin/env bash
# Step 0: get the reference.
#   - HISAT2 index of the human genome GRCh38 (prebuilt; ~4 GB download, needs ~8 GB RAM to use)
#   - Ensembl gene annotation (GTF): where each gene and exon is on the genome
#   - Ensembl transcript sequences (cDNA) + a Salmon index, for the alignment-free route
# Everything goes to data/reference/ and is downloaded only once.
source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
need curl tar gzip hisat2-inspect
mkdir -p "$REF_DIR"

# download <url> <destination>: resumable, only if missing
download() {
  local url=$1 dest=$2
  if [[ -s $dest ]]; then log "  $(basename "$dest") already present"; return; fi
  log "  downloading $(basename "$dest") ..."
  curl -fL --retry 3 -C - -sS -o "$dest.part" "$url" || die "download failed: $url (rerun to resume)"
  mv "$dest.part" "$dest"
}

log "Genome index (HISAT2, GRCh38)"
if [[ ! -s $HISAT2_INDEX.1.ht2 ]]; then
  download "$HISAT2_INDEX_URL" "$REF_DIR/grch38_genome.tar.gz"
  log "  unpacking ..."
  tar -xzf "$REF_DIR/grch38_genome.tar.gz" -C "$REF_DIR"
  [[ -s $HISAT2_INDEX.1.ht2 ]] || die "unpacked the index but $HISAT2_INDEX.1.ht2 is missing; check HISAT2_INDEX_URL"
  rm "$REF_DIR/grch38_genome.tar.gz"
else
  log "  already present"
fi

log "Gene annotation (Ensembl GTF, release $ENSEMBL_RELEASE)"
download "$GTF_URL" "$GTF"

# The genome index and the GTF must name chromosomes the same way ("1" vs "chr1"),
# otherwise no read is ever assigned to a gene. Check once, here, rather than find out later.
idx_chr=$(hisat2-inspect -n "$HISAT2_INDEX" | awk '{print $1}' | sort -u)
gtf_chr=$(gzip -dc "$GTF" | awk -F'\t' '!/^#/ {print $1}' | sort -u)
shared=$(comm -12 <(echo "$idx_chr") <(echo "$gtf_chr") | wc -l)
(( shared >= 20 )) || die "genome index and GTF share only $shared chromosome names (index: $(echo $idx_chr | head -c 60)...; GTF: $(echo $gtf_chr | head -c 60)...). Use matching sources."
log "  genome index and GTF share $shared chromosome names: OK"

if [[ $RUN_SALMON == yes ]]; then
  need salmon
  log "Transcript sequences (Ensembl cDNA) and Salmon index"
  download "$CDNA_URL" "$CDNA"
  if [[ ! -s $TX2GENE ]]; then
    # Header lines look like: >ENST00000415118.1 cdna chromosome:GRCh38:14:... gene:ENSG00000223997.1 ...
    gzip -dc "$CDNA" | awk '/^>/ { tx = substr($1, 2); for (i = 2; i <= NF; i++) if ($i ~ /^gene:/) { g = substr($i, 6); sub(/\.[0-9]+$/, "", g); print tx "\t" g } }' > "$TX2GENE.part"
    mv "$TX2GENE.part" "$TX2GENE"
  fi
  log "  transcript-to-gene table: $(wc -l < "$TX2GENE") transcripts"
  if ! salmon_index_ready; then
    log "  building Salmon index (a few minutes) ..."
    rm -rf "$SALMON_INDEX"
    salmon index -t "$CDNA" -i "$SALMON_INDEX" -p "$THREADS" > "$REF_DIR/salmon_index.log" 2>&1 ||
      die "salmon index failed; see $REF_DIR/salmon_index.log"
  else
    log "  Salmon index already present"
  fi
fi

log "Reference ready in $REF_DIR"
