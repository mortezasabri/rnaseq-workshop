# Settings for run_all.sh and every step script.
#
# Each line has the form  NAME=${NAME:-default}  so any value can be overridden for one run
# from the command line, e.g.:   THREADS=4 N_READS=500000 ./run_all.sh
# The R scripts read this same file, so keep that one-line form for every setting.

# --- Input ------------------------------------------------------------------------------
SAMPLES=${SAMPLES:-config/samples.tsv}  # sample sheet; see config/samples.template.tsv
N_READS=${N_READS:-2000000}             # reads (pairs) kept from the start of each FASTQ; 0 = all
THREADS=${THREADS:-8}                   # CPU threads for trimming, alignment, counting

# --- Reference: human genome GRCh38 with Ensembl gene annotation -----------------------
ENSEMBL_RELEASE=${ENSEMBL_RELEASE:-116}
HISAT2_INDEX_URL=${HISAT2_INDEX_URL:-https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz}
GTF_URL=${GTF_URL:-https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz}
CDNA_URL=${CDNA_URL:-https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz}

# --- Library and trimming ---------------------------------------------------------------
STRANDED=${STRANDED:-0}                 # featureCounts -s: 0 unstranded, 1 stranded, 2 reverse
TRIM_WINDOW=${TRIM_WINDOW:-4:20}        # Trimmomatic SLIDINGWINDOW: cut when 4-base mean quality < 20
TRIM_MINLEN=${TRIM_MINLEN:-25}          # Trimmomatic MINLEN: drop reads shorter than this after trimming
RUN_SALMON=${RUN_SALMON:-yes}           # yes = also quantify with Salmon (alignment-free route)

# --- Downstream (R) ---------------------------------------------------------------------
COUNTS_SOURCE=${COUNTS_SOURCE:-pipeline}  # pipeline = our featureCounts table
                                          # airway   = published full-depth table (no upstream needed)
REFERENCE_LEVEL=${REFERENCE_LEVEL:-untreated}  # the condition everything is compared against
PADJ_CUTOFF=${PADJ_CUTOFF:-0.05}        # adjusted p-value (FDR) cutoff for "significant"
LFC_CUTOFF=${LFC_CUTOFF:-1}             # |log2 fold change| cutoff for gene lists (1 = 2-fold)
SHOWCASE_GENE=${SHOWCASE_GENE:-CRISPLD2}  # gene shown as read coverage, from raw reads to result
SEED=${SEED:-42}                        # random seed (GSEA permutations), for reproducible results
