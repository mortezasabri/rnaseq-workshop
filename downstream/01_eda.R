# 1_eda: exploratory data analysis. Before testing anything: do the samples look sane, and do
# they group by treatment?
#   01_library_sizes          reads per sample
#   02_normalization          size factors: each sample before and after dividing by its own
#   03_pca                    samples in two dimensions (variance-stabilized data)
#   04_sample_distances       how similar is each sample to every other one
# Input: the count matrix (COUNTS_SOURCE in config/config.sh) + the sample sheet.
# Output for the next steps: results/downstream/objects/dds.rds, vsd.rds
source("downstream/_common.R")
suppressPackageStartupMessages({
  library(DESeq2)
  library(org.Hs.eg.db)
  library(pheatmap)
  library(dplyr)
  library(tidyr)
})
step_header("1_eda: exploratory data analysis")
out <- out_dir("1_eda")
samples <- read_samples()

# --- The count matrix: one row per gene, one column per sample ------------------------------
if (cfg$COUNTS_SOURCE == "pipeline") {
  counts <- read.delim(need_file(file.path(UP, "featureCounts", "gene_counts.tsv"),
                                 "the upstream part (./run_all.sh upstream), or set COUNTS_SOURCE=airway"),
                       row.names = 1, check.names = FALSE)
  source_label <- "this pipeline (featureCounts)"
} else if (cfg$COUNTS_SOURCE == "airway") {
  # The published counts for the same eight runs, at full sequencing depth.
  suppressPackageStartupMessages(library(SummarizedExperiment))
  data("airway", package = "airway", envir = environment())
  if (!all(samples$sample %in% colnames(airway))) {
    stop("COUNTS_SOURCE=airway only works with the airway sample sheet (config/samples.tsv).", call. = FALSE)
  }
  counts <- assay(airway)
  source_label <- "published airway counts (full depth)"
} else {
  stop("COUNTS_SOURCE must be 'pipeline' or 'airway', not '", cfg$COUNTS_SOURCE, "'.", call. = FALSE)
}
missing <- setdiff(samples$sample, colnames(counts))
if (length(missing) > 0) stop("Samples missing from the count matrix: ", paste(missing, collapse = ", "), call. = FALSE)
counts <- as.matrix(counts[, samples$sample])
message("  counts from ", source_label, ": ", nrow(counts), " genes x ", ncol(counts), " samples")

# --- DESeqDataSet: counts + sample information + the model ----------------------------------
# With a batch column (here: the donor cell line) the model compares treated vs. untreated
# *within* each donor, which removes donor-to-donor differences from the comparison.
design <- if (has_batch(samples)) ~ batch + condition else ~ condition
dds <- DESeqDataSetFromMatrix(countData = counts, colData = samples, design = design)

# Drop genes with almost no reads: at least 10 reads in at least as many samples as the
# smallest group. They carry no information and only add to the multiple-testing burden.
min_samples <- min(table(samples$condition))
keep <- rowSums(counts(dds) >= 10) >= min_samples
message(sprintf("  keeping %d of %d genes with >= 10 reads in >= %d samples", sum(keep), nrow(dds), min_samples))
dds <- dds[keep, ]

# Gene symbols and Entrez IDs, used for labels and for the enrichment analyses later on
ids <- rownames(dds)
mcols(dds)$symbol <- suppressMessages(mapIds(org.Hs.eg.db, ids, "SYMBOL", "ENSEMBL", multiVals = "first"))
mcols(dds)$entrez <- suppressMessages(mapIds(org.Hs.eg.db, ids, "ENTREZID", "ENSEMBL", multiVals = "first"))

dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)   # variance-stabilized log2-like values, for plots and distances
save_object(dds, "dds")
save_object(vsd, "vsd")
cols <- condition_colors(samples)

# --- 1. Library sizes -----------------------------------------------------------------------
lib <- data.frame(sample = colnames(dds), condition = samples$condition,
                  reads = colSums(counts(dds)) / 1e6, size_factor = sizeFactors(dds))
write_table(lib, out, "01_library_sizes.tsv")
p <- ggplot(lib, aes(reorder(sample, as.integer(condition)), reads, fill = condition)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Reads counted per sample", subtitle = source_label,
       x = NULL, y = "million reads", fill = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
save_plot(p, out, "01_library_sizes.png", width = 6.5, height = 4.5)

# --- 2. Normalization: dividing by size factors lines the samples up ------------------------
# Median-of-ratios, drawn as what it is: every gene's count divided by that gene's average
# (geometric mean) across samples. The median of these ratios is the sample's size factor;
# after dividing by it, every sample's median ratio is 1 (log2 = 0).
# Our pipeline takes the same number of reads from every sample, so its depths are nearly equal
# and there is little to correct. For the airway samples this plot therefore uses the published
# full-depth counts (15-31 M reads per sample), where the correction is visible.
data("airway", package = "airway", envir = environment())
if (all(samples$sample %in% colnames(airway))) {
  suppressPackageStartupMessages(library(SummarizedExperiment))
  norm_counts <- assay(airway)[, samples$sample]
  norm_counts <- norm_counts[rowSums(norm_counts >= 10) >= min_samples, ]
  norm_source <- "published airway counts at full depth"
} else {
  norm_counts <- counts(dds)
  norm_source <- source_label
}
norm_sf <- estimateSizeFactorsForMatrix(norm_counts)
norm_info <- data.frame(sample = colnames(norm_counts), million_reads = colSums(norm_counts) / 1e6,
                        size_factor = norm_sf, counts = norm_source)
write_table(norm_info, out, "02_normalization.tsv")

all_positive <- rowSums(norm_counts == 0) == 0       # the genes the size factors are computed from
log_ratio <- log2(norm_counts[all_positive, ]) - rowMeans(log2(norm_counts[all_positive, ]))
x_labels <- sprintf("%s\n%.0f M reads", norm_info$sample, norm_info$million_reads)
norm_long <- bind_rows(
  data.frame(log_ratio, check.names = FALSE) |> mutate(type = "before: raw counts"),
  data.frame(sweep(log_ratio, 2, log2(norm_sf)), check.names = FALSE) |> mutate(type = "after: divided by size factor")
) |>
  pivot_longer(-type, names_to = "sample", values_to = "log2ratio") |>
  mutate(type = factor(type, c("before: raw counts", "after: divided by size factor")),
         condition = samples$condition[match(sample, samples$sample)],
         sample = factor(sample, norm_info$sample, x_labels))
whiskers <- sapply(split(norm_long$log2ratio, interaction(norm_long$type, norm_long$sample)),
                   function(v) boxplot.stats(v)$stats[c(1, 5)])
y_lim <- range(whiskers) + c(-0.05, 0.35)
sf_labels <- data.frame(type = factor("before: raw counts", levels(norm_long$type)),
                        sample = factor(x_labels, x_labels), y = y_lim[2] - 0.12,
                        label = sprintf("SF\n%.2f", norm_sf))
p <- ggplot(norm_long, aes(sample, log2ratio, fill = condition)) +
  geom_hline(yintercept = 0, color = COL_MUTED, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3, width = 0.7) +
  geom_text(data = sf_labels, aes(sample, y, label = label), inherit.aes = FALSE, size = 3.2, lineheight = 0.9) +
  facet_wrap(~type) +
  coord_cartesian(ylim = y_lim) +
  scale_fill_manual(values = cols) +
  labs(title = "Normalization: divide each sample by its size factor (SF)",
       subtitle = paste0("Each gene's count relative to its average across samples; ", norm_source),
       x = NULL, y = "log2(count / gene's average)", fill = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "top")
save_plot(p, out, "02_normalization.png", width = 9, height = 4.8)

# --- 3. PCA -------------------------------------------------------------------------------------
pca <- plotPCA(vsd, intgroup = c("condition", "batch"), returnData = TRUE, ntop = 500)
pct <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, color = condition, shape = batch)) +
  geom_point(size = 4) +
  scale_color_manual(values = cols) +
  scale_shape_manual(values = c(16, 17, 15, 18, 3, 4, 8)[seq_len(nlevels(samples$batch))]) +
  labs(title = "PCA: samples separate by treatment",
       subtitle = "500 most variable genes, variance-stabilized counts",
       x = paste0("PC1: ", pct[1], "% of variance"), y = paste0("PC2: ", pct[2], "% of variance"),
       color = "condition", shape = if (has_batch(samples)) "cell line / batch" else NULL)
if (!has_batch(samples)) p <- p + guides(shape = "none")
save_plot(p, out, "03_pca.png", width = 7, height = 5)

# --- 4. Sample-to-sample distances ---------------------------------------------------------
d <- dist(t(assay(vsd)))
m <- as.matrix(d)
rownames(m) <- colnames(m) <- paste(vsd$sample, vsd$condition, vsd$batch, sep = " | ")
annotation <- data.frame(condition = vsd$condition, row.names = rownames(m))
save_png(out, "04_sample_distances.png", width = 8, height = 5.5,
  pheatmap(m, clustering_distance_rows = d, clustering_distance_cols = d,
           color = sequential_ramp, border_color = NA, show_colnames = FALSE,
           annotation_row = annotation, annotation_colors = list(condition = cols),
           main = "Sample-to-sample distances", silent = FALSE))
