# 0_upstream: what happened between the FASTQ files and the count matrix.
#   01_read_funnel        how many reads survive each step (trim -> align -> assigned to a gene)
#   02_coverage_<gene>    the reads on one gene, untreated vs treated: the result is visible already
#   03_featurecounts_vs_salmon   two quantification routes, same answer?
#   04_pipeline_vs_published     our 2-million-read subset vs the published full-depth counts
# Needs the upstream results (./run_all.sh upstream). Skipped when they are not there.
source("downstream/_common.R")
step_header("0_upstream: from reads to counts")

fc_file <- file.path(UP, "featureCounts", "gene_counts.tsv")
if (!file.exists(fc_file)) {
  message("  No upstream results (", fc_file, "); skipping. Run ./run_all.sh upstream to create them.")
  quit(save = "no", status = 0)
}
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(patchwork)
})
samples <- read_samples()
out <- out_dir("0_upstream")

# --- 1. Read funnel ----------------------------------------------------------------------
first_number <- function(pattern, lines) {   # the number right after the first match of pattern
  hit <- grep(pattern, lines, value = TRUE)[1]
  as.numeric(regmatches(hit, regexec(paste0("(?:", pattern, ")\\s*([0-9.]+)"), hit, perl = TRUE))[[1]][2])
}
fc_summary <- read.delim(file.path(UP, "featureCounts", "counts.txt.summary"), check.names = FALSE)
funnel <- bind_rows(lapply(samples$sample, function(s) {
  trim <- readLines(file.path(UP, "trimmomatic", paste0(s, ".trimmomatic.log")))
  input <- first_number("Input Read Pairs:|Input Reads:", trim)
  kept  <- first_number("Both Surviving:|Surviving:", trim)
  hisat <- readLines(file.path(UP, "hisat2", paste0(s, ".hisat2.txt")))
  rate  <- as.numeric(sub(".*?([0-9.]+)%.*", "\\1", grep("overall alignment rate", hisat, value = TRUE, ignore.case = TRUE),
                          perl = TRUE))   # "98.8% overall alignment rate" or "Overall alignment rate: 98.8%"
  assigned <- fc_summary[fc_summary$Status == "Assigned", grep(paste0("/", s, "\\.bam$"), names(fc_summary))]
  data.frame(sample = s, sequenced = as.integer(input), after_trimming = as.integer(kept),
             aligned = as.integer(round(kept * rate / 100)), assigned_to_gene = as.integer(assigned))
}))
write_table(funnel, out, "01_read_funnel.tsv")

funnel_long <- funnel |>
  pivot_longer(-sample, names_to = "step", values_to = "reads") |>
  mutate(step = factor(gsub("_", " ", step), levels = gsub("_", " ", names(funnel)[-1])),
         percent = 100 * reads / rep(funnel$sequenced, each = 4))
funnel_mid <- funnel_long |> group_by(step) |> summarize(percent = median(percent), .groups = "drop")
p <- ggplot(funnel_long, aes(step, percent)) +
  geom_col(data = funnel_mid, fill = COL_REF, width = 0.6) +
  geom_point(position = position_nudge(x = 0.4), size = 1.8, color = COL_MUTED) +
  geom_text(data = funnel_mid, aes(label = sprintf("%.0f%%", percent)), vjust = -0.5, size = 4) +
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(title = "Where do the reads go?",
       subtitle = sprintf("Median over %d samples (gray dots: each sample); %s reads sequenced per sample",
                          nrow(funnel), format(median(funnel$sequenced), big.mark = ",")),
       x = NULL, y = "% of sequenced reads")
save_plot(p, out, "01_read_funnel.png", width = 7, height = 4.5)

# --- 2. Coverage of one gene -------------------------------------------------------------
gtf <- file.path("data/reference", basename(cfg$GTF_URL))
bams <- file.path("data/bam", paste0(samples$sample, ".bam"))
gene <- cfg$SHOWCASE_GENE
if (file.exists(gtf) && all(file.exists(bams)) && nzchar(Sys.which("samtools"))) {
  # Exons of the gene's canonical transcript (Ensembl tag), read straight from the GTF
  cmd <- sprintf("gzip -dc %s | awk -F'\\t' '$3 == \"exon\" && $9 ~ /gene_name \"%s\";/ && $9 ~ /Ensembl_canonical/'",
                 shQuote(gtf), gene)
  exons <- read.delim(pipe(cmd), header = FALSE)[, c(1, 4, 5)]
  names(exons) <- c("chr", "start", "end")
  if (nrow(exons) == 0) stop("Gene ", gene, " (SHOWCASE_GENE) not found in ", gtf)
  region <- sprintf("%s:%d-%d", exons$chr[1], min(exons$start) - 500, max(exons$end) + 500)

  depth <- read.delim(pipe(paste("samtools depth -a -r", region, paste(shQuote(bams), collapse = " "))),
                      header = FALSE)
  names(depth) <- c("chr", "pos", samples$sample)
  # CPM: scale each sample to reads per million assigned reads, so deeper samples do not look "higher"
  scale <- funnel$assigned_to_gene / 1e6
  bin_width <- 50
  cov <- depth |>
    pivot_longer(all_of(samples$sample), names_to = "sample", values_to = "depth") |>
    mutate(depth = depth / scale[match(sample, funnel$sample)],
           bin = floor(pos / bin_width) * bin_width,
           condition = samples$condition[match(sample, samples$sample)]) |>
    group_by(condition, bin) |>
    summarize(depth = mean(depth), .groups = "drop")

  cols <- condition_colors(samples)
  p_cov <- ggplot(cov, aes(bin / 1e6, depth, fill = condition)) +
    geom_area(color = NA) +
    facet_grid(condition ~ .) +
    scale_fill_manual(values = cols, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = paste0(gene, ": the reads themselves already show the effect"),
         subtitle = "Reads covering each position, in CPM (per million reads assigned to genes); mean of each condition",
         x = NULL, y = "read depth (CPM)") +
    theme(strip.text.y = element_text(angle = 0), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p_gene <- ggplot(exons) +
    annotate("segment", x = min(exons$start) / 1e6, xend = max(exons$end) / 1e6, y = 0, yend = 0, color = COL_MUTED) +
    geom_rect(aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -1, ymax = 1), fill = "#3b3a36") +
    labs(x = paste0("position on chromosome ", exons$chr[1], " (Mb)"), y = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank())
  save_plot(p_cov / p_gene + plot_layout(heights = c(5, 1)), out, paste0("02_coverage_", gene, ".png"),
            width = 8, height = 4.5)
} else {
  message("  BAM files or GTF not found; skipping the coverage plot.")
}

# --- 3. featureCounts vs Salmon ----------------------------------------------------------
fc <- read.delim(fc_file, row.names = 1, check.names = FALSE)
salmon_files <- file.path(UP, "salmon", samples$sample, "quant.sf")
if (all(file.exists(salmon_files))) {
  suppressPackageStartupMessages(library(tximport))
  tx2gene <- read.delim("data/reference/tx2gene.tsv", header = FALSE)
  txi <- suppressMessages(tximport(setNames(salmon_files, samples$sample), type = "salmon", tx2gene = tx2gene))
  salmon_counts <- round(txi$counts)
  write.table(data.frame(gene_id = rownames(salmon_counts), salmon_counts, check.names = FALSE),
              file.path(UP, "salmon", "gene_counts_salmon.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  shared <- intersect(rownames(fc), rownames(salmon_counts))
  both <- data.frame(featureCounts = rowSums(fc[shared, ]), Salmon = rowSums(salmon_counts[shared, ])) |>
    filter(featureCounts + Salmon > 0)
  r <- cor(log10(both$featureCounts + 1), log10(both$Salmon + 1))
  strand <- vapply(samples$sample, function(s) {
    j <- readLines(file.path(UP, "salmon", s, "lib_format_counts.json"), warn = FALSE)
    sub('.*"expected_format": "([A-Z]+)".*', "\\1", grep('"expected_format"', j, value = TRUE))
  }, "")
  p <- ggplot(both, aes(featureCounts + 1, Salmon + 1)) +
    geom_abline(color = COL_MUTED, linetype = 2) +
    geom_point(size = 0.5, alpha = 0.25, color = COL_REF) +
    scale_x_log10(labels = scales::label_comma()) + scale_y_log10(labels = scales::label_comma()) +
    annotate("text", x = 1.5, y = max(both$Salmon) * 0.6, hjust = 0, size = 4.2,
             label = sprintf("r = %.3f (log scale)\n%s genes", r, format(nrow(both), big.mark = ","))) +
    labs(title = "Two routes, one answer: alignment + counting vs. Salmon",
         subtitle = sprintf("Reads per gene, summed over samples. Salmon detected library type: %s",
                            paste(unique(strand), collapse = ", ")),
         x = "featureCounts (HISAT2 alignment) + 1", y = "Salmon (tximport, gene level) + 1")
  save_plot(p, out, "03_featurecounts_vs_salmon.png", width = 6.5, height = 5.5)
} else {
  message("  Salmon results not found; skipping the featureCounts vs Salmon comparison.")
}

# --- 4. Our subset vs the published full-depth counts (airway study only) -----------------
suppressPackageStartupMessages(library(SummarizedExperiment))
data("airway", package = "airway", envir = environment())
if (all(samples$sample %in% colnames(airway))) {
  published <- assay(airway)[, samples$sample]
  shared <- intersect(rownames(fc), rownames(published))
  both <- data.frame(pipeline = rowSums(fc[shared, ]), published = rowSums(published[shared, ])) |>
    filter(pipeline + published > 0)
  r <- cor(log10(both$pipeline + 1), log10(both$published + 1))
  depth_ratio <- sum(published) / sum(fc[shared, ])
  p <- ggplot(both, aes(published + 1, pipeline + 1)) +
    geom_point(size = 0.5, alpha = 0.25, color = COL_REF) +
    scale_x_log10(labels = scales::label_comma()) + scale_y_log10(labels = scales::label_comma()) +
    annotate("text", x = 1.5, y = max(both$pipeline) * 0.6, hjust = 0, size = 4.2,
             label = sprintf("r = %.3f (log scale)\npublished data is %.0fx deeper", r, depth_ratio)) +
    labs(title = "Our subset vs. the published counts",
         subtitle = "Reads per gene, summed over samples; published = full runs (airway Bioconductor package)",
         x = "published counts + 1", y = "this pipeline + 1")
  save_plot(p, out, "04_pipeline_vs_published.png", width = 6.5, height = 5.5)
}
