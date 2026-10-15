# 2_dea: differential expression analysis with DESeq2.
# For every gene: how much does it change with the treatment (log2 fold change), and how sure
# are we (adjusted p-value)?
#   deseq2_all_genes.tsv      every tested gene
#   deseq2_significant.tsv    padj < PADJ_CUTOFF and |log2FC| >= LFC_CUTOFF, sorted by padj
#   summary.txt               how many genes go up / down
#   01_deseq2_vs_edger        a second method (edgeR) on the same counts: do they agree?
#   02_depth_check            (pipeline counts only) which of the full-depth data's genes the subset finds
# Output for the next steps: results/downstream/objects/dds_fitted.rds, de_results.rds
source("downstream/_common.R")
suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(dplyr)
})
step_header("2_dea: differential expression (DESeq2)")
out <- out_dir("2_dea")
dds <- load_object("dds", "downstream/01_eda.R")
ref <- levels(dds$condition)[1]
trt <- levels(dds$condition)[2]

# --- DESeq2: normalization, dispersion estimation, negative binomial GLM, Wald test ----------
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("condition", trt, ref), alpha = PADJ_CUTOFF)
# Shrunken fold changes: genes with few reads get noisy, extreme fold changes; apeglm pulls
# those toward zero. Use them for ranking and plotting; the p-values stay the same.
coef <- paste0("condition_", make.names(trt), "_vs_", make.names(ref))
shrunk <- lfcShrink(dds, coef = coef, type = "apeglm", quiet = TRUE)
save_object(dds, "dds_fitted")

de <- data.frame(
  gene_id = rownames(res),
  symbol = mcols(dds)$symbol,
  entrez = mcols(dds)$entrez,
  baseMean = res$baseMean,
  log2FC = res$log2FoldChange,
  log2FC_shrunk = shrunk$log2FoldChange,
  stat = res$stat,
  pvalue = res$pvalue,
  padj = res$padj
) |>
  mutate(direction = case_when(
    !is.na(padj) & padj < PADJ_CUTOFF & log2FC_shrunk >= LFC_CUTOFF ~ "up",
    !is.na(padj) & padj < PADJ_CUTOFF & log2FC_shrunk <= -LFC_CUTOFF ~ "down",
    TRUE ~ "not significant")) |>
  arrange(padj)
save_object(de, "de_results")   # full precision: GSEA ranks on stat, rounding would create ties
rounded <- mutate(de, across(c(baseMean), ~ round(.x, 1)), across(c(log2FC, log2FC_shrunk, stat), ~ round(.x, 3)),
                  across(c(pvalue, padj), ~ signif(.x, 3)))
write_table(rounded, out, "deseq2_all_genes.tsv")
write_table(filter(rounded, direction != "not significant"), out, "deseq2_significant.tsv")

n_sig <- sum(de$padj < PADJ_CUTOFF, na.rm = TRUE)
summary_lines <- c(
  sprintf("Comparison: %s vs. %s (design: %s)", trt, ref, paste(deparse(design(dds)), collapse = "")),
  sprintf("Genes tested: %d", nrow(de)),
  sprintf("Genes with padj < %g: %d (up %d, down %d)", PADJ_CUTOFF, n_sig,
          sum(de$padj < PADJ_CUTOFF & de$log2FC > 0, na.rm = TRUE),
          sum(de$padj < PADJ_CUTOFF & de$log2FC < 0, na.rm = TRUE)),
  sprintf("... and |shrunken log2FC| >= %g: up %d, down %d", LFC_CUTOFF,
          sum(de$direction == "up"), sum(de$direction == "down")),
  sprintf("Genes without padj (too few reads or outlier sample): %d", sum(is.na(de$padj))),
  "",
  "Top 10 genes by adjusted p-value:",
  capture.output(print(head(rounded[, c("symbol", "baseMean", "log2FC_shrunk", "padj")], 10), row.names = FALSE))
)
writeLines(summary_lines, file.path(out, "summary.txt"))
message(paste0("  ", summary_lines[1:4], collapse = "\n"))

# --- edgeR on the same counts (quasi-likelihood F-test) ------------------------------------
y <- DGEList(counts = counts(dds), samples = as.data.frame(colData(dds)))
y <- normLibSizes(y)
design_matrix <- model.matrix(design(dds), data = y$samples)
y <- estimateDisp(y, design_matrix)
fit <- glmQLFit(y, design_matrix)
qlf <- glmQLFTest(fit, coef = ncol(design_matrix))
edger <- topTags(qlf, n = Inf, sort.by = "none")$table
both <- data.frame(gene_id = rownames(edger), edgeR_logFC = edger$logFC, edgeR_FDR = edger$FDR) |>
  inner_join(de, by = "gene_id") |>
  mutate(called_by = case_when(
    padj < PADJ_CUTOFF & edgeR_FDR < PADJ_CUTOFF ~ "both",
    !is.na(padj) & padj < PADJ_CUTOFF ~ "DESeq2 only",
    edgeR_FDR < PADJ_CUTOFF ~ "edgeR only",
    TRUE ~ "neither"))
categories <- c("both", "DESeq2 only", "edgeR only", "neither")
overlap <- count(both, called_by = factor(called_by, categories), name = "genes", .drop = FALSE)
write_table(overlap, out, "01_deseq2_vs_edger_overlap.tsv")
n <- setNames(overlap$genes, as.character(overlap$called_by))
lfc_r <- cor(both$log2FC, both$edgeR_logFC, use = "complete.obs")
both <- mutate(both, shown = case_when(called_by == "both" ~ "significant in both",
                                       called_by == "neither" ~ "neither",
                                       TRUE ~ "significant in one method only")) |>
  arrange(match(shown, c("neither", "significant in one method only", "significant in both")))
p <- ggplot(both, aes(log2FC, edgeR_logFC, color = shown)) +
  geom_abline(color = COL_MUTED, linetype = 2) +
  geom_point(size = 0.6, alpha = 0.5) +
  scale_color_manual(values = c(`significant in both` = COL_TREAT, `significant in one method only` = "#3b3a36",
                                neither = COL_NS), breaks = c("significant in both", "significant in one method only", "neither"),
                     name = NULL) +
  labs(title = "DESeq2 and edgeR: same fold changes",
       subtitle = sprintf("r = %.3f. Significant (FDR < %g): both %d, DESeq2 only %d, edgeR only %d",
                          lfc_r, PADJ_CUTOFF, n[["both"]], n[["DESeq2 only"]], n[["edgeR only"]]),
       x = "DESeq2 log2 fold change", y = "edgeR log2 fold change") +
  theme(legend.position = "top")
save_plot(p, out, "01_deseq2_vs_edger.png", width = 6, height = 5.5)

# --- Depth check: is a 2-million-read subset enough? (airway study, pipeline counts) ----------
suppressPackageStartupMessages(library(SummarizedExperiment))
data("airway", package = "airway", envir = environment())
if (cfg$COUNTS_SOURCE == "pipeline" && all(dds$sample %in% colnames(airway))) {
  full <- DESeqDataSetFromMatrix(assay(airway)[, dds$sample], colData = as.data.frame(colData(dds))[, c("sample", "condition", "batch")],
                                 design = design(dds))
  full <- full[rowSums(counts(full) >= 10) >= min(table(full$condition)), ]
  full <- DESeq(full, quiet = TRUE)
  full_res <- results(full, contrast = c("condition", trt, ref), alpha = PADJ_CUTOFF)
  cmp <- data.frame(gene_id = rownames(full_res), full_log2FC = full_res$log2FoldChange, full_padj = full_res$padj) |>
    inner_join(de, by = "gene_id") |>
    filter(!is.na(full_padj))
  sig_full <- sum(cmp$full_padj < PADJ_CUTOFF, na.rm = TRUE)
  sig_both <- sum(cmp$full_padj < PADJ_CUTOFF & cmp$padj < PADJ_CUTOFF, na.rm = TRUE)
  depth <- round(median(colSums(assay(airway)[, dds$sample])) / median(colSums(counts(dds))))
  writeLines(c(sprintf("Genes tested in both analyses: %d", nrow(cmp)),
               sprintf("Of these, significant (padj < %g) in the full-depth data (about %dx more reads): %d",
                       PADJ_CUTOFF, depth, sig_full),
               sprintf("Found again with this pipeline's subset: %d (%.0f%%)", sig_both, 100 * sig_both / sig_full),
               sprintf("Correlation of log2 fold changes (genes tested in both): %.3f",
                       cor(cmp$full_log2FC, cmp$log2FC, use = "complete.obs"))),
             file.path(out, "02_depth_check.txt"))
  # Which of the full data's significant genes does the subset find? Split by how much they change.
  bins <- c("less than 1.5-fold", "1.5- to 2-fold", "2- to 4-fold", "more than 4-fold")
  recovery <- cmp |>
    filter(full_padj < PADJ_CUTOFF) |>
    mutate(change = cut(2^abs(full_log2FC), c(1, 1.5, 2, 4, Inf), labels = bins, right = FALSE),
           found = !is.na(padj) & padj < PADJ_CUTOFF) |>
    group_by(change) |>
    summarize(genes = n(), found = sum(found), .groups = "drop") |>
    mutate(percent_found = 100 * found / genes)
  write_table(recovery, out, "02_depth_check.tsv")
  # percent labels sit just above the dark bar, or inside it (white) when the bar is nearly full
  labels <- mutate(recovery, inside = percent_found > 85,
                   label = sprintf("%.0f%%", percent_found),
                   counts = sprintf("%s of %s genes", scales::comma(found), scales::comma(genes)))
  p <- ggplot(recovery, aes(change, percent_found)) +
    geom_col(aes(y = 100), fill = "#e1e0d9", width = 0.65) +
    geom_col(fill = "#0d366b", width = 0.65) +
    geom_text(data = labels, aes(label = label, vjust = ifelse(inside, 1.6, -0.5), color = inside),
              size = 5, fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "#0b0b0b")) +
    geom_text(data = labels, aes(y = 100, label = counts), vjust = -0.6, size = 3.4, color = "#52514e") +
    scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 25), labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0))) +
    labs(title = "With 2 M reads we find the big changes and miss the small ones",
         subtitle = sprintf("The %s genes significant in the full data (%d times more reads), by how much they change",
                            format(sig_full, big.mark = ","), depth),
         x = "change with dexamethasone, in the full data", y = "also found with 2 M reads") +
    theme(panel.grid.major.x = element_blank())
  save_plot(p, out, "02_depth_check.png", width = 8, height = 5)
}
