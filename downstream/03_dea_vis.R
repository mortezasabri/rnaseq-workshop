# 3_dea_vis: looking at the differential expression results.
#   01_ma_plot            fold change vs. expression level: every gene is a dot
#   02_volcano            fold change vs. significance, top genes labeled
#   03_heatmap_top30      the 30 most significant genes across all samples
#   04_top_genes          normalized counts of the top genes, one line per donor
#   05_pvalue_histogram   a quick health check of the test
#   06_ma_plain_vs_shrunken  the MA plot before and after fold-change shrinkage
source("downstream/_common.R")
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
  library(pheatmap)
})
step_header("3_dea_vis: visualizing the results")
out <- out_dir("3_dea_vis")
de  <- load_object("de_results", "downstream/02_dea.R")
dds <- load_object("dds_fitted", "downstream/02_dea.R")
vsd <- load_object("vsd", "downstream/01_eda.R")
samples <- as.data.frame(colData(dds))[, c("sample", "condition", "batch")]
cols <- condition_colors(samples)
trt <- levels(dds$condition)[2]
de$label <- ifelse(is.na(de$symbol), de$gene_id, de$symbol)
de$direction <- factor(de$direction, c("up", "down", "not significant"))
n_up <- sum(de$direction == "up"); n_down <- sum(de$direction == "down")
direction_labels <- c(up = paste0("up (", n_up, ")"), down = paste0("down (", n_down, ")"),
                      `not significant` = "not significant")
cutoff_text <- sprintf("colored: padj < %g and |log2FC| >= %g", PADJ_CUTOFF, LFC_CUTOFF)

# --- 1. MA plot -------------------------------------------------------------------------------
ma <- filter(de, baseMean > 0) |> arrange(desc(direction))   # significant genes drawn on top
p <- ggplot(ma, aes(baseMean, log2FC_shrunk, color = direction)) +
  geom_hline(yintercept = 0, color = COL_MUTED) +
  geom_point(size = 0.7, alpha = 0.6) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_color_manual(values = direction_colors, labels = direction_labels, name = NULL) +
  labs(title = paste("MA plot:", trt, "vs.", levels(dds$condition)[1]),
       subtitle = paste("Shrunken log2 fold changes;", cutoff_text),
       x = "mean normalized count", y = "log2 fold change") +
  theme(legend.position = "top")
save_plot(p, out, "01_ma_plot.png", width = 7, height = 5)

# --- 2. Volcano plot ---------------------------------------------------------------------------
vol <- filter(de, !is.na(padj)) |> mutate(neglog10 = -log10(padj)) |> arrange(desc(direction))
top <- vol |> filter(direction != "not significant") |> arrange(padj) |> head(15)
p <- ggplot(vol, aes(log2FC_shrunk, neglog10, color = direction)) +
  geom_point(size = 0.7, alpha = 0.6) +
  geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = 2, color = COL_MUTED) +
  geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = 2, color = COL_MUTED) +
  geom_text_repel(data = top, aes(label = label), color = "#0b0b0b", size = 3.3,
                  max.overlaps = Inf, box.padding = 0.4, min.segment.length = 0, seed = SEED) +
  scale_color_manual(values = direction_colors, labels = direction_labels, name = NULL) +
  labs(title = "Volcano plot", subtitle = paste("15 most significant genes labeled;", cutoff_text),
       x = "log2 fold change (shrunken)", y = "-log10 adjusted p-value") +
  theme(legend.position = "top")
save_plot(p, out, "02_volcano.png", width = 7, height = 5.5)

# --- 3. Heatmap of the top 30 genes ------------------------------------------------------------
top30 <- de |> filter(!is.na(padj)) |> arrange(padj) |> head(30)
m <- assay(vsd)[top30$gene_id, ]
m <- m - rowMeans(m)                         # each gene relative to its own average
rownames(m) <- top30$label
colnames(m) <- vsd$sample
ann <- data.frame(condition = vsd$condition, batch = vsd$batch, row.names = vsd$sample)
batch_cols <- setNames(c("#3b3a36", "#898781", "#c3c2b7", "#e1e0d9", "#52514e", "#a8a69d")[seq_len(nlevels(vsd$batch))],
                       levels(vsd$batch))
lim <- max(abs(m))
save_png(out, "03_heatmap_top30.png", width = 7, height = 8,
  pheatmap(m, color = diverging_ramp, breaks = seq(-lim, lim, length.out = 102), border_color = NA,
           annotation_col = ann, annotation_colors = list(condition = cols, batch = batch_cols),
           main = "Top 30 genes (variance-stabilized, centered per gene)", fontsize_row = 9, silent = FALSE))

# --- 4. Counts of the top genes, paired by donor ----------------------------------------------
genes <- unique(c(head(top30$gene_id, 5), de$gene_id[which(de$symbol == cfg$SHOWCASE_GENE)]))
norm <- counts(dds, normalized = TRUE)[genes, , drop = FALSE]
long <- data.frame(gene = de$label[match(genes, de$gene_id)], norm, check.names = FALSE) |>
  pivot_longer(-gene, names_to = "sample", values_to = "count") |>
  left_join(samples, by = "sample") |>
  mutate(gene = factor(gene, unique(de$label[match(genes, de$gene_id)])))
p <- ggplot(long, aes(condition, count + 1)) +
  geom_line(aes(group = batch), color = COL_MUTED, linewidth = 0.4) +
  geom_point(aes(color = condition), size = 2.5) +
  facet_wrap(~gene, scales = "free_y", nrow = 2) +
  scale_y_log10(labels = scales::label_comma()) +
  scale_color_manual(values = cols, guide = "none") +
  labs(title = "Top genes, sample by sample",
       subtitle = "Normalized counts; gray lines connect the two samples of the same donor",
       x = NULL, y = "normalized count + 1 (log scale)") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p, out, "04_top_genes.png", width = 9, height = 5.5)

# --- 5. p-value histogram ----------------------------------------------------------------------
p <- ggplot(filter(de, !is.na(pvalue), baseMean > 1), aes(pvalue)) +
  geom_histogram(binwidth = 0.025, boundary = 0, fill = COL_REF, color = "white", linewidth = 0.2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "p-value histogram",
       subtitle = "Flat = genes with no change; the peak near 0 = real differences",
       x = "p-value", y = "genes")
save_plot(p, out, "05_pvalue_histogram.png", width = 6, height = 4)

# --- 6. Why shrink? The same MA plot with plain and with shrunken fold changes ---------------
# Plain fold changes of genes with few reads spread wide by chance; apeglm pulls them toward 0
# unless the data support them. p-values are not affected, only the fold changes.
moved <- ma |> filter(!is.na(symbol)) |> mutate(change = abs(log2FC - log2FC_shrunk)) |>
  arrange(desc(change)) |> head(3)
shrink_long <- bind_rows(
  mutate(ma, type = "plain fold change", lfc = log2FC),
  mutate(ma, type = "shrunken fold change (apeglm)", lfc = log2FC_shrunk)
) |> mutate(type = factor(type, c("plain fold change", "shrunken fold change (apeglm)")))
moved_long <- filter(shrink_long, gene_id %in% moved$gene_id)
p <- ggplot(shrink_long, aes(baseMean, lfc, color = direction)) +
  geom_hline(yintercept = 0, color = COL_MUTED) +
  geom_point(size = 0.7, alpha = 0.6) +
  geom_point(data = moved_long, shape = 21, size = 3, stroke = 0.8, color = "#0b0b0b", fill = NA) +
  geom_text_repel(data = moved_long, aes(label = label), color = "#0b0b0b", size = 3.2,
                  box.padding = 0.5, min.segment.length = 0, seed = SEED) +
  facet_wrap(~type) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_color_manual(values = direction_colors, labels = direction_labels, name = NULL) +
  labs(title = "Why shrink fold changes?",
       subtitle = "Same genes, same p-values; circled: the 3 genes that shrinkage moves most",
       x = "mean normalized count", y = "log2 fold change") +
  theme(legend.position = "top")
save_plot(p, out, "06_ma_plain_vs_shrunken.png", width = 10, height = 5)
