# Shared settings and helpers for the downstream R scripts.
# Every script starts with  source("downstream/_common.R")  and is run from the repository
# folder, either by run_all.sh or by hand (Rscript downstream/04_go.R, or line by line in RStudio).

if (!file.exists("config/config.sh")) {
  stop("Run this from the repository folder (the one that contains config/config.sh).", call. = FALSE)
}

suppressPackageStartupMessages({
  library(ggplot2)
})

# --- Settings: read config/config.sh, so shell and R steps share one set of values -------
# Lines look like  NAME=${NAME:-value}; an environment variable of the same name wins.
read_config <- function(path = "config/config.sh") {
  lines <- readLines(path)
  m <- regmatches(lines, regexec("^([A-Z0-9_]+)=\\$\\{\\1:-([^}]*)\\}", lines, perl = TRUE))
  m <- m[lengths(m) == 3]
  cfg <- setNames(vapply(m, `[`, "", 3), vapply(m, `[`, "", 2))
  env <- Sys.getenv(names(cfg), unset = "")
  cfg[nzchar(env)] <- env[nzchar(env)]
  as.list(cfg)
}
cfg <- read_config()
PADJ_CUTOFF <- as.numeric(cfg$PADJ_CUTOFF)
LFC_CUTOFF  <- as.numeric(cfg$LFC_CUTOFF)
SEED        <- as.integer(cfg$SEED)

# --- Folders --------------------------------------------------------------------------
UP   <- "results/upstream"
DOWN <- "results/downstream"
RDS  <- file.path(DOWN, "objects")   # R objects handed from one step to the next

# Each step owns one folder and empties it first, so no file from an earlier run survives.
out_dir <- function(name) {
  d <- file.path(DOWN, name)
  unlink(d, recursive = TRUE)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# need_file(path, step): stop with a pointer to the step that makes the file
need_file <- function(path, step) {
  if (!file.exists(path)) {
    stop("Missing ", path, ". Run ", step, " first.", call. = FALSE)
  }
  path
}

load_object <- function(name, step) readRDS(need_file(file.path(RDS, paste0(name, ".rds")), step))
save_object <- function(x, name) {
  dir.create(RDS, recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, file.path(RDS, paste0(name, ".rds")))
}

# --- Samples --------------------------------------------------------------------------
read_samples <- function(path = cfg$SAMPLES) {
  s <- read.delim(need_file(path, "nothing (check SAMPLES in config/config.sh)"),
                  comment.char = "#", colClasses = "character", check.names = FALSE)
  expected <- c("sample", "condition", "batch", "fastq_1", "fastq_2")
  if (!identical(names(s), expected)) {
    stop(path, ": columns must be ", paste(expected, collapse = ", "), call. = FALSE)
  }
  levels <- unique(s$condition)
  if (length(levels) != 2 || !cfg$REFERENCE_LEVEL %in% levels) {
    stop(path, ": needs exactly two conditions, one of them REFERENCE_LEVEL ('",
         cfg$REFERENCE_LEVEL, "'). Found: ", paste(levels, collapse = ", "), call. = FALSE)
  }
  s$condition <- factor(s$condition, levels = c(cfg$REFERENCE_LEVEL, setdiff(levels, cfg$REFERENCE_LEVEL)))
  s$batch <- factor(s$batch)
  rownames(s) <- s$sample
  s[, c("sample", "condition", "batch")]
}

has_batch <- function(samples) nlevels(samples$batch) > 1

# --- Looks: one palette for every figure ----------------------------------------------
# Condition colors: blue = reference (untreated), orange = treated. Genes that go up with the
# treatment are drawn orange, genes that go down blue, everything else gray.
COL_REF   <- "#2a78d6"
COL_TREAT <- "#eb6834"
COL_NS    <- "#bdbcb6"
COL_MUTED <- "#898781"
condition_colors <- function(samples) setNames(c(COL_REF, COL_TREAT), levels(samples$condition))
direction_colors <- c(up = COL_TREAT, down = COL_REF, `not significant` = COL_NS)
diverging_ramp   <- colorRampPalette(c(COL_REF, "#f0efec", "#e34948"))(101)
sequential_ramp  <- colorRampPalette(c("#0d366b", "#2a78d6", "#cde2fb", "#fcfcfb"))(100)

theme_set(theme_bw(base_size = 12) +
            theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_line(color = "#e1e0d9", linewidth = 0.3),
                  panel.border = element_rect(color = COL_MUTED, linewidth = 0.4),
                  plot.title = element_text(face = "bold"),
                  plot.title.position = "plot",
                  legend.key = element_blank()))

# Enrichment plots (clusterProfiler/enrichplot) color by adjusted p-value; use one blue ramp,
# dark = most significant.
options(enrichplot.colours = c("#0d366b", "#9ec5f4"))

# Standard dot plot for an enrichment result, with wrapped term names
enrich_dotplot <- function(x, title, n = 15, ...) {
  enrichplot::dotplot(x, showCategory = n, label_format = 45, ...) +
    labs(title = title) + theme(panel.grid.major.y = element_blank())
}

# Gene-concept network: terms and their genes; genes colored by fold change (blue down, orange up)
enrich_cnetplot <- function(x, fold_change, title, n = 4) {
  suppressMessages(
    enrichplot::cnetplot(x, showCategory = n, foldChange = fold_change, node_label = "all") +
      scale_color_gradient2(low = COL_REF, mid = "#d9d8d2", high = COL_TREAT, midpoint = 0,
                            name = "log2 fold change") +
      labs(title = title, subtitle = "genes colored by fold change: orange = up, blue = down") +
      theme(plot.margin = margin(8, 8, 8, 16))
  )
}

# --- Saving ---------------------------------------------------------------------------
save_plot <- function(plot, dir, file, width = 7, height = 5) {
  ggsave(file.path(dir, file), plot, width = width, height = height, dpi = 150, bg = "white")
  message("  wrote ", file.path(dir, file))
}

# For base-graphics or grid plots (pheatmap, DESeq2::plotMA): save_png(dir, file, { code })
save_png <- function(dir, file, expr, width = 7, height = 5) {
  png(file.path(dir, file), width = width, height = height, units = "in", res = 150, bg = "white")
  on.exit(dev.off())
  force(expr)
  message("  wrote ", file.path(dir, file))
}

write_table <- function(df, dir, file) {
  write.table(df, file.path(dir, file), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  message("  wrote ", file.path(dir, file))
}

# Enrichment results: save the table, and say plainly when nothing is enriched, because an empty
# result is a valid answer and the plotting functions fail on it.
has_terms <- function(x) !is.null(x) && nrow(as.data.frame(x)) > 0
save_enrichment <- function(x, dir, file) {
  if (has_terms(x)) {
    write_table(as.data.frame(x), dir, file)
  } else {
    writeLines("No term passed the significance cutoff.", file.path(dir, sub("\\.tsv$", "_NONE.txt", file)))
    message("  no significant terms for ", file)
  }
  invisible(has_terms(x))
}

step_header <- function(title) message("\n== ", title, " ==")

# Term maps and trees (emapplot, treeplot) need enough terms to cluster and can fail on small
# results. Draw them when possible; otherwise leave a note saying why, and carry on.
#   save_term_plot(x, n_min, dir, file, width, height, plot_fun)
save_term_plot <- function(x, n_min, dir, file, width, height, plot_fun) {
  n <- nrow(as.data.frame(x))
  note <- file.path(dir, sub("\\.png$", "_NOT_DRAWN.txt", file))
  if (n < n_min) {
    writeLines(sprintf("Not drawn: %d significant terms, this plot needs at least %d.", n, n_min), note)
    message("  ", file, ": only ", n, " terms, not drawn")
    return(invisible(FALSE))
  }
  ok <- tryCatch({ save_plot(plot_fun(x), dir, file, width, height); TRUE },
                 error = function(e) {
                   writeLines(paste("Not drawn:", conditionMessage(e)), note)
                   message("  ", file, " not drawn: ", conditionMessage(e))
                   FALSE
                 })
  invisible(ok)
}

# --- Gene lists for the enrichment steps (4-12) ---------------------------------------------
# Built from the DESeq2 results (step 2). Enrichment tools use Entrez gene IDs.
#   up, down, all   significant genes (padj < PADJ_CUTOFF and |shrunken log2FC| >= LFC_CUTOFF)
#   universe        every gene that was tested: the background for over-representation tests
#   ranked          all tested genes ranked by the DESeq2 Wald statistic, for GSEA
#   fold_change     shrunken log2FC of the significant genes, to color network plots
#   fold_change_all shrunken log2FC of every tested gene
gene_lists <- function() {
  de <- load_object("de_results", "downstream/02_dea.R")
  de <- de[!is.na(de$entrez) & !is.na(de$padj), ]
  de <- de[!duplicated(de$entrez), ]          # one row per Entrez ID (sorted by padj, best kept)
  sig <- de[de$direction != "not significant", ]
  ranked <- sort(setNames(de$stat, de$entrez), decreasing = TRUE)
  list(up = sig$entrez[sig$direction == "up"],
       down = sig$entrez[sig$direction == "down"],
       all = sig$entrez,
       universe = de$entrez,
       ranked = ranked,
       fold_change = setNames(sig$log2FC_shrunk, sig$entrez),
       fold_change_symbol = setNames(sig$log2FC_shrunk, sig$symbol),
       fold_change_all = setNames(de$log2FC_shrunk, de$entrez))
}

# Wrap online lookups (KEGG, WikiPathways, Europe PMC): retry once, then fail with a clear message.
with_retry <- function(expr, what) {
  call <- substitute(expr)
  env <- parent.frame()
  for (attempt in 1:2) {
    result <- tryCatch(eval(call, env), error = function(e) e)
    if (!inherits(result, "error")) return(result)
    message("  ", what, " failed (attempt ", attempt, "): ", conditionMessage(result))
    if (attempt == 1) Sys.sleep(5)
  }
  stop(what, " failed twice. It needs internet access; check the connection and rerun this step.", call. = FALSE)
}
