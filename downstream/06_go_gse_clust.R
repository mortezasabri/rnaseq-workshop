# 6_go_gse_clust: Gene Set Enrichment Analysis (GSEA) on GO, plus grouping of the terms.
# ORA needs a cutoff ("significant or not"). GSEA uses every gene: rank all genes from most up to
# most down (here by the DESeq2 Wald statistic) and ask whether a term's genes sit at the top or
# bottom of the ranking. It catches terms where many genes move a little.
#   01_dotplot          top terms, split into activated (up) and suppressed (down)
#   02_gseaplot         the classic running-score plot for the top terms
#   03_ridgeplot        fold-change distribution of the core genes of each term
#   04_emapplot, 05_treeplot   terms grouped as in 5_go_ora_clust
source("downstream/_common.R")
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
})
step_header("6_go_gse_clust: GO gene set enrichment analysis")
out <- out_dir("6_go_gse_clust")
g <- gene_lists()

set.seed(SEED)
gse <- gseGO(geneList = g$ranked, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
             minGSSize = 15, maxGSSize = 500, pvalueCutoff = 0.05, seed = TRUE, verbose = FALSE)
gse <- setReadable(gse, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
save_object(gse, "gse_go")
if (save_enrichment(gse, out, "go_gsea_bp.tsv")) {
  p <- enrich_dotplot(gse, "GSEA, GO Biological Process", n = 10, split = ".sign") +
    facet_grid(. ~ .sign)
  save_plot(p, out, "01_dotplot.png", width = 10, height = 7)

  res <- as.data.frame(gse)
  top_up <- head(res$ID[res$NES > 0], 2); top_down <- head(res$ID[res$NES < 0], 1)
  p <- gseaplot2(gse, geneSetID = c(top_up, top_down), pvalue_table = FALSE,
                 title = "Running enrichment score: where do the term's genes fall in the ranking?",
                 color = c(COL_TREAT, "#eda100", COL_REF)[seq_along(c(top_up, top_down))])
  save_plot(p, out, "02_gseaplot.png", width = 9, height = 6.5)

  p <- ridgeplot(gse, showCategory = 15, label_format = 45) +
    labs(title = "Fold changes of the core genes per term", x = "Wald statistic")
  save_plot(p, out, "03_ridgeplot.png", width = 9, height = 7.5)

  sim <- pairwise_termsim(gse)
  save_term_plot(sim, 3, out, "04_emapplot.png", 10, 8, function(x)
    emapplot(x, showCategory = 30) + labs(title = "GSEA enrichment map"))
  save_term_plot(sim, 10, out, "05_treeplot.png", 13, 9, function(x)
    treeplot(x, showCategory = 20, label_format = 28) + labs(title = "GSEA terms clustered into themes"))
}
