# 12_ReactomePA: Reactome, an expert-curated database of reactions and pathways (ReactomePA).
#   01_ora_dotplot     over-represented pathways
#   02_ora_emapplot    pathways that share genes, linked
#   03_gsea_dotplot    pathways shifted up or down (GSEA)
#   04_gseaplot        running score of the top pathways
#   05_pathway_network the genes of the top pathway as a network, colored by fold change
source("downstream/_common.R")
suppressPackageStartupMessages({
  library(ReactomePA)
  library(clusterProfiler)
  library(enrichplot)
})
step_header("12_ReactomePA: Reactome pathways")
out <- out_dir("12_ReactomePA")
g <- gene_lists()

rp <- enrichPathway(gene = g$all, universe = g$universe, organism = "human", pvalueCutoff = 0.05, readable = TRUE)
if (save_enrichment(rp, out, "reactome_ora.tsv")) {
  save_plot(enrich_dotplot(rp, "Reactome: over-represented among significant genes"),
            out, "01_ora_dotplot.png", width = 8, height = 6.5)
  save_term_plot(pairwise_termsim(rp), 3, out, "02_ora_emapplot.png", 10, 8, function(x)
    emapplot(x, showCategory = 25) + labs(title = "Reactome pathways sharing genes"))
}

set.seed(SEED)
grp <- gsePathway(geneList = g$ranked, organism = "human", minGSSize = 15, maxGSSize = 500,
                  pvalueCutoff = 0.05, seed = TRUE, verbose = FALSE)
if (has_terms(grp)) grp <- setReadable(grp, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
if (save_enrichment(grp, out, "reactome_gsea.tsv")) {
  save_plot(enrich_dotplot(grp, "Reactome: GSEA", n = 10, split = ".sign") + facet_grid(. ~ .sign),
            out, "03_gsea_dotplot.png", width = 10, height = 6.5)
  top <- head(as.data.frame(grp)$ID, 3)
  save_plot(gseaplot2(grp, geneSetID = top, pvalue_table = FALSE, title = "Reactome GSEA: top pathways",
                      color = c(COL_TREAT, "#eda100", COL_REF)[seq_along(top)]),
            out, "04_gseaplot.png", width = 9, height = 6.5)
}

# Network of the most significant pathway small enough to read (<= 80 genes in our universe)
if (has_terms(rp)) {
  res <- as.data.frame(rp)
  size <- as.numeric(sub("/.*", "", res$BgRatio))
  if (any(size <= 80)) {
    top_pathway <- res$Description[size <= 80][1]
    message("  network of: ", top_pathway)
    p <- suppressMessages(viewPathway(top_pathway, organism = "human", readable = TRUE, foldChange = g$fold_change_all) +
      scale_color_gradient2(low = COL_REF, mid = "#d9d8d2", high = COL_TREAT, midpoint = 0, name = "log2 fold change") +
      labs(title = paste("Reactome:", top_pathway),
           subtitle = "every gene colored by fold change: orange = up, blue = down, dark gray = not detected"))
    save_plot(p, out, "05_pathway_network.png", width = 10, height = 8)
  }
}
