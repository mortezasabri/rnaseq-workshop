# 9_do_gse: Disease Ontology GSEA on the full gene ranking (see 6_go_gse_clust for the idea).
#   01_dotplot     top disease terms, activated vs. suppressed
#   02_gseaplot    running score of the top terms
#   03_ridgeplot   statistic distribution of the core genes per term
source("downstream/_common.R")
suppressPackageStartupMessages({
  library(DOSE)
  library(enrichplot)
  library(org.Hs.eg.db)
})
step_header("9_do_gse: Disease Ontology GSEA")
out <- out_dir("9_do_gse")
g <- gene_lists()

set.seed(SEED)
gdo <- gseDO(geneList = g$ranked, ont = "HDO", organism = "hsa", minGSSize = 15, maxGSSize = 500,
             pvalueCutoff = 0.05, seed = TRUE, verbose = FALSE)
if (has_terms(gdo)) gdo <- setReadable(gdo, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
if (save_enrichment(gdo, out, "do_gsea.tsv")) {
  save_plot(enrich_dotplot(gdo, "Disease Ontology: GSEA", n = 10, split = ".sign") + facet_grid(. ~ .sign),
            out, "01_dotplot.png", width = 10, height = 6.5)
  top <- head(as.data.frame(gdo)$ID, 3)
  save_plot(gseaplot2(gdo, geneSetID = top, pvalue_table = FALSE, title = "Disease Ontology GSEA: top terms",
                      color = c(COL_TREAT, "#eda100", COL_REF)[seq_along(top)]),
            out, "02_gseaplot.png", width = 9, height = 6.5)
  save_plot(ridgeplot(gdo, showCategory = 15, label_format = 45) +
              labs(title = "Statistic of the core genes per disease term", x = "Wald statistic"),
            out, "03_ridgeplot.png", width = 9, height = 7.5)
}
