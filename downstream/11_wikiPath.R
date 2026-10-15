# 11_wikiPath: WikiPathways, a community-curated, openly licensed pathway collection.
#   01_ora_dotplot    over-represented pathways
#   02_gsea_dotplot   pathways shifted up or down (GSEA)
#   03_cnetplot       genes behind the top over-represented pathways
# Needs internet: the current WikiPathways release is downloaded on each run.
source("downstream/_common.R")
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
})
step_header("11_wikiPath: WikiPathways")
out <- out_dir("11_wikiPath")
g <- gene_lists()

wp <- with_retry(enrichWP(gene = g$all, universe = g$universe, organism = "Homo sapiens", pvalueCutoff = 0.05),
                 "WikiPathways over-representation")
if (has_terms(wp)) wp <- setReadable(wp, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
if (save_enrichment(wp, out, "wikipathways_ora.tsv")) {
  save_plot(enrich_dotplot(wp, "WikiPathways: over-represented among significant genes"),
            out, "01_ora_dotplot.png", width = 8, height = 6.5)
  save_plot(enrich_cnetplot(wp, g$fold_change_symbol, "Top WikiPathways and their genes"),
            out, "03_cnetplot.png", width = 10, height = 8)
}

set.seed(SEED)
gwp <- with_retry(gseWP(geneList = g$ranked, organism = "Homo sapiens", minGSSize = 15, maxGSSize = 500,
                        pvalueCutoff = 0.05, seed = TRUE, verbose = FALSE), "WikiPathways GSEA")
if (has_terms(gwp)) gwp <- setReadable(gwp, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
if (save_enrichment(gwp, out, "wikipathways_gsea.tsv")) {
  save_plot(enrich_dotplot(gwp, "WikiPathways: GSEA", n = 10, split = ".sign") + facet_grid(. ~ .sign),
            out, "02_gsea_dotplot.png", width = 10, height = 6.5)
}
