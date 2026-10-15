# 4_go: Gene Ontology over-representation analysis (ORA).
# Question: among the significant genes, are some GO terms more frequent than expected by chance,
# given all the genes we could have detected (the "universe")? Fisher's / hypergeometric test.
#   01_dotplot_all_ontologies   top terms in Biological Process, Molecular Function, Cellular Component
#   02_dotplot_bp               Biological Process in detail
#   03_barplot_bp               the same as bars (number of genes per term)
#   04_cnetplot_bp              which genes drive the top terms, colored by fold change
# Output for later steps: results/downstream/objects/ego_bp.rds
source("downstream/_common.R")
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
})
step_header("4_go: GO over-representation")
out <- out_dir("4_go")
g <- gene_lists()
message(sprintf("  %d significant genes (up %d, down %d) against a universe of %d tested genes",
                length(g$all), length(g$up), length(g$down), length(g$universe)))

ego_all <- enrichGO(gene = g$all, universe = g$universe, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                    ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
save_enrichment(ego_all, out, "go_ora_all_ontologies.tsv")
if (has_terms(ego_all)) {
  p <- enrich_dotplot(ego_all, "GO over-representation: significant genes", n = 6, split = "ONTOLOGY") +
    facet_grid(ONTOLOGY ~ ., scales = "free", space = "free")
  save_plot(p, out, "01_dotplot_all_ontologies.png", width = 8, height = 8)
}

ego_bp <- enrichGO(gene = g$all, universe = g$universe, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                   ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
save_object(ego_bp, "ego_bp")
if (save_enrichment(ego_bp, out, "go_ora_bp.tsv")) {
  save_plot(enrich_dotplot(ego_bp, "GO Biological Process: over-represented terms"),
            out, "02_dotplot_bp.png", width = 8, height = 6.5)
  p <- barplot(ego_bp, showCategory = 15, label_format = 45) +
    labs(title = "GO Biological Process: genes per term", x = "significant genes in the term")
  save_plot(p, out, "03_barplot_bp.png", width = 8, height = 6.5)
  save_plot(enrich_cnetplot(ego_bp, g$fold_change_symbol, "Top GO terms and the genes behind them"),
            out, "04_cnetplot_bp.png", width = 10, height = 8)
}
