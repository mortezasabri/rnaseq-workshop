# 5_go_ora_clust: making a long GO list readable. GO terms overlap heavily ("blood circulation"
# and "circulatory system process" share most genes), so group similar terms.
#   01_compare_up_down     ORA separately for up- and down-regulated genes, side by side
#   02_dotplot_simplified  redundant terms removed (semantic similarity > 0.7)
#   03_emapplot            enrichment map: terms as nodes, linked when they share genes
#   04_treeplot            terms clustered into themes, each theme named by its frequent words
source("downstream/_common.R")
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
})
step_header("5_go_ora_clust: grouping GO over-representation results")
out <- out_dir("5_go_ora_clust")
g <- gene_lists()
ego_bp <- load_object("ego_bp", "downstream/04_go.R")

cc <- compareCluster(geneClusters = list(up = g$up, down = g$down), fun = "enrichGO",
                     universe = g$universe, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                     pvalueCutoff = 0.05, readable = TRUE)
if (save_enrichment(cc, out, "go_ora_up_vs_down.tsv")) {
  save_plot(enrich_dotplot(cc, "GO Biological Process: up- vs. down-regulated genes", n = 8),
            out, "01_compare_up_down.png", width = 8, height = 7.5)
}

if (has_terms(ego_bp)) {
  simple <- simplify(ego_bp, cutoff = 0.7, by = "p.adjust", select_fun = min)
  save_enrichment(simple, out, "go_ora_bp_simplified.tsv")
  message(sprintf("  simplify: %d terms -> %d", nrow(as.data.frame(ego_bp)), nrow(as.data.frame(simple))))
  save_plot(enrich_dotplot(simple, "GO Biological Process, redundant terms removed"),
            out, "02_dotplot_simplified.png", width = 8, height = 6.5)

  sim <- pairwise_termsim(ego_bp)
  save_term_plot(sim, 3, out, "03_emapplot.png", 10, 8, function(x)
    emapplot(x, showCategory = 30) + labs(title = "Enrichment map: terms sharing genes are linked"))
  save_term_plot(sim, 10, out, "04_treeplot.png", 13, 9, function(x)
    treeplot(x, showCategory = 20, label_format = 28) + labs(title = "Terms clustered into themes"))
}
