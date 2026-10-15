# 8_do_ora: Disease Ontology (DO) over-representation. Which diseases are the significant genes
# known to be linked to? (DOSE package; Human Disease Ontology)
#   01_dotplot     top disease terms
#   02_cnetplot    the genes behind the top terms, colored by fold change
#   03_upsetplot   how the top terms share genes
# Output for later steps: results/downstream/objects/do_ora.rds (used by 10_pubmed)
source("downstream/_common.R")
suppressPackageStartupMessages({
  library(DOSE)
  library(enrichplot)
})
step_header("8_do_ora: Disease Ontology over-representation")
out <- out_dir("8_do_ora")
g <- gene_lists()

do <- enrichDO(gene = g$all, universe = g$universe, ont = "HDO", organism = "hsa",
               pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
save_object(do, "do_ora")
if (save_enrichment(do, out, "do_ora.tsv")) {
  save_plot(enrich_dotplot(do, "Disease Ontology: over-represented among significant genes"),
            out, "01_dotplot.png", width = 8, height = 6.5)
  save_plot(enrich_cnetplot(do, g$fold_change_symbol, "Top disease terms and their genes"),
            out, "02_cnetplot.png", width = 10, height = 8)
  save_plot(upsetplot(do, n = 8) + labs(title = "Genes shared between the top disease terms"),
            out, "03_upsetplot.png", width = 9, height = 5.5)
}
