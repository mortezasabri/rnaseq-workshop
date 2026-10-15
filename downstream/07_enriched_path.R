# 7_enriched_path: KEGG pathways. GO terms describe functions; pathways describe how genes work
# together (signaling cascades, metabolism). Same two questions as for GO: ORA and GSEA.
#   01_kegg_ora_dotplot     pathways over-represented among the significant genes
#   02_kegg_gsea_dotplot    pathways shifted up or down in the full ranking
#   03_kegg_gseaplot        running score of the top pathways
#   hsaXXXXX.pathview.png   the top pathway map, each gene box colored by its fold change
# Needs internet: KEGG is queried online (rest.kegg.jp).
source("downstream/_common.R")
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
})
step_header("7_enriched_path: KEGG pathways")
out <- out_dir("7_enriched_path")
g <- gene_lists()

kk <- with_retry(enrichKEGG(gene = g$all, universe = g$universe, organism = "hsa", pvalueCutoff = 0.05),
                 "KEGG over-representation")
if (has_terms(kk)) kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
if (save_enrichment(kk, out, "kegg_ora.tsv")) {
  save_plot(enrich_dotplot(kk, "KEGG pathways: over-represented among significant genes"),
            out, "01_kegg_ora_dotplot.png", width = 8, height = 6.5)
}

set.seed(SEED)
gk <- with_retry(gseKEGG(geneList = g$ranked, organism = "hsa", minGSSize = 15, maxGSSize = 500,
                         pvalueCutoff = 0.05, seed = TRUE, verbose = FALSE), "KEGG GSEA")
if (has_terms(gk)) gk <- setReadable(gk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
if (save_enrichment(gk, out, "kegg_gsea.tsv")) {
  save_plot(enrich_dotplot(gk, "KEGG pathways: GSEA", n = 10, split = ".sign") + facet_grid(. ~ .sign),
            out, "02_kegg_gsea_dotplot.png", width = 10, height = 6.5)
  top <- head(as.data.frame(gk)$ID, 3)
  save_plot(gseaplot2(gk, geneSetID = top, pvalue_table = FALSE, title = "KEGG GSEA: top pathways",
                      color = c(COL_TREAT, "#eda100", COL_REF)[seq_along(top)]),
            out, "03_kegg_gseaplot.png", width = 9, height = 6.5)
}

# Pathway map of the most significant real pathway. Skipped: disease maps (hsa05...), which
# collect genes from many pathways, and the ligand-receptor catalogs (e.g. "Cytokine-cytokine
# receptor interaction"), which are lists of receptors rather than pathways. Both make poor pictures.
pick_from <- if (has_terms(kk)) as.data.frame(kk) else if (has_terms(gk)) as.data.frame(gk) else NULL
catalogs <- c("hsa04060", "hsa04061", "hsa04080", "hsa04081", "hsa04082", "hsa04512", "hsa04514")
if (!is.null(pick_from)) pick_from <- pick_from[!grepl("^hsa05", pick_from$ID) & !pick_from$ID %in% catalogs, ]
if (!is.null(pick_from) && nrow(pick_from) > 0) {
  pathway <- pick_from$ID[1]
  message("  pathview map for ", pathway, " (", pick_from$Description[1], ")")
  suppressPackageStartupMessages(library(pathview))
  de <- load_object("de_results", "downstream/02_dea.R")
  de <- de[!is.na(de$entrez) & !duplicated(de$entrez), ]
  fc <- setNames(de$log2FC_shrunk, de$entrez)
  kegg_files <- file.path(tempdir(), "kegg")
  dir.create(kegg_files, showWarnings = FALSE)
  old <- setwd(out)                    # pathview writes its picture into the working folder
  with_retry(suppressMessages(pathview(gene.data = fc, pathway.id = pathway, species = "hsa",
                                       kegg.dir = kegg_files, limit = list(gene = 2, cpd = 1),
                                       low = list(gene = COL_REF), mid = list(gene = "#f0efec"),
                                       high = list(gene = COL_TREAT))), "KEGG pathway map")
  setwd(old)
  message("  wrote ", file.path(out, paste0(pathway, ".pathview.png")))
}
