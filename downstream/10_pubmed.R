# 10_pubmed: how much has been published about our findings, year by year? Counts articles in
# Europe PMC (which includes PubMed) and plots them as a share of all articles that year.
#   01_trend_diseases   the top Disease Ontology terms from 8_do_ora
#   02_trend_genes      the top up-regulated genes together with "glucocorticoid OR dexamethasone"
# Needs internet (Europe PMC REST API).
source("downstream/_common.R")
suppressPackageStartupMessages(library(enrichplot))
step_header("10_pubmed: publication trends")
out <- out_dir("10_pubmed")
last_year <- as.integer(format(Sys.Date(), "%Y")) - 1      # last complete year
period <- 2005:last_year

trend_plot <- function(queries, labels, title) {
  p <- with_retry(pmcplot(queries, period, proportion = TRUE), "Europe PMC query")
  # one color per query, in the fixed order of the palette
  p + scale_color_manual(values = c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4")[seq_along(queries)],
                         breaks = queries, labels = labels, name = NULL) +
    scale_y_continuous(labels = function(x) x * 1e6) +
    labs(title = title, x = NULL, y = "articles per million published") +
    theme(legend.position = "right")
}

do <- tryCatch(load_object("do_ora", "downstream/08_do_ora.R"), error = function(e) NULL)
if (has_terms(do)) {
  terms <- head(as.data.frame(do)$Description, 5)
  write_table(data.frame(query = terms), out, "01_trend_diseases_queries.tsv")
  save_plot(trend_plot(terms, terms, "Publications on the top disease terms"),
            out, "01_trend_diseases.png", width = 9, height = 4.5)
} else {
  message("  no Disease Ontology terms (8_do_ora); skipping the disease trend")
}

de <- load_object("de_results", "downstream/02_dea.R")
genes <- head(unique(c(cfg$SHOWCASE_GENE, de$symbol[de$direction == "up" & !is.na(de$symbol)])), 5)
queries <- paste0(genes, " AND (glucocorticoid OR dexamethasone)")
write_table(data.frame(query = queries), out, "02_trend_genes_queries.tsv")
save_plot(trend_plot(queries, genes, "Publications linking top genes to glucocorticoids"),
          out, "02_trend_genes.png", width = 9, height = 4.5)
