#!/usr/bin/env Rscript
# Install the R packages the downstream scripts and slides need. Installed packages are skipped.
# Not needed if you use the conda environment (environment.yml), which already has them.
#
# Usage:  Rscript install.R           install what is missing (10-30 minutes the first time)
#         Rscript install.R --check   only report; exit code 1 if anything is missing
# Needs R >= 4.4 and an internet connection.

required_packages <- function() {
  list(
    cran = c("ggplot2", "ggrepel", "ggbeeswarm", "patchwork", "pheatmap", "RColorBrewer",
             "dplyr", "tidyr", "readr", "ggupset", "ggnewscale", "ggridges", "europepmc",
             "rmarkdown", "knitr"),
    bioc = c("DESeq2", "edgeR", "apeglm", "airway", "tximport", "org.Hs.eg.db",
             "clusterProfiler", "enrichplot", "DOSE", "ReactomePA", "pathview")
  )
}

missing_packages <- function() {
  pkgs <- unlist(required_packages(), use.names = FALSE)
  pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
}

if (identical(commandArgs(trailingOnly = TRUE), "--check")) {
  missing <- missing_packages()
  if (length(missing) > 0) {
    message("Missing R packages: ", paste(missing, collapse = ", "),
            "\nInstall them with:  Rscript install.R   (or use the conda environment, see README.md)")
    quit(status = 1)
  }
  message("All required R packages are installed.")
  quit(status = 0)
}

missing <- missing_packages()
if (length(missing) == 0) {
  message("All required R packages are already installed.")
} else {
  message("Installing: ", paste(missing, collapse = ", "))
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install(missing, update = FALSE, ask = FALSE)
  still_missing <- missing_packages()
  if (length(still_missing) > 0) {
    stop("These packages failed to install: ", paste(still_missing, collapse = ", "),
         ". Scroll up for the compiler or download error. On Linux this is usually a missing ",
         "system library (e.g. libcurl, libxml2, libssl development packages).")
  }
  message("All required packages are installed.")
}
