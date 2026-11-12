## 00_setup.R : install dependencies into a project-local library (binary packages only)
RB  <- Sys.getenv("GLP2R_RB", "W:/Re_build")   # repository root (holds 05_Analysis/ and _work/Rlib)
LIB <- file.path(RB, "_work", "Rlib")
dir.create(LIB, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(LIB, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), pkgType = "binary",
        install.packages.compile.from.source = "never", timeout = 600)

cran <- c("BiocManager", "data.table", "dplyr", "tidyr", "readr", "stringr", "jsonlite", "httr", "xml2",
          "ggplot2", "ggrepel", "patchwork", "cowplot", "scales", "RColorBrewer", "gridExtra",
          "survminer", "glmnet", "timeROC", "survivalROC", "mice", "maxstat", "pROC", "rms", "Hmisc",
          "officer", "flextable", "UCSCXenaTools", "msigdbr", "pheatmap", "ggpubr", "ggsci", "reshape2",
          "tidyestimate", "R.utils", "ragg", "systemfonts")
need <- setdiff(cran, rownames(installed.packages()))
if (length(need)) {
  for (p in need) tryCatch(install.packages(p, lib = LIB, dependencies = c("Depends", "Imports")),
                           error = function(e) message("FAILED ", p, ": ", conditionMessage(e)))
}

bioc <- c("S4Vectors", "IRanges", "GenomicRanges", "SummarizedExperiment", "DESeq2", "limma", "edgeR",
          "fgsea", "GSVA", "GEOquery", "sva", "impute", "preprocessCore", "GO.db", "AnnotationDbi",
          "org.Hs.eg.db", "WGCNA")
options(BioC_mirror = "https://bioconductor.org")
need <- setdiff(bioc, rownames(installed.packages()))
for (p in need) {
  tryCatch({
    if (p == "WGCNA") install.packages(p, lib = LIB, dependencies = c("Depends", "Imports"))
    else BiocManager::install(p, lib = LIB, update = FALSE, ask = FALSE, type = "binary", dependencies = c("Depends", "Imports"))
  }, error = function(e) message("FAILED ", p, ": ", conditionMessage(e)))
}

ip <- rownames(installed.packages(lib.loc = LIB))
cat("\nInstalled in project lib:", length(ip), "\n")
missing <- setdiff(c(cran, bioc), rownames(installed.packages()))
cat("Still missing:", if (length(missing)) paste(missing, collapse = ", ") else "none", "\n")
