## 00b_install_extra.R : install packages that the first pass missed (can be re-run safely)
LIB <- file.path(Sys.getenv("GLP2R_RB", "W:/Re_build"), "_work", "Rlib")
.libPaths(c(LIB, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), BioC_mirror = "https://bioconductor.org", timeout = 900,
        install.packages.compile.from.source = "never")
extra_bioc <- c("GenomeInfoDbData", "hgu133plus2.db", "org.Hs.eg.db", "GSVA", "fgsea", "GEOquery", "WGCNA", "impute", "preprocessCore", "GO.db")
for (p in extra_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat("installing", p, "\n")
    tryCatch(BiocManager::install(p, lib = LIB, update = FALSE, ask = FALSE, type = "binary"),
             error = function(e) cat("FAILED", p, conditionMessage(e), "\n"))
  }
}
extra_cran <- c("jomo", "mice", "tidyestimate", "UCSCXenaTools", "msigdbr", "officer", "flextable", "ragg", "systemfonts", "WGCNA", "pheatmap")
for (p in extra_cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat("installing", p, "\n")
    tryCatch(install.packages(p, lib = LIB, type = "binary"), error = function(e) cat("FAILED", p, conditionMessage(e), "\n"))
  }
}
chk <- c("SummarizedExperiment", "DESeq2", "GenomeInfoDbData", "GSVA", "fgsea", "limma", "edgeR", "mice", "tidyestimate", "msigdbr", "WGCNA", "officer", "flextable", "survminer", "timeROC", "glmnet", "maxstat", "ggplot2", "patchwork", "ggrepel", "hgu133plus2.db", "org.Hs.eg.db", "GEOquery", "ragg")
print(sapply(chk, requireNamespace, quietly = TRUE))
