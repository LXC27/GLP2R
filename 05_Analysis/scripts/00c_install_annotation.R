## 00c: annotation data packages are only distributed as source (no compilation needed)
LIB <- file.path(Sys.getenv("GLP2R_RB", "W:/Re_build"), "_work", "Rlib")
.libPaths(c(LIB, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), BioC_mirror = "https://bioconductor.org", timeout = 1800)
for (p in c("GenomeInfoDbData", "GO.db", "org.Hs.eg.db", "hgu133plus2.db")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat("installing (source)", p, "\n")
    tryCatch(BiocManager::install(p, lib = LIB, update = FALSE, ask = FALSE, type = "source", INSTALL_opts = "--no-multiarch"),
             error = function(e) cat("FAILED", p, conditionMessage(e), "\n"))
  }
}
chk <- c("GenomeInfoDbData", "GenomeInfoDb", "SummarizedExperiment", "DESeq2", "GSVA", "GEOquery", "GO.db", "org.Hs.eg.db", "hgu133plus2.db", "WGCNA")
for (p in chk) {
  r <- tryCatch({ suppressPackageStartupMessages(loadNamespace(p)); "OK" }, error = function(e) paste("ERR:", conditionMessage(e)))
  cat(sprintf("%-22s %s\n", p, r))
}
