RB  <- Sys.getenv("GLP2R_RB", "W:/Re_build")
LIB <- if (getRversion() < "4.5") file.path(RB, "_work", "Rlib") else Sys.getenv("R_LIBS_USER")
if (!nzchar(LIB)) LIB <- file.path(RB, "_work", "Rlib")
dir.create(LIB, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(LIB, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), BioC_mirror = "https://bioconductor.org", timeout = 1800)
for (p in c("DOSE", "clusterProfiler", "enrichplot")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat("installing", p, "\n")
    tryCatch(BiocManager::install(p, lib = LIB, update = FALSE, ask = FALSE, type = "binary"), error = function(e) cat("FAILED", p, conditionMessage(e), "\n"))
  }
}
for (p in c("DOSE", "clusterProfiler")) cat(p, tryCatch({ suppressPackageStartupMessages(loadNamespace(p)); "OK" }, error = function(e) conditionMessage(e)), "\n")
