## 01a_rdata_inventory.R : inventory of the author's original R objects (V:/ = Z:\文章\Rdata)
## Purpose: record class / dimensions / key fields of every .rds so that usable vs. unusable objects are documented.
options(width = 200, warn = -1)
RB <- Sys.getenv("GLP2R_RB", "W:/Re_build")
if (getRversion() < "4.5") { .libPaths(c(file.path(RB, "_work", "Rlib"), .libPaths())) } else { ul <- Sys.getenv("R_LIBS_USER"); if (nzchar(ul) && dir.exists(ul)) .libPaths(c(ul, .libPaths())) }
suppressWarnings(for (p in c("SummarizedExperiment", "DESeq2")) try(suppressMessages(library(p, character.only = TRUE)), silent = TRUE))
files <- list.files(Sys.getenv("GLP2R_ORIG", "V:/"), pattern = "\\.rds$", full.names = TRUE)
sink(file.path(RB, "05_Analysis", "logs", "01a_rdata_inventory.txt"), split = TRUE)
cat("R:", R.version.string, "\nDate:", format(Sys.time()), "\n")
show <- function(x, nm, depth = 0) {
  ind <- strrep("  ", depth)
  d <- tryCatch(if (!is.null(dim(x))) paste(dim(x), collapse = "x") else length(x), error = function(e) "?")
  cat(ind, "- ", nm, ": ", paste(class(x), collapse = ","), " [", d, "]\n", sep = "")
  if (depth < 1 && is.list(x) && !is.data.frame(x) && !isS4(x)) {
    for (n in head(names(x), 30)) show(x[[n]], n, depth + 1)
  }
}
for (f in files) {
  cat("\n=== ", basename(f), " (", round(file.size(f) / 1024), " KB)\n", sep = "")
  x <- tryCatch(readRDS(f), error = function(e) { cat("READ ERROR: ", conditionMessage(e), "\n"); NULL })
  if (is.null(x)) next
  tryCatch(show(x, basename(f)), error = function(e) cat("  (show failed: ", conditionMessage(e), "; class=", paste(class(x), collapse = ","), ")\n", sep = ""))
  if (is.data.frame(x)) { cat("colnames: ", paste(head(colnames(x), 40), collapse = ", "), "\n") }
}
sink()
