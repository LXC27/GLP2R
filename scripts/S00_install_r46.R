## S00_install_r46.R : the bundled Rlib was built for R 4.4.1 and cannot be loaded by R >= 4.5.
## Installs the CRAN packages needed by the figure scripts (95/97/98) and S01 into the user library of the running R.
if (getRversion() >= "4.5") {
  ul <- Sys.getenv("R_LIBS_USER"); dir.create(ul, recursive = TRUE, showWarnings = FALSE); .libPaths(c(ul, .libPaths()))
  need <- c("ggplot2", "dplyr", "tidyr", "tibble", "patchwork", "ggrepel", "cowplot", "ragg", "scales", "glmnet", "survival")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  cat("R", as.character(getRversion()), "| user library:", ul, "\nmissing:", paste(miss, collapse = ", "), "\n")
  for (repo in c("https://mirrors.tuna.tsinghua.edu.cn/CRAN", "https://cloud.r-project.org")) {
    miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]; if (!length(miss)) break
    try(install.packages(miss, lib = ul, repos = repo, type = "binary"))
  }
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  cat("still missing:", if (length(miss)) paste(miss, collapse = ", ") else "none", "\n")
  for (p in need) cat(p, as.character(tryCatch(packageVersion(p), error = function(e) NA)), "\n")
} else cat("R < 4.5: bundled library is used; nothing to install\n")
