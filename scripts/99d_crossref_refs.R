## 99d_crossref_refs.R : retrieve bibliographic metadata for every verified DOI from Crossref (title, authors, journal, volume, pages) -> 06_References/crossref_metadata.csv
if (getRversion() < "4.5") { .libPaths(c(file.path(Sys.getenv("GLP2R_RB", "W:/Re_build"), "_work", "Rlib"), .libPaths())) } else { ul <- Sys.getenv("R_LIBS_USER"); if (nzchar(ul) && dir.exists(ul)) .libPaths(c(ul, .libPaths())) }
suppressPackageStartupMessages(library(jsonlite))
RB <- Sys.getenv("GLP2R_RB", "W:/Re_build")
ref <- read.csv(file.path(RB, "06_References", "reference_verification.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
skip <- c("4", "8", "9", "10", "34", "35")
ref <- ref[!ref$orig_no %in% skip, ]
doi <- sub("^doi:", "", ref$verified_identifier); doi <- trimws(sub("\\s+\\(.*$", "", doi))
cache <- file.path(RB, "05_Analysis", "data_cache", "crossref_cache"); dir.create(cache, showWarnings = FALSE)
getmeta <- function(d) {
  f <- file.path(cache, paste0(gsub("[^A-Za-z0-9]", "_", d), ".json"))
  if (!file.exists(f) || file.size(f) < 50) {
    u <- paste0("https://api.crossref.org/works/", utils::URLencode(d, reserved = TRUE), "?mailto=glp2r.revision@example.org")
    for (k in 1:4) { system2("curl.exe", c("-s", "-L", "--max-time", "40", "-o", shQuote(f), shQuote(u))); if (file.exists(f) && file.size(f) > 50 && grepl("^\\{", readLines(f, n = 1, warn = FALSE))) break; Sys.sleep(3) }
    Sys.sleep(1.5)
  }
  tryCatch(fromJSON(f, simplifyVector = FALSE)$message, error = function(e) NULL)
}
out <- list()
for (i in seq_len(nrow(ref))) {
  m <- getmeta(doi[i]); if (is.null(m)) { out[[i]] <- data.frame(orig_no = ref$orig_no[i], doi = doi[i], found = FALSE, stringsAsFactors = FALSE); cat("NOT RETRIEVED:", ref$orig_no[i], doi[i], "\n"); next }
  au <- vapply(m$author, function(a) { fam <- if (!is.null(a$family)) a$family else if (!is.null(a$name)) a$name else ""; gv <- if (!is.null(a$given)) a$given else ""; ini <- paste(substring(strsplit(gsub("[-.]", " ", gv), " +")[[1]], 1, 1), collapse = ""); trimws(paste(fam, ini)) }, "")
  yr <- tryCatch(m$issued$`date-parts`[[1]][[1]], error = function(e) NA)
  ct <- if (length(m$`container-title`)) m$`container-title`[[1]] else if (!is.null(m$institution)) "bioRxiv" else ""
  out[[i]] <- data.frame(orig_no = ref$orig_no[i], doi = doi[i], found = TRUE, authors = paste(au, collapse = "; "), n_authors = length(au), title = gsub("<[^>]+>", "", if (length(m$title)) m$title[[1]] else ""), journal = ct, year = yr, volume = if (!is.null(m$volume)) m$volume else "", issue = if (!is.null(m$issue)) m$issue else "", page = if (!is.null(m$page)) m$page else if (!is.null(m$`article-number`)) m$`article-number` else "", stringsAsFactors = FALSE)
  cat(sprintf("%-4s %-10s %s (%s) %s\n", ref$orig_no[i], substr(au[1], 1, 10), substr(out[[i]]$title, 1, 60), yr, substr(ct, 1, 25)))
}
res <- do.call(rbind, lapply(out, function(z) { cols <- c("orig_no", "doi", "found", "authors", "n_authors", "title", "journal", "year", "volume", "issue", "page"); for (c in setdiff(cols, names(z))) z[[c]] <- NA; z[, cols] }))
write.csv(res, file.path(RB, "06_References", "crossref_metadata.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("retrieved:", sum(res$found), "of", nrow(res), "\n")
