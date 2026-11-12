## 02_extract_external.R : extract the pieces of the public data needed for the analyses
##  - UCSC Xena legacy TCGA-COAD/READ: HM450 beta values for GLP2R (+ control genes), GISTIC thresholded CNV, clinicalMatrix (MSI/CIMP/BRAF)
##  - GEO GSE39582 / GSE17536: expression (gene level) + phenotype
options(warn = 1, stringsAsFactors = FALSE)
RB <- Sys.getenv("GLP2R_RB", "W:/Re_build"); DC <- file.path(RB, "05_Analysis", "data_cache")
.libPaths(c(file.path(RB, "_work", "Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(data.table) })
log <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

## ---------------------------------------------------------------- HM450 (stream; the files are too large to read whole)
pm <- fread(file.path(DC, "illuminaMethyl450_probeMap.tsv"), sep = "\t")
setnames(pm, c("id", "gene", "chrom", "chromStart", "chromEnd", "strand"))
ctrl <- c("GLP2R", "MLH1", "CACNA1G", "IGF2", "NEUROG1", "RUNX3", "SOCS1", "CDKN2A")
has_gene <- function(s, g) vapply(strsplit(s, ",", fixed = TRUE), function(z) any(z %in% g), logical(1))
pm_sub <- pm[has_gene(gene, ctrl)]
log("HM450 probes kept:", nrow(pm_sub), "(GLP2R:", sum(has_gene(pm_sub$gene, "GLP2R")), ")")
fwrite(pm_sub, file.path(DC, "hm450_probe_annotation_subset.csv"))

read_hm450 <- function(f, keep) {
  con <- gzfile(f, "rt"); on.exit(close(con))
  hdr <- readLines(con, n = 1)
  out <- list(); n <- 0
  repeat {
    ch <- readLines(con, n = 20000)
    if (!length(ch)) break
    n <- n + length(ch)
    id <- sub("\t.*", "", ch)
    k <- ch[id %in% keep]
    if (length(k)) out[[length(out) + 1]] <- k
  }
  log("  rows scanned:", n)
  txt <- c(hdr, unlist(out))
  as.data.frame(fread(text = paste(txt, collapse = "\n"), sep = "\t", header = TRUE), check.names = FALSE)
}
hm <- lapply(c(COAD = "COAD_HM450.tsv.gz", READ = "READ_HM450.tsv.gz"), function(f) {
  log("reading", f); d <- read_hm450(file.path(DC, f), pm_sub$id)
  rownames(d) <- d[[1]]; d[[1]] <- NULL; as.matrix(d)
})
common <- intersect(rownames(hm$COAD), rownames(hm$READ))
hm_all <- cbind(hm$COAD[common, , drop = FALSE], hm$READ[common, , drop = FALSE])
log("HM450 subset matrix:", paste(dim(hm_all), collapse = " x "))
saveRDS(list(beta = hm_all, probes = as.data.frame(pm_sub)), file.path(DC, "hm450_subset.rds"))

## ---------------------------------------------------------------- GISTIC thresholded CNV (gene x sample, -2..2)
rd_gz <- function(f) { d <- as.data.frame(fread(cmd = paste("gzip -dc", shQuote(file.path(DC, f))), sep = "\t", header = TRUE), check.names = FALSE); rownames(d) <- d[[1]]; d[[1]] <- NULL; as.matrix(d) }
gi <- tryCatch({
  a <- rd_gz("COAD_gistic_thresholded.tsv.gz"); b <- rd_gz("READ_gistic_thresholded.tsv.gz")
  cg <- intersect(rownames(a), rownames(b)); cbind(a[cg, ], b[cg, ])
}, error = function(e) { log("GISTIC read via gzip failed:", conditionMessage(e)); NULL })
if (is.null(gi)) {
  rd2 <- function(f) { con <- gzfile(file.path(DC, f), "rt"); d <- read.delim(con, check.names = FALSE, row.names = 1); as.matrix(d) }
  a <- rd2("COAD_gistic_thresholded.tsv.gz"); b <- rd2("READ_gistic_thresholded.tsv.gz"); cg <- intersect(rownames(a), rownames(b)); gi <- cbind(a[cg, ], b[cg, ])
}
log("GISTIC matrix:", paste(dim(gi), collapse = " x "), "; GLP2R present:", "GLP2R" %in% rownames(gi))
saveRDS(gi, file.path(DC, "gistic_thresholded_coadread.rds"))

## ---------------------------------------------------------------- Xena clinical matrices
cl <- rbindlist(list(fread(file.path(DC, "COAD_clinicalMatrix.tsv"), sep = "\t")[, cohort := "COAD"],
                     fread(file.path(DC, "READ_clinicalMatrix.tsv"), sep = "\t")[, cohort := "READ"]), fill = TRUE)
saveRDS(as.data.frame(cl), file.path(DC, "xena_clinical_coadread.rds"))
log("Xena clinical rows:", nrow(cl))
log("done")
