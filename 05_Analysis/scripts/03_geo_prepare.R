## 03_geo_prepare.R : GSE39582 / GSE17536 expression (gene level) and phenotype from the GEO series-matrix files
options(warn = 1, stringsAsFactors = FALSE)
RB <- Sys.getenv("GLP2R_RB", "W:/Re_build"); DC <- file.path(RB, "05_Analysis", "data_cache")
.libPaths(c(file.path(RB, "_work", "Rlib"), .libPaths()))
suppressPackageStartupMessages(library(data.table))
log <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

parse_series <- function(f) {
  con <- gzfile(f, "rt"); on.exit(close(con))
  lines <- readLines(con)
  b <- which(lines == "!series_matrix_table_begin"); e <- which(lines == "!series_matrix_table_end")
  meta <- lines[seq_len(b - 1)]
  tab <- fread(text = paste(lines[(b + 1):(e - 1)], collapse = "\n"), sep = "\t", header = TRUE, quote = "\"")
  gsm <- gsub("\"", "", strsplit(meta[grep("^!Sample_geo_accession", meta)], "\t")[[1]][-1])
  ch <- meta[grep("^!Sample_characteristics_ch1", meta)]
  ph <- data.frame(geo_accession = gsm, stringsAsFactors = FALSE)
  for (ln in ch) {
    v <- gsub("\"", "", strsplit(ln, "\t")[[1]][-1])
    k <- sub(":.*", "", v[nzchar(v)][1])
    ph[[k]] <- ifelse(grepl(":", v), trimws(sub("^[^:]*:", "", v)), NA)
  }
  ph$title <- gsub("\"", "", strsplit(meta[grep("^!Sample_title", meta)], "\t")[[1]][-1])
  expr <- as.matrix(tab[, -1]); rownames(expr) <- gsub("\"", "", tab[[1]])
  colnames(expr) <- gsub("\"", "", colnames(expr))
  list(expr = expr, pheno = ph, meta = meta[grep("^!Series_(title|summary|overall_design|platform_id|submission_date|last_update_date)|^!Sample_(data_processing)", meta)][1:12])
}

## GPL570 platform table -> gene symbol per probe
gp <- fread(file.path(DC, "GPL570_table.txt"), sep = "\t", skip = "ID\t", quote = "")
gp <- gp[!grepl("^!", ID)]
sym_col <- grep("^Gene Symbol$", names(gp), value = TRUE)
gp$sym <- trimws(sub(" ///.*", "", gp[[sym_col]]))
gp <- gp[nzchar(sym) & sym != "---"]
log("GPL570 probes with symbols:", nrow(gp))

collapse_genes <- function(expr) {
  keep <- intersect(rownames(expr), gp$ID)
  m <- expr[keep, , drop = FALSE]
  sym <- gp$sym[match(keep, gp$ID)]
  mu <- rowMeans(m, na.rm = TRUE)
  ord <- order(sym, -mu)
  first <- !duplicated(sym[ord])
  sel <- ord[first]
  g <- m[sel, , drop = FALSE]; rownames(g) <- sym[sel]
  list(gene = g, probe_of_gene = setNames(keep[sel], sym[sel]))
}

out <- list()
for (acc in c("GSE39582", "GSE17536")) {
  log("parsing", acc)
  s <- parse_series(file.path(DC, paste0(acc, "_series_matrix.txt.gz")))
  log(acc, "probes x samples:", paste(dim(s$expr), collapse = " x "), "; expr range:", paste(round(range(s$expr, na.rm = TRUE), 2), collapse = " .. "))
  log(acc, "pheno columns:", paste(names(s$pheno), collapse = " | "))
  cg <- collapse_genes(s$expr)
  log(acc, "genes:", nrow(cg$gene), "; GLP2R present:", "GLP2R" %in% rownames(cg$gene), "; GLP2R probe:", cg$probe_of_gene["GLP2R"])
  out[[acc]] <- list(expr = cg$gene, pheno = s$pheno, probe_of_gene = cg$probe_of_gene,
                     glp2r_probes = intersect(rownames(s$expr), gp$ID[gp$sym == "GLP2R"]),
                     glp2r_probe_expr = s$expr[intersect(rownames(s$expr), gp$ID[gp$sym == "GLP2R"]), , drop = FALSE])
}
saveRDS(out, file.path(DC, "geo_cohorts.rds"))
log("done")
