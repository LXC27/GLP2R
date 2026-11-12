## 22a_gdc_paired_download.R : the supplied matrix has no tumour/normal pairs; fetch the matched primary-tumour STAR-Counts files of the 33 normal-tissue patients from GDC (same GENCODE v36 / STAR workflow as the supplied normals)
source(file.path(Sys.getenv("GLP2R_SCRIPTS", "W:/Re_build/05_Analysis/scripts"), "00_common.R"))
suppressPackageStartupMessages({library(jsonlite)})
sink(file.path(LOGD, "22a_gdc_paired_download.txt"), split = TRUE)
co <- readRDS(file.path(DC, "cohort.rds")); cl <- co$clin
npat <- sort(unique(cl$patient[cl$sample_type == "Solid Tissue Normal"])); cat("normal-tissue patients:", length(npat), "\n")
flt <- list(op = "and", content = list(
  list(op = "in", content = list(field = "cases.submitter_id", value = as.list(npat))),
  list(op = "=", content = list(field = "data_type", value = "Gene Expression Quantification")),
  list(op = "=", content = list(field = "analysis.workflow_type", value = "STAR - Counts")),
  list(op = "=", content = list(field = "access", value = "open"))))
body <- toJSON(list(filters = flt, fields = "file_id,file_name,md5sum,cases.submitter_id,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.portions.analytes.aliquots.submitter_id", format = "JSON", size = 500), auto_unbox = TRUE)
tf <- tempfile(); writeLines(body, tf)
out <- system2("curl", c("-s", "-S", "-X", "POST", "-H", shQuote("Content-Type: application/json"), "--data-binary", shQuote(paste0("@", tf)), "https://api.gdc.cancer.gov/files"), stdout = TRUE)
js <- fromJSON(paste(out, collapse = ""), simplifyVector = FALSE)
if (length(js$data$hits) == 0) { cat("RAW RESPONSE:", substr(paste(out, collapse = ""), 1, 800), "\n") }
hits <- js$data$hits; cat("files returned:", length(hits), "\n")
tab <- do.call(rbind, lapply(hits, function(h) { cs <- h$cases[[1]]; sm <- cs$samples[[1]]
  data.frame(file_id = h$file_id, file_name = h$file_name, md5 = h$md5sum, patient = cs$submitter_id, sample_type = sm$sample_type, sample_id = sm$submitter_id,
             aliquot = tryCatch(sm$portions[[1]]$analytes[[1]]$aliquots[[1]]$submitter_id, error = function(e) NA), stringsAsFactors = FALSE) }))
print(table(tab$sample_type))
tum <- tab[tab$sample_type == "Primary Tumor", ]
tum <- tum[order(tum$patient, tum$sample_id, tum$aliquot), ]; dup <- tum$patient[duplicated(tum$patient)]; cat("patients with >1 primary tumour file:", length(unique(dup)), "\n")
## keep one file per patient (lowest sample vial/aliquot barcode, prespecified)
tum1 <- tum[!duplicated(tum$patient), ]; cat("patients with matched primary tumour:", nrow(tum1), " of ", length(npat), "\n"); cat("missing:", paste(setdiff(npat, tum1$patient), collapse = ", "), "\n")
dir.create(file.path(DC, "gdc_paired"), showWarnings = FALSE)
for (i in seq_len(nrow(tum1))) { dest <- file.path(DC, "gdc_paired", paste0(tum1$file_id[i], ".tsv"))
  if (!file.exists(dest) || file.size(dest) < 1000) download.file(paste0("https://api.gdc.cancer.gov/data/", tum1$file_id[i]), dest, mode = "wb", quiet = TRUE)
  tum1$bytes[i] <- file.size(dest); tum1$md5_ok[i] <- as.character(tools::md5sum(dest)) == tum1$md5[i] }
print(tum1[, c("patient", "sample_id", "bytes", "md5_ok")]); write_csv_utf8(tum1, file.path(DC, "gdc_paired_manifest.csv"))
## validation: the supplied normal counts must equal the GDC STAR unstranded counts of the same aliquot
nrm <- tab[tab$sample_type == "Solid Tissue Normal", ][1:2, ]
for (i in seq_len(nrow(nrm))) { dest <- file.path(DC, "gdc_paired", paste0("normal_check_", nrm$file_id[i], ".tsv")); if (!file.exists(dest)) download.file(paste0("https://api.gdc.cancer.gov/data/", nrm$file_id[i]), dest, mode = "wb", quiet = TRUE)
  x <- read.delim(dest, comment.char = "#", stringsAsFactors = FALSE); x <- x[grepl("^ENSG", x$gene_id), ]
  bc <- grep(paste0("^", nrm$aliquot[i]), colnames(co$counts), value = TRUE); if (!length(bc)) bc <- grep(paste0("^", nrm$patient[i], "-11"), colnames(co$counts), value = TRUE)
  cat(nrm$patient[i], "->", bc, " identical counts:", identical(as.numeric(co$counts[x$gene_id, bc[1]]), as.numeric(x$unstranded)), " max abs diff:", max(abs(co$counts[x$gene_id, bc[1]] - x$unstranded)), "\n") }
write_csv_utf8(tab, file.path(DC, "gdc_normal_tumour_file_listing.csv"))
log_msg("done"); sink()
