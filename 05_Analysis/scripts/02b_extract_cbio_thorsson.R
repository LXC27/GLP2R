## 02b_extract_cbio_thorsson.R : parse cBioPortal PanCanAtlas clinical (MSI scores, subtype, survival) and Thorsson 2018 immune subtypes
options(width = 200, warn = -1, timeout = 900, repos = c(CRAN = "https://cloud.r-project.org"))
RB <- Sys.getenv("GLP2R_RB", "W:/Re_build"); LIB <- file.path(RB, "_work", "Rlib")
.libPaths(c(LIB, .libPaths()))
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl", lib = LIB, type = "binary")
suppressPackageStartupMessages({library(jsonlite); library(readxl)})
DC <- file.path(RB, "05_Analysis", "data_cache")
sink(file.path(RB, "05_Analysis", "logs", "02b_extract_cbio_thorsson.txt"), split = TRUE)
## ---- cBioPortal sample + patient clinical -> wide
rd <- function(f) { d <- fromJSON(file.path(DC, f)); d[, c("sampleId", "patientId", "clinicalAttributeId", "value")[c("sampleId", "patientId", "clinicalAttributeId", "value") %in% names(d)]] }
s <- rd("cbio_pancan_sample.json"); p <- rd("cbio_pancan_patient.json")
cat("sample rows:", nrow(s), " patient rows:", nrow(p), "\n")
sw <- reshape(s[, c("sampleId", "patientId", "clinicalAttributeId", "value")], idvar = c("sampleId", "patientId"), timevar = "clinicalAttributeId", direction = "wide")
colnames(sw) <- sub("^value\\.", "", colnames(sw))
pw <- reshape(p[, c("patientId", "clinicalAttributeId", "value")], idvar = "patientId", timevar = "clinicalAttributeId", direction = "wide")
colnames(pw) <- sub("^value\\.", "", colnames(pw))
cat("sample-wide:", dim(sw), " patient-wide:", dim(pw), "\n")
cbio <- list(sample = sw, patient = pw)
saveRDS(cbio, file.path(DC, "cbio_pancan_coadread.rds"))
print(table(sw$SUBTYPE, useNA = "ifany"))
print(summary(as.numeric(sw$MSI_SCORE_MANTIS))); print(summary(as.numeric(sw$MSI_SENSOR_SCORE)))
print(table(sw$SAMPLE_TYPE, useNA = "ifany"))
cat("OS_STATUS:\n"); print(table(pw$OS_STATUS, useNA = "ifany"))
## ---- Thorsson xlsx
zf <- file.path(DC, "Thorsson_download.zip"); xf <- file.path(DC, "Thorsson_2018_immune_landscape.xlsx")
file.copy(zf, xf, overwrite = TRUE)
sh <- excel_sheets(xf); print(sh)
for (i in seq_along(sh)) {
  d <- tryCatch(read_excel(xf, sheet = i, n_max = 5), error = function(e) NULL)
  if (!is.null(d)) { cat("\n== sheet", i, sh[i], ": ", ncol(d), "cols\n"); print(head(colnames(d), 30)) }
}
sink()
