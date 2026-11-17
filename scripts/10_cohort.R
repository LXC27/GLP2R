## 10_cohort.R : build the analysis cohort (TCGA-COAD/READ), harmonise clinical covariates, merge external annotations, REMARK flow
source(file.path(Sys.getenv("GLP2R_SCRIPTS", "W:/Re_build/05_Analysis/scripts"), "00_common.R"))
suppressPackageStartupMessages({library(SummarizedExperiment); library(org.Hs.eg.db); library(readxl)})
sink(file.path(LOGD, "10_cohort.txt"), split = TRUE)
log_msg("load TCGA_CRC_combined")
x <- readRDS(file.path(Sys.getenv("GLP2R_ORIG", "V:"), "TCGA_CRC_combined.rds"))
counts <- x$counts; cl <- as.data.frame(x$clinical); stopifnot(identical(colnames(counts), rownames(cl)))
## ---- gene map (Ensembl -> HGNC symbol)
ens <- sub("\\..*$", "", rownames(counts))
sym <- suppressMessages(mapIds(org.Hs.eg.db, keys = ens, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
gmap <- data.frame(ensembl_full = rownames(counts), ensembl = ens, symbol = unname(sym), stringsAsFactors = FALSE)
cat("genes:", nrow(gmap), " mapped symbols:", sum(!is.na(gmap$symbol)), " unique symbols:", length(unique(na.omit(gmap$symbol))), "\n")
print(gmap[gmap$symbol %in% c("GLP2R", "GLP1R", "GIP", "GIPR", "GCG", "DPP4", "DNMT1", "DNMT3A", "DNMT3B", "TET2", "CD274"), ])
## ---- clinical harmonisation
tum <- cl$tissue_type == "Tumor"; prim <- cl$sample_type == "Primary Tumor"
cl$barcode_sample <- rownames(cl); cl$patient <- substr(rownames(cl), 1, 12)
cl$condition <- ifelse(cl$tissue_type == "Normal", "Normal", "Tumor")
cl$stage_grp <- factor(dplyr::case_when(grepl("^Stage IV", cl$ajcc_pathologic_stage) ~ "IV", grepl("^Stage III", cl$ajcc_pathologic_stage) ~ "III",
                                        grepl("^Stage II", cl$ajcc_pathologic_stage) ~ "II", grepl("^Stage I", cl$ajcc_pathologic_stage) ~ "I", TRUE ~ NA_character_), levels = c("I", "II", "III", "IV"))
cl$age_y <- suppressWarnings(as.numeric(cl$age_at_index)); i <- is.na(cl$age_y); cl$age_y[i] <- suppressWarnings(-as.numeric(cl$days_to_birth[i]) / 365.25)
cl$sex <- factor(ifelse(cl$gender %in% c("female", "male"), cl$gender, NA), levels = c("female", "male"))
cl$site <- factor(cl$cancer_type, levels = c("COAD", "READ"))
cl$vital <- as.character(cl$vital_status)
## original-style OS from raw GDC fields (used for the reconciliation only)
cl$os_gdc_time <- ifelse(cl$vital == "Dead", as.numeric(cl$days_to_death), as.numeric(cl$days_to_last_follow_up))
cl$os_gdc_event <- ifelse(cl$vital == "Dead", 1L, ifelse(cl$vital == "Alive", 0L, NA_integer_))
## ---- external annotations
pat <- cl$patient
cb <- readRDS(file.path(DC, "cbio_pancan_coadread.rds"))
sw <- cb$sample; sw$patient <- sw$patientId
cat("cBio patient-level columns:", paste(colnames(cb$patient), collapse = ", "), "\n")
m <- match(pat, sw$patient)
cl$mantis <- as.numeric(sw$MSI_SCORE_MANTIS[m]); cl$msisensor <- as.numeric(sw$MSI_SENSOR_SCORE[m])
cl$tmb <- as.numeric(sw$TMB_NONSYNONYMOUS[m]); cl$mut_count <- as.numeric(sw$MUTATION_COUNT[m]); cl$fga <- as.numeric(sw$FRACTION_GENOME_ALTERED[m]); cl$aneuploidy <- as.numeric(sw$ANEUPLOIDY_SCORE[m])
pw <- cb$patient; mp <- match(pat, pw$patientId)
cat("cBio SUBTYPE:\n"); print(table(pw$SUBTYPE, useNA = "ifany"))
cl$pancan_subtype <- pw$SUBTYPE[mp]
## TCGA-CDR (Liu et al. 2018) curated endpoints
xf <- file.path(DC, "Thorsson_2018_immune_landscape.xlsx")
cdr <- as.data.frame(read_excel(xf, sheet = "TCGA-CDR")); ex <- as.data.frame(read_excel(xf, sheet = "ExtraEndpoints"))
cdr <- cdr[cdr$type %in% c("COAD", "READ"), c("bcr_patient_barcode", "OS", "OS.time", "DSS", "DSS.time", "DFI", "DFI.time", "PFI", "PFI.time")]
cat("CDR COAD/READ rows:", nrow(cdr), "\n")
mc <- match(pat, cdr$bcr_patient_barcode)
for (v in c("OS", "OS.time", "DSS", "DSS.time", "DFI", "DFI.time", "PFI", "PFI.time")) cl[[paste0("cdr_", sub("\\.", "_", v))]] <- suppressWarnings(as.numeric(cdr[[v]][mc]))
## ABSOLUTE purity (Taylor et al. 2018 / PanCanAtlas)
ab <- read.delim(file.path(DC, "ABSOLUTE_mastercalls.txt"), check.names = FALSE, stringsAsFactors = FALSE)
ab$patient <- substr(ab$array, 1, 12); ab <- ab[grepl("-01$", ab$array), ]
ma <- match(pat, ab$patient); cl$purity_absolute <- suppressWarnings(as.numeric(ab$purity[ma])); cl$ploidy_absolute <- suppressWarnings(as.numeric(ab$ploidy[ma]))
## Thorsson immune subtype
th <- read.delim(gzfile(file.path(DC, "Thorsson_Subtype_Immune_Model_Based.txt.gz")), check.names = FALSE, stringsAsFactors = FALSE)
cat("Thorsson columns:", paste(colnames(th), collapse = ", "), "\n"); print(head(th, 3))
idc <- colnames(th)[1]; stc <- grep("Immune", colnames(th), value = TRUE)[1]
mt <- match(pat, substr(th[[idc]], 1, 12)); cl$immune_subtype <- th[[stc]][mt]
## TCGA 2012 paper fields (already in colData) + Xena
xc <- readRDS(file.path(DC, "xena_clinical_coadread.rds"))
cat("Xena MSI-like columns:\n"); print(grep("msi|MSI|CIMP|MLH1|braf|kras|hyper", colnames(xc), value = TRUE))
mx <- match(pat, substr(xc$sampleID, 1, 12))
for (v in c("CIMP", "MSI_updated_Oct62011", "microsatellite_instability", "hypermutation", "AWG_MLH1_silencing", "braf_gene_analysis_result", "kras_mutation_found")) if (v %in% colnames(xc)) cl[[paste0("xena_", v)]] <- xc[[v]][mx]
cl$paper_msi <- as.character(cl$paper_MSI_status)
## ---- MSI definition: concordance of independent sources
tt <- cl[prim, ]
cat("\npaper_MSI_status vs MANTIS>0.6 (primary tumours):\n"); print(table(paper = tt$paper_msi, MANTIS_H = tt$mantis > 0.6, useNA = "ifany"))
cat("MANTIS by paper status:\n"); print(tapply(tt$mantis, tt$paper_msi, summary))
cat("MSIsensor by paper status:\n"); print(tapply(tt$msisensor, tt$paper_msi, median, na.rm = TRUE))
ev <- tt[tt$paper_msi %in% c("MSI-H", "MSI-L", "MSS") & !is.na(tt$mantis), ]
a <- ev$paper_msi == "MSI-H"; b <- ev$mantis > 0.6
po <- mean(a == b); pe <- mean(a) * mean(b) + (1 - mean(a)) * (1 - mean(b)); cat(sprintf("n=%d agreement=%.3f Cohen kappa=%.3f\n", nrow(ev), po, (po - pe) / (1 - pe)))
cl$msi_mantis <- factor(ifelse(is.na(cl$mantis), NA, ifelse(cl$mantis > 0.6, "MSI-H", "MSS/MSI-L")), levels = c("MSS/MSI-L", "MSI-H"))
cl$msi_paper <- factor(ifelse(cl$paper_msi == "MSI-H", "MSI-H", ifelse(cl$paper_msi %in% c("MSS", "MSI-L"), "MSS/MSI-L", NA)), levels = c("MSS/MSI-L", "MSI-H"))
## primary MSI definition: TCGA-2012 expert call where available, otherwise MANTIS (cross-validated above)
cl$msi <- cl$msi_paper; j <- is.na(cl$msi) & !is.na(cl$msi_mantis); cl$msi[j] <- cl$msi_mantis[j]
cl$msi_source <- ifelse(!is.na(cl$msi_paper), "TCGA-2012 call", ifelse(!is.na(cl$msi_mantis), "MANTIS (PanCanAtlas)", NA))
cat("\nMSI final (primary tumours):\n"); print(table(cl$msi[prim], cl$msi_source[prim], useNA = "ifany"))
## ---- survival endpoints (primary = TCGA-CDR OS)
cl$os_time <- cl$cdr_OS_time; cl$os_event <- cl$cdr_OS
cl$os_time_m <- cl$os_time / 30.4375
cat("\nOS availability among primary tumours (CDR):", sum(prim & !is.na(cl$os_time) & cl$os_time > 0 & !is.na(cl$os_event)), "\n")
cat("OS availability (GDC raw, original definition):", sum(prim & !is.na(cl$os_gdc_time) & cl$os_gdc_time > 0 & !is.na(cl$os_gdc_event)), "\n")
cat("events CDR:", sum(cl$os_event[prim] == 1, na.rm = TRUE), " events GDC:", sum(cl$os_gdc_event[prim] == 1, na.rm = TRUE), "\n")
## ---- REMARK-style flow
fl <- list()
fl$total_rnaseq <- ncol(counts); fl$normal <- sum(cl$tissue_type == "Normal"); fl$tumor_all <- sum(tum); fl$metastatic_excl <- sum(tum & !prim); fl$primary <- sum(prim)
ok_os <- prim & !is.na(cl$os_time) & cl$os_time > 0 & !is.na(cl$os_event)
fl$with_os <- sum(ok_os); fl$os_events <- sum(cl$os_event[ok_os] == 1)
fl$no_os <- fl$primary - fl$with_os
fl$cc_age_sex_stage_site <- sum(ok_os & !is.na(cl$age_y) & !is.na(cl$sex) & !is.na(cl$stage_grp))
fl$cc_plus_msi <- sum(ok_os & !is.na(cl$age_y) & !is.na(cl$sex) & !is.na(cl$stage_grp) & !is.na(cl$msi))
fl$cc_original_4var <- sum(!is.na(x$surv_df$age) & !is.na(x$surv_df$stage))
fl$with_wes <- sum(!is.na(cl$tmb[prim])); fl$with_purity_abs <- sum(prim & !is.na(cl$purity_absolute)); fl$with_msi <- sum(prim & !is.na(cl$msi)); fl$with_stage <- sum(prim & !is.na(cl$stage_grp))
print(unlist(fl))
## original surv_df re-check (n=503)
sv <- x$surv_df; cat("\noriginal surv_df: n=", nrow(sv), " events=", sum(sv$status), " in recomputed OS set: ", sum(sv$sample_id %in% rownames(cl)[ok_os]), "\n")
## ---- save
saveRDS(list(counts = counts, gmap = gmap, clin = cl, prim = prim, flow = fl, ok_os = ok_os), file.path(DC, "cohort.rds"))
write_csv_utf8(cbind(sample = rownames(cl), cl[, c("patient", "condition", "sample_type", "site", "sex", "age_y", "stage_grp", "msi", "msi_source", "mantis", "msisensor", "tmb", "purity_absolute", "immune_subtype", "os_time", "os_event", "os_gdc_time", "os_gdc_event", "cdr_PFI", "cdr_PFI_time")]), file.path(RES, "cohort_clinical_harmonised.csv"))
log_msg("done"); sink()
