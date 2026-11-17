## 01b2_explore_fields.R : key categorical fields + Xena barcode matching
options(width = 220, warn = -1)
RB <- Sys.getenv("GLP2R_RB", "W:/Re_build")
if (getRversion() < "4.5") { .libPaths(c(file.path(RB, "_work", "Rlib"), .libPaths())) } else { ul <- Sys.getenv("R_LIBS_USER"); if (nzchar(ul) && dir.exists(ul)) .libPaths(c(ul, .libPaths())) }
suppressPackageStartupMessages({library(SummarizedExperiment)})
sink(file.path(RB, "05_Analysis", "logs", "01b2_explore_fields.txt"), split = TRUE)
x <- readRDS(file.path(Sys.getenv("GLP2R_ORIG", "V:"), "TCGA_CRC_combined.rds"))
cl <- x$clinical
cat("class of each column (non-atomic only):\n")
print(names(Filter(function(z) !is.atomic(z), as.list(cl))))
tb <- function(v) { cat("\n#", v, "\n"); r <- tryCatch(print(table(as.character(cl[[v]]), useNA = "ifany")), error = function(e) cat("ERR", conditionMessage(e), "\n")) }
for (v in c("primary_site", "disease_type", "tumor_descriptor", "paper_MSI_status", "paper_methylation_subtype", "paper_MLH1_silencing", "paper_expression_subtype",
            "paper_hypermutated", "paper_tumor_site", "paper_anatomic_organ_subdivision", "paper_tumor_stage", "paper_vital_status", "paper_cancer", "paper_histological_type",
            "sample_group", "cancer_type", "paper_sequenced", "paper_RNAseq", "paper_methylation", "paper_SNP6", "paper_WGS", "paper_complete",
            "ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_pathologic_m", "classification_of_tumor", "prior_treatment", "prior_malignancy", "synchronous_malignancy", "gender", "paper_gender")) tb(v)
cat("\nage_at_index summary:\n"); print(summary(as.numeric(cl$age_at_index)))
cat("age_at_diagnosis summary (days):\n"); print(summary(as.numeric(cl$age_at_diagnosis)))
cat("days_to_death summary:\n"); print(summary(as.numeric(cl$days_to_death)))
cat("days_to_last_follow_up summary:\n"); print(summary(as.numeric(cl$days_to_last_follow_up)))
cat("\nunique patients (tumors):", length(unique(cl$patient[cl$tissue_type == "Tumor"])), " tumors:", sum(cl$tissue_type == "Tumor"), "\n")
cat("duplicate patients among tumors:\n"); print(table(table(cl$patient[cl$tissue_type == "Tumor"])))
cat("\nsample_group vs project:\n"); print(table(cl$sample_group, cl$project_id))
sv <- x$surv_df
cat("\nsurv_df stage counts:\n"); print(table(sv$stage, useNA = "ifany"))
cat("surv_df MSI counts:\n"); print(table(sv$MSI, useNA = "ifany"))
cat("surv_df status:\n"); print(table(sv$status))
cat("surv_df time summary:\n"); print(summary(sv$time))
cat("surv_df age summary:\n"); print(summary(sv$age))
## Xena
xc <- readRDS(file.path(RB, "05_Analysis", "data_cache", "xena_clinical_coadread.rds"))
cat("\nXena clinical class/dim:", class(xc), dim(xc), "\n")
print(head(colnames(xc), 5)); print(head(xc[, 1:3], 3))
cat("Xena sampleID example:", head(xc$sampleID, 3), "\n")
pat_xena <- substr(xc$sampleID, 1, 12)
pat_loc <- unique(cl$patient)
cat("Xena patients:", length(unique(pat_xena)), " matches to local patients:", sum(unique(pat_xena) %in% pat_loc), " of ", length(pat_loc), "\n")
sink()
