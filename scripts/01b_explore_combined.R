## 01b_explore_combined.R : explore TCGA_CRC_combined clinical fields needed for the cohort build
options(width = 220, warn = -1)
RB <- Sys.getenv("GLP2R_RB", "W:/Re_build")
if (getRversion() < "4.5") { .libPaths(c(file.path(RB, "_work", "Rlib"), .libPaths())) } else { ul <- Sys.getenv("R_LIBS_USER"); if (nzchar(ul) && dir.exists(ul)) .libPaths(c(ul, .libPaths())) }
suppressPackageStartupMessages({library(SummarizedExperiment); library(DESeq2)})
sink(file.path(RB, "05_Analysis", "logs", "01b_explore_combined.txt"), split = TRUE)
x <- readRDS(file.path(Sys.getenv("GLP2R_ORIG", "V:"), "TCGA_CRC_combined.rds"))
cl <- x$clinical
cat("clinical dim:", dim(cl), "\n")
cat("colnames:\n"); print(colnames(cl))
cat("\nrownames head:\n"); print(head(rownames(cl)))
cat("\ncolnames(counts) head:\n"); print(head(colnames(x$counts)))
cat("rownames(counts) head:\n"); print(head(rownames(x$counts)))
cat("identical(colnames(counts), rownames(clinical))? ", identical(colnames(x$counts), rownames(cl)), "\n")
cat("\ndds class / dim: ", class(x$dds), dim(x$dds), "\n")
cat("dds design: "); print(design(x$dds))
cat("dds colData names:\n"); print(colnames(colData(x$dds)))
cat("rowData names:\n"); print(colnames(rowData(x$dds)))
cat("rowData head:\n"); print(head(as.data.frame(rowData(x$dds))))
cat("\nvsd class: ", class(x$vsd), dim(x$vsd), "\n")
cat("expr_mat head rownames: "); print(head(rownames(x$expr_mat)))
cat("\n--- surv_df ---\n"); print(dim(x$surv_df)); print(colnames(x$surv_df)); print(head(x$surv_df, 3))
for (v in c("sample_type", "shortLetterCode", "definition", "tissue_type", "ajcc_pathologic_stage", "vital_status", "gender", "project_id", "primary_site", "disease_type",
            "MSI_status", "msi_status", "msi", "cohort", "dataset", "tumor_site", "site")) {
  if (v %in% colnames(cl)) { cat("\n#", v, "\n"); print(table(cl[[v]], useNA = "ifany")) }
}
cat("\nColumns matching MSI/MMR/stage/age/gender/purity/site/project:\n")
print(grep("msi|mmr|micro|stage|age|gender|sex|purity|project|site|coad|read|cohort|group|time|days|status|event|subtype|race|ethn", colnames(cl), ignore.case = TRUE, value = TRUE))
cat("\nlast 15 columns head:\n"); print(head(cl[, tail(colnames(cl), 15)], 4))
sink()
