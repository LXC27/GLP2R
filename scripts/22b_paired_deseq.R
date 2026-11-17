## 22b_paired_deseq.R : paired tumour-vs-normal sensitivity analysis (32 patients with matched primary tumour + normal), DESeq2 ~ patient + condition
source(file.path(Sys.getenv("GLP2R_SCRIPTS", "W:/Re_build/05_Analysis/scripts"), "00_common.R"))
suppressPackageStartupMessages({library(DESeq2); library(SummarizedExperiment)})
sink(file.path(LOGD, "22b_paired_deseq.txt"), split = TRUE)
co <- readRDS(file.path(DC, "cohort.rds")); cl <- co$clin; gmap <- co$gmap
man <- read.csv(file.path(DC, "gdc_paired_manifest.csv"), stringsAsFactors = FALSE); man <- man[man$md5_ok %in% c(TRUE, "TRUE"), ]
rd <- function(f) { x <- read.delim(file.path(DC, "gdc_paired", paste0(f, ".tsv")), comment.char = "#", stringsAsFactors = FALSE); x <- x[grepl("^ENSG", x$gene_id), ]; setNames(x$unstranded, x$gene_id) }
tc <- sapply(man$file_id, rd); colnames(tc) <- man$sample_id; stopifnot(all(rownames(tc) %in% rownames(co$counts)))
nm <- rownames(cl)[cl$sample_type == "Solid Tissue Normal" & cl$patient %in% man$patient]; nc <- co$counts[rownames(tc), nm]
cnt <- cbind(nc, tc); cnt <- round(cnt); storage.mode(cnt) <- "integer"
cd <- data.frame(patient = factor(c(cl[nm, "patient"], man$patient)), condition = factor(rep(c("Normal", "Tumor"), c(length(nm), nrow(man))), levels = c("Normal", "Tumor")), row.names = colnames(cnt))
cat("pairs:", nlevels(cd$patient), " samples:", ncol(cnt), "\n"); print(table(cd$condition)); print(table(cl[nm, "site"]))
keep <- rowSums(cnt >= 10) >= 10; cat("genes tested:", sum(keep), "\n")
dd <- DESeq(DESeqDataSetFromMatrix(cnt[keep, ], cd, ~ patient + condition), quiet = TRUE)
r2 <- as.data.frame(results(dd, contrast = c("condition", "Tumor", "Normal"), alpha = 0.05)); r2$ensembl_full <- rownames(r2); r2$symbol <- gmap$symbol[match(r2$ensembl_full, gmap$ensembl_full)]; r2 <- r2[order(r2$pvalue), ]
saveRDS(r2, file.path(DC, "de_paired.rds"))
vs <- assay(varianceStabilizingTransformation(dd, blind = FALSE)); gid <- gmap$ensembl_full[which(gmap$symbol == "GLP2R")][1]
pg <- data.frame(patient = cd$patient, condition = cd$condition, vst = vs[gid, ], row.names = NULL)
pw <- reshape(pg, idvar = "patient", timevar = "condition", direction = "wide"); pw$delta <- pw$vst.Tumor - pw$vst.Normal; pw$patient <- as.character(pw$patient)
wt <- wilcox.test(pw$vst.Tumor, pw$vst.Normal, paired = TRUE, exact = FALSE)
cat(sprintf("GLP2R paired: median VST normal %.2f, tumour %.2f; n down in tumour = %d/%d; Wilcoxon signed-rank p=%s\n", median(pw$vst.Normal), median(pw$vst.Tumor), sum(pw$delta < 0), nrow(pw), fmt_p(wt$p.value)))
write_csv_utf8(pw, file.path(RES, "paired_glp2r_vst.csv"))
g <- function(d, s) d[which(d$symbol == s), c("symbol", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]
cat("\nPaired model, key genes:\n"); print(do.call(rbind, lapply(c("GLP2R", "GLP1R", "GIP", "GIPR", "GCG", "DPP4", "DNMT1", "DNMT3A", "DNMT3B", "CD274"), function(s) g(r2, s))))
sig <- r2$padj < 0.05 & abs(r2$log2FoldChange) > 1 & !is.na(r2$padj); cat(sprintf("Paired DEGs |log2FC|>1 & FDR<0.05: %d (%d up, %d down)\n", sum(sig), sum(sig & r2$log2FoldChange > 0), sum(sig & r2$log2FoldChange < 0)))
log_msg("done"); sink()
