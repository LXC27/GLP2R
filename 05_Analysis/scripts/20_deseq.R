## 20_deseq.R : tumour-vs-normal differential expression (DESeq2) with COAD/READ covariate + paired sensitivity; VST matrices
source(file.path(Sys.getenv("GLP2R_SCRIPTS", "W:/Re_build/05_Analysis/scripts"), "00_common.R"))
suppressPackageStartupMessages({library(DESeq2); library(SummarizedExperiment)})
sink(file.path(LOGD, "20_deseq.txt"), split = TRUE)
co <- readRDS(file.path(DC, "cohort.rds")); cl <- co$clin; counts <- co$counts; gmap <- co$gmap
sel <- cl$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")
cl1 <- cl[sel, ]; cnt <- round(counts[, rownames(cl1)]); storage.mode(cnt) <- "integer"
cl1$condition <- factor(cl1$condition, levels = c("Normal", "Tumor")); cl1$site <- droplevels(cl1$site)
cat("samples:", ncol(cnt), " tumour:", sum(cl1$condition == "Tumor"), " normal:", sum(cl1$condition == "Normal"), "\n"); print(table(cl1$condition, cl1$site))
keep <- rowSums(cnt >= 10) >= 10
cat("genes tested after filter (>=10 counts in >=10 samples):", sum(keep), " of ", nrow(cnt), "\n")
log_msg("DESeq2 main model: ~ site + condition")
dds <- DESeqDataSetFromMatrix(cnt[keep, ], colData = cl1[, c("condition", "site", "patient")], design = ~ site + condition)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("condition", "Tumor", "Normal"), alpha = 0.05)
rdf <- as.data.frame(res); rdf$ensembl_full <- rownames(rdf); rdf$symbol <- gmap$symbol[match(rdf$ensembl_full, gmap$ensembl_full)]
rdf <- rdf[order(rdf$pvalue), ]
saveRDS(rdf, file.path(DC, "de_tumor_normal.rds"))
log_msg("VST (blind = FALSE)")
vs <- assay(varianceStabilizingTransformation(dds, blind = FALSE))
## symbol-level matrices (keep the Ensembl gene with highest mean VST for duplicated symbols)
sym <- gmap$symbol[match(rownames(vs), gmap$ensembl_full)]
ok <- !is.na(sym); vs <- vs[ok, ]; sym <- sym[ok]
o <- order(rowMeans(vs), decreasing = TRUE); vs <- vs[o, ]; sym <- sym[o]; vs <- vs[!duplicated(sym), ]; rownames(vs) <- sym[!duplicated(sym)]
saveRDS(vs, file.path(DC, "vst_all_624.rds"))
## the supplied matrix contains no patient with both a primary tumour and a normal sample (checked here); matched pairs are retrieved from GDC in 22_paired_gdc.R
pp <- intersect(cl1$patient[cl1$condition == "Normal"], cl1$patient[cl1$condition == "Tumor"]); cat("patients with both tumour and normal in the supplied matrix:", length(pp), "\n")
## summary
g <- function(d, s) d[which(d$symbol == s), c("symbol", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]
cat("\nMain model, key genes:\n"); print(do.call(rbind, lapply(c("GLP2R", "GLP1R", "GIP", "GIPR", "GCG", "DPP4", "DNMT1", "DNMT3A", "DNMT3B", "TET2", "CD274"), function(s) g(rdf, s))))
sig <- rdf$padj < 0.05 & abs(rdf$log2FoldChange) > 1 & !is.na(rdf$padj)
cat(sprintf("\nDEGs |log2FC|>1 & FDR<0.05: %d total (%d up, %d down)\n", sum(sig), sum(sig & rdf$log2FoldChange > 0), sum(sig & rdf$log2FoldChange < 0)))
cat("Rank of GLP2R by padj in main model:", which(rdf$symbol == "GLP2R")[1], " rank by |stat|: ", rank(-abs(rdf$stat), na.last = "keep")[which(rdf$symbol == "GLP2R")[1]], "\n")
log_msg("done"); sink()
