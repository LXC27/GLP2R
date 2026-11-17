## 21_expr_tests.R : GLP2R expression statistics (tumour vs normal, stage, MSI, site, GLP1R) + tumour-only VST
source(file.path(Sys.getenv("GLP2R_SCRIPTS", "W:/Re_build/05_Analysis/scripts"), "00_common.R"))
suppressPackageStartupMessages({library(DESeq2)})
sink(file.path(LOGD, "21_expr_tests.txt"), split = TRUE)
co <- readRDS(file.path(DC, "cohort.rds")); cl <- co$clin; counts <- co$counts; gmap <- co$gmap; prim <- co$prim
rdf <- readRDS(file.path(DC, "de_tumor_normal.rds")); r2 <- readRDS(file.path(DC, "de_paired.rds")); vsa <- readRDS(file.path(DC, "vst_all_624.rds"))
## ---- tumour-only VST (blind) for all tumour-internal analyses
tn <- rownames(cl)[prim]; cn <- round(counts[, tn]); storage.mode(cn) <- "integer"; cn <- cn[rowSums(cn >= 10) >= 10, ]
log_msg("tumour-only VST on", nrow(cn), "genes x", ncol(cn), "samples")
vt <- to_symbol_matrix(DESeq2::vst(cn, blind = TRUE, nsub = 2000), gmap); saveRDS(vt, file.path(DC, "vst_tumor_591.rds"))
## also log2(CPM+1) tumour matrix (sensitivity measure for GLP2R)
lib <- colSums(counts[, tn]); lcpm <- to_symbol_matrix(log2(sweep(counts[, tn], 2, lib / 1e6, "/") + 1), gmap); saveRDS(lcpm, file.path(DC, "log2cpm_tumor_591.rds"))
RM <- list(); add <- function(k, v, d = "") RM[[length(RM) + 1]] <<- data.frame(key = k, value = as.character(v), description = d)
## ---- 1. tumour vs normal
g <- "GLP2R"
tum <- rownames(cl)[cl$sample_type == "Primary Tumor"]; nor <- rownames(cl)[cl$sample_type == "Solid Tissue Normal"]
xt <- vsa[g, tum]; xn <- vsa[g, nor]
w <- wilcox.test(xt, xn); auc <- auc_mw(xn, xt)
set.seed(SEED); bs <- replicate(2000, auc_mw(sample(xn, replace = TRUE), sample(xt, replace = TRUE))); aci <- quantile(bs, c(.025, .975))
cat(sprintf("GLP2R VST median tumour %.2f vs normal %.2f; Wilcoxon p=%s; AUC(normal>tumour)=%.3f (%.3f-%.3f)\n", median(xt), median(xn), fmt_p(w$p.value), auc, aci[1], aci[2]))
pw <- read.csv(file.path(RES, "paired_glp2r_vst.csv")); pp <- pw$patient   # 32 matched pairs retrieved from GDC (22a/22b); the supplied matrix has none
wp <- wilcox.test(pw$vst.Tumor, pw$vst.Normal, paired = TRUE, exact = FALSE); cat("paired n =", length(pp), " paired Wilcoxon p =", fmt_p(wp$p.value), " median paired diff (VST):", round(median(pw$delta), 2), " fraction pairs T<N:", mean(pw$delta < 0), "\n")
rawc <- counts[gmap$ensembl_full[match(g, gmap$symbol)], ]
det_t <- mean(rawc[tum] >= 10); det_n <- mean(rawc[nor] >= 10)
cat(sprintf("Detection (>=10 raw counts): tumour %.1f%% normal %.1f%%; median raw counts tumour %.0f normal %.0f\n", 100 * det_t, 100 * det_n, median(rawc[tum]), median(rawc[nor])))
gm <- rdf[which(rdf$symbol == g)[1], ]; gp <- r2[which(r2$symbol == g)[1], ]
add("glp2r_log2fc", signif(gm$log2FoldChange, 4), "DESeq2 (~site+condition) tumour vs normal"); add("glp2r_lfcSE", signif(gm$lfcSE, 4)); add("glp2r_padj", signif(gm$padj, 4)); add("glp2r_pvalue", signif(gm$pvalue, 4)); add("glp2r_baseMean", signif(gm$baseMean, 4))
add("glp2r_log2fc_paired", signif(gp$log2FoldChange, 4), "paired DESeq2 (~patient+condition)"); add("glp2r_padj_paired", signif(gp$padj, 4)); add("glp2r_pvalue_paired", signif(gp$pvalue, 4)); add("n_pairs", length(pp))
add("glp2r_auc_tumor_vs_normal", sprintf("%.3f (%.3f-%.3f)", auc, aci[1], aci[2]), "AUC normal>tumour on VST"); add("glp2r_wilcox_p", signif(w$p.value, 4))
add("glp2r_paired_wilcox_p", signif(wp$p.value, 4)); add("glp2r_paired_n_down", sum(pw$delta < 0), "pairs in which tumour GLP2R VST is lower than matched normal"); add("glp2r_detect_tumour_pct", round(100 * det_t, 1)); add("glp2r_detect_normal_pct", round(100 * det_n, 1))
add("glp2r_median_vst_tumour", round(median(xt), 2)); add("glp2r_median_vst_normal", round(median(xn), 2))
for (s in c("GLP1R", "GIP", "GIPR", "GCG", "DPP4")) { d <- rdf[which(rdf$symbol == s)[1], ]; dp <- r2[which(r2$symbol == s)[1], ]
  add(paste0(tolower(s), "_baseMean"), signif(d$baseMean, 4)); add(paste0(tolower(s), "_log2fc"), signif(d$log2FoldChange, 4)); add(paste0(tolower(s), "_pvalue"), signif(d$pvalue, 4)); add(paste0(tolower(s), "_padj"), signif(d$padj, 4))
  add(paste0(tolower(s), "_padj_paired"), signif(dp$padj, 4)) }
sig <- !is.na(rdf$padj) & rdf$padj < 0.05 & abs(rdf$log2FoldChange) > 1
add("deg_total", sum(sig)); add("deg_up", sum(sig & rdf$log2FoldChange > 0)); add("deg_down", sum(sig & rdf$log2FoldChange < 0)); add("genes_tested", nrow(rdf))
sig2 <- !is.na(r2$padj) & r2$padj < 0.05 & abs(r2$log2FoldChange) > 1; add("deg_paired_total", sum(sig2)); add("deg_paired_up", sum(sig2 & r2$log2FoldChange > 0)); add("deg_paired_down", sum(sig2 & r2$log2FoldChange < 0))
add("glp2r_rank_by_abs_stat", rank(-abs(rdf$stat), na.last = "keep")[which(rdf$symbol == g)[1]], "rank among all tested genes"); add("glp2r_in_top_down_logfc_pct", round(100 * mean(rdf$log2FoldChange[!is.na(rdf$padj)] > gm$log2FoldChange), 2), "% of tested genes with LFC higher than GLP2R (i.e. GLP2R percentile from bottom)")
## ---- 2. tumour-internal: stage, MSI, site, covariates   (tumour-only VST `vt`)
d <- cl[tn, ]; d$glp2r <- vt[g, tn]; d$glp1r <- vt["GLP1R", tn]; d$glp2r_lcpm <- lcpm[g, tn]
d$stage_num <- as.numeric(d$stage_grp)
ds <- d[!is.na(d$stage_grp), ]
kw <- kruskal.test(glp2r ~ stage_grp, ds); sp <- cor.test(ds$glp2r, ds$stage_num, method = "spearman", exact = FALSE)
cat(sprintf("\nStage: n=%d KW p=%s; Spearman rho=%.3f p=%s\n", nrow(ds), fmt_p(kw$p.value), sp$estimate, fmt_p(sp$p.value))); print(table(ds$stage_grp))
print(tapply(ds$glp2r, ds$stage_grp, median))
lm_st <- lm(glp2r ~ stage_num + site + msi, ds); cat("lm GLP2R ~ stage(trend)+site+MSI:\n"); print(round(summary(lm_st)$coefficients, 4))
add("stage_n", nrow(ds)); add("stage_kw_p", signif(kw$p.value, 4)); add("stage_spearman_rho", round(sp$estimate, 3)); add("stage_spearman_p", signif(sp$p.value, 4)); add("stage_adj_trend_p", signif(summary(lm_st)$coefficients["stage_num", 4], 4), "adjusted for site and MSI")
dm <- d[!is.na(d$msi), ]; wm <- wilcox.test(glp2r ~ msi, dm); aucm <- auc_mw(dm$glp2r[dm$msi == "MSI-H"], dm$glp2r[dm$msi == "MSS/MSI-L"])
cat(sprintf("\nMSI: n=%d (MSI-H %d); Wilcoxon p=%s; AUC(MSI-H > MSS)=%.3f; medians %s\n", nrow(dm), sum(dm$msi == "MSI-H"), fmt_p(wm$p.value), aucm, paste(round(tapply(dm$glp2r, dm$msi, median), 2), collapse = " / ")))
lm_m <- lm(glp2r ~ msi + site + age_y + sex + stage_grp, dm); cat("lm GLP2R ~ MSI + site + age + sex + stage:\n"); print(round(summary(lm_m)$coefficients, 4))
for (s in c("COAD", "READ")) { z <- dm[dm$site == s, ]; if (length(unique(z$msi)) == 2) cat(sprintf("  %s: n=%d MSI-H=%d p=%s\n", s, nrow(z), sum(z$msi == "MSI-H"), fmt_p(wilcox.test(glp2r ~ msi, z)$p.value))) }
add("msi_n", nrow(dm)); add("msi_h_n", sum(dm$msi == "MSI-H")); add("msi_wilcox_p", signif(wm$p.value, 4)); add("msi_auc", round(aucm, 3), "AUC MSI-H > MSS/MSI-L"); add("msi_median_diff_vst", round(median(dm$glp2r[dm$msi == "MSI-H"]) - median(dm$glp2r[dm$msi == "MSS/MSI-L"]), 2))
add("msi_adj_beta", round(summary(lm_m)$coefficients["msiMSI-H", 1], 3), "VST units, adjusted"); add("msi_adj_p", signif(summary(lm_m)$coefficients["msiMSI-H", 4], 4))
## MSS-only stage trend (MSI confounding check)
dmss <- ds[which(ds$msi == "MSS/MSI-L"), ]; spm <- cor.test(dmss$glp2r, dmss$stage_num, method = "spearman", exact = FALSE); cat(sprintf("Stage trend within MSS/MSI-L: n=%d rho=%.3f p=%s\n", nrow(dmss), spm$estimate, fmt_p(spm$p.value)))
add("stage_trend_mss_rho", round(spm$estimate, 3)); add("stage_trend_mss_p", signif(spm$p.value, 4))
wsite <- wilcox.test(glp2r ~ site, d); cat("Site COAD vs READ: p=", fmt_p(wsite$p.value), " medians ", paste(round(tapply(d$glp2r, d$site, median), 2), collapse = "/"), "\n"); add("site_wilcox_p", signif(wsite$p.value, 4))
add("age_spearman_rho", round(cor(d$glp2r, d$age_y, method = "spearman", use = "pair"), 3)); add("sex_wilcox_p", signif(wilcox.test(glp2r ~ sex, d)$p.value, 4))
## GLP1R vs GLP2R
cc <- cor.test(d$glp1r, d$glp2r, method = "spearman", exact = FALSE); set.seed(SEED); bb <- replicate(2000, { i <- sample.int(nrow(d), replace = TRUE); cor(d$glp1r[i], d$glp2r[i], method = "spearman") })
cat(sprintf("GLP1R~GLP2R tumour Spearman rho=%.3f (%.3f-%.3f) p=%s\n", cc$estimate, quantile(bb, .025), quantile(bb, .975), fmt_p(cc$p.value)))
add("glp1r_glp2r_rho", round(cc$estimate, 3)); add("glp1r_glp2r_rho_lo", round(quantile(bb, .025), 3)); add("glp1r_glp2r_rho_hi", round(quantile(bb, .975), 3)); add("glp1r_glp2r_p", signif(cc$p.value, 4))
add("glp1r_median_raw_count_tumour", median(counts[gmap$ensembl_full[match("GLP1R", gmap$symbol)], tn]), "very low expression")
add("glp1r_detect_pct_tumour", round(100 * mean(counts[gmap$ensembl_full[match("GLP1R", gmap$symbol)], tn] >= 10), 1)); add("glp1r_detect_pct_normal", round(100 * mean(counts[gmap$ensembl_full[match("GLP1R", gmap$symbol)], nor] >= 10), 1))
add("gip_detect_pct_tumour", round(100 * mean(counts[gmap$ensembl_full[match("GIP", gmap$symbol)], tn] >= 10), 1))
## ---- 3. Table S1 (50 genes) with raw P and FDR from both models
s1 <- read.csv(file.path(Sys.getenv("GLP2R_INPUT", "W:/Supplementary"), "Table_S1_GLP_genes.csv"), fileEncoding = "UTF-8"); genes <- s1$Gene
tabS1 <- do.call(rbind, lapply(genes, function(s) { a <- rdf[which(rdf$symbol == s)[1], ]; b <- r2[which(r2$symbol == s)[1], ]
  data.frame(Gene = s, baseMean = a$baseMean, log2FC = a$log2FoldChange, lfcSE = a$lfcSE, P = a$pvalue, FDR = a$padj, log2FC_paired = b$log2FoldChange, FDR_paired = b$padj, stringsAsFactors = FALSE) }))
tabS1$Expressed <- ifelse(is.na(tabS1$baseMean), "not tested (<10 counts in >=10 samples not met)", ifelse(tabS1$baseMean >= 10, "yes", "low (baseMean<10)"))
tabS1$mean_raw_count <- sapply(tabS1$Gene, function(s) { i <- gmap$ensembl_full[which(gmap$symbol == s)[1]]; if (is.na(i)) NA else mean(counts[i, ]) }); tabS1 <- tabS1[order(tabS1$P), ]
write_csv_utf8(tabS1, file.path(RES, "TableS1_candidate_genes_DE.csv")); cat("\nTable S1 genes:", nrow(tabS1), " lowly expressed (baseMean<10):", paste(tabS1$Gene[which(tabS1$baseMean < 10)], collapse = ", "), "; not tested:", paste(tabS1$Gene[is.na(tabS1$baseMean)], collapse = ", "), "\n")
add("tableS1_n_genes", nrow(tabS1)); add("tableS1_n_FDR05", sum(tabS1$FDR < 0.05, na.rm = TRUE))
rank1 <- which(tabS1$Gene == "GLP2R"); add("glp2r_rank_in_tableS1_by_P", rank1)
saveRDS(list(d = d, tabS1 = tabS1), file.path(RES, "expr_tests.rds"))
rm_write(do.call(rbind, RM), "expr"); log_msg("done"); sink()
