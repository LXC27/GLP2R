## 40_immune.R : GLP2R and the tumour immune microenvironment (purity-aware; BH-FDR; Charoentong 28 cell types, ESTIMATE, MCP-counter, checkpoints, Thorsson subtypes)
source(file.path(Sys.getenv("GLP2R_SCRIPTS", "W:/Re_build/05_Analysis/scripts"), "00_common.R"))
suppressPackageStartupMessages({library(GSVA); library(tidyestimate)})
sink(file.path(LOGD, "40_immune.txt"), split = TRUE)
ex <- readRDS(file.path(RES, "expr_tests.rds")); d <- ex$d; vt <- readRDS(file.path(DC, "vst_tumor_591.rds")); gmap <- readRDS(file.path(DC, "cohort.rds"))$gmap
stopifnot(all(rownames(d) %in% colnames(vt))); vt <- vt[, rownames(d)]
RM <- list(); add <- function(k, v, desc = "") RM[[length(RM) + 1]] <<- data.frame(key = k, value = as.character(v), description = desc)
## ---- ESTIMATE (expression) and ABSOLUTE (genomic) purity
est <- tryCatch({ df <- as.data.frame(vt); df <- filter_common_genes(df, id = "hgnc_symbol", tell_missing = FALSE, find_alias = TRUE); estimate_score(df, is_affymetrix = FALSE) }, error = function(e) { cat("ESTIMATE failed:", conditionMessage(e), "\n"); NULL })
if (!is.null(est)) { print(head(est, 3)); est <- as.data.frame(est); rownames(est) <- est[[1]]; d$stromal <- est[rownames(d), "stromal"]; d$immune <- est[rownames(d), "immune"]; d$estimate <- est[rownames(d), "estimate"]
  d$purity_est <- cos(0.6049872018 + 0.0001467884 * d$estimate) }
cat("purity: ABSOLUTE n =", sum(!is.na(d$purity_absolute)), " ESTIMATE-derived n =", sum(!is.na(d$purity_est)), "; Spearman between them =", round(cor(d$purity_absolute, d$purity_est, method = "spearman", use = "pair"), 3), "\n")
add("purity_abs_vs_est_rho", round(cor(d$purity_absolute, d$purity_est, method = "spearman", use = "pair"), 3)); add("purity_abs_n", sum(!is.na(d$purity_absolute)))
add("glp2r_vs_purity_abs_rho", round(cor(d$glp2r, d$purity_absolute, method = "spearman", use = "pair"), 3), "GLP2R expression vs ABSOLUTE purity")
add("glp2r_vs_purity_est_rho", round(cor(d$glp2r, d$purity_est, method = "spearman", use = "pair"), 3))
## ---- partial Spearman: rank-residual correlation
pcor_rank <- function(x, y, Z) { ok <- complete.cases(x, y, Z); x <- rank(x[ok]); y <- rank(y[ok]); Z <- Z[ok, , drop = FALSE]
  Zm <- model.matrix(~ ., data = as.data.frame(Z))[, -1, drop = FALSE]; rx <- resid(lm(x ~ Zm)); ry <- resid(lm(y ~ Zm)); r <- cor(rx, ry); n <- sum(ok); k <- ncol(Zm)
  tt <- r * sqrt((n - 2 - k) / (1 - r^2)); c(rho = r, p = 2 * pt(-abs(tt), n - 2 - k), n = n) }
Zadj <- d[, c("purity_absolute", "msi", "site")]
## ---- gene sets
gmt <- readLines(file.path(DC, "Charoentong2017.gmt")); gs <- lapply(strsplit(gmt, "\t"), function(f) unique(f[-(1:2)])); names(gs) <- trimws(gsub("\\s+", " ", sub("memeory", "memory", sapply(strsplit(gmt, "\t"), `[`, 1))))
gs <- lapply(gs, function(g) intersect(g, rownames(vt))); cat("Charoentong sets: ", length(gs), "; mean genes found: ", round(mean(sapply(gs, length)), 1), "\n")
ss <- gsva(ssgseaParam(as.matrix(vt), gs, alpha = 0.25, normalize = TRUE), verbose = FALSE); cat("ssGSEA matrix:", dim(ss), "\n")
## ---- GLP2R vs 28 cell types: Spearman, partial (purity+MSI+site), high/low
g <- d$glp2r; hi <- g > median(g)
resI <- do.call(rbind, lapply(rownames(ss), function(ct) { y <- ss[ct, rownames(d)]; sc <- suppressWarnings(cor.test(g, y, method = "spearman", exact = FALSE)); pc <- pcor_rank(g, y, Zadj)
  pe <- if ("purity_est" %in% colnames(d)) pcor_rank(g, y, d[, c("purity_est", "msi", "site")]) else c(rho = NA, p = NA, n = NA)
  data.frame(cell_type = ct, rho = sc$estimate, p = sc$p.value, partial_rho = pc["rho"], partial_p = pc["p"], partial_n = pc["n"], partial_rho_est = pe["rho"], partial_p_est = pe["p"],
             median_high = median(y[hi]), median_low = median(y[!hi]), auc_high_vs_low = auc_mw(y[hi], y[!hi]), p_wilcox = wilcox.test(y[hi], y[!hi])$p.value, row.names = NULL) }))
for (v in c("p", "partial_p", "partial_p_est", "p_wilcox")) resI[[paste0("FDR_", v)]] <- p.adjust(resI[[v]], "BH")
resI <- resI[order(resI$p), ]; print(resI[, c("cell_type", "rho", "FDR_p", "partial_rho", "FDR_partial_p", "partial_rho_est", "FDR_partial_p_est", "auc_high_vs_low", "FDR_p_wilcox")], digits = 3)
write_csv_utf8(resI, file.path(RES, "immune_celltypes_charoentong.csv"))
add("n_celltypes", nrow(resI)); add("n_sig_unadj_rho_FDR05", sum(resI$FDR_p < 0.05)); add("n_sig_partial_FDR05", sum(resI$FDR_partial_p < 0.05)); add("n_sig_partial_est_FDR05", sum(resI$FDR_partial_p_est < 0.05, na.rm = TRUE)); add("n_sig_highlow_FDR05", sum(resI$FDR_p_wilcox < 0.05))
add("celltypes_sig_partial", paste(resI$cell_type[resI$FDR_partial_p < 0.05], collapse = "; ")); add("celltypes_sig_unadj", paste(resI$cell_type[resI$FDR_p < 0.05], collapse = "; "))
add("max_abs_rho", round(max(abs(resI$rho)), 3)); add("max_abs_partial_rho", round(max(abs(resI$partial_rho)), 3))
## ---- stromal / immune / ESTIMATE scores and MCP-counter populations
mcp <- read.delim(file.path(DC, "MCPcounter_genes.txt"), check.names = FALSE, stringsAsFactors = FALSE); colnames(mcp)[1:2] <- c("gene", "pop")
mcps <- sapply(split(mcp$gene, mcp$pop), function(gg) { gg <- intersect(gg, rownames(vt)); if (length(gg) < 2) rep(NA_real_, ncol(vt)) else colMeans(vt[gg, , drop = FALSE]) }); rownames(mcps) <- colnames(vt); cat("MCP-counter populations scored:", paste(colnames(mcps)[colSums(!is.na(mcps)) > 0], collapse = ", "), "\n")
sc_cols <- cbind(d[, intersect(c("stromal", "immune", "estimate"), colnames(d))], mcps[rownames(d), , drop = FALSE]); sc_res <- do.call(rbind, lapply(colnames(sc_cols), function(v) { y <- sc_cols[[v]]; if (sum(is.finite(y)) < 30) return(NULL); a <- cor.test(g, y, method = "spearman", exact = FALSE); pc <- pcor_rank(g, y, Zadj)
  data.frame(score = v, rho = a$estimate, p = a$p.value, partial_rho = pc["rho"], partial_p = pc["p"], row.names = NULL) }))
sc_res$FDR_p <- p.adjust(sc_res$p, "BH"); sc_res$FDR_partial_p <- p.adjust(sc_res$partial_p, "BH"); print(sc_res, digits = 3); write_csv_utf8(sc_res, file.path(RES, "immune_scores_mcp_estimate.csv"))
for (i in seq_len(nrow(sc_res))) if (sc_res$score[i] %in% c("stromal", "immune", "estimate", "Fibroblasts", "Endothelial cells", "CD8 T cells", "T cells", "Myeloid dendritic cells", "NK cells", "Cytotoxic lymphocytes"))
  add(paste0("rho_", gsub(" ", "_", sc_res$score[i])), sprintf("rho=%.3f FDR=%s; partial rho=%.3f FDR=%s", sc_res$rho[i], fmt_p(sc_res$FDR_p[i]), sc_res$partial_rho[i], fmt_p(sc_res$FDR_partial_p[i])))
## ---- immune checkpoints (Pearson + Spearman, BH)
ck <- c(PD1 = "PDCD1", `PD-L1` = "CD274", `PD-L2` = "PDCD1LG2", `CTLA-4` = "CTLA4", `LAG-3` = "LAG3", `TIM-3` = "HAVCR2", TIGIT = "TIGIT", IDO1 = "IDO1")
ckr <- do.call(rbind, lapply(names(ck), function(n) { y <- vt[ck[n], rownames(d)]; a <- cor.test(g, y, method = "pearson"); b <- cor.test(g, y, method = "spearman", exact = FALSE); pc <- pcor_rank(g, y, Zadj)
  data.frame(Gene = ck[n], Name = n, Pearson_r = a$estimate, Pearson_p = a$p.value, Spearman_rho = b$estimate, Spearman_p = b$p.value, partial_rho = pc["rho"], partial_p = pc["p"], row.names = NULL) }))
ckr$Pearson_FDR <- p.adjust(ckr$Pearson_p, "BH"); ckr$Spearman_FDR <- p.adjust(ckr$Spearman_p, "BH"); ckr$partial_FDR <- p.adjust(ckr$partial_p, "BH"); print(ckr, digits = 3)
write_csv_utf8(ckr, file.path(RES, "TableS5_checkpoints_recomputed.csv")); add("checkpoint_n_sig_spearman_FDR05", sum(ckr$Spearman_FDR < 0.05)); add("checkpoint_n_sig_pearson_FDR05", sum(ckr$Pearson_FDR < 0.05)); add("checkpoint_n_sig_partial_FDR05", sum(ckr$partial_FDR < 0.05))
add("checkpoint_max_abs_spearman", round(max(abs(ckr$Spearman_rho)), 3))
## ---- T-cell inflamed GEP and IFN-gamma signature (Ayers 2017), CD40LG
gep <- intersect(c("CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6", "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT"), rownames(vt))
ifn <- intersect(c("IDO1", "CXCL10", "CXCL9", "HLA-DRA", "STAT1", "IFNG"), rownames(vt)); d$gep <- colMeans(vt[gep, rownames(d)]); d$ifng <- colMeans(vt[ifn, rownames(d)]); d$cd40lg <- vt["CD40LG", rownames(d)]
for (v in c("gep", "ifng", "cd40lg")) { a <- cor.test(g, d[[v]], method = "spearman", exact = FALSE); pc <- pcor_rank(g, d[[v]], Zadj); add(paste0("rho_", v), sprintf("rho=%.3f p=%s; partial rho=%.3f p=%s", a$estimate, fmt_p(a$p.value), pc["rho"], fmt_p(pc["p"]))) }
add("gep_genes_found", length(gep)); add("gep_hilo_p", signif(wilcox.test(d$gep ~ hi)$p.value, 4)); cat("GEP: genes", length(gep), " rho vs GLP2R:", round(cor(g, d$gep, method = "spearman"), 3), "\n")
## ---- Thorsson immune subtypes
th <- d[!is.na(d$immune_subtype), ]; th$sub <- sub(".*\\((Immune )?(C[0-9]).*", "\\2", th$immune_subtype); cat("\nThorsson subtype counts:\n"); print(table(th$sub))
kw <- kruskal.test(glp2r ~ sub, th); cat("Kruskal-Wallis across all subtypes present: p =", fmt_p(kw$p.value), "\n"); print(tapply(th$glp2r, th$sub, median))
add("thorsson_n", nrow(th)); add("thorsson_counts", paste(names(table(th$sub)), table(th$sub), sep = "=", collapse = "; ")); add("thorsson_kw_p", signif(kw$p.value, 4))
pw <- pairwise.wilcox.test(th$glp2r, th$sub, p.adjust.method = "BH"); print(pw$p.value, digits = 3)
th$msi_s <- th$msi; mm <- th[!is.na(th$msi), ]; fit_s <- lm(rank(glp2r) ~ msi + site + sub, mm); fit_0 <- lm(rank(glp2r) ~ msi + site, mm); add("thorsson_adj_msi_site_anova_p", signif(anova(fit_0, fit_s)$`Pr(>F)`[2], 4), "rank-based partial F-test for subtype after MSI and site (subtype entered last)")
mss <- th[which(th$msi == "MSS/MSI-L"), ]; if (length(unique(mss$sub)) > 1) add("thorsson_kw_p_MSS_only", signif(kruskal.test(glp2r ~ sub, mss)$p.value, 4))
## ---- MSS-only stratified correlations (key cell types)
mss_d <- d[which(d$msi == "MSS/MSI-L"), ]; strat <- do.call(rbind, lapply(c("Activated dendritic cell", "Immature dendritic cell", "Macrophage", "Activated CD8 T cell", "Regulatory T cell", "Natural killer cell"), function(ct) {
  y <- ss[ct, rownames(mss_d)]; a <- cor.test(mss_d$glp2r, y, method = "spearman", exact = FALSE); y2 <- ss[ct, rownames(d)[which(d$msi == "MSI-H")]]; b <- cor.test(d$glp2r[which(d$msi == "MSI-H")], y2, method = "spearman", exact = FALSE)
  data.frame(cell_type = ct, MSS_rho = a$estimate, MSS_p = a$p.value, MSS_n = nrow(mss_d), MSIH_rho = b$estimate, MSIH_p = b$p.value, MSIH_n = length(y2), row.names = NULL) })); print(strat, digits = 3)
write_csv_utf8(strat, file.path(RES, "immune_msi_stratified.csv"))
## ---- TMB
tm <- d[!is.na(d$tmb), ]; a <- cor.test(tm$glp2r, log10(tm$tmb + 0.1), method = "spearman", exact = FALSE); mt <- tm[tm$msi == "MSS/MSI-L" & !is.na(tm$msi), ]; b <- cor.test(mt$glp2r, log10(mt$tmb + 0.1), method = "spearman", exact = FALSE)
cat(sprintf("TMB: n=%d rho=%.3f p=%s; MSS-only n=%d rho=%.3f p=%s\n", nrow(tm), a$estimate, fmt_p(a$p.value), nrow(mt), b$estimate, fmt_p(b$p.value)))
add("tmb_rho_all", sprintf("rho=%.3f p=%s n=%d", a$estimate, fmt_p(a$p.value), nrow(tm))); add("tmb_rho_mss", sprintf("rho=%.3f p=%s n=%d", b$estimate, fmt_p(b$p.value), nrow(mt)))
saveRDS(list(d = d, ss = ss, resI = resI, sc_res = sc_res, ckr = ckr, th = th, strat = strat), file.path(RES, "immune_objects.rds"))
rm_write(do.call(rbind, RM), "immune"); log_msg("done"); sink()
