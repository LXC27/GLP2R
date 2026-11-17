## 60_epi_cnv.R : GLP2R promoter methylation (HM450, TCGA-COADREAD) and copy number (GISTIC2 thresholded) vs expression; pre-specified decision rule for the "epigenetic silencing" wording
source(file.path(Sys.getenv("GLP2R_SCRIPTS", "W:/Re_build/05_Analysis/scripts"), "00_common.R"))
sink(file.path(LOGD, "60_epi_cnv.txt"), split = TRUE)
ex <- readRDS(file.path(RES, "expr_tests.rds")); d <- ex$d; cl <- readRDS(file.path(DC, "cohort.rds"))$clin; vsa <- readRDS(file.path(DC, "vst_all_624.rds"))
h <- readRDS(file.path(DC, "hm450_subset.rds")); B <- h$beta; pr <- h$probes; gi <- readRDS(file.path(DC, "gistic_thresholded_coadread.rds"))
RM <- list(); add <- function(k, v, desc = "") RM[[length(RM) + 1]] <<- data.frame(key = k, value = as.character(v), description = desc)
cat("HM450 subset:", dim(B), "\n"); print(table(substr(colnames(B), 14, 15)))
d$pid <- substr(rownames(d), 1, 12); typ <- substr(colnames(B), 14, 15); pid_b <- substr(colnames(B), 1, 12)
bt <- B[, typ == "01", drop = FALSE]; colnames(bt) <- substr(colnames(bt), 1, 12); bt <- bt[, !duplicated(colnames(bt)), drop = FALSE]
bn <- B[, typ == "11", drop = FALSE]; colnames(bn) <- substr(colnames(bn), 1, 12)
cat("methylation tumours:", ncol(bt), " normals:", ncol(bn), " tumours also in RNA-seq primary set:", sum(colnames(bt) %in% d$pid), "\n"); add("hm450_n_tumour", ncol(bt)); add("hm450_n_normal", ncol(bn)); add("hm450_n_tumour_with_rnaseq", sum(colnames(bt) %in% d$pid))
## ---- GLP2R probes (pre-specified: TSS-proximal window chr17:9,727,368-9,729,368, hg19)
g2 <- pr[pr$chrom == "chr17" & pr$chromStart >= 9727368 & pr$chromStart <= 9729368, ]; cat("\nGLP2R-window probes in the subset:\n"); print(g2)
g2b <- pr[grepl("GLP2R", pr$gene), ]; cat("probes annotated to GLP2R:", paste(g2b$id, collapse = ","), "\n")
pids <- intersect(g2$id, rownames(B)); cat("probes used (n =", length(pids), "):", paste(pids, collapse = ", "), "\n"); add("glp2r_window_probes", paste(pids, collapse = ","))
## ---- tumour vs normal
tn <- do.call(rbind, lapply(pids, function(p) { a <- bt[p, ]; b <- bn[p, ]; w <- wilcox.test(a, b); data.frame(probe = p, pos = pr$chromStart[match(p, pr$id)], median_tumour = median(a, na.rm = TRUE), median_normal = median(b, na.rm = TRUE), delta = median(a, na.rm = TRUE) - median(b, na.rm = TRUE), p = w$p.value, n_t = sum(!is.na(a)), n_n = sum(!is.na(b))) }))
tn$FDR <- p.adjust(tn$p, "BH"); tn <- tn[order(tn$pos), ]; print(tn, digits = 3); write_csv_utf8(tn, file.path(RES, "glp2r_methylation_tumour_vs_normal.csv"))
mt <- colMeans(bt[pids, , drop = FALSE], na.rm = TRUE); mn <- colMeans(bn[pids, , drop = FALSE], na.rm = TRUE); wm <- wilcox.test(mt, mn)
cat(sprintf("Window mean beta: tumour median %.3f vs normal %.3f; delta %.3f; Wilcoxon p=%s\n", median(mt), median(mn), median(mt) - median(mn), fmt_p(wm$p.value)))
add("glp2r_meth_tumour_median", round(median(mt), 3)); add("glp2r_meth_normal_median", round(median(mn), 3)); add("glp2r_meth_delta", round(median(mt) - median(mn), 3)); add("glp2r_meth_tn_p", signif(wm$p.value, 4)); add("glp2r_meth_n_probes_FDR05_delta010", sum(tn$FDR < 0.05 & tn$delta > 0.10)); add("glp2r_meth_max_probe_delta", round(max(tn$delta), 3))
## ---- tumour-internal: beta vs expression
common <- intersect(colnames(bt), d$pid); dm <- d[match(common, d$pid), ]; dm$meth <- mt[common]
rho_all <- cor.test(dm$meth, dm$glp2r, method = "spearman", exact = FALSE); ok <- complete.cases(dm[, c("meth", "glp2r", "purity_absolute", "msi", "site")])
pc <- function(x, y, Z) { ok <- complete.cases(x, y, Z); Zm <- model.matrix(~ ., as.data.frame(Z[ok, , drop = FALSE]))[, -1, drop = FALSE]; rx <- resid(lm(rank(x[ok]) ~ Zm)); ry <- resid(lm(rank(y[ok]) ~ Zm)); r <- cor(rx, ry); n <- sum(ok); tt <- r * sqrt((n - 2 - ncol(Zm)) / (1 - r^2)); c(rho = r, p = 2 * pt(-abs(tt), n - 2 - ncol(Zm)), n = n) }
pa <- pc(dm$meth, dm$glp2r, dm[, c("purity_absolute", "msi", "site")]); cat(sprintf("\nWindow beta vs GLP2R expression: Spearman rho=%.3f (p=%s, n=%d); adjusted (purity, MSI, site) rho=%.3f (p=%s, n=%d)\n", rho_all$estimate, fmt_p(rho_all$p.value), nrow(dm), pa["rho"], fmt_p(pa["p"]), pa["n"]))
add("glp2r_meth_expr_rho", sprintf("%.3f (p=%s, n=%d)", rho_all$estimate, fmt_p(rho_all$p.value), nrow(dm))); add("glp2r_meth_expr_rho_adj", sprintf("%.3f (p=%s, n=%d)", pa["rho"], fmt_p(pa["p"]), pa["n"]))
pp <- do.call(rbind, lapply(pids, function(p) { x <- bt[p, common]; a <- cor.test(x, dm$glp2r, method = "spearman", exact = FALSE); q <- pc(x, dm$glp2r, dm[, c("purity_absolute", "msi", "site")]); data.frame(probe = p, rho = a$estimate, p = a$p.value, rho_adj = q["rho"], p_adj = q["p"], row.names = NULL) })); pp$FDR_adj <- p.adjust(pp$p_adj, "BH"); print(pp, digits = 3); write_csv_utf8(pp, file.path(RES, "glp2r_meth_vs_expr_probes.csv"))
add("glp2r_meth_probe_min_rho_adj", round(min(pp$rho_adj), 3)); add("glp2r_meth_probe_max_rho_adj", round(max(pp$rho_adj), 3))
if (any(dm$msi == "MSI-H", na.rm = TRUE)) { wmsi <- wilcox.test(meth ~ msi, dm); cat(sprintf("GLP2R window beta in MSI-H (median %.3f, n=%d) vs MSS/MSI-L (%.3f, n=%d): p=%s\n", median(dm$meth[dm$msi == "MSI-H"], na.rm = TRUE), sum(dm$msi == "MSI-H", na.rm = TRUE), median(dm$meth[dm$msi == "MSS/MSI-L"], na.rm = TRUE), sum(dm$msi == "MSS/MSI-L", na.rm = TRUE), fmt_p(wmsi$p.value))); add("glp2r_meth_msih_p", signif(wmsi$p.value, 4)) }
## ---- positive control: MLH1 promoter methylation (CIMP/MSI-H) vs MLH1 expression
ml <- pr[grepl("(^|,)MLH1(,|$)", pr$gene), ]; cat("\nMLH1 probes in subset:", nrow(ml), "\n")
if (nrow(ml) > 0) { dtum <- d[match(common, d$pid), ]; dtum$mlh1 <- vsa["MLH1", rownames(d)][match(common, d$pid)]; mm <- colMeans(bt[ml$id, common, drop = FALSE], na.rm = TRUE); dtum$mlh1_meth <- mm
  a <- cor.test(dtum$mlh1_meth, dtum$mlh1, method = "spearman", exact = FALSE); wmsi <- wilcox.test(mlh1_meth ~ msi, dtum); cat(sprintf("MLH1 promoter beta vs MLH1 expression: rho=%.3f p=%s; MSI-H vs MSS beta p=%s (median %.2f vs %.2f)\n", a$estimate, fmt_p(a$p.value), fmt_p(wmsi$p.value), median(dtum$mlh1_meth[which(dtum$msi == "MSI-H")], na.rm = TRUE), median(dtum$mlh1_meth[which(dtum$msi == "MSS/MSI-L")], na.rm = TRUE)))
  add("mlh1_control_rho", round(a$estimate, 3), "positive control: promoter beta vs expression"); add("mlh1_control_p", signif(a$p.value, 4)); add("mlh1_meth_msih_p", signif(wmsi$p.value, 4)) }
## ---- copy number (GISTIC2 thresholded: -2 deep del, -1 shallow del, 0, 1 gain, 2 amp)
cat("\nGISTIC: GLP2R present:", "GLP2R" %in% rownames(gi), "\n"); colnames(gi) <- substr(colnames(gi), 1, 12); cn_ok <- intersect(colnames(gi), d$pid)
if ("GLP2R" %in% rownames(gi)) { cnv <- gi["GLP2R", cn_ok]; dc <- d[match(cn_ok, d$pid), ]; dc$cn <- factor(cnv, levels = -2:2, labels = c("deep del", "shallow del", "diploid", "gain", "amp")); print(table(dc$cn)); add("cnv_n", nrow(dc)); add("cnv_table", paste(names(table(dc$cn)), table(dc$cn), sep = "=", collapse = "; "))
  kw <- kruskal.test(glp2r ~ cn, dc); sp <- cor.test(as.numeric(as.character(cnv)), dc$glp2r, method = "spearman", exact = FALSE); cat(sprintf("GLP2R expression across CN states: KW p=%s; Spearman(CN, expr) rho=%.3f p=%s\n", fmt_p(kw$p.value), sp$estimate, fmt_p(sp$p.value))); print(tapply(dc$glp2r, dc$cn, median))
  add("cnv_kw_p", signif(kw$p.value, 4)); add("cnv_spearman_rho", round(sp$estimate, 3)); add("cnv_loss_pct", round(100 * mean(cnv < 0), 1)); add("cnv_gain_pct", round(100 * mean(cnv > 0), 1)); add("cnv_deep_del_pct", round(100 * mean(cnv == -2), 2))
  wl <- wilcox.test(dc$glp2r[cnv < 0], dc$glp2r[cnv >= 0]); add("cnv_loss_vs_nonloss_p", signif(wl$p.value, 4)); cat("loss vs non-loss p=", fmt_p(wl$p.value), " medians ", round(median(dc$glp2r[cnv < 0]), 2), "vs", round(median(dc$glp2r[cnv >= 0]), 2), "\n")
  pcn <- pc(as.numeric(as.character(cnv)), dc$glp2r, dc[, c("purity_absolute", "msi", "site")]); add("cnv_rho_adj", sprintf("%.3f (p=%s)", pcn["rho"], fmt_p(pcn["p"]))) }
## ---- decision rule (pre-specified in the analysis plan)
c1 <- any(tn$FDR < 0.05 & tn$delta > 0.10) && (median(mt) - median(mn) > 0.10 && wm$p.value < 0.05); c2 <- !is.na(pa["rho"]) && pa["rho"] <= -0.30 && pa["p"] < 0.05
cat(sprintf("\nDECISION RULE: (i) tumour-vs-normal promoter delta-beta>0.10 at FDR<0.05 = %s; (ii) adjusted beta-expression rho <= -0.30 = %s => 'epigenetic silencing' wording %s\n", c1, c2, ifelse(c1 && c2, "RETAINED", "DOWNGRADED")))
add("decision_i_delta_beta", c1); add("decision_ii_rho", c2); add("decision_epigenetic_silencing_retained", c1 && c2)
saveRDS(list(tn = tn, pp = pp, dm = dm, mt = mt, mn = mn, c1 = c1, c2 = c2), file.path(RES, "epi_objects.rds")); rm_write(do.call(rbind, RM), "epi"); log_msg("done"); sink()
