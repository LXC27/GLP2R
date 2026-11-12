## 98_figS.R : Supplementary Figures S1-S6
source(file.path(Sys.getenv("GLP2R_SCRIPTS", "W:/Re_build/05_Analysis/scripts"), "90_fig_helpers.R"))
suppressPackageStartupMessages(library(survival))
rd <- function(f) readRDS(file.path(RES, f)); rdc <- function(f) readRDS(file.path(DC, f))
leg_in <- function(x = 0.98, y = 0.98, hj = 1, vj = 1) theme(legend.position = "inside", legend.position.inside = c(x, y), legend.justification = c(hj, vj), legend.key.height = unit(2.5, "mm"))

## ================= S6: LASSO risk scores (waterfall) =================
lo <- rd("lasso_objects.rds"); ttl <- c(TCGA_train = "TCGA training", TCGA_test = "TCGA held-out test", GSE39582 = "GSE39582", GSE17536 = "GSE17536")
wf <- lapply(names(ttl), function(n) { x <- lo$cohorts[[n]]; x <- x[order(x$s), ]; x$r <- seq_len(nrow(x)); x$st <- factor(ifelse(x$event == 1, "Death", "Censored"), c("Censored", "Death")); cut <- if (n == "TCGA_test") median(lo$cohorts$TCGA_train$s) else median(x$s)
  ggplot(x, aes(r, s, colour = st)) + geom_point(size = 0.6, stroke = 0) + geom_hline(yintercept = cut, linetype = 2, linewidth = 0.3, colour = "grey40") + scale_colour_manual(values = c(Censored = "grey65", Death = "#C8472B"), name = NULL) +
    labs(x = "Patients ranked by signature score", y = "Signature score", title = sprintf("%s (n = %d, %d events)", ttl[[n]], nrow(x), sum(x$event))) + leg_in(0.02, 0.98, 0, 1) })
fS1 <- wrap_plots(wf, ncol = 2) + plot_annotation(tag_levels = "A")
save_fig(fS1, "FigureS6_LASSO_risk_scores", w = 180, h = 110, dir = FIGS)

## ================= S5: GEO external cohorts (probe signal & KM) =================
gc <- rdc("geo_cohorts.rds"); go <- rd("geo_objects.rds")
dens <- function(nm, ttl) { e <- gc[[nm]]$expr; pm <- rowMeans(e); gp <- mean(gc[[nm]]$glp2r_probe_expr); pr <- 100 * mean(pm <= gp)
  ggplot(data.frame(m = pm), aes(m)) + geom_density(fill = "grey85", colour = "grey40", linewidth = 0.3) + geom_vline(xintercept = gp, colour = "#C8472B", linewidth = 0.5) +
    annotate("text", x = max(pm) - 0.1, y = Inf, label = sprintf("Red line: GLP2R probe 221312_at\nmean log\u2082 = %.2f\n%.0fth percentile of %s probesets", gp, pr, format(nrow(e), big.mark = ",")), hjust = 1, vjust = 1.3, size = TXT, lineheight = 0.95) +
    labs(x = "Mean log\u2082 expression per probeset", y = "Density", title = ttl) }
pA <- dens("GSE39582", "GSE39582: probe signal vs array background"); pB <- dens("GSE17536", "GSE17536: probe signal vs array background")
kmg <- function(x, ttl) { x <- x[!is.na(x$time) & !is.na(x$event) & !is.na(x$glp2r), ]; x$grp <- factor(ifelse(x$glp2r > median(x$glp2r), "High", "Low"), c("Low", "High")); lr <- survdiff(Surv(time, event) ~ grp, x); p <- 1 - pchisq(lr$chisq, 1); h <- cox_tab(coxph(Surv(time, event) ~ grp, x))[1, ]
  wrap_elements(full = km_panel(x, "grp", col_grp, xmax = 120, step = 24, title = sprintf("%s (n = %d, %d events)", ttl, nrow(x), sum(x$event)), annot = sprintf("HR (high vs low) = %s\nLog-rank %s", fmt_ci(h$HR, h$lo, h$hi), plab(p)), legend_pos = c(0.03, 0.03), leg_br = TRUE)) }
t39 <- go$t39; t17 <- go$t17
pC <- kmg(t39, "GSE39582"); pD <- kmg(t17, "GSE17536")
p39 <- gc$GSE39582$pheno; x39 <- as.numeric(gc$GSE39582$glp2r_probe_expr[1, ]); names(x39) <- colnames(gc$GSE39582$glp2r_probe_expr); p39$glp2r <- x39[rownames(p39)]; if (all(is.na(p39$glp2r))) p39$glp2r <- x39[p39$geo_accession]
bn <- data.frame(g = factor(ifelse(p39$dataset == "Non Tumoral", "Non-tumoral mucosa", "Tumour"), c("Non-tumoral mucosa", "Tumour")), y = p39$glp2r); wn <- wilcox.test(y ~ g, bn)
pE <- box_panel(bn, "g", "y", c("Non-tumoral mucosa" = col_tn[["Normal"]], Tumour = col_tn[["Tumor"]]), "Probe 221312_at (log\u2082)", "GSE39582: tumour vs mucosa", sprintf("Wilcoxon P = %.4f\nAUC = %.2f", wn$p.value, auc_mw(bn$y[bn$g != "Tumour"], bn$y[bn$g == "Tumour"])))
pN <- ggplot() + xlim(0, 1) + ylim(0, 1) + theme_void() + labs(tag = " ") + annotate("text", x = 0.02, y = 0.62, hjust = 0, vjust = 1, size = TXT * 1.05, lineheight = 1.1, label = paste(strwrap("Interpretation: probe 221312_at lies in the lowest decile of probeset signal in both cohorts (red lines in A and B), i.e. close to array background. Differences between groups of a few tenths of a log\u2082 unit are therefore not reliably distinguishable from noise, and these microarray data cannot be taken as independent validation of the RNA-seq prognostic association.", 60), collapse = "\n"))
fS2 <- (pA | pB) / (pC | pD) / (pE | plot_spacer()) + plot_layout(heights = c(0.8, 1.15, 0.8)) + plot_annotation(tag_levels = "A")
save_fig(fS2, "FigureS5_GEO_cohorts", w = 180, h = 205, dir = FIGS)

## ================= S3: WGCNA details & cell-of-origin signal =================
w <- rd("wgcna_objects.rds"); mt <- w$modtab; vt <- rdc("vst_tumor_591.rds"); va <- rdc("vst_all_624.rds")
fi <- w$fi; pA <- ggplot(fi, aes(Power, SFT.R.sq)) + geom_hline(yintercept = 0.85, linetype = 2, linewidth = 0.3, colour = "grey50") + geom_line(linewidth = 0.4) + geom_point(size = 1) + geom_point(data = fi[fi$Power == w$power, ], colour = "#C8472B", size = 2) +
  annotate("text", x = w$power + 1, y = 0.1, label = sprintf("Selected power = %d\nsigned R\u00b2 = %.2f", w$power, fi$SFT.R.sq[fi$Power == w$power]), hjust = 0, size = TXT, lineheight = 0.95) + labs(x = "Soft-thresholding power", y = "Scale-free fit (signed R\u00b2)", title = "Soft-threshold selection")
ms <- mt[order(-mt$size), ]; ms$module <- factor(ms$module, rev(ms$module))
pB <- ggplot(ms, aes(size, module, fill = module)) + geom_col(colour = "grey30", linewidth = 0.2, width = 0.75) + scale_fill_manual(values = setNames(as.character(ms$module), as.character(ms$module)), guide = "none") + geom_text(aes(label = size), hjust = -0.15, size = TXT) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) + labs(x = "Genes per module (grey: 1,999 unassigned)", y = NULL, title = "Module sizes")
ht <- rbind(data.frame(module = mt$module, trait = "GLP2R", rho = mt$rho_GLP2R), data.frame(module = mt$module, trait = "GLP2R\n(partial)", rho = mt$rho_partial), data.frame(module = mt$module, trait = "MSI-H", rho = mt$rho_MSI_H), data.frame(module = mt$module, trait = "Purity", rho = mt$rho_purity))
ht$module <- factor(ht$module, rev(mt$module[order(-abs(mt$rho_GLP2R))])); ht$trait <- factor(ht$trait, c("GLP2R", "GLP2R\n(partial)", "MSI-H", "Purity"))
pC <- ggplot(ht, aes(trait, module, fill = rho)) + geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.2f", rho)), size = TXT * 0.9) + scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#C8472B", limits = c(-0.6, 0.6), name = "\u03c1") +
  labs(x = NULL, y = NULL, title = "Module\u2013trait correlations (tumours)") + theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.key.height = unit(6, "mm"), legend.key.width = unit(2.5, "mm"))
tq <- names(w$mc)[w$mc == "turquoise"]; me <- w$ME[, "turquoise"]; kme <- as.numeric(cor(t(vt[tq, rownames(w$ME)]), me)); names(kme) <- tq; kg <- kme[["GLP2R"]]
pD <- ggplot(data.frame(k = kme), aes(k)) + geom_histogram(bins = 40, fill = "grey80", colour = "grey40", linewidth = 0.2) + geom_vline(xintercept = kg, colour = "#C8472B", linewidth = 0.5) +
  annotate("text", x = kg - 0.02, y = Inf, label = sprintf("GLP2R\nkME %.2f\nrank %s of %s", kg, format(sum(kme >= kg), big.mark = ","), format(length(kme), big.mark = ",")), hjust = 1, vjust = 1.3, size = TXT, lineheight = 0.95) + scale_x_continuous(limits = c(0, 1)) + labs(x = "Intramodular connectivity (kME)", y = "Genes", title = "GLP2R in the turquoise module")
## E: hub genes vs GLP2R in tumours and normals
hubs <- names(sort(kme, decreasing = TRUE))[1:100]; hubs <- setdiff(hubs, "GLP2R"); ntm <- grep("-11[AB]-", colnames(va), value = TRUE); tum <- colnames(vt)
rho_t <- vapply(hubs, function(g) cor(vt["GLP2R", tum], vt[g, tum], method = "spearman"), 0); rho_n <- vapply(hubs, function(g) cor(va["GLP2R", ntm], va[g, ntm], method = "spearman"), 0)
he <- data.frame(g = factor(rep(c("Tumours (n = 591)", sprintf("Normal mucosa (n = %d)", length(ntm))), each = length(hubs)), c(sprintf("Normal mucosa (n = %d)", length(ntm)), "Tumours (n = 591)")), r = c(rho_t, rho_n))
cat(sprintf("S3E hub-gene median rho tumours %.3f normals %.3f (n normal = %d)\n", median(rho_t), median(rho_n), length(ntm)))
rm_write(data.frame(key = c("top100_hub_median_rho_tumour", "top100_hub_median_rho_normal", "top100_hub_n_normal"), value = c(round(median(rho_t), 3), round(median(rho_n), 3), length(ntm)),
                    description = c("Median Spearman rho between GLP2R and the top-100 turquoise kME hub genes (excl. GLP2R), 591 tumours", "same, normal mucosa", "normals used")), "figS")
pE <- ggplot(he, aes(g, r, fill = g)) + geom_hline(yintercept = 0, linewidth = 0.3) + geom_boxplot(outlier.size = 0.4, linewidth = 0.3, width = 0.55, alpha = 0.7) + scale_fill_manual(values = setNames(c(col_tn[["Normal"]], col_tn[["Tumor"]]), levels(he$g)), guide = "none") +
  annotate("text", x = 1.5, y = 0.8, label = sprintf("Median \u03c1: normal %.2f; tumour %.2f", median(rho_n), median(rho_t)), size = TXT, hjust = 0.5) + labs(x = NULL, y = "Spearman \u03c1 with GLP2R", title = "Top-100 turquoise hub genes vs GLP2R")
pN3 <- ggplot() + xlim(0, 1) + ylim(0, 1) + theme_void() + labs(tag = " ") + annotate("text", x = 0.02, y = 0.62, hjust = 0, vjust = 1, size = TXT * 1.05, lineheight = 1.1, label = paste(strwrap("Interpretation: the turquoise module (1,747 genes) is enriched for epithelial\u2013mesenchymal transition, myogenesis, blood-vessel morphogenesis and stromal/endothelial signatures, and is strongly anti-correlated with tumour purity (\u03c1 = \u22120.56). GLP2R is a peripheral member (kME rank 1,516 of 1,747). Its module-membership pattern therefore reflects the non-epithelial (stromal/vascular/neuro-muscular) compartment rather than a tumour-cell programme.", 60), collapse = "\n"))
fS3 <- (pA | pB) / (pC | pD) / (pE | plot_spacer()) + plot_layout(heights = c(1, 1.15, 0.9)) + plot_annotation(tag_levels = "A")
save_fig(fS3, "FigureS3_WGCNA_and_cell_origin", w = 180, h = 205, dir = FIGS)

## ================= S4: checkpoint and compartment-score heatmaps =================
io <- rd("immune_objects.rds"); ck <- io$ckr
hm <- function(df, ttl, xl, w_) { ggplot(df, aes(m, g, fill = v)) + geom_tile(colour = "white") + geom_text(aes(label = lab), size = TXT * 0.9) + scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#C8472B", limits = c(-w_, w_), name = "\u03c1 / r", oob = scales::squish) + labs(x = NULL, y = NULL, title = ttl) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.key.height = unit(6, "mm"), legend.key.width = unit(2.5, "mm")) }
star <- function(q) ifelse(q < 0.001, "***", ifelse(q < 0.01, "**", ifelse(q < 0.05, "*", "")))
ck$lab_g <- ifelse(ck$Name == ck$Gene, ck$Gene, paste0(ck$Name, " (", ck$Gene, ")")); gord <- rev(ck$lab_g)
dfc <- rbind(data.frame(g = ck$lab_g, m = "Pearson r", v = ck$Pearson_r, q = ck$Pearson_FDR), data.frame(g = ck$lab_g, m = "Spearman \u03c1", v = ck$Spearman_rho, q = ck$Spearman_FDR), data.frame(g = ck$lab_g, m = "Partial\u1d43 \u03c1", v = ck$partial_rho, q = ck$partial_FDR))
dfc$g <- factor(dfc$g, gord); dfc$m <- factor(dfc$m, c("Pearson r", "Spearman \u03c1", "Partial\u1d43 \u03c1")); dfc$lab <- paste0(sprintf("%.2f", dfc$v), star(dfc$q))
pA <- hm(dfc, "Immune-checkpoint genes", "", 0.3)
sr <- io$sc_res; lm_ <- c(stromal = "ESTIMATE stromal", immune = "ESTIMATE immune", estimate = "ESTIMATE score"); sr$score <- ifelse(sr$score %in% names(lm_), lm_[sr$score], sr$score); sr$g <- factor(sr$score, rev(sr$score)); dfs <- rbind(data.frame(g = sr$g, m = "Spearman \u03c1", v = sr$rho, q = sr$FDR_p), data.frame(g = sr$g, m = "Partial\u1d43 \u03c1", v = sr$partial_rho, q = sr$FDR_partial_p)); dfs$m <- factor(dfs$m, c("Spearman \u03c1", "Partial\u1d43 \u03c1")); dfs$lab <- paste0(sprintf("%.2f", dfs$v), star(dfs$q))
pB <- hm(dfs, "MCP-counter and ESTIMATE scores", "", 0.5)
fS4 <- (pA | pB) + plot_annotation(tag_levels = "A")
save_fig(fS4, "FigureS4_checkpoints_and_compartment_scores", w = 180, h = 85, dir = FIGS)

## ================= S1: cohort flow (REMARK), PCA, detection =================
co <- rdc("cohort.rds"); fl <- co$flow
bx <- function(x, y, w, h, lab, fill = "white") list(annotate("rect", xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2, fill = fill, colour = "grey30", linewidth = 0.3), annotate("text", x = x, y = y, label = lab, size = TXT, lineheight = 0.95))
ar <- function(x, y, xe, ye) annotate("segment", x = x, y = y, xend = xe, yend = ye, linewidth = 0.3, arrow = arrow(length = unit(1.4, "mm"), type = "closed"))
ft <- function(n) format(n, big.mark = ",")
pF <- ggplot() + xlim(0, 100) + ylim(0, 100) + theme_void() +
  bx(26, 88, 44, 10, sprintf("TCGA-COAD/READ RNA-seq samples\nn = %d", fl$total_rnaseq), "grey95") + ar(26, 83, 26, 71) +
  bx(26, 66, 44, 10, sprintf("Primary colorectal tumours\nn = %d (431 COAD, 160 READ)", fl$primary), "#DCE9F7") + ar(26, 61, 26, 51) +
  bx(26, 46, 44, 10, sprintf("With overall survival (TCGA-CDR)\nn = %d; %d events (excluded: 4 no OS record, 21 zero follow-up)", fl$with_os, fl$os_events), "#DCE9F7") + ar(26, 41, 26, 33) +
  bx(26, 28, 44, 10, sprintf("Complete cases: age, sex, stage, site\nn = %d (stage missing 107); with MSI n = %d (77 events)", fl$cc_age_sex_stage_site, fl$cc_plus_msi), "#DCE9F7") +
  ar(48, 44, 55, 39) + bx(77, 37, 44, 12, "Signature development: event-stratified ~70/30 split\nof the 566 patients with OS\ntraining n = 396 (80 events); test n = 170 (34)", "#DCE9F7") + ar(77, 31, 77, 22) +
  ar(48, 88, 55, 88) + bx(77, 88, 44, 20, sprintf("Not in the tumour analyses\nNormal mucosa n = %d (tumour\u2013normal comparison)\nMetastatic sample n = %d (excluded)\nPaired design: normal samples of 32 of these patients reused;\ntheir 32 primary tumours from GDC (not among the 591)", fl$normal, fl$metastatic_excl), "white") +
  ar(48, 66, 55, 66) + bx(77, 62, 44, 26, "Tumour-only analyses (primary tumours)\nABSOLUTE purity n = 574; MSI status n = 575\nStage n = 480; WES/MC3 mutations n = 504\nHM450 methylation n = 395 (370 with RNA-seq)\nGISTIC2 copy number n = 581\nImmune/pathway/WGCNA: all 591", "white") +
  bx(77, 14, 44, 16, "External cohorts (GEO, HG-U133 Plus 2)\nGSE39582: 566 tumours (556 with OS)\nGSE17536: 177 tumours\nGLP2R probe 221312_at", "white")
## PCA
tt <- substr(colnames(va), 14, 15); cond <- ifelse(tt == "11", "Normal", "Tumour"); site <- co$clin$project_id[match(colnames(va), co$clin$barcode)]; v <- apply(va, 1, var); top <- names(sort(v, decreasing = TRUE))[1:2000]; pc <- prcomp(t(va[top, ]), scale. = FALSE); ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
pcd <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], g = ifelse(cond == "Normal", "Normal", site)); pcd$g <- factor(pcd$g, c("Normal", "TCGA-COAD", "TCGA-READ"))
pP <- ggplot(pcd, aes(PC1, PC2, colour = g)) + geom_point(size = 0.6, alpha = 0.7, stroke = 0) + scale_colour_manual(values = c(Normal = "#4C78A8", `TCGA-COAD` = "#D95F02", `TCGA-READ` = "#7B3294"), name = NULL) + labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]), title = "PCA, top-2000 variable genes (VST)") + leg_in(0.98, 0.98)
## detection
dg <- rdc("de_tumor_normal.rds"); eg <- dg$ensembl_full[which(dg$symbol == "GLP2R")]; cn <- co$counts[eg, ]; cn <- cn[match(colnames(va), colnames(co$counts))]; names(cn) <- colnames(va)
dtc <- data.frame(g = factor(ifelse(cond == "Normal", sprintf("Normal (n = %d)", sum(cond == "Normal")), sprintf("Tumour (n = %d)", sum(cond == "Tumour"))), c(sprintf("Normal (n = %d)", sum(cond == "Normal")), sprintf("Tumour (n = %d)", sum(cond == "Tumour")))), y = log10(as.numeric(cn) + 1)); dtc$det <- cn >= 10
pr <- tapply(dtc$det, dtc$g, mean) * 100
pQ <- ggplot(dtc, aes(g, y, fill = g)) + geom_hline(yintercept = log10(11), linetype = 2, linewidth = 0.3, colour = "grey40") + geom_boxplot(outlier.shape = NA, linewidth = 0.3, width = 0.55, alpha = 0.7) + geom_jitter(width = 0.15, size = 0.3, alpha = 0.4, stroke = 0) +
  scale_fill_manual(values = setNames(c(col_tn[["Normal"]], col_tn[["Tumor"]]), levels(dtc$g)), guide = "none") + annotate("text", x = 1:2, y = max(dtc$y) * 1.05, label = sprintf("\u226510 reads: %.1f%%", pr), size = TXT) + annotate("text", x = 0.45, y = log10(11) + 0.08, label = "10 reads", hjust = 0, size = TXT, colour = "grey30") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) + labs(x = NULL, y = expression(log[10](raw~reads + 1)), title = "GLP2R read counts")
fS5 <- pF / (pP | pQ) + plot_layout(heights = c(1.4, 1)) + plot_annotation(tag_levels = "A")
save_fig(fS5, "FigureS1_cohort_flow_and_QC", w = 180, h = 175, dir = FIGS)

## ================= S2: survival sensitivity =================
so <- rd("survival_objects.rds"); dd <- so$dd; rs <- read.csv(file.path(RES, "rm_survival.csv")); R <- setNames(rs$value, rs$key)
dd$grp <- factor(ifelse(dd$glp2r > median(dd$glp2r), "High", "Low"), c("Low", "High"))
kms <- function(x, ttl, tm = "time", ev = "event", xmax = 120, step = 24, ylab = "Overall survival") { x$time <- x[[tm]]; x$event <- x[[ev]]; x <- x[!is.na(x$time) & !is.na(x$event), ]; lr <- survdiff(Surv(time, event) ~ grp, x); p <- 1 - pchisq(lr$chisq, 1); h <- cox_tab(coxph(Surv(time, event) ~ grp, x))[1, ]
  wrap_elements(full = km_panel(x, "grp", col_grp, xmax = xmax, step = step, title = sprintf("%s (n = %d, %d events)", ttl, nrow(x), sum(x$event)), annot = sprintf("HR (high vs low) = %s\nLog-rank %s", fmt_ci(h$HR, h$lo, h$hi), plab(p)), ylab = ylab, legend_pos = c(0.03, 0.03), leg_br = TRUE)) }
pA <- kms(dd[dd$site == "COAD", ], "OS, colon (COAD)"); pB <- kms(dd[dd$site == "READ", ], "OS, rectum (READ)")
dp <- dd[!is.na(dd$cdr_PFI_time) & dd$cdr_PFI_time > 0 & !is.na(dd$cdr_PFI), ]; dp$pfi_m <- dp$cdr_PFI_time / 30.4375; dp$grp <- factor(ifelse(dp$glp2r > median(dp$glp2r), "High", "Low"), c("Low", "High"))
pC <- kms(dp, "Progression-free interval", "pfi_m", "cdr_PFI", ylab = "Progression-free interval")
pN6 <- ggplot() + xlim(0, 1) + ylim(0, 1) + theme_void() + labs(tag = " ") + annotate("text", x = 0.02, y = 0.92, hjust = 0, vjust = 1, size = TXT * 1.05, lineheight = 1.1, label = paste(strwrap("Notes: groups are defined by the median of GLP2R in the full OS cohort (A, B) or of the PFI cohort (C). Subgroup curves rest on few events (rectum: 24), so apparent separation is imprecise and is not interpreted as site-specific effect modification (interaction P = 0.88). GLP2R has no association with progression-free interval (adjusted HR 0.95, P = 0.70).", 58), collapse = "\n"))
pr_ <- function(s) { m <- regmatches(s, regexec("^([0-9.]+) \\(([0-9.]+)-([0-9.]+)\\); p=([0-9.]+)", s))[[1]]; as.numeric(m[2:5]) }
sens <- list(c("Primary model M2 (OS)", R[["cox_glp2r_hr_M2_MSI_primary_model_"]]), c("Follow-up truncated at 60 months", R[["cox_m2_trunc60_glp2r"]]), c("Stratified by stage", R[["cox_m2_stage_strata_glp2r"]]), c("log\u2082CPM instead of VST", R[["cox_m2_log2cpm_glp2r"]]), c("MSS/MSI-L only (n = 393)", R[["cox_glp2r_in_MSS"]]), c("Multiple imputation (n = 566)", R[["mi_glp2r_hr"]]), c("PFI, adjusted (n = 451)", R[["pfi_adj_glp2r"]]), c("High vs low GLP2R (binary)", R[["cox_m2_median_binary_glp2r"]]))
fs <- do.call(rbind, lapply(sens, function(z) { v <- pr_(z[2]); data.frame(label = z[1], HR = v[1], lo = v[2], hi = v[3], p = v[4]) })); fs$text <- paste0(fmt_ci(fs$HR, fs$lo, fs$hi), "; P = ", sprintf("%.2f", fs$p)); fs$col <- ifelse(seq_len(nrow(fs)) == 1, "b", "a")
pD <- forest_panel(fs[, c("label", "HR", "lo", "hi", "text", "col")], xlab = "Hazard ratio for GLP2R (per SD unless stated; 95% CI)", xlim = c(0.3, 1.5), breaks = c(0.5, 0.75, 1, 1.25), title = "Robustness of the adjusted association", wid = c(1, 0.5), txt_head = "HR (95% CI); P")
fS6 <- (pA | pB) / (pC | plot_spacer()) / wrap_elements(full = pD) + plot_layout(heights = c(1, 1, 0.9)) + plot_annotation(tag_levels = "A")
save_fig(fS6, "FigureS2_survival_sensitivity", w = 180, h = 215, dir = FIGS); cat("Supplementary figures saved\n")
