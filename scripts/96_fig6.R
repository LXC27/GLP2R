## 96_fig6.R : Figure 6 - epigenetic, copy-number and somatic context of GLP2R loss
source(file.path(Sys.getenv("GLP2R_SCRIPTS", "W:/Re_build/05_Analysis/scripts"), "90_fig_helpers.R"))
e <- readRDS(file.path(RES, "epi_objects.rds")); x <- readRDS(file.path(RES, "epi_extra_objects.rds")); rme <- setNames(read.csv(file.path(RES, "rm_epi.csv"))$value, read.csv(file.path(RES, "rm_epi.csv"))$key)
cT <- col_tn[["Tumor"]]; cN <- col_tn[["Normal"]]; cl_sel <- "#C8472B"
## A: probe-level delta-beta
tn <- e$tn[order(e$tn$pos), ]; tn$probe <- factor(tn$probe, tn$probe); tn$sig <- ifelse(tn$FDR < 0.05, "FDR < 0.05", "n.s.")
pA <- ggplot(tn, aes(delta, factor(probe, rev(levels(probe))), colour = sig)) + geom_vline(xintercept = 0, linewidth = 0.3) + geom_vline(xintercept = c(-0.1, 0.1), linetype = 2, linewidth = 0.25, colour = "grey50") + geom_segment(aes(xend = 0, yend = factor(probe, rev(levels(probe)))), linewidth = 0.4) + geom_point(size = 1.5) +
  scale_colour_manual(values = col_sig, name = NULL) + scale_x_continuous(limits = c(-0.2, 0.12), breaks = seq(-0.2, 0.1, 0.1)) + labs(y = NULL, x = "\u0394\u03b2 (tumour \u2212 normal; dashed: \u00b10.10)", title = "Probe-level methylation (5\u2032\u21923\u2032)") +
  theme(axis.text.y = element_text(size = 5.2), legend.position = "inside", legend.position.inside = c(0.02, 0.02), legend.justification = c(0, 0), legend.key.height = unit(2.5, "mm"))
## B: window mean beta tumour vs normal
bw <- data.frame(beta = c(e$mt, e$mn), g = factor(rep(c("Tumour", "Normal"), c(length(e$mt), length(e$mn))), c("Normal", "Tumour"))); wm <- wilcox.test(beta ~ g, bw)
pB <- box_panel(bw, "g", "beta", c(Normal = cN, Tumour = cT), "Mean \u03b2 of 8 promoter probes", "Promoter window", sprintf("\u0394 = %+.3f\nWilcoxon %s", median(e$mt) - median(e$mn), plab(wm$p.value)))
## C: beta vs expression
dm <- e$dm; ct <- cor.test(dm$meth, dm$glp2r, method = "spearman", exact = FALSE)
pC <- ggplot(dm, aes(meth, glp2r, colour = msi)) + geom_point(size = 0.5, alpha = 0.6, stroke = 0) + scale_colour_manual(values = col_msi, breaks = names(col_msi), na.value = "grey85", name = NULL) + geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "black", linewidth = 0.4, formula = y ~ x) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.3))) + labs(x = "Promoter-window mean \u03b2", y = "GLP2R (VST)", title = "Methylation vs expression (tumours)") + theme(legend.position = "inside", legend.position.inside = c(0.98, 0.98), legend.justification = c(1, 1), legend.key.height = unit(2.5, "mm")) +
  annotate("text", x = min(dm$meth, na.rm = TRUE), y = Inf, label = sprintf("Spearman \u03c1 = %.2f\n(%s)\nadjusted\u1d43 \u03c1 = %s", ct$estimate, plab(ct$p.value), sub(" \\(p=.*", "", rme[["glp2r_meth_expr_rho_adj"]])), hjust = 0, vjust = 1.3, size = TXT, lineheight = 0.95)
## D: MLH1 positive control
ml <- x$mlh1[complete.cases(x$mlh1[, c("mlh1", "mlh1_meth")]), ]; cm <- cor.test(ml$mlh1_meth, ml$mlh1, method = "spearman", exact = FALSE)
pD <- ggplot(ml, aes(mlh1_meth, mlh1, colour = msi)) + geom_point(size = 0.5, alpha = 0.6, stroke = 0) + scale_colour_manual(values = col_msi, breaks = names(col_msi), na.value = "grey85", name = NULL) + geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "black", linewidth = 0.4, formula = y ~ x) +
  labs(x = "MLH1 promoter mean \u03b2", y = "MLH1 (VST)", title = "Positive control: MLH1") + theme(legend.position = "inside", legend.position.inside = c(0.02, 0.02), legend.justification = c(0, 0), legend.key.height = unit(2.5, "mm")) +
  annotate("text", x = max(ml$mlh1_meth), y = max(ml$mlh1), label = sprintf("Spearman \u03c1 = %.2f\n%s", cm$estimate, plab(cm$p.value)), hjust = 1, vjust = 1, size = TXT, lineheight = 0.95)
## E: CNV
dc <- x$cnv; dc$cn <- droplevels(dc$cn); kw <- kruskal.test(glp2r ~ cn, dc); sp <- cor.test(dc$cnv, dc$glp2r, method = "spearman", exact = FALSE)
pE <- box_panel(dc, "cn", "glp2r", c("Deep loss" = "#2C7BB6", "Shallow loss" = "#ABD9E9", Diploid = "#BDBDBD", Gain = "#FDAE61", Amp = "#D7191C"), "GLP2R (VST)", "GISTIC2 copy-number state",
  sprintf("Spearman \u03c1 = %.2f (%s)\nadjusted\u1d43 \u03c1 = %s", sp$estimate, plab(sp$p.value), sub(" \\(p=.*", "", rme[["cnv_rho_adj"]]))) + theme(axis.text.x = element_text(size = 5.2))
## F: TMB
tm <- x$tm; ta <- cor.test(tm$glp2r, tm$ltmb, method = "spearman", exact = FALSE); ms <- tm[which(tm$msi == "MSS/MSI-L"), ]; tb <- cor.test(ms$glp2r, ms$ltmb, method = "spearman", exact = FALSE)
pF <- ggplot(tm, aes(glp2r, ltmb, colour = msi)) + geom_point(size = 0.5, alpha = 0.6, stroke = 0) + scale_colour_manual(values = col_msi, breaks = names(col_msi), na.value = "grey85", name = NULL) + geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "black", linewidth = 0.4, formula = y ~ x) +
  labs(x = "GLP2R (VST)", y = expression(log[10](TMB + 0.1)), title = "Tumour mutational burden") + theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98), legend.justification = c(0, 1), legend.key.height = unit(2.5, "mm")) +
  annotate("text", x = max(tm$glp2r), y = min(tm$ltmb), label = sprintf("All: \u03c1 = %.2f (%s)\nMSS only: \u03c1 = %.2f (%s)", ta$estimate, plab(ta$p.value), tb$estimate, plab(tb$p.value)), hjust = 1, vjust = 0, size = TXT, lineheight = 0.95)
## G: DNMT/TET tumour vs normal
rg <- x$reg; rg$gene <- factor(rg$gene, rev(rg$gene))
g1 <- rbind(data.frame(gene = rg$gene, lfc = rg$log2FC_unpaired, q = rg$padj_unpaired, des = "Unpaired (n = 33 vs 591)"), data.frame(gene = rg$gene, lfc = rg$log2FC_paired, q = rg$padj_paired, des = "Paired (32 pairs)"))
g1$sig <- g1$q < 0.05
pG <- ggplot(g1, aes(lfc, gene, shape = des, fill = sig, colour = des)) + geom_vline(xintercept = 0, linewidth = 0.3) + geom_point(size = 1.8, position = position_dodge(width = 0.55), stroke = 0.5) +
  scale_shape_manual(values = c(`Unpaired (n = 33 vs 591)` = 21, `Paired (32 pairs)` = 23), name = NULL) + scale_colour_manual(values = c(`Unpaired (n = 33 vs 591)` = "grey20", `Paired (32 pairs)` = cl_sel), name = NULL) + scale_fill_manual(values = c(`TRUE` = "black", `FALSE` = "white"), guide = "none") +
  labs(x = "log\u2082 fold change, tumour vs normal (filled: FDR < 0.05)", y = NULL, title = "Methylation regulators: expression") + theme(legend.position = "inside", legend.position.inside = c(0.98, 0.02), legend.justification = c(1, 0), legend.key.height = unit(2.5, "mm")) + guides(colour = guide_legend(override.aes = list(fill = "grey60")))
## H: DNMT/TET correlation with GLP2R
h1 <- rbind(data.frame(gene = rg$gene, rho = rg$rho, q = rg$FDR_rho, ty = "Spearman"), data.frame(gene = rg$gene, rho = rg$rho_adj, q = rg$FDR_rho_adj, ty = "Partial\u1d43")); h1$sig <- h1$q < 0.05
pH <- ggplot(h1, aes(rho, gene, shape = ty, fill = sig, colour = ty)) + geom_vline(xintercept = 0, linewidth = 0.3) + geom_point(size = 1.8, position = position_dodge(width = 0.55), stroke = 0.5) +
  scale_shape_manual(values = setNames(c(21, 23), c("Spearman", "Partial\u1d43")), name = NULL) + scale_colour_manual(values = setNames(c("grey20", cl_sel), c("Spearman", "Partial\u1d43")), name = NULL) + scale_fill_manual(values = c(`TRUE` = "black", `FALSE` = "white"), guide = "none") +
  labs(x = "\u03c1 with GLP2R in tumours (filled: FDR < 0.05)", y = NULL, title = "Regulators vs GLP2R (tumours, n = 561)") + theme(legend.position = "inside", legend.position.inside = c(0.02, 0.02), legend.justification = c(0, 0), legend.key.height = unit(2.5, "mm")) + guides(colour = guide_legend(override.aes = list(fill = "grey60")))
## I: mutation frequency
fq <- x$freq; fq$gene <- factor(fq$gene, fq$gene[order(fq$pct)]); fq$hl <- fq$gene == "GLP2R"
pI <- ggplot(fq, aes(pct, gene, fill = hl)) + geom_col(width = 0.7) + geom_text(aes(label = sprintf("%.1f%%", pct)), hjust = -0.1, size = TXT) + scale_fill_manual(values = c(`TRUE` = cl_sel, `FALSE` = "grey65"), guide = "none") + scale_x_continuous(limits = c(0, 90), breaks = seq(0, 80, 20), expand = c(0, 0)) +
  labs(x = sprintf("Tumours with non-silent mutation (%%; n = %d)", x$n_seq), y = NULL, title = "Somatic mutation frequency")
## J: mutation status vs GLP2R (partial rho with approximate CI)
mo <- x$mo; mo$se <- 1 / sqrt(mo$n_mut + mo$n_wt - 6); mo$lo <- tanh(atanh(mo$rho_adj) - 1.96 * mo$se); mo$hi <- tanh(atanh(mo$rho_adj) + 1.96 * mo$se); mo$gene <- factor(mo$gene, rev(mo$gene)); mo$sig <- ifelse(mo$FDR_adj < 0.05, "FDR < 0.05", "n.s.")
mo$lab <- sprintf("n = %d; MSS-only %s", mo$n_mut, plab(mo$p_mss, "P"))
pJ <- ggplot(mo, aes(rho_adj, gene, colour = sig)) + geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey40") + geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.35, orientation = "y") + geom_point(size = 1.4) + geom_text(aes(x = 0.30, label = lab), hjust = 0, size = TXT * 0.9, colour = "black") +
  scale_colour_manual(values = col_sig, name = NULL) + scale_x_continuous(limits = c(-0.35, 0.75), breaks = c(-0.3, -0.15, 0, 0.15)) + labs(x = "Partial \u03c1 of GLP2R with mutation status (95% CI)", y = NULL, title = "GLP2R vs driver-gene mutation (tumours)") +
  theme(legend.position = "inside", legend.position.inside = c(0.02, 0.02), legend.justification = c(0, 0), legend.key.height = unit(2.5, "mm"))
fig <- (pA | pB | pC) / (pD | pE | pF) / (pG | pH) / (pI | pJ) + plot_layout(heights = c(1, 1, 0.9, 1)) + plot_annotation(tag_levels = "A")
save_fig(fig, "Figure6_GLP2R_epigenetic_genomic", w = 180, h = 215); cat("Figure 6 saved\n")
