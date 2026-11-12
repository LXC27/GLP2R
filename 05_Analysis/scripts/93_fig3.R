## 93_fig3.R : Figure 3 - immune microenvironment associations of GLP2R (tumour-only; adjusted for purity, MSI and site)
source(file.path(Sys.getenv("GLP2R_SCRIPTS", "W:/Re_build/05_Analysis/scripts"), "90_fig_helpers.R"))
sink(file.path(LOGD, "93_fig3.txt"), split = TRUE)
io <- readRDS(file.path(RES, "immune_objects.rds")); resI <- io$resI; sc <- io$sc_res; ck <- io$ckr; th <- io$th; st <- io$strat
cl_raw <- "grey55"; cl_par <- "#C8472B"
dumbbell <- function(df, lab, title, xlab = "Spearman \u03c1 with GLP2R (filled: BH-FDR < 0.05)", sz = 5.4) {
  df <- df[order(df$rho), ]; df$lab <- factor(df[[lab]], df[[lab]])
  lg <- rbind(data.frame(lab = df$lab, v = df$rho, type = "Unadjusted", sig = df$FDR_p < 0.05), data.frame(lab = df$lab, v = df$partial_rho, type = "Partial (purity, MSI, site)", sig = df$FDR_partial_p < 0.05))
  ggplot(lg, aes(v, lab)) + geom_vline(xintercept = 0, linewidth = 0.3) + geom_line(aes(group = lab), colour = "grey80", linewidth = 0.4) + geom_point(aes(colour = type, fill = ifelse(sig, type, NA), shape = type), size = 1.3, stroke = 0.5) +
    scale_colour_manual(values = c(Unadjusted = cl_raw, `Partial (purity, MSI, site)` = cl_par), name = NULL) + scale_fill_manual(values = c(Unadjusted = cl_raw, `Partial (purity, MSI, site)` = cl_par), na.value = "white", guide = "none") + scale_shape_manual(values = c(Unadjusted = 21, `Partial (purity, MSI, site)` = 21), guide = "none") +
    labs(x = xlab, y = NULL, title = title) + theme(axis.text.y = element_text(size = sz), legend.position = "bottom", legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-4, 0, 0, 0), panel.grid.major.y = element_line(colour = "grey94", linewidth = 0.25))
}
nsig_par <- sum(resI$FDR_partial_p < 0.05); cat("cell types:", nrow(resI), " partial FDR<0.05:", nsig_par, " max|rho|:", max(abs(resI$rho)), " max|partial|:", max(abs(resI$partial_rho)), "\n")
pA <- dumbbell(resI, "cell_type", sprintf("Immune cell infiltration (ssGSEA, %d signatures)", nrow(resI)))
sc$score2 <- sc$score; pB <- dumbbell(sc, "score2", "MCP-counter / ESTIMATE scores", sz = 5.8)
## checkpoints
ckl <- do.call(rbind, list(data.frame(g = ck$Gene, v = ck$Pearson_r, m = "Pearson r", fdr = ck$Pearson_FDR), data.frame(g = ck$Gene, v = ck$Spearman_rho, m = "Spearman \u03c1", fdr = ck$Spearman_FDR), data.frame(g = ck$Gene, v = ck$partial_rho, m = "Partial \u03c1", fdr = ck$partial_FDR)))
ck_cols <- setNames(c("#3B6FB6", "grey35", cl_par), c("Pearson r", "Spearman \u03c1", "Partial \u03c1")); ckl$g <- factor(ckl$g, rev(ck$Gene)); ckl$m <- factor(ckl$m, c("Pearson r", "Spearman \u03c1", "Partial \u03c1"))
pC <- ggplot(ckl, aes(v, g, colour = m, shape = m, fill = ifelse(fdr < 0.05, as.character(m), NA))) + geom_vline(xintercept = 0, linewidth = 0.3) + geom_point(size = 1.5, stroke = 0.5, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = ck_cols, name = NULL) + scale_fill_manual(values = ck_cols, na.value = "white", guide = "none") + scale_shape_manual(values = c(21, 22, 24), name = NULL) +
  scale_x_continuous(limits = c(-0.2, 0.2)) + labs(x = "Correlation with GLP2R (filled: BH-FDR < 0.05)", y = NULL, title = "Immune-checkpoint genes") + theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-4, 0, 0, 0))
## Thorsson subtypes
th$sub <- factor(th$sub); kw <- kruskal.test(glp2r ~ sub, th); mm <- th[!is.na(th$msi_s), ]; p_an <- anova(lm(rank(glp2r) ~ msi_s + site, mm), lm(rank(glp2r) ~ msi_s + site + sub, mm))$`Pr(>F)`[2]; mss <- th[th$msi_s == "MSS/MSI-L", ]; kw2 <- kruskal.test(glp2r ~ sub, droplevels(mss))
cat(sprintf("Thorsson: n=%d KW=%.3g adjANOVA=%.3g MSS KW=%.3g\n", nrow(th), kw$p.value, p_an, kw2$p.value))
cols_sub <- setNames(c("#E41A1C", "#FF7F00", "#4DAF4A", "#377EB8", "#984EA3"), levels(th$sub)[1:5][seq_len(nlevels(th$sub))])
names(cols_sub) <- levels(th$sub)
pD <- box_panel(th, "sub", "glp2r", cols_sub, "GLP2R expression (VST)", "Thorsson pan-cancer immune subtypes", sprintf("Kruskal\u2013Wallis %s\nMSI/site-adjusted (rank-based F-test) %s\nMSS/MSI-L only: %s", plab(kw$p.value), plab(p_an), plab(kw2$p.value)))
## MSI-stratified
sl <- rbind(data.frame(ct = st$cell_type, rho = st$MSS_rho, g = sprintf("MSS/MSI-L (n=%d)", st$MSS_n[1]), p = st$MSS_p), data.frame(ct = st$cell_type, rho = st$MSIH_rho, g = sprintf("MSI-H (n=%d)", st$MSIH_n[1]), p = st$MSIH_p)); sl$ct <- factor(sl$ct, rev(st$cell_type))
sl$g <- factor(sl$g, unique(sl$g)); pE <- ggplot(sl, aes(rho, ct, colour = g, shape = g, fill = ifelse(p < 0.05, as.character(g), NA))) + geom_vline(xintercept = 0, linewidth = 0.3) + geom_point(size = 1.6, stroke = 0.5, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = setNames(c("grey35", "#C8472B"), levels(sl$g)), name = NULL) + scale_fill_manual(values = setNames(c("grey35", "#C8472B"), levels(sl$g)), na.value = "white", guide = "none") + scale_shape_manual(values = c(21, 24), name = NULL) +
  labs(x = "Spearman \u03c1 with GLP2R (filled: P < 0.05, unadjusted)", y = NULL, title = "MSI-stratified analysis") + theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-4, 0, 0, 0))
fig <- wrap_plots(pA, pB, pC, pD, pE, design = "AB\nAB\nAC\nDE") + plot_layout(heights = c(1.1, 0.8, 0.8, 1.2), widths = c(1.05, 1)) + plot_annotation(tag_levels = "A")
save_fig(fig, "Figure3_GLP2R_immune", w = 180, h = 215)
cat("Figure 3 saved\n"); sink()
