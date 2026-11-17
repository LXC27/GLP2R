## 97_fig7.R : Figure 7 - immune-functional readouts and evidence summary
source(file.path(Sys.getenv("GLP2R_SCRIPTS", "W:/Re_build/05_Analysis/scripts"), "90_fig_helpers.R"))
o <- readRDS(file.path(RES, "immune_objects.rds")); d <- o$d
rmi <- read.csv(file.path(RES, "rm_immune.csv")); rmi <- setNames(rmi$value, rmi$key)
sc <- function(v, ylab, title, key, pos = "tl") {
  dd <- d[!is.na(d[[v]]), ]; ct <- cor.test(dd$glp2r, dd[[v]], method = "spearman", exact = FALSE)
  pr <- regmatches(rmi[[key]], regexec("partial rho=(-?[0-9.]+) p=([0-9.e+-]+)", rmi[[key]]))[[1]]
  pp <- as.numeric(pr[3]); ptxt_adj <- plab(pp)
  fr <- function(r) ifelse(abs(r) < 0.01, sprintf("%.3f", r), sprintf("%.2f", r)); lab <- sprintf("Spearman \u03c1 = %s (%s)\npartial\u1d43 \u03c1 = %s (%s)", fr(ct$estimate), plab(ct$p.value), fr(as.numeric(pr[2])), ptxt_adj)
  ggplot(dd, aes(glp2r, .data[[v]], colour = msi)) + geom_point(size = 0.45, alpha = 0.55, stroke = 0) + geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "black", linewidth = 0.4, formula = y ~ x) +
    scale_colour_manual(values = col_msi, breaks = names(col_msi), na.value = "grey85", name = NULL) + scale_y_continuous(expand = expansion(mult = c(0.03, 0.3))) +
    annotate("text", x = min(dd$glp2r), y = Inf, label = lab, hjust = 0, vjust = 1.3, size = TXT, lineheight = 0.95) + labs(x = "GLP2R (VST)", y = ylab, title = title) +
    theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-4, 0, 0, 0), legend.key.height = unit(2.5, "mm"))
}
pA <- sc("cd40lg", "CD40LG (VST)", "CD40LG", "rho_cd40lg")
pB <- sc("ifng", "IFN-\u03b3 signature (mean VST)", "IFN-\u03b3 response signature", "rho_ifng")
pC <- sc("gep", "T-cell-inflamed GEP (mean VST)", "T-cell-inflamed GEP (18 genes)", "rho_gep")

## D: evidence summary
V <- c(S = "Supported", W = "Weak / descriptive", N = "Not supported")
vcol <- c(S = "#1B7837", W = "#E0A100", N = "#B2182B")
ev <- data.frame(
  dom = c("Expression", "Expression", "Expression", "Prognosis", "Prognosis", "Validation", "Immune", "Immune", "Co-expression", "Epigenetics", "Copy number / mutation", "TMB", "Signature"),
  txt = c("Lower in tumour vs normal: log\u2082FC \u22124.07 (FDR 2\u00d710\u207b\u2076\u2076); paired 32/32 pairs lower, log\u2082FC \u22123.82",
          "Lower in MSI-H than MSS/MSI-L (AUC 0.19; P 3.5\u00d710\u207b\u00b9\u2077)",
          "No association with stage (Spearman \u03c1 0.05; P 0.31)",
          "Unadjusted OS: HR per SD 0.82 (0.68\u20130.98), P 0.032",
          "After age, sex, stage, site, MSI: HR 0.86 (0.69\u20131.07), P 0.18",
          "Adjusted HR 0.83 (TCGA), 1.00 (GSE39582), 0.71 (GSE17536); 3-cohort summary 0.86 (0.70\u20131.04), I\u00b2 65%",
          "ssGSEA: adjusted |\u03c1| \u2264 0.24 (NK 0.14); MCP-counter: endothelial \u03c1 0.42, cytotoxic and NK scores null",
          "Checkpoint genes: all Spearman FDR > 0.4; none significant after adjustment",
          "In the stromal/EMT turquoise module; non-hub gene (kME rank 1,516 of 1,747 module genes)",
          "8 promoter CpGs: window \u0394\u03b2 \u22120.03; adjusted \u03c1 \u22120.06 (P 0.29); decision rule not met",
          "Deep deletion 1.6% (any loss 58.9%); GLP2R mutated in 2.2%; adjusted CNV \u03c1 0.00",
          "\u03c1 \u22120.29 falls to \u22120.13 after MSI/site adjustment (MSS-only \u22120.14)",
          "8-gene LASSO: held-out test C-index 0.518; GSE39582 0.548; GSE17536 0.545"),
  lab = c("Tumour loss of expression", "MSI-H association", "Stage association", "Prognostic association", "Independent prognostic value", "External validation (microarray)", "Immune-cell infiltration", "Immune-checkpoint genes", "WGCNA context", "Promoter hypermethylation", "Copy-number loss / somatic mutation as cause", "Mutational-burden link", "Prognostic signature"),
  v = c("S", "S", "N", "W", "N", "N", "W", "N", "W", "N", "N", "W", "N"), stringsAsFactors = FALSE)
ev$txt <- vapply(strwrap(ev$txt, 86, simplify = FALSE), function(z) paste(z, collapse = "\n"), "")
ev$lab <- vapply(strwrap(ev$lab, 26, simplify = FALSE), function(z) paste(z, collapse = "\n"), "")
ev$y <- rev(seq_len(nrow(ev))); ev$stripe <- ev$y %% 2 == 0
xl <- c(0, 40, 140); xv <- 142
pD <- ggplot(ev) + geom_rect(data = ev[ev$stripe, ], aes(xmin = -1, xmax = 175, ymin = y - 0.5, ymax = y + 0.5), fill = "grey95", colour = NA) +
  geom_text(aes(x = xl[1], y = y, label = lab), hjust = 0, size = TXT, lineheight = 0.9, fontface = "bold") + geom_text(aes(x = xl[2], y = y, label = txt), hjust = 0, size = TXT, lineheight = 0.9) +
  geom_tile(aes(x = xv, y = y, fill = v), width = 2.2, height = 0.45) + geom_text(aes(x = xv + 3, y = y, label = V[v]), hjust = 0, size = TXT) + scale_fill_manual(values = vcol, guide = "none") +
  annotate("text", x = c(xl[1], xl[2], xv), y = nrow(ev) + 0.95, label = c("Question", "Key result (TCGA-COADREAD unless stated)", "Verdict"), hjust = 0, size = TXT, fontface = "bold") +
  annotate("segment", x = -1, xend = 175, y = nrow(ev) + 0.55, yend = nrow(ev) + 0.55, linewidth = 0.3) +
  scale_x_continuous(limits = c(-1, 175), expand = c(0, 0)) + scale_y_continuous(limits = c(0.4, nrow(ev) + 1.3), expand = c(0, 0)) + labs(title = "Evidence summary") + theme_void(base_size = 7) + theme(plot.title = element_text(size = 6.8, face = "bold", margin = margin(0, 0, 2, 0)))
fig <- (pA | pB | pC) / pD + plot_layout(heights = c(0.9, 2)) + plot_annotation(tag_levels = "A")
save_fig(fig, "Figure7_GLP2R_immune_function_evidence", w = 180, h = 175); cat("Figure 7 saved\n")
