## 91_fig1.R : Figure 1 - GLP2R loss in CRC (unpaired + paired), genome-wide DE, clinicopathological associations, candidate-gene panel
source(file.path(Sys.getenv("GLP2R_SCRIPTS", "W:/Re_build/05_Analysis/scripts"), "90_fig_helpers.R"))
sink(file.path(LOGD, "91_fig1.txt"), split = TRUE)
ex <- readRDS(file.path(RES, "expr_tests.rds")); d <- ex$d; tab1 <- ex$tabS1; vsa <- readRDS(file.path(DC, "vst_all_624.rds")); de <- readRDS(file.path(DC, "de_tumor_normal.rds")); dp <- readRDS(file.path(DC, "de_paired.rds"))
pr <- read.csv(file.path(RES, "paired_glp2r_vst.csv")); rmx <- setNames(read.csv(file.path(RES, "rm_expr.csv"))$value, read.csv(file.path(RES, "rm_expr.csv"))$key)
nrm <- colnames(vsa)[substr(colnames(vsa), 14, 15) == "11"]; tum <- rownames(d)
dtn <- data.frame(group = c(rep("Normal", length(nrm)), rep("Tumor", length(tum))), vst = c(vsa["GLP2R", nrm], vsa["GLP2R", tum])); dtn$group <- factor(dtn$group, c("Normal", "Tumor"))
g2 <- de[which(de$symbol == "GLP2R"), ]; w <- wilcox.test(vst ~ group, dtn); lab_a <- sprintf("DESeq2 log\u2082FC = %.2f\nFDR = %s\nWilcoxon %s", g2$log2FoldChange, ptxt(g2$padj), plab(w$p.value))
pA <- box_panel(dtn, "group", "vst", col_tn, "GLP2R expression (VST)", "Tumour vs normal (unpaired)", lab_a)
## paired
pr$id <- seq_len(nrow(pr)); pl <- rbind(data.frame(id = pr$id, g = "Normal", v = pr$vst.Normal), data.frame(id = pr$id, g = "Tumor", v = pr$vst.Tumor)); pl$g <- factor(pl$g, c("Normal", "Tumor"))
pg <- dp[which(dp$symbol == "GLP2R"), ]; wp <- wilcox.test(pr$vst.Normal, pr$vst.Tumor, paired = TRUE); nlow <- sum(pr$vst.Tumor < pr$vst.Normal)
pB <- ggplot(pl, aes(g, v)) + geom_line(aes(group = id), colour = "grey70", linewidth = 0.3) + geom_point(aes(colour = g), size = 0.9) + scale_colour_manual(values = col_tn) + scale_x_discrete(labels = function(x) paste0(x, "\n(n=", nrow(pr), ")")) +
  labs(x = NULL, y = "GLP2R expression (VST)", title = "Matched pairs (sensitivity)") + theme(legend.position = "none") + scale_y_continuous(expand = expansion(mult = c(0.03, 0.45))) +
  annotate("text", x = 1.5, y = max(pl$v) + 0.2, label = sprintf("lower in tumour: %d/%d pairs\npaired DESeq2 log\u2082FC = %.2f\nFDR = %s", nlow, nrow(pr), pg$log2FoldChange, ptxt(pg$padj)), size = TXT, vjust = 0, lineheight = 0.95)
## volcano
vo <- de[!is.na(de$padj), ]; vo$cls <- ifelse(vo$padj < 0.05 & vo$log2FoldChange > 1, "Up", ifelse(vo$padj < 0.05 & vo$log2FoldChange < -1, "Down", "NS")); vo$y <- pmin(-log10(vo$padj), 300); nu <- sum(vo$cls == "Up"); nd <- sum(vo$cls == "Down")
cat("DEGs:", nu + nd, "up", nu, "down", nd, "tested", nrow(vo), "\n"); lab <- vo[vo$symbol %in% c("GLP2R", "GLP1R", "GCG", "PRKG2"), ]
pC <- ggplot(vo, aes(log2FoldChange, y, colour = cls)) + geom_point(size = 0.25, alpha = 0.5, stroke = 0) + scale_colour_manual(values = col_dir, guide = "none") + geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = 0.25) + geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.25) +
  geom_point(data = lab, shape = 21, fill = "black", colour = "black", size = 1.2) + geom_text_repel(data = lab, aes(label = symbol), colour = "black", size = TXT, min.segment.length = 0, box.padding = 0.4, max.overlaps = 50, seed = 1) +
  labs(x = "log\u2082 fold change (tumour / normal)", y = expression(-log[10]~FDR), title = "Genome-wide differential expression") + annotate("text", x = -9.8, y = 292, label = sprintf("Down %s", format(nd, big.mark = ",")), colour = col_dir["Down"], hjust = 0, size = TXT) + annotate("text", x = 9.8, y = 292, label = sprintf("Up %s", format(nu, big.mark = ",")), colour = col_dir["Up"], hjust = 1, size = TXT) + coord_cartesian(xlim = c(-10, 10))
## stage / MSI / site
d$stage <- factor(d$stage_grp, c("I", "II", "III", "IV")); sk <- kruskal.test(glp2r ~ stage, d); sp <- cor.test(as.numeric(d$stage), d$glp2r, method = "spearman", exact = FALSE)
pD <- box_panel(d[!is.na(d$stage), ], "stage", "glp2r", setNames(c("#9ECAE1", "#6BAED6", "#3182BD", "#08519C"), c("I", "II", "III", "IV")), "GLP2R expression (VST)", "AJCC stage", sprintf("Kruskal\u2013Wallis %s\nSpearman \u03c1 = %.2f (%s)", plab(sk$p.value), sp$estimate, plab(sp$p.value)))
d$msi2 <- factor(d$msi, c("MSS/MSI-L", "MSI-H")); wm <- wilcox.test(glp2r ~ msi2, d); am <- auc_mw(d$glp2r[d$msi2 == "MSI-H"], d$glp2r[d$msi2 == "MSS/MSI-L"])
pE <- box_panel(d[!is.na(d$msi2), ], "msi2", "glp2r", col_msi, "GLP2R expression (VST)", "MSI status", sprintf("Wilcoxon %s\nAUC (MSI-H > MSS) = %.2f", plab(wm$p.value), am))
d$site2 <- factor(d$site, c("COAD", "READ")); ws <- wilcox.test(glp2r ~ site2, d)
pF <- box_panel(d, "site2", "glp2r", col_site, "GLP2R expression (VST)", "Tumour site", sprintf("Wilcoxon %s", plab(ws$p.value)))
cs <- cor.test(d$glp1r, d$glp2r, method = "spearman", exact = FALSE); b1 <- tab1[tab1$Gene == "GLP1R", ]
pG <- ggplot(d, aes(glp2r, glp1r)) + geom_point(size = 0.4, alpha = 0.5, stroke = 0, colour = "grey30") + geom_smooth(method = "lm", se = FALSE, colour = "#C8472B", linewidth = 0.5, formula = y ~ x) + labs(x = "GLP2R (VST)", y = "GLP1R (VST)", title = "GLP1R vs GLP2R (tumours)") +
  annotate("text", x = min(d$glp2r), y = max(d$glp1r), label = sprintf("Spearman \u03c1 = %.2f, %s\nGLP1R baseMean = %.1f\n(near detection limit)", cs$estimate, plab(cs$p.value), b1$baseMean), hjust = 0, vjust = 1, size = TXT, lineheight = 0.95)
## candidate genes (Table S1)
tc <- tab1; tc$lo <- tc$log2FC - 1.96 * tc$lfcSE; tc$hi <- tc$log2FC + 1.96 * tc$lfcSE; tc$cls <- ifelse(tc$Expressed == "yes", ifelse(!is.na(tc$FDR) & tc$FDR < 0.05, "Expressed, FDR < 0.05", "Expressed, n.s."), "Low expression (baseMean < 10 or filtered)")
tc$cls[is.na(tc$FDR)] <- "Low expression (baseMean < 10 or filtered)"; tc <- tc[!is.na(tc$log2FC), ]; tc <- tc[order(tc$log2FC), ]; tc$Gene <- factor(tc$Gene, tc$Gene)
pH <- ggplot(tc, aes(Gene, log2FC, colour = cls)) + geom_hline(yintercept = 0, linewidth = 0.3) + geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.3) + geom_point(size = 1) +
  scale_colour_manual(values = c("Expressed, FDR < 0.05" = "#C8472B", "Expressed, n.s." = "grey40", "Low expression (baseMean < 10 or filtered)" = "#E6A117"), name = NULL) +
  labs(x = NULL, y = "log\u2082FC (tumour / normal)", title = sprintf("Candidate GLP-signalling genes (n = %d tested; %d not testable owing to filtering)", nrow(tc), nrow(tab1) - nrow(tc))) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5.5), legend.position = "top", legend.key.height = unit(2.5, "mm"), legend.margin = margin(0, 0, 0, 0))
fig <- (pA | pB | pC) / (pD + theme(axis.text.x = element_text(size = 5.5)) | pE | pF | pG) / pH + plot_layout(heights = c(1.15, 1, 0.9)) + plot_annotation(tag_levels = "A")
save_fig(fig, "Figure1_GLP2R_expression", w = 180, h = 200)
cat("Figure 1 saved\n"); cat("stage KW p:", sk$p.value, " MSI p:", wm$p.value, " site p:", ws$p.value, " GLP1R-GLP2R rho:", cs$estimate, "\n"); sink()
