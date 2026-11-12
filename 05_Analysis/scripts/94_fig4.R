## 94_fig4.R : Figure 4 - pathway, compartment and co-expression-module context of GLP2R
source(file.path(Sys.getenv("GLP2R_SCRIPTS", "W:/Re_build/05_Analysis/scripts"), "90_fig_helpers.R"))
sink(file.path(LOGD, "94_fig4.txt"), split = TRUE)
pw <- readRDS(file.path(RES, "pathway_objects.rds")); w <- readRDS(file.path(RES, "wgcna_objects.rds"))
cl_pos <- "#C8472B"; cl_neg <- "#3B6FB6"
clean <- function(s, n = 46) { s <- sub("^(HALLMARK|KEGG|GOBP|REACTOME)_", "", s); s <- tolower(gsub("_", " ", s)); rep <- c("\\buv\\b" = "UV", "\\bdna\\b" = "DNA", "\\bil6\\b" = "IL6", "\\bjak\\b" = "JAK", "\\bstat3\\b" = "STAT3", "\\bkras\\b" = "KRAS", "\\bmyc\\b" = "MYC", "\\be2f\\b" = "E2F", "\\bg2m\\b" = "G2M", "\\bmtorc1\\b" = "mTORC1", "\\bdn\\b" = "down", "\\bv1\\b" = "V1", "\\bv2\\b" = "V2", "\\bgao\\b" = "GAO", "\\bcj\\b" = "CJ", "\\bhu\\b" = "HU", "\\bfan\\b" = "FAN", "\\bmuraro\\b" = "MURARO", "\\bjones\\b" = "JONES", "\\bctx\\b" = "CTX"); for (k in names(rep)) s <- gsub(k, rep[[k]], s); s <- paste0(toupper(substr(s, 1, 1)), substring(s, 2)); ifelse(nchar(s) > n, paste0(substr(s, 1, n - 1), "\u2026"), s) }
## A  Hallmark (partial-ranking)
hm <- pw$allf[pw$allf$collection == "Hallmark" & pw$allf$ranking == "adjusted", ]; hm <- hm[hm$padj < 0.05, ]
hs <- rbind(head(hm[order(-hm$NES), ], 8), head(hm[order(hm$NES), ], 8)); hs <- hs[order(hs$NES), ]; hs$lab <- factor(clean(hs$pathway), clean(hs$pathway))
cat("Hallmark FDR<0.05 (adjusted ranking):", nrow(hm), " of 50\n")
pA <- ggplot(hs, aes(NES, lab, fill = NES > 0)) + geom_col(width = 0.7) + geom_vline(xintercept = 0, linewidth = 0.3) + scale_fill_manual(values = c(`TRUE` = cl_pos, `FALSE` = cl_neg), guide = "none") +
  labs(x = "Normalised enrichment score (GSEA)", y = NULL, title = "MSigDB Hallmark (ranked by partial \u03c1; FDR < 0.05)") + theme(axis.text.y = element_text(size = 5.6))
## B  compartment markers
mr <- pw$mr; mr$compartment <- factor(mr$compartment); ord <- tapply(mr$rho_tumour, mr$compartment, mean); mr$compartment <- factor(mr$compartment, names(sort(ord)))
ml <- rbind(data.frame(comp = mr$compartment, gene = mr$gene, v = mr$rho_tumour, g = "Tumour (n = 591)"), data.frame(comp = mr$compartment, gene = mr$gene, v = mr$rho_normal, g = "Normal (n = 33)"))
mm <- aggregate(v ~ comp + g, ml, mean)
pB <- ggplot(ml, aes(v, comp, colour = g)) + geom_vline(xintercept = 0, linewidth = 0.3) + geom_point(alpha = 0.45, size = 0.8, shape = 16, position = position_dodge(width = 0.55)) + geom_point(data = mm, shape = 18, size = 2.6, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = c(`Tumour (n = 591)` = cl_pos, `Normal (n = 33)` = cl_neg), name = NULL) + labs(x = "Spearman \u03c1 with GLP2R", y = NULL, title = "Compartment marker genes") +
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-4, 0, 0, 0), axis.text.y = element_text(size = 5.8))
print(mm)
## C  top genes
cg <- pw$cg; top <- rbind(head(cg[order(-cg$rho_partial), ], 12), head(cg[order(cg$rho_partial), ], 12)); top <- top[order(top$rho_partial), ]; top$gene <- factor(top$gene, top$gene)
nsig <- sum(cg$FDR_partial < 0.05); n3 <- sum(abs(cg$rho_partial) >= 0.3); cat("genes:", nrow(cg), " FDR<0.05:", nsig, " |rho|>=0.3:", n3, "\n")
pC <- ggplot(top, aes(rho_partial, gene, fill = rho_partial > 0)) + geom_col(width = 0.7) + geom_vline(xintercept = 0, linewidth = 0.3) + scale_fill_manual(values = c(`TRUE` = cl_pos, `FALSE` = cl_neg), guide = "none") +
  labs(x = "Partial Spearman \u03c1 with GLP2R (tumours, n = 561)", y = NULL, title = sprintf("Strongest co-expressed genes (%s tested)", format(nrow(cg), big.mark = ","))) + theme(axis.text.y = element_text(size = 5.6))
## D  candidate genes
cc <- pw$cc; cc <- cc[order(cc$partial_rho), ]; cc$Gene <- factor(cc$Gene, cc$Gene)
dl <- rbind(data.frame(g = cc$Gene, v = cc$rho, t = "Unadjusted", sig = cc$FDR < 0.05), data.frame(g = cc$Gene, v = cc$partial_rho, t = "Partial (purity, MSI, site)", sig = cc$partial_FDR < 0.05))
cat("candidates:", nrow(cc), " FDR<0.05 raw:", sum(cc$FDR < 0.05), " partial:", sum(cc$partial_FDR < 0.05), "\n")
pD <- ggplot(dl, aes(v, g)) + geom_vline(xintercept = 0, linewidth = 0.3) + geom_line(aes(group = g), colour = "grey82", linewidth = 0.35) + geom_point(aes(colour = t, fill = ifelse(sig, t, NA)), shape = 21, size = 1.2, stroke = 0.45) +
  scale_colour_manual(values = c(Unadjusted = "grey55", `Partial (purity, MSI, site)` = cl_pos), name = NULL) + scale_fill_manual(values = c(Unadjusted = "grey55", `Partial (purity, MSI, site)` = cl_pos), na.value = "white", guide = "none") +
  labs(x = "Spearman \u03c1 (filled: FDR < 0.05)", y = NULL, title = "Candidate genes (47 tested)") + theme(axis.text.y = element_text(size = 5), legend.position = "bottom", legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-4, 0, 0, 0), panel.grid.major.y = element_line(colour = "grey95", linewidth = 0.2))
## E  WGCNA modules
mt <- w$modtab; gm <- w$gl_mod; mt$lab <- paste0(mt$module, " (", prettyNum(mt$size, big.mark = ","), ")", ifelse(mt$module == gm, " \u2190 GLP2R", "")); mt <- mt[order(mt$rho_GLP2R), ]; mt$lab <- factor(mt$lab, mt$lab)
cat("GLP2R module:", gm, "\n")
pE <- ggplot(mt, aes(rho_GLP2R, lab, fill = module)) + geom_col(width = 0.72, colour = "grey30", linewidth = 0.2) + geom_point(aes(x = rho_partial), shape = 21, fill = "white", size = 1.2, stroke = 0.4) + geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_fill_identity() + labs(x = "Eigengene vs GLP2R, \u03c1 (bar: raw; circle: partial)", y = NULL, title = "WGCNA modules (591 tumours)") + theme(axis.text.y = element_text(size = 5.6))
## F  turquoise enrichment
en <- w$en; en <- en[en$module == gm & en$collection %in% c("Hallmark", "GO_BP", "C8"), ]
tp <- do.call(rbind, lapply(split(en, en$collection), function(x) head(x[order(x$FDR), ], 4))); tp$coll <- c(Hallmark = "Hallmark", GO_BP = "GO-BP", C8 = "Cell type (C8)")[tp$collection]
tp$lab <- paste0(clean(tp$set, 38), " [", tp$overlap, "/", tp$set_size, "]"); tp$frac <- 100 * tp$overlap / tp$set_size; tp$fl <- ifelse(tp$FDR < 1e-300, "FDR < 1e-300", sprintf("FDR %.0e", tp$FDR))
tp <- tp[order(factor(tp$coll, c("Hallmark", "GO-BP", "Cell type (C8)")), tp$FDR, decreasing = c(TRUE, TRUE), method = "radix"), ]; tp$lab <- factor(tp$lab, tp$lab)
pF <- ggplot(tp, aes(frac, lab, fill = coll)) + geom_col(width = 0.7) + geom_text(aes(label = fl), hjust = -0.08, size = 1.8) + scale_fill_manual(values = c(Hallmark = "#7A4EAB", `GO-BP` = "#2E8B7A", `Cell type (C8)` = "#D98C1F"), name = NULL) + scale_x_continuous(limits = c(0, 125), breaks = c(0, 25, 50, 75, 100), expand = c(0, 0)) +
  labs(x = "Gene-set members in module (%)", y = NULL, title = sprintf("Enrichment of the %s module", gm)) + theme(axis.text.y = element_text(size = 5.2), legend.position = "bottom", legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-4, 0, 0, 0), legend.key.size = unit(2.5, "mm"))
des <- c(patchwork::area(1, 1, 4, 1), patchwork::area(1, 2, 3, 2), patchwork::area(5, 1, 9, 1), patchwork::area(4, 2, 9, 2), patchwork::area(10, 1, 12, 1), patchwork::area(10, 2, 12, 2))
fig <- wrap_plots(pA, pB, pC, pD, pE, pF, design = des) + plot_annotation(tag_levels = "A")
save_fig(fig, "Figure4_GLP2R_pathways", w = 180, h = 210)
cat("Figure 4 saved\n"); sink()
