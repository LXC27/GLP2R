## 50_pathway.R : GLP2R-correlated transcriptome in tumours (purity/MSI/site-adjusted partial Spearman ranking) -> fgsea (Hallmark, KEGG legacy, GO BP, C8 cell-type signatures); candidate pathway genes; compartment-marker correlations
source(file.path(Sys.getenv("GLP2R_SCRIPTS", "W:/Re_build/05_Analysis/scripts"), "00_common.R"))
suppressPackageStartupMessages({library(fgsea); library(msigdbr)})
sink(file.path(LOGD, "50_pathway.txt"), split = TRUE)
ex <- readRDS(file.path(RES, "expr_tests.rds")); d <- ex$d; vt <- readRDS(file.path(DC, "vst_tumor_591.rds"))[, rownames(d)]
vsa <- readRDS(file.path(DC, "vst_all_624.rds")); cl <- readRDS(file.path(DC, "cohort.rds"))$clin
RM <- list(); add <- function(k, v, desc = "") RM[[length(RM) + 1]] <<- data.frame(key = k, value = as.character(v), description = desc)
ok <- complete.cases(d[, c("purity_absolute", "msi", "site")]); cat("tumours with purity+MSI+site:", sum(ok), "\n"); add("pathway_adj_n", sum(ok))
## ---- genome-wide partial Spearman (matrix form)
R <- apply(vt[, ok], 1, rank); y <- rank(d$glp2r[ok]); Z <- model.matrix(~ purity_absolute + msi + site, d[ok, ])[, -1, drop = FALSE]
Q <- qr.Q(qr(cbind(1, Z))); res_mat <- function(M) M - Q %*% (t(Q) %*% M)
Rr <- res_mat(R); yr <- as.numeric(res_mat(matrix(y, ncol = 1)))
Rc <- sweep(Rr, 2, colMeans(Rr)); yc <- yr - mean(yr); r_adj <- as.numeric(crossprod(Rc, yc)) / (sqrt(colSums(Rc^2)) * sqrt(sum(yc^2))); names(r_adj) <- rownames(vt)
r_raw <- apply(vt, 1, function(x) cor(x, d$glp2r, method = "spearman")); r_raw <- r_raw[names(r_adj)]
dfree <- sum(ok) - 2 - ncol(Z); tt <- r_adj * sqrt(dfree / (1 - r_adj^2)); p_adj <- 2 * pt(-abs(tt), dfree); cg <- data.frame(gene = names(r_adj), rho_partial = r_adj, rho_raw = r_raw[names(r_adj)], p_partial = p_adj, FDR_partial = p.adjust(p_adj, "BH"))
cg <- cg[order(-abs(cg$rho_partial)), ]; cg <- cg[cg$gene != "GLP2R", ]; write_csv_utf8(cg, file.path(RES, "glp2r_correlated_genes.csv"))
cat("genes tested:", nrow(cg), "; |partial rho|>=0.3 & FDR<0.05:", sum(abs(cg$rho_partial) >= 0.3 & cg$FDR_partial < 0.05), "; FDR<0.05:", sum(cg$FDR_partial < 0.05), "\n")
cat("Top positive:\n"); print(head(cg[order(-cg$rho_partial), ], 15), digits = 3); cat("Top negative:\n"); print(head(cg[order(cg$rho_partial), ], 10), digits = 3)
add("n_genes_corr_tested", nrow(cg)); add("n_genes_partial_FDR05", sum(cg$FDR_partial < 0.05)); add("n_genes_abs_rho_ge_0.3_FDR05", sum(abs(cg$rho_partial) >= 0.3 & cg$FDR_partial < 0.05)); add("max_abs_partial_rho", round(max(abs(cg$rho_partial)), 3))
add("top_positive_genes", paste(head(cg$gene[order(-cg$rho_partial)], 10), collapse = ", ")); add("top_negative_genes", paste(head(cg$gene[order(cg$rho_partial)], 10), collapse = ", "))
## ---- fgsea
stats_adj <- sort(setNames(cg$rho_partial, cg$gene), decreasing = TRUE); stats_raw <- sort(setNames(cg$rho_raw, cg$gene), decreasing = TRUE)
getsets <- function(coll, sub = NULL) { m <- if (is.null(sub)) msigdbr(species = "Homo sapiens", collection = coll) else msigdbr(species = "Homo sapiens", collection = coll, subcollection = sub); split(m$gene_symbol, m$gs_name) }
cols <- list(Hallmark = getsets("H"), KEGG = getsets("C2", "CP:KEGG_LEGACY"), GO_BP = getsets("C5", "GO:BP"), C8_celltype = getsets("C8"))
out <- list(); set.seed(SEED)
for (nm in names(cols)) { pw <- cols[[nm]]; for (st in c("adjusted", "unadjusted")) { s <- if (st == "adjusted") stats_adj else stats_raw
    f <- fgseaMultilevel(pw, s, minSize = 15, maxSize = 500, eps = 0, nproc = 1); f$collection <- nm; f$ranking <- st; f$leadingEdge <- sapply(f$leadingEdge, function(z) paste(head(z, 12), collapse = ";")); f <- f[order(f$pval), ]
    f$FDR_collection <- p.adjust(f$pval, "BH"); out[[paste(nm, st)]] <- as.data.frame(f); cat(sprintf("%s (%s): %d sets, %d at FDR<0.05\n", nm, st, nrow(f), sum(f$FDR_collection < 0.05, na.rm = TRUE))) } }
allf <- do.call(rbind, out); rownames(allf) <- NULL; write_csv_utf8(allf, file.path(RES, "fgsea_all.csv"))
for (nm in c("Hallmark", "KEGG", "C8_celltype", "GO_BP")) { z <- allf[allf$collection == nm & allf$ranking == "adjusted", ]; z <- z[order(z$pval), ]; cat("\n---", nm, "(adjusted ranking), top 12 by P ---\n"); print(head(z[, c("pathway", "NES", "size", "pval", "FDR_collection")], 12), digits = 3)
  add(paste0("fgsea_", nm, "_n_FDR05_adj"), sum(z$FDR_collection < 0.05, na.rm = TRUE)); add(paste0("fgsea_", nm, "_top_pos"), paste(head(z$pathway[z$FDR_collection < 0.05 & z$NES > 0][order(-z$NES[z$FDR_collection < 0.05 & z$NES > 0])], 5), collapse = "; "))
  add(paste0("fgsea_", nm, "_top_neg"), paste(head(z$pathway[z$FDR_collection < 0.05 & z$NES < 0][order(z$NES[z$FDR_collection < 0.05 & z$NES < 0])], 5), collapse = "; ")) }
## ---- candidate pathway genes (the 50-gene set of Table S1) within tumours
s1 <- read.csv(file.path(RES, "TableS1_candidate_genes_DE.csv"), fileEncoding = "UTF-8"); gs1 <- intersect(setdiff(s1$Gene, "GLP2R"), rownames(vt))
cc <- do.call(rbind, lapply(gs1, function(s) { a <- cor.test(d$glp2r, vt[s, rownames(d)], method = "spearman", exact = FALSE); data.frame(Gene = s, rho = a$estimate, p = a$p.value, partial_rho = cg$rho_partial[match(s, cg$gene)], partial_p = cg$p_partial[match(s, cg$gene)], row.names = NULL) }))
cc$FDR <- p.adjust(cc$p, "BH"); cc$partial_FDR <- p.adjust(cc$partial_p, "BH"); cc <- cc[order(cc$p), ]; write_csv_utf8(cc, file.path(RES, "candidate_genes_correlation.csv")); cat("\nTable S1 genes correlated with GLP2R in tumours (n tested =", nrow(cc), "), FDR<0.05:", sum(cc$FDR < 0.05), " partial FDR<0.05:", sum(cc$partial_FDR < 0.05), "\n"); print(head(cc, 10), digits = 3)
add("candidate_genes_n_tested", nrow(cc)); add("candidate_genes_n_FDR05", sum(cc$FDR < 0.05)); add("candidate_genes_n_partial_FDR05", sum(cc$partial_FDR < 0.05))
## ---- compartment markers: where does tumour GLP2R signal come from?
mk <- list(Epithelial = c("EPCAM", "CDX2", "KRT20", "VIL1"), `Enteroendocrine/L-cell` = c("CHGA", "PYY", "GCG", "NEUROD1"), `Enteric neuron/glia` = c("ELAVL4", "PHOX2B", "HAND2", "UCHL1", "S100B", "PLP1", "GFAP"), Fibroblast = c("COL1A1", "COL1A2", "PDGFRB", "DCN", "LUM"),
           `Smooth muscle` = c("ACTA2", "DES", "MYH11", "CNN1"), Endothelial = c("PECAM1", "VWF", "CDH5", "CLDN5"), `T cell` = c("CD3D", "CD3E", "CD8A"), Myeloid = c("CD68", "CD163", "ITGAM"))
mr <- do.call(rbind, lapply(names(mk), function(cn) do.call(rbind, lapply(intersect(mk[[cn]], rownames(vt)), function(s) { a <- cor.test(d$glp2r, vt[s, rownames(d)], method = "spearman", exact = FALSE); data.frame(compartment = cn, gene = s, rho_tumour = a$estimate, p_tumour = a$p.value, partial_rho = cg$rho_partial[match(s, cg$gene)], row.names = NULL) }))))
nm_ <- rownames(cl)[cl$sample_type == "Solid Tissue Normal"]; mr$rho_normal <- sapply(mr$gene, function(s) if (s %in% rownames(vsa)) cor(vsa["GLP2R", nm_], vsa[s, nm_], method = "spearman") else NA); mr$FDR_tumour <- p.adjust(mr$p_tumour, "BH")
write_csv_utf8(mr, file.path(RES, "glp2r_compartment_markers.csv")); cat("\nCompartment markers (rho with GLP2R):\n"); print(aggregate(cbind(rho_tumour, partial_rho, rho_normal) ~ compartment, mr, mean), digits = 3)
ag <- aggregate(rho_tumour ~ compartment, mr, mean); for (i in seq_len(nrow(ag))) add(paste0("compartment_mean_rho_", gsub("[^A-Za-z]+", "_", ag$compartment[i])), round(ag$rho_tumour[i], 3))
saveRDS(list(cg = cg, allf = allf, cc = cc, mr = mr), file.path(RES, "pathway_objects.rds")); rm_write(do.call(rbind, RM), "pathway"); log_msg("done"); sink()
