## 95_fig5.R : Figure 5 - exploratory LASSO-Cox signature: derivation, stability and (lack of) external validation
source(file.path(Sys.getenv("GLP2R_SCRIPTS", "W:/Re_build/05_Analysis/scripts"), "90_fig_helpers.R"))
suppressPackageStartupMessages(library(glmnet))
sink(file.path(LOGD, "95_fig5.txt"), split = TRUE)
o <- readRDS(file.path(RES, "lasso_objects.rds")); rm <- read.csv(file.path(RES, "rm_lasso.csv")); rmv <- setNames(rm$value, rm$key)
cl_sel <- "#C8472B"; sel <- names(o$sel_min)[abs(o$sel_min) > 1e-6]; cat("selected (non-zero):", paste(sel, collapse = ","), "\n")
## A coefficient path
B <- as.matrix(o$fit$beta); lam <- o$fit$lambda; pd <- data.frame(gene = rep(rownames(B), ncol(B)), ll = rep(log(lam), each = nrow(B)), b = as.vector(B)); pd$s <- pd$gene %in% sel
lab <- pd[pd$s & pd$ll == min(pd$ll), ]; lab <- lab[order(lab$b), ]; lab$yl <- lab$b; for (k in seq_len(nrow(lab))[-1]) lab$yl[k] <- max(lab$yl[k], lab$yl[k - 1] + 0.075); lmin_l <- log(o$lmin); l1_l <- log(o$l1se)
pA <- ggplot(pd, aes(ll, b, group = gene)) + geom_hline(yintercept = 0, linewidth = 0.3) + geom_line(data = pd[!pd$s, ], colour = "grey75", linewidth = 0.3) + geom_line(data = pd[pd$s, ], aes(colour = gene), linewidth = 0.5, show.legend = FALSE) +
  geom_vline(xintercept = c(lmin_l, l1_l), linetype = c(2, 3), linewidth = 0.3) + geom_segment(data = lab, aes(x = ll, xend = ll - 0.12, y = b, yend = yl, colour = gene), linewidth = 0.2, show.legend = FALSE) + geom_text(data = lab, aes(x = ll - 0.14, y = yl, label = gene, colour = gene), hjust = 1, size = 1.9, show.legend = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0.3, 0.02))) + annotate("text", x = lmin_l, y = max(pd$b) * 0.98, label = "\u03bb(min)", hjust = 1.1, size = 1.9) + annotate("text", x = l1_l, y = max(pd$b) * 0.98, label = "\u03bb(1-SE)", hjust = 1.1, size = 1.9) +
  labs(x = expression(log(lambda)), y = "Coefficient (standardised)", title = "LASSO-Cox path (45 genes; TCGA train, n = 396)")
## B cross-validation curve
cv <- o$cvs; cd <- data.frame(ll = log(cv$lambda), m = cv$cvm, lo = cv$cvlo, hi = cv$cvup, nz = cv$nzero)
pB <- ggplot(cd, aes(ll, m)) + geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey85") + geom_line(colour = cl_sel, linewidth = 0.5) + geom_point(colour = cl_sel, size = 0.7) + geom_vline(xintercept = c(lmin_l, l1_l), linetype = c(2, 3), linewidth = 0.3) +
  annotate("text", x = lmin_l, y = max(cd$hi), label = "\u03bb(min) = median of 20 repeated\n10-fold CVs", hjust = 1.05, vjust = 1, size = 1.9, lineheight = 0.95) + labs(x = expression(log(lambda)), y = "Partial-likelihood deviance", title = "Cross-validation (single 10-fold run shown)")
## C bootstrap selection frequency
sf <- o$sf_tab[order(-o$sf_tab$selection_freq), ]; sf <- head(sf, 15); sf <- sf[order(sf$selection_freq), ]; sf$gene <- factor(sf$gene, sf$gene)
pC <- ggplot(sf, aes(selection_freq * 100, gene, fill = selected_min)) + geom_col(width = 0.7) + geom_vline(xintercept = 80, linetype = 2, linewidth = 0.3) + scale_fill_manual(values = c(`TRUE` = cl_sel, `FALSE` = "grey65"), labels = c(`TRUE` = "selected at \u03bb(min)", `FALSE` = "not selected"), name = NULL) +
  labs(x = "Bootstrap selection frequency (%; B = 200)", y = NULL, title = "Selection stability (top 15 of 45 genes)") + theme(axis.text.y = element_text(size = 5.4), legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank(), legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-4, 0, 0, 0))
## D coefficients
cf <- data.frame(gene = names(o$sel_min), b = as.numeric(o$sel_min)); cf <- cf[abs(cf$b) > 1e-6, ]; cf <- cf[order(cf$b), ]; cf$gene <- factor(cf$gene, cf$gene)
pD <- ggplot(cf, aes(b, gene, fill = b > 0)) + geom_col(width = 0.7) + geom_vline(xintercept = 0, linewidth = 0.3) + scale_fill_manual(values = c(`TRUE` = cl_sel, `FALSE` = "#3B6FB6"), guide = "none") + labs(x = "Coefficient at \u03bb(min)", y = NULL, title = "Eight-gene signature (exploratory)") + theme(axis.text.y = element_text(size = 5.8))
## E-H Kaplan-Meier (train median cut for TCGA; within-cohort median for GEO)
cr <- c(`Low risk` = "#3B6FB6", `High risk` = "#C8472B"); ttl <- c(TCGA_train = "TCGA training", TCGA_test = "TCGA hold-out test", GSE39582 = "GSE39582", GSE17536 = "GSE17536")
kmp <- lapply(names(o$cohorts), function(n) { x <- o$cohorts[[n]]; lr <- survdiff(Surv(time, event) ~ grp, x); p <- 1 - pchisq(lr$chisq, 1); h <- cox_tab(coxph(Surv(time, event) ~ grp, x))[1, ]
  km_panel(x, "grp", cr, xmax = 120, step = 24, title = sprintf("%s (n = %d, %d events)", ttl[[n]], nrow(x), sum(x$event)), annot = sprintf("HR %s\nLog-rank %s", fmt_ci(h$HR, h$lo, h$hi), plab(p)), legend_pos = c(0.03, 0.03), leg_top = TRUE, parts = TRUE) })
## KM curves and risk tables enter the main grid unwrapped, so their x-axes align with panels A-D and I-J
## I C-index forest
pf <- o$perf; pf$lab <- sprintf("%s (n = %d, %d ev.)", c("TCGA training", "TCGA hold-out test", "GSE39582", "GSE17536"), pf$n, pf$events)
fd <- data.frame(label = pf$lab, HR = pf$C, lo = pf$C_lo, hi = pf$C_hi, col = c("b", "a", "a", "a"), stringsAsFactors = FALSE)
fd$text <- fmt_ci(fd$HR, fd$lo, fd$hi, 3)
pI <- ggplot(fd, aes(HR, factor(label, rev(label)), colour = col)) + geom_vline(xintercept = 0.5, linetype = 2, linewidth = 0.3, colour = "grey40") + geom_segment(aes(x = lo, xend = hi, yend = factor(label, rev(label))), linewidth = 0.45) + geom_point(shape = 15, size = 1.4) + geom_text(aes(x = 1.07, label = text), hjust = 1, size = TXT * 0.95, colour = "black") +
  scale_colour_manual(values = c(a = "grey25", b = cl_sel, c = "#3B6FB6"), guide = "none") + scale_x_continuous(limits = c(0.35, 1.08), breaks = seq(0.4, 0.8, 0.1)) + labs(x = "Harrell's C-index (95% CI)", y = NULL, title = "Discrimination") + theme(axis.text.y = element_text(size = 5.4))
## J permutation null
pm <- data.frame(c = o$perm); obs <- pf$C[pf$cohort == "TCGA_test"]; pp <- (1 + sum(pm$c >= obs)) / (1 + nrow(pm)); cat(sprintf("perm p = %.3f, empty = %.2f\n", pp, mean(pm$c == 0.5)))
pJ <- ggplot(pm, aes(c)) + geom_histogram(binwidth = 0.01, fill = "grey70", colour = "white", linewidth = 0.2) + geom_vline(xintercept = obs, colour = cl_sel, linewidth = 0.6) + annotate("text", x = obs + 0.012, y = Inf, vjust = 1.6, hjust = 0, size = TXT, label = sprintf("Observed C = %.3f\nPermutation P = %.2f\n(%.0f%% of null fits empty)", obs, pp, 100 * mean(pm$c == 0.5)), lineheight = 0.95) +
  labs(x = "Hold-out C-index under permuted survival labels", y = "Permutations (n = 100)", title = "Permutation null (test set)")
## K adjusted HR and incremental C
cinc <- c(TCGA_test = rmv[["lasso_Cindex_clin_TCGA_test"]], GSE39582 = rmv[["lasso_Cindex_clin_GSE39582"]], GSE17536 = rmv[["lasso_Cindex_clin_GSE17536"]])
ha <- do.call(rbind, lapply(names(o$hr_adj), function(n) cbind(o$hr_adj[[n]], cohort = n))); rmn <- c(TCGA_test = rmv[["lasso_adjHR_TCGA_test"]], GSE39582 = rmv[["lasso_adjHR_GSE39582"]], GSE17536 = rmv[["lasso_adjHR_GSE17536"]])
nn <- as.integer(sub(".*n=", "", rmn)); ha$label <- sprintf("%s (n = %d)", c("TCGA hold-out test", "GSE39582", "GSE17536"), nn)
ha$text <- sprintf("%s; P = %.3f;  C %s", fmt_ci(ha$HR, ha$lo, ha$hi), ha$p, gsub(" -> ", "\u2192", cinc[ha$cohort])); ha$col <- "a"
pK <- forest_panel(ha[, c("label", "HR", "lo", "hi", "text", "col")], xlab = "Adjusted HR per SD of signature score (age, sex, stage)", xlim = c(0.5, 2), breaks = c(0.5, 0.75, 1, 1.5, 2), txt_head = "HR (95% CI); P; C-index clinical\u2192+signature", wid = c(1, 1.1), title = "Incremental value beyond clinical covariates")
fig <- wrap_plots(pA, pB, pC, pD, kmp[[1]]$km, kmp[[2]]$km, kmp[[1]]$risk, kmp[[2]]$risk, kmp[[3]]$km, kmp[[4]]$km, kmp[[3]]$risk, kmp[[4]]$risk, pI, pJ, wrap_elements(full = pK),
  design = "AB\nCD\nEF\nGH\nIJ\nKL\nMN\nOO", heights = c(3, 3, 2.5, 0.6, 2.5, 0.6, 2.2, 3)) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E", "F", "", "", "G", "H", "", "", "I", "J", "K")))
save_fig(fig, "Figure5_GLP2R_signature", w = 210, h = 240)
cat("Figure 5 saved\n"); sink()
