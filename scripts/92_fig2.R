## 92_fig2.R : Figure 2 - prognostic association of GLP2R (TCGA OS), sensitivity analyses, subgroups, external cohorts
source(file.path(Sys.getenv("GLP2R_SCRIPTS", "W:/Re_build/05_Analysis/scripts"), "90_fig_helpers.R"))
suppressPackageStartupMessages(library(splines))
sink(file.path(LOGD, "92_fig2.txt"), split = TRUE)
so <- readRDS(file.path(RES, "survival_objects.rds")); dd <- so$dd; cc5 <- so$cc5
rs <- read.csv(file.path(RES, "rm_survival.csv")); R <- setNames(rs$value, rs$key)
dd$grp <- factor(ifelse(dd$glp2r > median(dd$glp2r), "High", "Low"), levels = c("Low", "High"))
lr <- survdiff(Surv(time, event) ~ grp, dd); p_lr <- 1 - pchisq(lr$chisq, 1); h0 <- cox_tab(coxph(Surv(time, event) ~ grp, dd))[1, ]
cat("KM:", p_lr, h0$HR, h0$lo, h0$hi, "\n")
pA <- km_panel(dd, "grp", col_grp, xmax = 120, step = 24, title = "GLP2R median split (TCGA-COADREAD, OS)",
  annot = sprintf("Log-rank %s\nHR (high vs low) = %s\nExploratory; maxstat-adjusted P = %.2f", plab(p_lr), fmt_ci(h0$HR, h0$lo, h0$hi), as.numeric(R["maxstat_corrected_p"])), legend_pos = c(0.03, 0.03), lab_map = NULL)
## B: models
cm <- read.csv(file.path(RES, "cox_glp2r_models.csv")); mi <- regmatches(R["mi_glp2r_hr"], regexec("^([0-9.]+) \\(([0-9.]+)-([0-9.]+)\\); p=([0-9.]+)", R["mi_glp2r_hr"]))[[1]]
fb <- data.frame(label = c("M0 unadjusted", "M0, same sample as M2", "M1 + age, sex, stage, site", "M2 + MSI (primary)", "M2, multiple imputation"),
  n = c(cm$n[1], cm$n[2], cm$n[3], cm$n[4], nrow(dd)), ev = c(cm$events[1:4], sum(dd$event)), HR = c(cm$HR[1:4], as.numeric(mi[2])), lo = c(cm$lo[1:4], as.numeric(mi[3])), hi = c(cm$hi[1:4], as.numeric(mi[4])), p = c(cm$p[1:4], as.numeric(mi[5])))
fb$label <- paste0(fb$label, "\n(n = ", fb$n, ", ", fb$ev, " events)"); fb$text <- paste0(fmt_ci(fb$HR, fb$lo, fb$hi), "; P = ", ifelse(fb$p < 0.001, "<0.001", sprintf("%.3f", fb$p))); fb$col <- ifelse(grepl("primary", fb$label), "b", "a")
pB <- wrap_elements(full = forest_panel(fb, xlab = "Hazard ratio per 1-SD increase in GLP2R (95% CI)", xlim = c(0.5, 1.4), breaks = c(0.5, 0.75, 1, 1.25), title = "GLP2R (per SD) and overall survival: effect of adjustment", wid = c(1, 0.42)))
## C: spline
pm <- coxph(Surv(time, event) ~ ns(glp2r_z, df = 3) + age10 + sex + stage_grp + site + msi, cc5); lin <- coxph(Surv(time, event) ~ glp2r_z + age10 + sex + stage_grp + site + msi, cc5); pnl <- anova(lin, pm)[2, "Pr(>|Chi|)"]
tm <- predict(pm, type = "terms", se.fit = TRUE, terms = 1); sp <- data.frame(z = cc5$glp2r_z, f = tm$fit[, 1], se = tm$se.fit[, 1]); sp <- sp[order(sp$z), ]
pC <- ggplot(sp, aes(z, exp(f))) + geom_ribbon(aes(ymin = exp(f - 1.96 * se), ymax = exp(f + 1.96 * se)), fill = "#9ECAE1", alpha = 0.5) + geom_line(colour = "#08519C", linewidth = 0.6) + geom_hline(yintercept = 1, linetype = 2, linewidth = 0.3) +
  geom_rug(aes(x = z), inherit.aes = FALSE, sides = "b", alpha = 0.25, linewidth = 0.2) + scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4)) + coord_cartesian(ylim = c(0.2, 5)) +
  labs(x = "GLP2R expression (z-score of VST)", y = "Adjusted HR (M2; vs sample mean)", title = "Dose-response (natural spline, 3 df)") + annotate("text", x = min(sp$z), y = 4.6, label = sprintf("Non-linearity: %s\nn = %d, %d events", plab(pnl), nrow(cc5), sum(cc5$event)), hjust = 0, vjust = 1, size = TXT, lineheight = 0.95)
## D: time-dependent AUC + incremental value
ta <- read.csv(file.path(RES, "glp2r_timeROC.csv"))
pD <- ggplot(ta, aes(factor(time_m), AUC)) + geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.3) + geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12, linewidth = 0.4, colour = "#08519C") + geom_point(size = 1.6, colour = "#08519C") +
  geom_text(aes(label = sprintf("%.2f", AUC)), nudge_x = 0.22, size = TXT) + scale_y_continuous(limits = c(0.3, 0.85), breaks = seq(0.3, 0.8, 0.1)) + scale_x_discrete(labels = c("12", "36", "60")) +
  labs(x = "Time (months)", y = "Time-dependent AUC (IPCW)", title = "Discrimination of GLP2R alone") + annotate("text", x = 0.55, y = 0.84, hjust = 0, vjust = 1, size = TXT, lineheight = 0.95,
  label = sprintf("Clinical model C-index %.3f\n+ GLP2R %.3f (\u0394 = %.4f)\nLRT %s", as.numeric(R["cindex_base"]), as.numeric(R["cindex_with_glp2r"]), 0.0045, plab(as.numeric(R["lrt_glp2r_added_p"]))))
## E: subgroup forest
sg <- read.csv(file.path(RES, "cox_subgroups.csv")); rows <- list()
for (g in unique(sg$Subgroup)) { x <- sg[sg$Subgroup == g, ]; rows[[length(rows) + 1]] <- data.frame(label = sprintf("%s  (P-int = %s; FDR = %.2f)", toupper(g), ifelse(x$p_interaction[1] < 0.001, "<0.001", sprintf("%.3f", x$p_interaction[1])), x$p_int_FDR[1]), HR = NA, lo = NA, hi = NA, text = "", hdr = TRUE, col = "a")
  rows[[length(rows) + 1]] <- data.frame(label = sprintf("   %s (n=%d, %d ev.)", x$Level, x$n, x$events), HR = x$HR, lo = x$lo, hi = x$hi, text = sprintf("%s; P = %.3f", fmt_ci(x$HR, x$lo, x$hi), x$p), hdr = FALSE, col = ifelse(g == "Stage" & x$Level == "I-II", "b", "a")) }
fe <- do.call(rbind, rows); fe$shp <- 15
pE <- wrap_elements(full = forest_panel(fe, xlim = c(0.25, 3.2), breaks = c(0.25, 0.5, 1, 2), title = "Pre-specified subgroups (M2-adjusted; exploratory)", wid = c(1, 0.5)))
## F: external cohorts
go <- readRDS(file.path(RES, "geo_objects.rds")); est <- go$est; ff <- data.frame(label = c("TCGA-COADREAD\n(RNA-seq)", "GSE39582*\n(Affymetrix)", "GSE17536*\n(Affymetrix)", "Random-effects\nsummary"), n = est$n, ev = est$events, HR = est$HR, lo = est$lo, hi = est$hi)
ff$label <- paste0(ff$label, "\nn = ", ff$n, ", ", ff$ev, " events"); ff$text <- fmt_ci(ff$HR, ff$lo, ff$hi); ff$shp <- c(15, 15, 15, 18); ff$col <- c("a", "a", "a", "b")
ff$text[4] <- sprintf("%s\nI\u00b2 = %.0f%%; P = %.3f", ff$text[4], go$I2, go$p_re)
pF <- forest_panel(ff, xlim = c(0.4, 1.6), breaks = c(0.5, 0.75, 1, 1.5), xlab = "Hazard ratio per SD (95% CI)\n*Probe 221312_at is near array background (see legend)", title = "Adjusted HR per SD in independent cohorts", wid = c(1, 0.42), txt_head = "HR (95% CI)")
pF <- wrap_elements(full = pF)
pA <- wrap_elements(full = pA)
fig <- wrap_plots(pA, pC, pB, pD, pF, pE, design = "AB\nCC\nDE\nFF") + plot_layout(heights = c(1.1, 0.9, 0.95, 1.2)) + plot_annotation(tag_levels = "A")
save_fig(fig, "Figure2_GLP2R_survival", w = 180, h = 225)
cat("Figure 2 saved; spline nonlinearity P =", pnl, "\n"); sink()
