## 90_fig_helpers.R : shared plotting helpers for all rebuilt figures (sourced after 00_common.R)
source(file.path(Sys.getenv("GLP2R_SCRIPTS", "W:/Re_build/05_Analysis/scripts"), "00_common.R"))
suppressPackageStartupMessages({ library(patchwork); library(ggrepel); library(cowplot) })
TXT <- 2.15   # annotation text size (mm), approx. 6 pt
theme_update(plot.title = element_text(size = 6.8, face = "bold", hjust = 0, margin = margin(b = 2)), plot.tag = element_text(size = 9, face = "bold"))
col_dir <- c(Up = "#C8472B", Down = "#3B6FB6", NS = "grey75")
col_sig <- c(`FDR < 0.05` = "#C8472B", `n.s.` = "grey60")
sup <- function(x) { m <- c("0" = "\u2070", "1" = "\u00b9", "2" = "\u00b2", "3" = "\u00b3", "4" = "\u2074", "5" = "\u2075", "6" = "\u2076", "7" = "\u2077", "8" = "\u2078", "9" = "\u2079", "-" = "\u207b"); paste(m[strsplit(as.character(x), "")[[1]]], collapse = "") }
ptxt <- function(p) vapply(p, function(z) if (is.na(z)) "NA" else if (z < 0.001) { e <- floor(log10(z)); sprintf("%.1f\u00d710%s", z / 10^e, sup(e)) } else sprintf("%.3f", z), "")
plab <- function(p, pre = "P") { s <- vapply(p, function(z) if (is.na(z)) "NA" else if (z < 0.001) { e <- floor(log10(z)); sprintf("%.1f\u00d710%s", z / 10^e, sup(e)) } else sprintf("%.3f", z), ""); paste0(pre, " = ", s) }
panel_tag <- function(p, tag) p + labs(tag = tag)

## ---- Kaplan-Meier curve with numbers-at-risk (returns a patchwork)
km_panel <- function(df, grp = "grp", cols = col_grp, xmax = 120, step = 24, title = NULL, annot = NULL, ylab = "Overall survival", xlab = "Months", ymin = 0, legend_pos = c(0.03, 0.03), lab_map = NULL, legend_h = FALSE, leg_br = FALSE, leg_top = FALSE, parts = FALSE) {
  df$.g <- factor(df[[grp]]); lv <- levels(df$.g); fit <- survfit(Surv(time, event) ~ .g, data = df)
  s <- summary(fit, censored = TRUE); st <- sub(".g=", "", as.character(s$strata)); cd <- data.frame(time = s$time, surv = s$surv, g = factor(st, levels = lv), nc = s$n.censor)
  cd <- rbind(data.frame(time = 0, surv = 1, g = factor(lv, levels = lv), nc = 0), cd); cd <- cd[cd$time <= xmax, ]
  nl <- table(df$.g); labs_leg <- if (is.null(lab_map)) setNames(paste0(lv, " (n=", nl[lv], ")"), lv) else lab_map
  p1 <- ggplot(cd, aes(time, surv, colour = g)) + geom_step(linewidth = 0.5) + geom_point(data = cd[cd$nc > 0, ], shape = 3, size = 0.7, stroke = 0.25, show.legend = FALSE) +
    scale_colour_manual(values = cols, labels = labs_leg[lv], name = NULL) + scale_x_continuous(limits = c(0, xmax), breaks = seq(0, xmax, step), expand = expansion(mult = c(0.05, 0.03))) + scale_y_continuous(limits = c(ymin, 1), breaks = seq(0, 1, 0.25), expand = c(0.01, 0)) +
    labs(x = NULL, y = ylab, title = title) + theme(legend.position = "inside", legend.position.inside = if (leg_br) c(0.99, 0.03) else c(0.98, 0.98), legend.justification = if (leg_br) c(1, 0) else c(1, 1), legend.key.height = unit(2.5, "mm"), legend.direction = if (legend_h) "horizontal" else "vertical", axis.text.x = element_blank())
  if (leg_top) p1 <- p1 + theme(legend.position = "top", legend.justification = c(1, 0), legend.direction = "horizontal", legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(2, 0, -2, 0))
  if (!is.null(annot)) p1 <- p1 + annotate("text", x = xmax * legend_pos[1], y = ymin + (1 - ymin) * legend_pos[2], label = annot, hjust = 0, vjust = 0, size = TXT, lineheight = 0.95)
  brk <- seq(0, xmax, step); ss <- summary(fit, times = brk, extend = TRUE); rk <- data.frame(time = ss$time, n = ss$n.risk, g = factor(sub(".g=", "", as.character(ss$strata)), levels = lv))
  p2 <- ggplot(rk, aes(time, g, label = n, colour = g)) + geom_text(size = TXT) + scale_colour_manual(values = cols, guide = "none") + scale_x_continuous(limits = c(0, xmax), breaks = brk, expand = expansion(mult = c(0.05, 0.03))) + scale_y_discrete(limits = rev(lv)) +
    labs(x = xlab, y = NULL) + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(size = 5.5), plot.margin = margin(0, 4, 2, 3))
  if (parts) return(list(km = p1 + theme(plot.margin = margin(5.5, 5.5, 0, 5.5)), risk = p2))
  p1 / p2 + plot_layout(heights = c(4.2, 1))
}
## ---- generic box + jitter panel
box_panel <- function(df, x, y, cols, ylab, title = NULL, annot = NULL, n_in_label = TRUE, ann_y = NULL, jit = 0.18, alpha_pt = 0.5, sz = 0.45) {
  df <- df[!is.na(df[[x]]) & !is.na(df[[y]]), ]; df$.x <- factor(df[[x]]); nl <- table(df$.x)
  if (n_in_label) levels(df$.x) <- paste0(levels(df$.x), "\n(n=", nl[levels(df$.x)], ")"); cn <- setNames(cols[names(nl)], levels(df$.x))
  p <- ggplot(df, aes(.x, .data[[y]], fill = .x, colour = .x)) + geom_jitter(width = jit, size = sz, alpha = alpha_pt, stroke = 0) + geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.3, alpha = 0.35, colour = "black", fatten = 1.4) +
    scale_fill_manual(values = cn) + scale_colour_manual(values = cn) + labs(x = NULL, y = ylab, title = title) + theme(legend.position = "none")
  if (!is.null(annot)) { yr <- range(df[[y]]); nl_ann <- length(strsplit(annot, "\n")[[1]]); hd <- 0.11 * nl_ann + 0.04
    p <- p + scale_y_continuous(expand = expansion(mult = c(0.03, hd))) + annotate("text", x = (length(nl) + 1) / 2, y = if (is.null(ann_y)) yr[2] + 0.03 * diff(yr) else ann_y, label = annot, size = TXT, vjust = 0, lineheight = 0.95) }
  p
}
## ---- forest plot (HR with CI) + text column; df needs label, HR, lo, hi, text (+ optional col, shp, hdr)
forest_panel <- function(df, xlab = "Hazard ratio (95% CI)", xlim = c(0.3, 3), ref = 1, breaks = c(0.25, 0.5, 1, 2), txt_head = "HR (95% CI); P", wid = c(1, 0.62), title = NULL) {
  n <- nrow(df); df$y <- rev(seq_len(n)); if (is.null(df$col)) df$col <- "a"; if (is.null(df$shp)) df$shp <- 15; if (is.null(df$hdr)) df$hdr <- FALSE
  d1 <- df[!df$hdr & !is.na(df$HR), ]; cl <- function(x) pmin(pmax(x, xlim[1]), xlim[2])
  p <- ggplot(d1, aes(HR, y, colour = col)) + geom_vline(xintercept = ref, linetype = 2, linewidth = 0.3, colour = "grey40") +
    geom_segment(aes(x = cl(lo), xend = cl(hi), y = y, yend = y), linewidth = 0.4) + geom_point(aes(shape = factor(shp), size = factor(shp))) + scale_shape_manual(values = c(`15` = 15, `18` = 18)) + scale_size_manual(values = c(`15` = 1.3, `18` = 2.6)) +
    scale_x_log10(limits = xlim, breaks = breaks, labels = breaks) + scale_y_continuous(breaks = df$y, labels = df$label, limits = c(0.4, n + 1.0), expand = c(0, 0)) + scale_colour_manual(values = c(a = "grey20", b = "#C8472B", c = "#3B6FB6")) +
    labs(x = xlab, y = NULL, title = title) + theme(axis.text.y = element_text(size = 5.6, hjust = 1), legend.position = "none")
  pt <- ggplot(df, aes(0, y, label = text)) + geom_text(hjust = 0, size = TXT * 0.95, lineheight = 0.9) + annotate("text", x = 0, y = n + 0.7, label = txt_head, hjust = 0, size = TXT * 0.95, fontface = "bold") +
    scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0.4, n + 1.0), expand = c(0, 0)) + theme_void()
  p + pt + plot_layout(widths = wid)
}
fmt_ci <- function(h, l, u, d = 2) sprintf(paste0("%.", d, "f (%.", d, "f\u2013%.", d, "f)"), h, l, u)
