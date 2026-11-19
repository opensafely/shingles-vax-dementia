# Figures from the released descriptive tables and model estimates. As a disclosure-control measure, scale the
# fitted birth-date densities by the rounded cohort count. Write placeholder figures when ggplot2 is not
# installed or there is nothing to draw, because the OpenSAFELY job runner fails an action whose declared
# outputs are missing.
#
# Author/contact: Samantha Ip (hyi20@cam.ac.uk)

PAL <- c(blue = "#2a78d6", orange = "#eb6834", ink = "#0b0b0b", ink2 = "#52514e", grid = "#e6e5e1")

save_plot <- function(p, file, width = 7.5, height = 4.5) {
  if (is.null(p) || !requireNamespace("ggplot2", quietly = TRUE)) {
    why <- if (requireNamespace("ggplot2", quietly = TRUE)) "No data to plot" else
      "Not drawn: the ggplot2 package is not installed"
    grDevices::png(file, width = width, height = height, units = "in", res = 150)
    graphics::plot.new()
    graphics::text(0.5, 0.5, why)
    grDevices::dev.off()
    log_line("  wrote ", basename(file), " (", why, ")")
  } else {
    ggplot2::ggsave(file, p, width = width, height = height, dpi = 150, bg = "white")
    log_line("  wrote ", basename(file))
  }
}

theme_rd <- function() {
  ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = PAL[["grid"]], linewidth = 0.3),
                   axis.text = ggplot2::element_text(colour = PAL[["ink2"]]),
                   axis.title = ggplot2::element_text(colour = PAL[["ink"]]),
                   plot.title = ggplot2::element_text(colour = PAL[["ink"]], face = "bold"),
                   plot.caption = ggplot2::element_text(colour = PAL[["ink2"]], hjust = 0),
                   plot.title.position = "plot", plot.caption.position = "plot",
                   legend.position = "top", legend.title = ggplot2::element_blank())
}

side_labels <- c("0" = "Born before 1 Sep 1958 (ineligible)", "1" = "Born on or after 1 Sep 1958 (eligible)")
side_scale <- function() ggplot2::scale_colour_manual(values = c(PAL[["blue"]], PAL[["orange"]]),
                                                      labels = side_labels, breaks = c("0", "1"))
cutoff_line <- function(cutoff) ggplot2::geom_vline(xintercept = cutoff, linetype = "dashed",
                                                    colour = PAL[["ink2"]], linewidth = 0.4)

# Crude birth-month incidence rates and model-based average rates over planned follow-up, from events_by_bin().
plot_rates <- function(ev, cutoff, title) {
  if (!nrow(ev) || !requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  ev$side <- as.character(ev$eligible)
  pts <- ev[!is.na(ev$rate_per_1000py), ]
  fit <- ev[!is.na(ev$model_rate_per_1000py), ]
  if (!nrow(pts) && !nrow(fit)) return(NULL)
  ggplot2::ggplot(mapping = ggplot2::aes(bin_middle, colour = side)) +
    cutoff_line(cutoff) +
    ggplot2::geom_point(data = pts, ggplot2::aes(y = rate_per_1000py), size = 1.6) +
    ggplot2::geom_line(data = fit, ggplot2::aes(y = model_rate_per_1000py, group = side), linewidth = 0.7) +
    side_scale() + theme_rd() +
    ggplot2::labs(x = "Birth month", y = "Events per 1000 person-years", title = title,
                  caption = paste("Points: crude rate in each birth month (counts rounded to 5 and withheld at 7 or fewer).",
                                  "Lines: the fitted Weibull RD model's rate for an adult followed from day 31 to the end of follow-up.", sep = "\n"))
}

# Birth-date density by recorded month or week, with the McCrary test's separate local linear fits on either
# side of the eligibility cutoff.
plot_density <- function(bins, mc, cutoff, title) {
  bins <- bins[!is.na(bins$n), ]
  if (!nrow(bins) || is.null(mc) || !requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  bins$side <- as.character(as.integer(bins$x >= 0))
  bins$date <- cutoff + bins$x * 365.25
  bins$per_day <- bins$n / bins$width_days
  h <- mc$bandwidth_days / 365.25
  line <- function(s, a, b) {
    x <- if (s == "1") seq(0, min(h, max(bins$x)), length.out = 50) else seq(-min(h, -min(bins$x)), 0, length.out = 50)
    data.frame(side = s, date = cutoff + x * 365.25, per_day = (a + b * x) * mc$n / 365.25)
  }
  fits <- rbind(line("0", mc$density_left, mc$slope_left), line("1", mc$density_right, mc$slope_right))
  ggplot2::ggplot(bins, ggplot2::aes(date, per_day, colour = side)) +
    cutoff_line(cutoff) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::geom_line(data = fits, linewidth = 0.7) +
    side_scale() + theme_rd() +
    ggplot2::labs(x = "Recorded birth month or week", y = "Adults per day of birth date", title = title,
                  caption = sprintf("McCrary test: log density ratio at the cutoff %.3f (SE %.3f), p = %.2f, bandwidth %.0f days. Counts rounded to 5.",
                                    mc$theta, mc$se, mc$p, mc$bandwidth_days))
}

# RZV uptake by birth month, from uptake_by_bin().
plot_uptake <- function(up, cutoff, title) {
  up <- up[!is.na(up$pct_vaccinated), ]
  if (!nrow(up) || !requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  up$side <- as.character(up$eligible)
  ggplot2::ggplot(up, ggplot2::aes(bin_middle, pct_vaccinated, colour = side)) +
    cutoff_line(cutoff) +
    ggplot2::geom_point(size = 1.6) +
    side_scale() + theme_rd() +
    ggplot2::labs(x = "Birth month", y = "Vaccinated (%)", title = title,
                  caption = paste("Adults with at least one RZV dose by the end of their follow-up.",
                                  "Counts rounded to 5 and withheld at 7 or fewer."))
}

# Cohort selection flowchart from one nation's selection_table() (cohort_selection.csv): the adults remaining
# after each SAP 5.2.1 criterion, in the SAP's order, the adults excluded at each criterion beside it, and the
# analysis cohort by eligibility at the foot. It uses the released counts, rounded to 5, and shows a withheld
# count as "7 or fewer". The table can also be read back from the released file and drawn outside the SDE.
plot_flowchart <- function(sel, title) {
  all <- sel[sel$group == "All", ]
  all <- all[order(all$step), ]
  if (nrow(all) < 2 || !requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  num <- function(x) ifelse(is.na(x), "7 or fewer", ifelse(x >= 10000, formatC(x, format = "d", big.mark = ","), as.character(x)))
  last <- function(g) { s <- sel[sel$group == g, ]; s$remaining[which.max(s$step)] }
  k <- nrow(all) - 1L                                            # criteria, after "Extracted"
  main <- data.frame(x = 0, y = -(0:k), label = paste0(c("Extracted", rep("Remaining", k - 1L), "Analysis cohort"), ": ", num(all$remaining)))
  excl <- data.frame(x = 0.62, y = -(seq_len(k) - 0.5), label = paste0(all$criterion[-1], "\n", num(all$excluded[-1]), " excluded"))
  foot <- data.frame(x = c(-0.55, 0.55), y = -(k + 1.15),
                     label = c(paste0("Eligible: ", num(last("Eligible"))), paste0("Ineligible: ", num(last("Ineligible")))))
  arrow <- grid::arrow(length = grid::unit(0.09, "inches"), type = "closed")
  down <- data.frame(x = 0, xend = 0, y = -(0:(k - 1)) - 0.22, yend = -(1:k) + 0.22)
  side <- data.frame(x = 0, xend = 0.6, y = excl$y, yend = excl$y)
  split <- data.frame(x = 0, xend = foot$x, y = -k - 0.22, yend = foot$y + 0.22)
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = rbind(down, split), ggplot2::aes(x, y, xend = xend, yend = yend),
                          colour = PAL[["ink2"]], linewidth = 0.4, arrow = arrow) +
    ggplot2::geom_segment(data = side, ggplot2::aes(x, y, xend = xend, yend = yend), colour = PAL[["ink2"]],
                          linewidth = 0.4, arrow = arrow) +
    ggplot2::geom_label(data = main, ggplot2::aes(x, y, label = label), size = 3.1, fill = "white", colour = PAL[["ink"]],
                        label.padding = grid::unit(0.35, "lines")) +
    ggplot2::geom_label(data = excl, ggplot2::aes(x, y, label = label), size = 2.8, hjust = 0, fill = "white",
                        colour = PAL[["ink2"]], lineheight = 0.95, label.padding = grid::unit(0.3, "lines")) +
    ggplot2::geom_label(data = foot, ggplot2::aes(x, y, label = label), size = 3.1, fill = "white",
                        colour = c(PAL[["orange"]], PAL[["blue"]]), label.padding = grid::unit(0.35, "lines")) +
    ggplot2::scale_x_continuous(limits = c(-1, 2.3)) +
    ggplot2::scale_y_continuous(limits = c(-(k + 1.6), 0.5)) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(colour = PAL[["ink"]], face = "bold"),
                   plot.caption = ggplot2::element_text(colour = PAL[["ink2"]], hjust = 0),
                   plot.title.position = "plot", plot.caption.position = "plot",
                   plot.margin = ggplot2::margin(8, 8, 8, 8)) +
    ggplot2::labs(title = title, caption = "Criteria of SAP 5.2.1, in its order. Counts rounded to 5 and withheld at 7 or fewer.")
}

# Forest plot with pooled estimates shown as diamonds. Preserve input row order within each panel, with panels
# in order of first appearance. Requires label, panel, est, lci, uci and pooled.
forest_plot <- function(df, xlab, title, log_x = TRUE, caption = NULL) {
  df <- df[is.finite(df$est), ]
  if (!nrow(df) || !requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  df$row <- factor(seq_len(nrow(df)), levels = rev(seq_len(nrow(df))))
  df$panel <- factor(df$panel, levels = unique(df$panel))
  df$text <- sprintf("%.2f (%.2f to %.2f)", df$est, df$lci, df$uci)
  df$kind <- ifelse(df$pooled, "Pooled", "Nation")
  lo <- min(df$lci, 1, na.rm = TRUE)
  hi <- max(df$uci, 1, na.rm = TRUE)
  x_scale <- if (log_x) ggplot2::scale_x_log10(limits = c(lo / 1.1, hi * 2.6)) else
    ggplot2::scale_x_continuous(limits = c(min(lo, 0), hi + 0.9 * (hi - min(lo, 0))))
  text_x <- if (log_x) hi * 1.25 else hi + 0.1 * (hi - min(lo, 0))
  ggplot2::ggplot(df, ggplot2::aes(y = row)) +
    ggplot2::geom_vline(xintercept = 1, colour = PAL[["ink2"]], linewidth = 0.4, linetype = "dashed") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lci, xmax = uci, colour = kind), orientation = "y", width = 0, linewidth = 0.7) +
    ggplot2::geom_point(ggplot2::aes(x = est, colour = kind, shape = kind), size = 2.6) +
    ggplot2::geom_text(ggplot2::aes(x = text_x, label = text), hjust = 0, size = 3, colour = PAL[["ink"]]) +
    ggplot2::scale_y_discrete(labels = stats::setNames(df$label, as.character(df$row))) +
    x_scale +
    ggplot2::scale_colour_manual(values = c(Nation = PAL[["blue"]], Pooled = PAL[["orange"]])) +
    ggplot2::scale_shape_manual(values = c(Nation = 16, Pooled = 18)) +
    ggplot2::facet_grid(panel ~ ., scales = "free_y", space = "free_y", switch = "y") +
    theme_rd() +
    ggplot2::theme(strip.placement = "outside", strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold"),
                   panel.grid.major.y = ggplot2::element_blank()) +
    ggplot2::labs(x = xlab, y = NULL, title = title, caption = caption)
}
