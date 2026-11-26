# McCrary (2008) density discontinuity test (SAP 6.4.1). Fit separate triangular-kernel local linear
# regressions on each side of the cutoff to birth-date densities estimated from recorded monthly (England) or
# weekly (Wales) bins. Report the log density ratio at the cutoff (theta, eligible over ineligible) and its
# asymptotic SE. Cap the rule-of-thumb bandwidth at max_h, which run_nation.R sets to the 2-year primary
# bandwidth, beyond which the cohort has no birth dates.
#
# Author/contact: Samantha Ip (hyi20@cam.ac.uk)

density_bins <- function(co, cutoff) {
  if (anyNA(co$bin_start)) stop("density_bins() needs the recorded birth month or week of every adult")
  agg <- rowsum(rep(1, nrow(co)), as.character(co$bin_start))
  start <- as.Date(rownames(agg))
  width <- co$bin_width[match(start, co$bin_start)]
  # x is the midpoint of the bin's first and last days, in years from the cutoff.
  data.frame(bin_start = start, width_days = width, n = agg[, 1],
             x = (as.numeric(start - cutoff) + (width - 1) / 2) / 365.25)
}

mccrary_bandwidth <- function(x, y, right) {
  side <- function(s) {
    xs <- x[right == s]
    ys <- y[right == s]
    if (length(xs) < 6) return(NA_real_)
    fit <- stats::lm(ys ~ poly(xs, 4, raw = TRUE))      # quartic pilot fit on each side
    b <- stats::coef(fit)
    f2 <- 2 * b[3] + 6 * b[4] * xs + 12 * b[5] * xs^2
    sigma2 <- sum(stats::resid(fit)^2) / (length(xs) - 5)
    h <- 3.348 * (sigma2 * max(abs(xs)) / sum(f2^2))^(1 / 5)
    if (is.finite(h) && h > 0) h else NA_real_          # NA if the fit has no curvature or no residual variance
  }
  mean(c(side(FALSE), side(TRUE)))
}

mccrary_test <- function(bins, max_h) {
  stopifnot(is.finite(max_h), max_h > 0)
  n <- sum(bins$n)
  x <- bins$x
  y <- bins$n / (n * bins$width_days / 365.25)          # proportion per year of birth date
  right <- x >= 0
  h <- mccrary_bandwidth(x, y, right)
  # Use max_h if the rule-of-thumb bandwidth exceeds it or either side gives none.
  if (!is.finite(h) || h <= 0 || h > max_h) h <- max_h
  side <- function(s) {
    w <- pmax(0, 1 - abs(x / h)) * (right == s)
    keep <- w > 0
    if (sum(keep) < 3) stop("fewer than 3 bins within the bandwidth on one side")
    stats::lm.wfit(cbind(1, x[keep]), y[keep], w[keep])$coefficients
  }
  bl <- side(FALSE)
  br <- side(TRUE)
  if (!isTRUE(bl[[1]] > 0 && br[[1]] > 0)) stop("a fitted density at the cutoff is not positive")
  theta <- log(br[[1]]) - log(bl[[1]])
  se <- sqrt((1 / (n * h)) * (24 / 5) * (1 / br[[1]] + 1 / bl[[1]]))
  data.frame(theta = theta, se = se, z = theta / se, p = 2 * stats::pnorm(-abs(theta / se)),
             bandwidth_days = h * 365.25, density_left = bl[[1]], density_right = br[[1]],
             slope_left = bl[[2]], slope_right = br[[2]], n = n,
             bins_left = sum(!right & abs(x) < h), bins_right = sum(right & abs(x) < h))
}
