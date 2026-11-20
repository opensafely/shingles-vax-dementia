# Approximate RZV receipt effects using fuzzy RD (SAP 6.5.1, secondary objective 3, and the sensitivity
# analyses of SAP 6.6). delta is the eligibility discontinuity in receipt, by the end of outcome-specific
# follow-up (derive.R), of at least one dose or, for SAP 6.6 (iii), two doses. DELTAS in config.R lists the
# estimates of delta.
#
# Scale the Weibull eligibility effects by 1 / delta: HR_CACE = HR_eligibility^(1 / delta) and
# TR_CACE = TR_eligibility^(1 / delta). inc_ineligible_hr_cace = 1 - (1 - F0)^HR_CACE applies HR_CACE to F0
# (inc_ineligible), the ineligible side's incidence at the cutoff. Reading the two as the compliers' risks with
# and without vaccination is a further assumption.
#
# A causal reading requires RD continuity, exclusion, monotonicity, a non-zero first stage and further
# hazard-model assumptions. These transformations do not exactly identify complier effects.
#
# delta's 95% CI is rdrobust's robust bias-corrected CI (SAP 6.5.1; Calonico, Cattaneo and Titiunik,
# Econometrica 2014), which allows for the smoothing bias of the local polynomial and for the uncertainty in
# estimating that bias. The other CIs are individual-level bootstrap percentile CIs (1000 samples in the SAP).
# The bootstrap SDs, of the logs for ratios, are the SEs for pooling across nations.
#
# Author/contact: Samantha Ip (hyi20@cam.ac.uk)

# Estimate delta by weighted least squares of receipt (column y) on 1, r, z and r*z, adding r^2 and r^2*z for
# p = 2, with triangular-kernel weights within the bandwidth h, in years. The coefficient of z is the
# discontinuity at r = 0. Frequency weights reproduce individual-level fitting, and the estimate is
# rdrobust's conventional estimate with the same p and h. Report the fitted uptake on each side at the cutoff.
# Predictions need not lie in [0, 1].
uptake_fit <- function(md, w, h, p = 1L, y = "vacc_by_end") {
  if (!is.numeric(h) || length(h) != 1 || !is.finite(h) || h <= 0) {
    stop("the bandwidth for delta is not a positive number")
  }
  w <- w * pmax(0, 1 - abs(md$r) / h)
  keep <- w > 0
  X <- cbind(1, md$r, md$z, md$rz)
  if (p == 2L) X <- cbind(X, md$r^2, md$r^2 * md$z)
  X <- X[keep, , drop = FALSE]
  if (qr(X)$rank < ncol(X)) stop("too few distinct birth dates on one side of the cutoff to estimate delta")
  b <- stats::lm.wfit(X, md[[y]][keep], w[keep])$coefficients
  c(uptake_ineligible = b[[1]], delta = b[[3]])
}
uptake_jump <- function(md, w, h, p = 1L, y = "vacc_by_end") uptake_fit(md, w, h, p, y)[["delta"]]

# rdrobust needs one row per adult. It warns of mass points whenever birth dates are binned, as they are by
# month or week here, and only that warning is silenced.
per_adult <- function(md) rep(seq_len(nrow(md)), md$n)
quiet_masspoints <- function(expr) withCallingHandlers(expr, warning = function(w) {
  if (grepl("mass points", conditionMessage(w), ignore.case = TRUE)) invokeRestart("muffleWarning")
})

# Choose the bandwidths, in years, for one of DELTAS: h for the estimate and b for the bias correction of its
# CI. "mserd" takes rdrobust's common MSE-optimal h and b for local linear regression with a triangular
# kernel, from the adults inside the analysis bandwidth, and caps both at that bandwidth, so that no kernel
# reaches beyond the analysis's adults. "analysis" uses the analysis bandwidth for both. Choose once on the
# full sample, and hold fixed in bootstrap refits (SAP 6.5.1).
delta_bandwidths <- function(md, analysis_bw, def) {
  if (def$bandwidth == "analysis") return(c(h = analysis_bw, b = analysis_bw))
  i <- per_adult(md)
  bw <- quiet_masspoints(rdrobust::rdbwselect(md[[def$receipt]][i], md$r[i], c = 0, p = def$p, kernel = "triangular",
                                              bwselect = "mserd"))$bws
  hb <- c(h = bw[1, 1], b = bw[1, 3])
  if (!all(is.finite(hb) & hb > 0)) stop("the bandwidth for delta is not a positive number")
  pmin(hb, analysis_bw)
}

# rdrobust's inference for one estimate of delta (a DELTAS entry with its bandwidths h and b): the
# conventional and robust SEs and the robust bias-corrected 95% CI, which uses a local polynomial one order
# higher, within b, to estimate the bias. The robust SE adds the variance of the estimated bias, so it is
# usually the larger. The bootstrap SD of delta, with h held fixed, estimates the conventional SE. Stop if
# rdrobust's conventional estimate is not the code's own.
delta_rbc <- function(md, spec) {
  i <- per_adult(md)
  f <- quiet_masspoints(rdrobust::rdrobust(md[[spec$receipt]][i], md$r[i], c = 0, p = spec$p, q = spec$p + 1L,
                                           h = spec$h, b = spec$b, kernel = "triangular"))
  if (abs(f$coef[1, 1] - uptake_jump(md, md$n, spec$h, spec$p, spec$receipt)) > 1e-8) {
    stop("rdrobust's estimate of delta differs from the code's own")
  }
  c(se_conventional = f$se[1, 1], se_robust = f$se[3, 1], lci = f$ci[3, 1], uci = f$ci[3, 2])
}

# The rows of results_receipt.csv for the estimates of delta named in defs, with missing values: one row per
# analysis and quantity, in the order receipt_estimates() gives them.
receipt_template <- function(defs, years = numeric(0)) {
  inc <- c("inc_ineligible_", "inc_ineligible_hr_cace_", "inc_difference_hr_cace_")
  one <- c("hr_eligibility", "delta", "uptake_ineligible", "uptake_eligible", "hr_cace", "tr_cace",
           as.vector(t(outer(inc, years, paste0))))
  q <- rep(one, length(defs))
  data.frame(analysis = rep(defs, each = length(one)), quantity = q,
             estimate = NA_real_, lci = NA_real_, uci = NA_real_,
             ci_method = ifelse(q == "delta", "robust bias-corrected", "bootstrap percentile"),
             boot_se = NA_real_, se_conventional = NA_real_, se_robust = NA_real_,
             se_scale = ifelse(q %in% c("hr_eligibility", "hr_cace", "tr_cace"), "log", "natural"),
             n_boot = NA_integer_, n_boot_failed = NA_integer_)
}

# Adults with and without receipt in each birth month that reaches into the bandwidth of an estimate of delta,
# counting the whole month. cutoff + r_days is the assumed birth date (README, convention 1), so for dementia's
# main delta these are the monthly counts that uptake_by_bin.csv releases rounded. In England a month is one
# birth bin, and in Wales several weeks.
month_counts <- function(md, spec) {
  m <- format(STUDY$cutoff + md$r_days, "%Y-%m")
  k <- m %in% m[abs(md$r) < spec$h]
  y <- md[[spec$receipt]][k]
  cbind(with = tapply(md$n[k] * y, m[k], sum), without = tapply(md$n[k] * (1 - y), m[k], sum))
}

# Calculate approximate receipt effects from a converged Weibull fit and each estimate of delta in specs
# (DELTAS entries with their bandwidths). Return missing complier effects if delta is zero or not finite, and
# missing values for all of an estimate of delta that cannot be made, without stopping the others. Derive
# incidence using the ineligible group's fitted cumulative hazard at the cutoff and the proportional hazards
# transformation. Name each estimate "analysis:quantity".
receipt_estimates <- function(md, w, f, specs, years = numeric(0), induction_days = 0) {
  if (!converged(f)) stop("the Weibull fit did not converge")
  b <- stats::coef(f)
  log_hr <- -b[["z"]] / f$scale
  unlist(lapply(names(specs), function(k) {
    s <- specs[[k]]
    u <- tryCatch(uptake_fit(md, w, s$h, s$p, s$receipt),
                  error = function(e) c(uptake_ineligible = NA_real_, delta = NA_real_))
    delta <- u[["delta"]]
    ok <- is.finite(delta) && delta != 0
    log_hr_cace <- if (ok) log_hr / delta else NA_real_
    v <- c(log_hr_eligibility = log_hr, delta = delta, uptake_ineligible = u[["uptake_ineligible"]],
           uptake_eligible = u[["uptake_ineligible"]] + delta, log_hr_cace = log_hr_cace,
           log_tr_cace = if (ok) b[["z"]] / delta else NA_real_)
    if (length(years)) {
      # ineligible-group cumulative hazard at the cutoff
      H0 <- exp((log(years - induction_days / 365.25) - b[["(Intercept)"]]) / f$scale)
      inc0 <- 1 - exp(-H0)
      inc1 <- 1 - exp(-H0 * exp(log_hr_cace))
      v <- c(v, stats::setNames(inc0, paste0("inc_ineligible_", years)),
             stats::setNames(inc1, paste0("inc_ineligible_hr_cace_", years)),
             stats::setNames(inc1 - inc0, paste0("inc_difference_hr_cace_", years)))
    }
    stats::setNames(v, paste0(k, ":", names(v)))
  }))
}

# Bootstrap individuals by multinomially resampling the collapsed-row counts, without stratification by
# eligibility. Refit the Weibull model and each delta on each resample, retaining their dependence. Initialise
# from the full-sample Weibull fit and hold the bandwidths for delta fixed.
#
# Each estimate of delta is judged on its own, so that one cannot remove the others. One that cannot be made on
# the full sample is left out, and the "delta_failures" attribute gives the reason. A bootstrap sample counts as
# failed for an analysis if any of its estimates is not finite there, and a failed Weibull refit fails the sample
# for every analysis. Record the reasons. Suppress an analysis's bootstrap CIs and SEs if its share of failed
# samples exceeds max_failed or fewer than two succeed. run_nation.R warns above 1%. The caller sets the
# original-sample event minimum. Bootstrap refits require at least one event per side.
#
# Return receipt_template()'s rows for the analyses made. Each delta row takes rdrobust's CI and SEs
# (delta_rbc()) and keeps the bootstrap SD in boot_se. If rdrobust fails for a delta, its CI is left missing and
# the "rbc_failures" attribute gives the reason. If a birth month reaching into the bandwidth for delta has fewer
# than min_count adults with receipt, or fewer without, the fitted uptakes at the cutoff could narrow down that
# withheld count, even when rounded, so both sides' are withheld and the "uptake_withheld" attribute names the
# analysis. delta is kept. min_count defaults to min_events_side, which run_nation.R sets from the nation's
# disclosure rule. round_receipt() rounds the results for release.
receipt_effect <- function(md, n_boot, seed, specs, min_events_side = 1L, years = numeric(0), induction_days = 0,
                           max_failed = STUDY$max_boot_failed, min_count = min_events_side) {
  stopifnot(n_boot >= 1, length(specs) >= 1)
  f <- fit_weibull_rd(md, min_events_side = min_events_side)
  why_not <- vapply(specs, function(s) tryCatch({ uptake_fit(md, md$n, s$h, s$p, s$receipt); "" },
                                                error = conditionMessage), "")
  if (all(nzchar(why_not))) stop(why_not[[1]])
  specs <- specs[!nzchar(why_not)]
  point <- receipt_estimates(md, md$n, f, specs, years, induction_days)
  rbc <- lapply(specs, function(s) tryCatch(delta_rbc(md, s), error = conditionMessage))
  init <- c(stats::coef(f), log(f$scale))
  n_total <- sum(md$n)
  p <- md$n / n_total
  set.seed(seed)
  boot <- matrix(NA_real_, n_boot, length(point), dimnames = list(NULL, names(point)))
  why <- rep("estimates not all finite", n_boot)
  for (b in seq_len(n_boot)) {
    w <- as.vector(stats::rmultinom(1, n_total, p))
    boot[b, ] <- tryCatch(receipt_estimates(md, w, fit_weibull_rd(md, w, init = init), specs, years, induction_days),
                          warning = function(e) { why[b] <<- paste("warning:", conditionMessage(e)); NA_real_ },
                          error = function(e) { why[b] <<- paste("error:", conditionMessage(e)); NA_real_ })
  }
  on_log <- grepl(":log_", names(point))                     # ratios are bootstrapped on the log scale
  back <- function(v, lg) ifelse(lg, exp(v), v)
  res <- receipt_template(names(specs), years)
  i <- match(sub(":log_", ":", names(point)), paste0(res$analysis, ":", res$quantity))
  res$estimate[i] <- back(unname(point), on_log)
  res$n_boot <- n_boot
  failed_any <- rep(FALSE, n_boot)
  withheld <- character(0)
  for (k in names(specs)) {
    cols <- startsWith(names(point), paste0(k, ":"))
    ok <- rowSums(!is.finite(boot[, cols, drop = FALSE])) == 0
    failed_any <- failed_any | !ok
    res$n_boot_failed[res$analysis == k] <- sum(!ok)
    if (mean(!ok) <= max_failed && sum(ok) >= 2L) {
      bk <- boot[ok, cols, drop = FALSE]
      q <- apply(bk, 2, stats::quantile, probs = c(0.025, 0.975), names = FALSE)
      res$lci[i[cols]] <- back(q[1, ], on_log[cols])         # percentile CIs also for HR_eligibility, unlike the
      res$uci[i[cols]] <- back(q[2, ], on_log[cols])         # Wald CIs of results_eligibility.csv
      res$boot_se[i[cols]] <- apply(bk, 2, stats::sd)
    }
    j <- res$analysis == k & res$quantity == "delta"
    res$lci[j] <- res$uci[j] <- NA_real_
    if (is.numeric(rbc[[k]])) res[j, names(rbc[[k]])] <- as.list(rbc[[k]])
    if (any(month_counts(md, specs[[k]]) < min_count)) {
      side <- res$analysis == k & res$quantity %in% c("uptake_ineligible", "uptake_eligible")
      res[side, c("estimate", "lci", "uci", "boot_se")] <- NA_real_
      withheld <- c(withheld, k)
    }
  }
  attr(res, "failures") <- table(why[failed_any])
  attr(res, "delta_failures") <- why_not[nzchar(why_not)]
  attr(res, "rbc_failures") <- unlist(Filter(is.character, rbc))
  attr(res, "uptake_withheld") <- withheld
  res
}

# Round the receipt results for release (results_receipt.csv). At full precision the fitted uptakes and delta
# are exact functions of the birth-bin counts. In England the bins are the months whose counts uptake_by_bin.csv
# releases rounded, and the two together give the exact counts back. HR_CACE, TR_CACE and the incidences with
# HR_CACE applied would give delta back, as results_eligibility.csv keeps the HR and TR at full precision. So
# estimates, CIs and the bandwidths for delta keep 3 decimal places, the incidences, which are proportions, 4,
# and SEs 3 significant figures.
round_receipt <- function(x) {
  d <- ifelse(startsWith(x$quantity, "inc_"), 4L, 3L)
  for (k in c("estimate", "lci", "uci")) x[[k]] <- round(x[[k]], d)
  for (k in c("boot_se", "se_conventional", "se_robust")) x[[k]] <- signif(x[[k]], 3)
  for (k in c("delta_bandwidth", "delta_bias_bandwidth")) x[[k]] <- round(x[[k]], 3)
  x
}
