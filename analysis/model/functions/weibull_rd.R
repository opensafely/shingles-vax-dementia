# Weibull proportional hazards (PH) RD model (SAP 6.4.1):
#   h(t | r, z) = p*t^(p-1)*exp(b0 + br*r + bz*z + brz*r*z),
# where t is years since the end of day 30 after the index date, r is years from the birth-date cutoff (the
# SAP's x_i^c, in years, to within a day), and z indicates eligibility. Thus bz is the log HR at the cutoff.
#
# Fit using survival::survreg's accelerated failure time (AFT) parameterisation, in which a_z, the
# coefficient of z, is the log time ratio (TR) at the cutoff, log HR = -a_z / scale and p = 1 / scale.
# Collapsing identical records into frequency-weighted rows preserves the likelihood, and makes an
# individual-level bootstrap sample a multinomial draw of the row frequencies.
#
# Author/contact: Samantha Ip (hyi20@cam.ac.uk)

collapse_rows <- function(x) {
  if (anyNA(x)) stop("collapse_rows() needs rows without missing values")   # NA would break the grouping
  n <- nrow(x)
  if (n) {
    x <- x[do.call(order, unname(as.list(x))), , drop = FALSE]
    new <- c(TRUE, rowSums(x[-1, , drop = FALSE] != x[-n, , drop = FALSE]) > 0)
    w <- tabulate(cumsum(new))
    x <- x[new, , drop = FALSE]
    x$n <- w
  } else {
    x$n <- integer(0)
  }
  rownames(x) <- NULL
  x$r <- x$r_days / 365.25                        # in years (the jump does not depend on the unit)
  x$z <- as.integer(x$r_days >= 0)                # eligibility from the recorded bin (derive_birth() checks)
  x$rz <- x$r * x$z
  x$time <- x$time_days / 365.25
  x
}

# Construct frequency-weighted model rows for the adults in `keep` who are at risk on day 31 for this
# outcome. Collapse by birth date, follow-up time, event status and the RZV receipt columns of fu named in
# `receipt` (by the end of follow-up), with n the number of adults in each row.
rd_data <- function(co, fu, keep = rep(TRUE, nrow(co)), receipt = "vacc_by_end") {
  k <- keep & fu$at_risk
  collapse_rows(cbind(data.frame(r_days = co$r_days[k], time_days = fu$time_days[k], event = fu$event[k]),
                      fu[k, receipt, drop = FALSE]))
}

# Fit the model with weights w (the counts n, or a bootstrap draw of them), adding the squared Newton
# decrement to the fit. Require events on both sides of the cutoff, without which the likelihood has no
# finite maximum, and a full-rank design, which needs at least two distinct birth dates on each side. The
# caller sets the minimum number of events per side, never below 1: run_nation.R passes the nation's
# min_events_side (8 in England, 10 in Wales) for every fit to the full data, and bootstrap refits keep the
# default of 1. The error message gives no count, as it goes to the released log.
fit_weibull_rd <- function(md, w = md$n, init = NULL, min_events_side = 1L) {
  md$w <- w
  md <- md[md$w > 0, , drop = FALSE]
  side <- vapply(0:1, function(s) sum((md$event * md$w)[md$z == s]), numeric(1))
  if (min(side) < max(1, min_events_side)) stop("too few events on one side of the cutoff to fit the model")
  if (qr(cbind(1, md$r, md$z, md$rz))$rank < 4) stop("too few distinct birth dates on one side of the cutoff to fit the model")
  if (is.null(init)) init <- exponential_start(md)
  # Right-censored Weibull regression with frequency weights and model-based standard errors. In Stata
  # terms: stset time [fweight=w], failure(event), then streg r z rz, distribution(weibull) time.
  f <- survival::survreg(survival::Surv(time, event) ~ r + z + rz, data = md, weights = md$w, dist = "weibull",
                         init = init, control = survival::survreg.control(maxiter = 100))
  f$newton_decrement <- newton_decrement(f, md)
  f
}

# Initialise the Weibull fit using an exponential PH model fitted by Poisson regression. In simulations,
# this improved stability under heavy censoring.
exponential_start <- function(md) {
  X <- cbind(1, md$r, md$z, md$rz)
  pg <- stats::glm.fit(X, md$event, weights = md$w, offset = log(md$time), family = stats::poisson())
  c(-pg$coefficients, 0)                  # with shape 1, the AFT coefficients are minus the PH ones
}

# Return the squared Newton decrement g'Vg, with g the score of the weighted log-likelihood and V survreg's
# covariance matrix (the inverse observed information), both at the estimate and in survreg's parameters
# (the AFT coefficients a and the log scale). A small value indicates approximate stationarity. The
# log-likelihood is concave in the log-hazard coefficients b and the shape p, as each adult contributes
# d log p, terms linear in (b, p) and -exp(p log t + x'b), each concave (d is the event indicator and
# x = (1, r, z, r*z)). survreg's parameters, a = -b / p and log scale = -log p, are a smooth function of
# (b, p) with a smooth inverse, so a stationary point in them is one in (b, p), and therefore the global
# maximum. A finite maximum need not exist, however, and the decrement can pass the convergence check
# without one, as when a side of the cutoff has no events (fit_weibull_rd() refuses such data).
newton_decrement <- function(f, md) {
  b <- stats::coef(f)
  V <- stats::vcov(f)
  if (!all(is.finite(b)) || !all(is.finite(V))) return(Inf)
  X <- cbind(1, md$r, md$z, md$rz)
  u <- (log(md$time) - drop(X %*% b)) / f$scale
  eu <- exp(u)
  g <- c(colSums(md$w * (eu - md$event) / f$scale * X), sum(md$w * (u * (eu - md$event) - md$event)))
  drop(t(g) %*% V %*% g)
}

# Flag approximate convergence when the squared Newton decrement is below 1e-6, so that a further Newton
# step would move each parameter by less than 0.001 of its SE. survreg's warnings, logged by fit_safely()
# in utils.R, do not set this flag, because survreg can stop short of the maximum without warning.
converged <- function(f) is.finite(f$newton_decrement) && f$newton_decrement < 1e-6

# For a fit that did not converge, keep the columns in `keep` (by default the convergence diagnostics) and
# leave the estimates blank (NA).
unconverged_na <- function(x, keep = c("iterations", "newton_decrement", "converged", "loglik")) {
  x[setdiff(names(x), keep)] <- NA
  x
}

# Return the eligibility log HR and log TR with model-based SEs (by the delta method for the log HR), the
# Weibull shape, the birth-date slopes either side of the cutoff and the fit diagnostics. The slopes are
# log HRs per year of later birth (per year younger), not per year of age.
weibull_effects <- function(f) {
  b <- stats::coef(f)
  V <- stats::vcov(f)
  s <- f$scale
  k <- c("z", "Log(scale)")
  g <- c(-1 / s, b[["z"]] / s)                    # gradient of -a_z / scale in (a_z, log scale), a_z = log TR
  data.frame(log_hr = -b[["z"]] / s,
             se_log_hr = sqrt(drop(t(g) %*% V[k, k] %*% g)),
             log_tr = b[["z"]],
             se_log_tr = sqrt(V["z", "z"]),
             shape = 1 / s,
             slope_ineligible = -b[["r"]] / s,          # log HR per year of birth date, below the cutoff
             slope_eligible = -(b[["r"]] + b[["rz"]]) / s,
             iterations = f$iter,
             newton_decrement = f$newton_decrement,
             converged = converged(f),
             loglik = f$loglik[2])
}

# Return model-based cumulative failure probabilities F(t) = 1 - S(t) at the cutoff (r = 0) for ineligible
# (z = 0) and eligible (z = 1) adults, and their difference. `years` counts from the index date, and the
# model's clock starts `induction_days` later (30 in run_nation.R), so each value is the probability of the
# outcome between the end of the induction period and `years` after the index date. Use
# log-cumulative-hazard CIs for each probability, and a delta-method SE and CI for the difference.
#
# Death and other exits are censored, so these are not competing-risks cumulative incidence functions.
# Estimates apply to adults free of the outcome and in follow-up at the start of day 31.
weibull_incidence <- function(f, years, induction_days = 0) {
  b <- stats::coef(f)
  V <- stats::vcov(f)
  s <- f$scale
  pars <- colnames(V)
  at <- function(t, z) {
    eta <- b[["(Intercept)"]] + z * b[["z"]]
    lh <- (log(t) - eta) / s                       # log cumulative hazard
    g <- stats::setNames(numeric(length(pars)), pars)
    g[c("(Intercept)", "z", "Log(scale)")] <- c(-1 / s, -z / s, -lh)
    H <- exp(lh)
    list(F = 1 - exp(-H), lh = lh, se = sqrt(drop(t(g) %*% V %*% g)), dF = exp(-H) * H * g)
  }
  do.call(rbind, lapply(years, function(y) {
    t <- y - induction_days / 365.25
    i0 <- at(t, 0)
    i1 <- at(t, 1)
    dd <- i1$dF - i0$dF
    se_diff <- sqrt(drop(t(dd) %*% V %*% dd))
    ci <- function(i, sgn) 1 - exp(-exp(i$lh + sgn * 1.96 * i$se))
    data.frame(years = y,
               inc_ineligible = i0$F, inc_ineligible_lci = ci(i0, -1), inc_ineligible_uci = ci(i0, 1),
               inc_eligible = i1$F, inc_eligible_lci = ci(i1, -1), inc_eligible_uci = ci(i1, 1),
               difference = i1$F - i0$F, se_difference = se_diff,
               difference_lci = i1$F - i0$F - 1.96 * se_diff, difference_uci = i1$F - i0$F + 1.96 * se_diff)
  }))
}

# Return the model-based average event rate per 1000 person-years, H(T) / T with H the cumulative hazard,
# over the planned follow-up of an adult born on each date of `dd` (from derive_birth()), from the end of
# day 30 after the index date to `end`. These rates are the fitted lines of the incidence-rate figures by
# birth month (fig_rates_*.png). Planned rather than observed follow-up is a disclosure control (README).
planned_rate <- function(f, dd, study, end) {
  b <- stats::coef(f)
  r <- dd$r_days / 365.25
  t <- as.numeric(end - dd$index_date - study$induction_days) / 365.25
  eta <- b[["(Intercept)"]] + b[["r"]] * r + b[["z"]] * dd$z + b[["rz"]] * r * dd$z
  ifelse(t > 0, 1000 * exp((log(pmax(t, 1e-9)) - eta) / f$scale) / t, NA_real_)
}
