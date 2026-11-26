# Inverse-variance fixed-effect meta-analysis of national estimates (SAP 6.4.1). Exclude non-finite estimates
# and those with non-finite or non-positive SEs. Report Cochran's Q. With at least two estimates, also
# calculate its p-value and, if Q is above 0, I-squared. Otherwise they are NA.
#
# Author/contact: Samantha Ip (hyi20@cam.ac.uk)

fe_pool <- function(est, se) {
  ok <- is.finite(est) & is.finite(se) & se > 0
  est <- est[ok]
  se <- se[ok]
  k <- length(est)
  if (!k) return(data.frame(estimate = NA_real_, se = NA_real_, k = 0L, q = NA_real_, p_het = NA_real_, i2 = NA_real_))
  w <- 1 / se^2
  b <- sum(w * est) / sum(w)
  q <- sum(w * (est - b)^2)
  data.frame(estimate = b, se = sqrt(1 / sum(w)), k = k, q = q,
             p_het = if (k > 1) stats::pchisq(q, k - 1, lower.tail = FALSE) else NA_real_,
             i2 = if (k > 1 && q > 0) max(0, (q - (k - 1)) / q) else NA_real_)
}
