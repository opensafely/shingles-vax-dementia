#!/usr/bin/env Rscript
# National regression discontinuity analysis with up to 3 years of follow-up (SAP 6.3, 6.4.1, 6.4.3, 6.5.1,
# 6.6 and 6.8). Run it separately in OpenSAFELY (England) and SAIL (Wales), and pool the released estimates
# with pool_nations.R.
#
# Author/contact: Samantha Ip (hyi20@cam.ac.uk)
#
# Usage:   Rscript run_nation.R --nation England|Wales --data FILE --out DIR [--boot N] [--data-end DATE]
# Options: --boot N          Number of bootstrap samples for the receipt effects (default 1000).
#          --data-end DATE   Last date of complete linked data, as YYYY-MM-DD (default: the nation's data_end in
#                            functions/config.R).
#
# Input:   --data, one row per adult with the variables in functions/required_variables.R, as a .csv, .csv.gz,
#          .rds, .parquet or .arrow file.
# Needs:   functions/*.R and the R packages survival, rdrobust (for delta), ggplot2 (for the figures) and arrow
#          (for a .arrow or .parquet extract).
#
# Outputs, in the --out folder (see README, "National analysis", "Outputs", for what each holds):
#   checks       data_checks.csv, records_by_calendar_month.csv
#   cohort       cohort_selection.csv, cohort_selection_by_bin.csv
#   descriptive  baseline_table.csv, baseline_smd.csv, events_summary.csv, events_by_bin.csv, uptake_by_bin.csv
#   density      density_test.csv, density_by_bin.csv
#   results      results_eligibility.csv, results_incidence.csv, results_receipt.csv
#   figures      fig_flowchart.png, fig_density.png, fig_rates_dementia.png, fig_rates_hz.png, fig_uptake.png
#   run record   settings.csv, log.txt, session_info.txt

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
this_file <- gsub("~+~", " ", this_file, fixed = TRUE)          # Rscript encodes spaces in its own path
here <- if (length(this_file)) dirname(normalizePath(this_file[1])) else getwd()
for (f in sort(list.files(file.path(here, "functions"), pattern = "\\.R$", full.names = TRUE))) source(f)

opt <- parse_args(commandArgs(TRUE), list(nation = NA, data = NA, out = NA, boot = STUDY$n_boot, "data-end" = NA))
nation <- opt$nation
if (!isTRUE(nation %in% names(NATIONS)) || is.na(opt$data) || is.na(opt$out)) {
  stop("Usage: Rscript run_nation.R --nation England|Wales --data FILE --out DIR [--boot N] [--data-end DATE]")
}
if (!requireNamespace("rdrobust", quietly = TRUE)) stop("delta needs the rdrobust package, which this R installation lacks")
cfg <- NATIONS[[nation]]
if (!file.exists(opt$data[1])) stop("--data: no such file: ", opt$data[1])
if (!is.na(opt[["data-end"]][1])) {
  cfg$data_end <- as.Date(opt[["data-end"]][1], format = "%Y-%m-%d")
  if (is.na(cfg$data_end)) stop("--data-end must be a date as YYYY-MM-DD")
}
n_boot <- suppressWarnings(as.numeric(opt$boot[1]))
if (!is.finite(n_boot) || n_boot < 1 || n_boot != round(n_boot)) stop("--boot must be a whole number of 1 or more")
n_boot <- as.integer(n_boot)
rule <- cfg$sdc
out <- opt$out
dir.create(out, recursive = TRUE, showWarnings = FALSE)
log_open(file.path(out, "log.txt"))
t0 <- Sys.time()
log_line("EPI-ZOSTER-110 regression discontinuity analysis: ", nation, " (", cfg$tre, ")")
log_line("Delta: local linear within the MSE-optimal bandwidth, with the robust bias-corrected CI, from rdrobust ",
         utils::packageVersion("rdrobust"))
log_line("Follow-up to ", format(min(STUDY$study_end, cfg$data_end)), ", with ", n_boot, " bootstrap samples")
log_line("Counts in this log are rounded to ", rule$round_to, ", and withheld at ", rule$redact_at_or_below, " or fewer")

# ---- Extract checks and derivations ----
# Check the extract and derive birth dates, eligibility, index dates, prior healthcare-use rates and, for each
# outcome, follow-up and vaccination receipt.
d <- prepare_extract(read_extract(opt$data), nation)
d <- derive_birth(d, cfg$bin, STUDY)
log_line("Extract: ", sdc_count(nrow(d), rule), " rows",
         if (attr(d, "dob_moved_to_middle") > 0) ". Birth dates on the first day of their month or week were moved to the middle" else "")
if (!"out_date_dementia_med" %in% names(d)) {
  log_line("out_date_dementia_med is not in the extract: dementia medication in the induction period is not excluded (SAP 5.2.1.2)")
}
d <- prior_rates(d)
fu <- lapply(STUDY$outcomes, function(v) follow_up(d, v, STUDY, cfg$data_end))
write_table(data_checks(d, nation, cfg$bin, STUDY, rule), out, "data_checks.csv")
# Doses recorded close together suggest that var_date_2vax_rzv is the next dose rather than the first at least
# 28 days after the first. A later second dose is then missing. Warn only when data_checks.csv releases the count.
if (sum(abs(d$var_date_2vax_rzv - d$var_date_1vax_rzv) < STUDY$min_dose_gap, na.rm = TRUE) > rule$redact_at_or_below) {
  log_line("WARNING: RZV doses recorded less than ", STUDY$min_dose_gap, " days apart (data_checks.csv). If the extract took the ",
           "next dose as the second, later second doses are missing (README, convention 3)")
}
# Monthly record counts help assess completeness near the study end.
by_month <- records_by_month(d, nation, STUDY, rule)
write_table(by_month, out, "records_by_calendar_month.csv")
for (lab in unique(by_month$record)) {
  seen <- by_month$month[by_month$record == lab & !is.na(by_month$n)]
  log_line("  last calendar month with more than ", rule$redact_at_or_below, " records, ", lab, ": ",
           if (length(seen)) max(seen) else "none")
}

# ---- Cohort selection (SAP 5.2.1, 6.3) ----
# Apply the selection criteria, which use the dementia follow-up, to define the cohort common to both outcomes.
crit <- cohort_criteria(d, nation, STUDY, fu[[STUDY$primary_outcome]])
labels <- names(crit)
ff <- first_failed(crit)
sel <- write_table(selection_table(d, ff, labels, nation, rule), out, "cohort_selection.csv")
write_table(selection_by_bin(d, ff, labels, nation, STUDY, rule), out, "cohort_selection_by_bin.csv")
# Report events on days 0 to 30 separately: dementia among adults who passed every earlier criterion, and HZ
# among adults in the analysis cohort.
in_induction <- function(o) (fu[[o]]$event == 1L & fu[[o]]$time_days <= 0L) %in% TRUE
early <- lapply(stats::setNames(names(fu), names(fu)), function(o) {
  among <- if (o == STUDY$primary_outcome) ff == match("No dementia in the 30-day induction period", labels) else ff == 0L
  d[among & in_induction(o), c("dob", "z"), drop = FALSE]
})
co <- d[ff == 0L, , drop = FALSE]
fu <- lapply(fu, function(x) x[ff == 0L, , drop = FALSE])
co$vacc_by_date <- fu[[STUDY$primary_outcome]]$vacc_by_date            # for the baseline table
co$two_doses_by_date <- fu[[STUDY$primary_outcome]]$two_doses_by_date
rm(d, crit)
log_line("Analysis cohort: ", sdc_count(nrow(co), rule), " adults")

# ---- Descriptive analyses (SAP 6.3, 6.5.1) ----
# The baseline table includes the available comorbidities, the number each adult has, and age on 31 August 2023
# and at the index date.
co$age_aug2023 <- age_years(co$dob, STUDY$age_date)
co$age_index <- age_years(co$dob, co$index_date)
present <- intersect(names(COMORBIDITIES), names(co))
table1 <- c(TABLE1_VARS, lapply(present, function(v) list(var = v, type = "bin", label = COMORBIDITIES[[v]])))
if (length(present)) {
  co$n_comorbidities <- comorbidity_count(co, present)
  table1 <- c(table1, list(list(var = "n_comorbidities", type = "cat", label = "Number of the comorbidities above")))
}
if (length(present) < length(COMORBIDITIES)) {
  log_line("Comorbidities not in the extract yet: ", paste(setdiff(names(COMORBIDITIES), present), collapse = ", "))
}
# The subgroups of SAP 6.6 (ii) count only the flags recorded for at least 0.1% of the cohort and more than 7
# adults (subgroup_flags() in describe.R).
sub_flags <- subgroup_flags(co, present, STUDY, rule)
if (length(sub_flags)) co$comorbidity_subgroup <- comorbidity_subgroup(co, sub_flags)
if (length(sub_flags) < length(present)) {
  log_line("Comorbidities recorded for under ", 100 * STUDY$min_flag_share, "% of the cohort, or for ", rule$redact_at_or_below,
           " or fewer adults, are not counted in the comorbidity subgroups: ", paste(setdiff(present, sub_flags), collapse = ", "))
}
bt <- write_table(baseline_table(co, table1, nation, rule), out, "baseline_table.csv")
write_table(baseline_smd(bt), out, "baseline_smd.csv")
write_table(do.call(rbind, lapply(names(fu), function(o) events_summary(co, fu[[o]], nation, o, STUDY, rule, early[[o]]))),
            out, "events_summary.csv")
# Birth-month uptake uses the dementia follow-up, as does the dementia delta.
uptake <- uptake_by_bin(co, fu[[STUDY$primary_outcome]], nation, rule)
write_table(uptake, out, "uptake_by_bin.csv")

# ---- Birth-date density test (SAP 6.4.1) ----
# Test the density in the analysis cohort, because cohort selection may itself introduce a discontinuity. Bins
# are the recorded birth months (England) or weeks (Wales).
bins <- density_bins(co, STUDY$cutoff)
mc <- if (nrow(bins)) fit_safely("McCrary test", mccrary_test(bins, max_h = STUDY$bandwidths[[STUDY$primary_bandwidth]])) else NULL
# Withhold the fitted densities and slopes, and round the test statistics, to limit reconstruction of the
# underlying birth-bin counts.
mc_keep <- c("theta", "se", "z", "p", "bandwidth_days", "n", "bins_left", "bins_right")
mc_round <- if (is.null(mc)) NULL else transform(mc, n = sdc_count(n, rule), theta = round(theta, 4), se = round(se, 4),
                                                 z = round(z, 2), p = round(p, 3), bandwidth_days = round(bandwidth_days))
mc_out <- if (is.null(mc)) as.data.frame(stats::setNames(rep(list(NA_real_), length(mc_keep)), mc_keep)) else mc_round[mc_keep]
write_table(cbind(nation = nation, mc_out), out, "density_test.csv")
bins$n <- sdc_count(bins$n, rule)
write_table(cbind(nation = rep(nation, nrow(bins)), bins), out, "density_by_bin.csv")

# ---- Eligibility effects (SAP 6.4.1, 6.4.3, 6.5.1) ----
effect_na <- data.frame(log_hr = NA_real_, se_log_hr = NA_real_, log_tr = NA_real_, se_log_tr = NA_real_,
                        shape = NA_real_, slope_ineligible = NA_real_, slope_eligible = NA_real_,
                        iterations = NA_integer_, newton_decrement = NA_real_, converged = FALSE, loglik = NA_real_)
elig <- list()
inc <- list()
fits <- list()
for (o in names(STUDY$outcomes)) {
  for (w in names(STUDY$bandwidths)) {
    inside <- in_bandwidth(co, STUDY, STUDY$bandwidths[[w]])
    for (s in STUDY$strata) {
      # Adults with identical model data share one row of md, and md$n counts them.
      md <- rd_data(co, fu[[o]], inside & (s == "All" | co$cov_cat_sex %in% s))
      f <- fit_safely(sprintf("%s, %s-year bandwidth, %s", o, w, s),
                      fit_weibull_rd(md, min_events_side = STUDY$min_events_side))
      key <- data.frame(nation = nation, outcome = o, bandwidth = w, stratum = s)
      counts <- data.frame(n = sdc_count(sum(md$n), rule), events = sdc_count(sum(md$n * md$event), rule),
                           person_years = sdc_stat(round(sum(md$n * md$time)), sum(md$n), rule))
      # Retain the fit diagnostics, but leave blank the estimates of a fit that did not converge.
      eff <- if (is.null(f)) effect_na else weibull_effects(f)
      if (!is.null(f) && !converged(f)) eff <- unconverged_na(eff)
      elig[[length(elig) + 1]] <- cbind(key, counts, eff,
                                        warnings = if (is.null(f)) "not fitted" else paste(attr(f, "warnings"), collapse = "; "))
      if (!is.null(f)) {
        ic <- weibull_incidence(f, STUDY$incidence_years, STUDY$induction_days)
        if (!converged(f)) ic <- unconverged_na(ic, keep = "years")
        inc[[length(inc) + 1]] <- cbind(key, ic, converged = converged(f))
        # The fitted line by birth month is drawn only from a primary fit that converged.
        if (w == STUDY$primary_bandwidth && s == "All" && converged(f)) fits[[o]] <- f
      }
    }
  }
}
elig <- do.call(rbind, elig)
elig$hr <- exp(elig$log_hr)
elig$hr_lci <- exp(elig$log_hr - 1.96 * elig$se_log_hr)
elig$hr_uci <- exp(elig$log_hr + 1.96 * elig$se_log_hr)
# Two-sided Wald test of log HR = 0 using its delta-method SE. survreg's z test of the AFT coefficient (log TR)
# tests the same null, but Wald tests depend on the parameterisation, so its p-value differs slightly.
elig$p <- 2 * stats::pnorm(-abs(elig$log_hr / elig$se_log_hr))
elig$tr <- exp(elig$log_tr)
elig$tr_lci <- exp(elig$log_tr - 1.96 * elig$se_log_tr)
elig$tr_uci <- exp(elig$log_tr + 1.96 * elig$se_log_tr)
write_table(elig, out, "results_eligibility.csv")
write_table(if (length(inc)) do.call(rbind, inc) else data.frame(nation = character(0)), out, "results_incidence.csv")

# Model-based rates over planned, not observed, follow-up for each birth month.
for (o in setdiff(names(STUDY$outcomes), names(fits))) {
  log_line("  No fitted line by birth month for ", o, ": the primary model was not fitted or did not converge")
}
months <- sort(unique(co$month_start))
at_month <- derive_birth(data.frame(var_date_birth = months), "month", STUDY)
ev <- lapply(names(STUDY$outcomes), function(o) {
  rate <- if (is.null(fits[[o]])) NULL else
    stats::setNames(planned_rate(fits[[o]], at_month, STUDY, min(STUDY$study_end, cfg$data_end)), as.character(months))
  events_by_bin(co, fu[[o]], nation, o, rule, rate)
})
names(ev) <- names(STUDY$outcomes)
write_table(do.call(rbind, ev), out, "events_by_bin.csv")

# ---- Receipt effects (SAP 6.5.1, 6.6) ----
# Estimate the receipt effects for each outcome and bandwidth with SAP 6.5.1's delta ("main"). SAP 6.6's
# sensitivity analyses concern secondary objective 3, the effect of receipt on dementia, and so are made for
# dementia at the primary bandwidth: the local quadratic delta (i) and two doses (iii) on the same bootstrap
# samples as main, and main within comorbidity subgroups (ii). HZ and the narrower bandwidths are SAP 6.4.3's
# own sensitivity analyses. Each analysis has its own fixed seed.
receipt_rows <- function(md, label, w, defs, seed, key, converged, p) {
  specs <- list()
  for (k in defs) {
    hb <- fit_safely(sprintf("%s: bandwidths for delta (%s)", label, k), delta_bandwidths(md, STUDY$bandwidths[[w]], DELTAS[[k]]))
    if (is.null(hb)) next
    specs[[k]] <- c(DELTAS[[k]], as.list(hb))
    if (DELTAS[[k]]$bandwidth == "mserd") {
      log_line("  ", label, ": MSE-optimal bandwidths for delta (", k, "), ", sprintf("%.3f", hb[["h"]]), " years and ",
               sprintf("%.3f", hb[["b"]]), " for the bias")
    }
  }
  t1 <- Sys.time()
  r <- if (length(specs)) fit_safely(label, receipt_effect(md, n_boot, seed, specs, STUDY$min_events_side,
                                                           STUDY$incidence_years, STUDY$induction_days)) else NULL
  if (!is.null(r)) {
    for (k in unique(r$analysis)) {
      nf <- r$n_boot_failed[r$analysis == k][1]
      if (nf > STUDY$max_boot_failed * n_boot) {
        log_line("  WARNING: ", label, " (", k, "): more than ", 100 * STUDY$max_boot_failed, "% of bootstrap samples failed, ",
                 "so the bootstrap intervals and SEs are left missing")
      } else if (nf > 0.01 * n_boot) {
        log_line("  WARNING: ", label, " (", k, "): more than 1% of bootstrap samples failed. ",
                 "The bootstrap intervals use the successful samples only and may not have their nominal coverage")
      }
      u <- r$estimate[r$analysis == k & r$quantity %in% c("uptake_ineligible", "uptake_eligible")]
      if (any(u < -1e-8 | u > 1 + 1e-8, na.rm = TRUE)) {
        log_line("  WARNING: ", label, ": the fitted uptake at the cutoff (", k, ") is outside 0 to 1, so the fit for delta extrapolates")
      }
      if (isTRUE(r$lci[r$analysis == k & r$quantity == "delta"] <= 0)) {
        log_line("  WARNING: ", label, ": the 95% CI for delta (", k, ") reaches 0 or below. ",
                 "HR_CACE = HR_eligibility^(1/delta), so the receipt effects are not reliable")
      }
    }
    for (k in names(attr(r, "delta_failures"))) {
      log_line("  ", label, ": no delta (", k, "): ", attr(r, "delta_failures")[[k]])
    }
    for (k in names(attr(r, "rbc_failures"))) {
      log_line("  ", label, ": no robust bias-corrected CI for delta (", k, "): ", attr(r, "rbc_failures")[[k]])
    }
    for (k in attr(r, "uptake_withheld")) {
      log_line("  ", label, ": fitted uptakes at the cutoff (", k, ") withheld, as one side has ", STUDY$min_events_side - 1L,
               " or fewer adults with, or without, receipt within the bandwidth for delta")
    }
    log_line("  ", label, ": ", round(as.numeric(Sys.time() - t1, units = "mins"), 1), " minutes, ",
             paste0(tapply(r$n_boot_failed, r$analysis, `[`, 1)[unique(r$analysis)], " of ", n_boot, " samples failed (",
                    unique(r$analysis), ")", collapse = ", "))
    fails <- attr(r, "failures")
    if (length(fails)) log_line("    failed samples by reason: ", paste0(names(fails), " (", as.integer(fails), ")", collapse = "; "))
  }
  rows <- receipt_template(defs, STUDY$incidence_years)
  rows$n_boot <- n_boot
  if (!is.null(r)) rows[match(paste(r$analysis, r$quantity), paste(rows$analysis, rows$quantity)), ] <- r
  bw_of <- function(k, x) if (is.null(specs[[k]])) NA_real_ else specs[[k]][[x]]
  cbind(key, rows,
        converged = converged,
        receipt_estimated = !is.null(r) & rows$analysis %in% r$analysis,
        p_eligibility = p,                            # also tests HR_CACE = 1, if delta is not 0
        delta_bandwidth = vapply(rows$analysis, bw_of, 0, x = "h", USE.NAMES = FALSE),
        delta_bias_bandwidth = vapply(rows$analysis, bw_of, 0, x = "b", USE.NAMES = FALSE))
}
rec <- list()
for (o in names(STUDY$outcomes)) {
  for (w in names(STUDY$bandwidths)) {
    defs <- if (o == STUDY$primary_outcome && w == STUDY$primary_bandwidth) names(DELTAS) else "main"
    md <- rd_data(co, fu[[o]], in_bandwidth(co, STUDY, STUDY$bandwidths[[w]]),
                  unique(vapply(DELTAS[defs], function(x) x$receipt, "")))
    # The receipt analysis uses the same sample and Weibull specification as the overall eligibility analysis
    # for this outcome and bandwidth. Inherit its convergence flag and p-value. receipt_estimated separately
    # records whether each receipt analysis ran without error.
    same <- elig$outcome == o & elig$bandwidth == w & elig$stratum == "All"
    rec[[length(rec) + 1]] <- receipt_rows(md, sprintf("Receipt effect, %s, %s-year bandwidth", o, w), w, defs,
                                           seed = STUDY$seed + 100L * match(o, names(STUDY$outcomes)) + match(w, names(STUDY$bandwidths)),
                                           key = data.frame(nation = nation, outcome = o, bandwidth = w, subgroup = "All"),
                                           converged = elig$converged[same], p = elig$p[same])
  }
}
# SAP 6.6 (ii): main within each comorbidity subgroup. Its Weibull fit is fitted here first, as the
# eligibility analysis has no such subgroups, to give the convergence flag and p-value.
if (length(sub_flags)) {
  o <- STUDY$primary_outcome
  w <- STUDY$primary_bandwidth
  inside <- in_bandwidth(co, STUDY, STUDY$bandwidths[[w]])
  for (g in COMORBIDITY_SUBGROUPS) {
    md <- rd_data(co, fu[[o]], inside & co$comorbidity_subgroup == g)
    label <- sprintf("Receipt effect, %s, %s-year bandwidth, %s", o, w, g)
    f <- fit_safely(sprintf("Eligibility effect, %s, %s-year bandwidth, %s", o, w, g),
                    fit_weibull_rd(md, min_events_side = STUDY$min_events_side))
    ok <- !is.null(f) && converged(f)
    eff <- if (ok) weibull_effects(f) else NULL
    rec[[length(rec) + 1]] <- receipt_rows(md, label, w, "main", seed = STUDY$seed + 1000L + match(g, COMORBIDITY_SUBGROUPS),
                                           key = data.frame(nation = nation, outcome = o, bandwidth = w, subgroup = g),
                                           converged = ok, p = if (ok) 2 * stats::pnorm(-abs(eff$log_hr / eff$se_log_hr)) else NA_real_)
  }
} else {
  log_line("No comorbidity flags to count, so no comorbidity subgroups (SAP 6.6 (ii))")
}
write_table(do.call(rbind, rec), out, "results_receipt.csv")

# ---- Figures ----
save_plot(plot_flowchart(sel, sprintf("%s: cohort selection", nation)), file.path(out, "fig_flowchart.png"),
          width = 8, height = 0.5 * (length(labels) + 4))
unit <- if (cfg$bin == "month") "month" else "week"
save_plot(plot_density(bins, mc_round, STUDY$cutoff, sprintf("%s: adults by recorded birth %s", nation, unit)),
          file.path(out, "fig_density.png"))
for (o in names(ev)) {
  save_plot(plot_rates(ev[[o]], STUDY$cutoff, sprintf("%s: %s incidence by birth month", nation, if (o == "hz") "HZ" else o)),
            file.path(out, sprintf("fig_rates_%s.png", o)))
}
save_plot(plot_uptake(uptake, STUDY$cutoff, sprintf("%s: RZV uptake by birth month", nation)),
          file.path(out, "fig_uptake.png"))

# Save the settings, which pool_nations.R checks are compatible across the nations.
code <- c(file.path(here, "run_nation.R"), sort(list.files(file.path(here, "functions"), pattern = "\\.R$", full.names = TRUE)))
settings <- c(nation = nation, cutoff = format(STUDY$cutoff), study_end = format(STUDY$study_end),
              follow_up_end = format(min(STUDY$study_end, cfg$data_end)), induction_days = STUDY$induction_days,
              bandwidths = paste(names(STUDY$bandwidths), collapse = " "), primary_bandwidth = STUDY$primary_bandwidth,
              incidence_years = paste(STUDY$incidence_years, collapse = " "), min_events_side = STUDY$min_events_side,
              min_dose_gap = STUDY$min_dose_gap, subgroup_flags = paste(sub_flags, collapse = " "),
              n_boot = n_boot, seed = STUDY$seed,
              code_md5 = code_checksum(code), r_version = R.version.string)
write_table(data.frame(setting = names(settings), value = unname(settings)), out, "settings.csv")
writeLines(utils::capture.output(utils::sessionInfo()), file.path(out, "session_info.txt"))
log_line("Done in ", round(as.numeric(Sys.time() - t0, units = "mins"), 1), " minutes")
