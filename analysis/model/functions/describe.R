# Extract checks and descriptive analyses (SAP 6.3, 6.5.1): baseline characteristics, SMDs, events,
# person-time, crude incidence rates and RZV uptake. Disclosure control (utils.R) is applied before release.
#
# Author/contact: Samantha Ip (hyi20@cam.ac.uk)

# ---- Data checks ----

# Count records violating extract conventions or plausibility checks, in the full extract after
# derive_birth(). Duplicate identifiers are counted but not removed.
data_checks <- function(d, nation, bin, study, rule) {
  idc <- intersect(c("patient_id", "alf_pe", "ALF_PE"), names(d))
  off <- as.integer(d$var_date_birth - d$bin_start)       # days into the recorded month or week
  index_gap <- abs(as.integer(d$var_date_index - d$index_date))
  earliest <- pmax(add_years(d$bin_start, 65L), study$nip_start)   # earliest possible 65th birthday, or 1 Sep 2023
  gap <- as.integer(d$var_date_2vax_rzv - d$var_date_1vax_rzv)     # days from the first RZV dose to the second
  n <- c(
    "Rows in the extract" = nrow(d),
    "Duplicate patient identifiers" = if (length(idc)) sum(duplicated(d[[idc[1]]])) else NA,
    "Birth date missing" = sum(is.na(d$var_date_birth)),
    "Birth date on the first day of its month or week (moved to the middle)" = sum(off %in% 0L),
    "Birth date at the middle of its month or week (the 15th or the Thursday)" = sum(off %in% if (bin == "month") 14L else 3L),
    "cov_int_age differs from the age derived from the birth date" = sum(d$cov_int_age != d$age_derived, na.rm = TRUE),
    "var_date_index differs from the derived index date" = sum(index_gap > 0, na.rm = TRUE),
    "var_date_index differs from the derived index date by more than 31 days" = sum(index_gap > 31, na.rm = TRUE),
    "Dementia dated before the index date" = sum(d$out_date_dementia < d$index_date, na.rm = TRUE),
    "HZ dated before the index date" = sum(d$out_date_hz < d$index_date, na.rm = TRUE),
    "First RZV dose before the index date (counted as vaccinated)" = sum(d$var_date_1vax_rzv < d$index_date, na.rm = TRUE),
    "First RZV dose before the earliest possible 65th birthday (1 September 2023 if ineligible)" =
      sum(d$var_date_1vax_rzv < earliest, na.rm = TRUE),
    "Second RZV dose with no first dose" = sum(is.na(d$var_date_1vax_rzv) & !is.na(d$var_date_2vax_rzv)),
    "Second RZV dose dated before the first (the two are taken in date order)" = sum(gap < 0, na.rm = TRUE),
    stats::setNames(sum(abs(gap) < study$min_dose_gap, na.rm = TRUE),
                    sprintf("RZV doses less than %d days apart (the later is not counted as a second dose)",
                            study$min_dose_gap)),
    "First RZV dose dated after death" = sum(d$var_date_1vax_rzv > d$var_date_death, na.rm = TRUE),
    "Registered on the index date but deregistered by then" =
      sum(d$var_bin_registration %in% TRUE & d$var_date_deregistered <= d$index_date, na.rm = TRUE),
    "Death dated before the birth date" = sum(d$var_date_death < d$dob, na.rm = TRUE),
    "var_date_lookback missing or not before the index date (prior rates missing)" =
      sum(!((d$var_date_lookback < d$index_date) %in% TRUE))
  )
  if ("out_date_dementia_med" %in% names(d)) {
    med <- d$out_date_dementia_med
    n <- c(n, "Dementia medication date recorded" = sum(!is.na(med)),
           "Dementia medication dated before the index date" = sum(med < d$index_date, na.rm = TRUE),
           "Dementia medication on the index date or in the 30 days after it" =
             sum(med >= d$index_date & med <= d$index_date + study$induction_days, na.rm = TRUE))
  }
  flags <- c("var_bin_residency", "var_bin_registration", "pre_bin_dementia",
             "pre_bin_ic", "pre_bin_vax_shingles")
  miss <- vapply(flags, function(v) sum(is.na(d[[v]])), numeric(1))
  names(miss) <- paste0(flags, " missing (", ifelse(startsWith(flags, "pre_"), "not excluded", "criterion not met"), ")")
  n <- c(n, miss)
  out <- data.frame(nation = nation, check = names(n), n = sdc_count(n, rule))
  out
}

# Count records by calendar month to assess data completeness (SAP Table 3, "end of data availability").
records_by_month <- function(d, nation, study, rule) {
  vars <- c("dementia (first record on or after the index date)" = "out_date_dementia",
            "HZ (first record on or after the index date)" = "out_date_hz", "death" = "var_date_death",
            "deregistration" = "var_date_deregistered", "first RZV dose" = "var_date_1vax_rzv")
  months <- seq(first_of_month(study$nip_start), first_of_month(study$study_end), by = "month")
  do.call(rbind, lapply(names(vars), function(lab) {
    x <- d[[vars[[lab]]]]
    x <- x[!is.na(x) & x >= study$nip_start & x <= study$study_end]
    n <- tabulate(match(format(x, "%Y-%m"), format(months, "%Y-%m")), length(months))
    data.frame(nation = nation, record = lab, month = format(months, "%Y-%m"), n = sdc_count(n, rule))
  }))
}

# ---- Baseline characteristics ----

# Define the baseline table's groups (SAP 6.3, Table 8): the whole cohort, by eligibility and, among eligible
# adults, by vaccination status.
baseline_groups <- function(co) list(
  "All"                      = rep(TRUE, nrow(co)),
  "Eligible"                 = co$z == 1L,
  "Ineligible"               = co$z == 0L,
  "Eligible, vaccinated"     = co$z == 1L & co$vacc_by_date == 1L,       # at least one dose by the administrative end date
  "Eligible, two doses"      = co$z == 1L & co$two_doses_by_date == 1L,  # a subset of the vaccinated
  "Eligible, not vaccinated" = co$z == 1L & co$vacc_by_date == 0L
)

# Count recorded comorbidities per adult among the flags given. Treat missing condition flags as absent, as
# SAP 6.4.2 does for a diagnosis with no record.
comorbidity_number <- function(co, flags) Reduce(`+`, lapply(flags, function(v) co[[v]] %in% TRUE), integer(nrow(co)))

# Group the number as SAP Table 6 does, for the baseline table: 0, 1, 2, 3, 4, 5 or more.
comorbidity_count <- function(co, flags) {
  n <- comorbidity_number(co, flags)
  ifelse(n >= 5, "5 or more", as.character(n))
}

# Group the number into COMORBIDITY_SUBGROUPS (config.R) for SAP 6.6 (ii), which estimates the effect of
# receipt by baseline health, as healthier adults may be more likely to be vaccinated. The SAP names no
# grouping. The number of conditions, which Table 6 also defines, puts baseline health in one variable, where
# a subgroup for each condition would give many small analyses. Three groups, 0, 1, and 2 or more, order
# adults by it while keeping each group large enough for its own Weibull model and delta.
comorbidity_subgroup <- function(co, flags) COMORBIDITY_SUBGROUPS[pmin(comorbidity_number(co, flags), 2L) + 1L]

# The flags that the subgroups count: those recorded for at least min_flag_share (0.1%) of the cohort and for
# more than 7 adults. An empty or unfinished codelist would otherwise put every adult in the 0 group, and a
# rarer condition, such as Huntington's disease at these ages, changes almost no adult's group. Leaving it
# out in both nations keeps their subgroups the same, as pooling requires.
subgroup_flags <- function(co, flags, study, rule) {
  n <- vapply(flags, function(v) sum(co[[v]] %in% TRUE), numeric(1))
  flags[n >= study$min_flag_share * nrow(co) & n > rule$redact_at_or_below]
}

# Recode categorical and yes/no variables, retaining missing values as a category.
# Return NULL for numeric variables.
as_levels <- function(x, spec) {
  if (spec$type == "num") return(NULL)
  v <- if (spec$type == "bin") ifelse(x, "Yes", "No") else as.character(x)
  v[is.na(v)] <- "Missing"
  v
}

# Omit variables that config.R restricts to other nations. Fix each variable's categories on the whole cohort,
# so that every group shows the same rows.
prepare_vars <- function(co, vars, nation) {
  vars <- Filter(function(s) is.null(s$nations) || nation %in% s$nations, vars)
  lapply(vars, function(s) {
    v <- as_levels(co[[s$var]], s)
    s$all_levels <- switch(s$type,
                           cat = c(if (is.null(CATEGORIES[[s$var]])) sort(unique(v[v != "Missing"])) else
                                     CATEGORIES[[s$var]], "Missing"),           # the SAP's categories, in its order
                           bin = c("Yes", "No", "Missing"),
                           num = c("Summary", "Missing"))
    s
  })
}

# Summarise numeric variables by non-missing count, mean, SD, median and quartiles (R type 7), with the
# missing count in a separate row. Tabulate categorical and yes/no variables, including missing values.
summarise_group <- function(co, keep, vars) {
  na <- NA_real_
  rows <- list(data.frame(variable = "N", label = "Adults", level = "", n = sum(keep),
                          mean = na, sd = na, median = na, q1 = na, q3 = na))
  for (s in vars) {
    x <- co[[s$var]][keep]
    if (s$type == "num") {
      ok <- !is.na(x)
      q <- if (any(ok)) stats::quantile(x[ok], c(0.5, 0.25, 0.75), names = FALSE) else rep(na, 3)
      rows[[length(rows) + 1]] <- data.frame(variable = s$var, label = s$label, level = s$all_levels,
                                             n = c(sum(ok), sum(!ok)),
                                             mean = c(if (any(ok)) mean(x[ok]) else na, na),
                                             sd = c(if (sum(ok) > 1) stats::sd(x[ok]) else na, na),
                                             median = c(q[1], na), q1 = c(q[2], na), q3 = c(q[3], na))
    } else {
      n <- as.numeric(table(factor(as_levels(x, s), levels = s$all_levels)))
      rows[[length(rows) + 1]] <- data.frame(variable = s$var, label = s$label, level = s$all_levels, n = n,
                                             mean = na, sd = na, median = na, q1 = na, q3 = na)
    }
  }
  do.call(rbind, rows)
}

# Apply count suppression and rounding before calculating percentages. Percentages are of the whole group,
# missing values included. Suppress a mean, SD, median or quartile that describes 7 or fewer adults, and
# round the rest to 3 decimal places to limit disclosure risk (README, "National analysis", "Outputs").
# Three places also keep the SMDs that baseline_smd() computes from the rounded means and SDs close to exact.
baseline_table <- function(co, vars, nation, rule) {
  vars <- prepare_vars(co, vars, nation)
  groups <- baseline_groups(co)
  parts <- list()
  for (g in names(groups)) {
    tab <- summarise_group(co, groups[[g]], vars)
    raw_n <- tab$n
    tab$n <- sdc_count(raw_n, rule)
    categorical <- !tab$level %in% c("", "Summary")
    tab$pct <- ifelse(categorical, 100 * tab$n / tab$n[1], NA)
    for (k in c("mean", "sd", "median", "q1", "q3")) tab[[k]] <- round(sdc_stat(tab[[k]], raw_n, rule), 3)
    parts[[length(parts) + 1]] <- cbind(nation = nation, group = g, tab)
  }
  do.call(rbind, parts)
}

smd_binary <- function(p1, p0) {
  s <- sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
  ifelse(s > 0, (p1 - p0) / s, ifelse(p1 == p0, 0, NA))
}

smd_numeric <- function(m1, s1, m0, s0) {
  s <- sqrt((s1^2 + s0^2) / 2)
  ifelse(s > 0, (m1 - m0) / s, ifelse(m1 == m0, 0, NA))
}

# Calculate SMDs (SAP 6.3, 6.5.1) from the released baseline table, not from the confidential counts, means
# and SDs, so that an SMD reveals nothing the table withholds. Compare eligible with ineligible adults and,
# among eligible adults, vaccinated with not vaccinated. Treat categories as binary indicators. Suppress
# SMDs when required inputs are withheld.
baseline_smd <- function(bt) {
  comps <- list("Eligible vs ineligible" = c("Eligible", "Ineligible"),
                "Vaccinated vs not vaccinated (eligible)" = c("Eligible, vaccinated", "Eligible, not vaccinated"))
  key <- function(x) paste(x$nation, x$variable, x$level, sep = "|")
  rows <- bt[bt$level != "", ]                                        # every row but the group sizes
  out <- do.call(rbind, lapply(names(comps), function(cname) {
    a <- rows[rows$group == comps[[cname]][1], ]
    b <- rows[rows$group == comps[[cname]][2], ]
    b <- b[match(key(a), key(b)), ]
    smd <- ifelse(a$level == "Summary", smd_numeric(a$mean, a$sd, b$mean, b$sd), smd_binary(a$pct / 100, b$pct / 100))
    data.frame(nation = a$nation, comparison = cname, variable = a$variable, level = a$level,
               smd = round(smd, 3))
  }))
  rownames(out) <- NULL
  out
}

# ---- Events, person-years and uptake by birth month (both nations) ----

month_middle <- function(months) months + 14L           # the 15th, where the figures place each birth month

# Tabulate adults, events, person-years and crude rates per 1000 person-years by birth month, for plotting.
# `model_rate` (optional) gives the fitted model's average rate over planned follow-up, named by each
# month's first day (planned_rate() in weibull_rd.R).
events_by_bin <- function(co, fu, nation, outcome, rule, model_rate = NULL) {
  k <- fu$at_risk
  if (!any(k)) return(data.frame(nation = character(0), outcome = character(0), bin_start = as.Date(character(0)),
                                 bin_middle = as.Date(character(0)), eligible = integer(0), n = numeric(0), events = numeric(0), person_years = numeric(0),
                                 rate_per_1000py = numeric(0), model_rate_per_1000py = numeric(0)))
  agg <- rowsum(cbind(n = 1, events = fu$event[k], days = fu$time_days[k]), as.character(co$month_start[k]))
  bins <- as.Date(rownames(agg))
  out <- data.frame(nation = nation, outcome = outcome, bin_start = bins, bin_middle = month_middle(bins),
                    eligible = as.integer(bins >= STUDY$cutoff),
                    n = sdc_count(agg[, "n"], rule), events = sdc_count(agg[, "events"], rule),
                    person_years = sdc_stat(round(agg[, "days"] / 365.25), agg[, "n"], rule))
  out$rate_per_1000py <- ifelse(out$person_years > 0, 1000 * out$events / out$person_years, NA_real_)   # from the released (rounded) counts
  out$model_rate_per_1000py <- if (is.null(model_rate)) NA_real_ else unname(model_rate[as.character(bins)])
  rownames(out) <- NULL
  out
}

# Tabulate adults, events from day 31 and person-years by eligibility group within each bandwidth, with
# crude rates per 1000 person-years. `early` holds the adults with the outcome on days 0 to 30, which the
# analysis leaves out. Count them in events_days_0_30.
events_summary <- function(co, fu, nation, outcome, study, rule, early = NULL) {
  parts <- list()
  for (w in names(study$bandwidths)) {
    inside <- in_bandwidth(co, study, study$bandwidths[[w]]) & fu$at_risk
    for (g in c("Eligible", "Ineligible")) {
      keep <- inside & co$z == (g == "Eligible")
      parts[[length(parts) + 1]] <- data.frame(
        nation = nation, outcome = outcome, bandwidth = w, group = g,
        n = sum(keep), events = sum(fu$event[keep]), person_years = sum(fu$time_days[keep]) / 365.25,
        events_days_0_30 = if (is.null(early)) NA_real_ else
          sum(in_bandwidth(early, study, study$bandwidths[[w]]) & early$z == (g == "Eligible")))
    }
  }
  out <- do.call(rbind, parts)
  out$person_years <- sdc_stat(round(out$person_years), out$n, rule)
  out$n <- sdc_count(out$n, rule)
  out$events <- sdc_count(out$events, rule)
  out$events_days_0_30 <- sdc_count(out$events_days_0_30, rule)
  out$rate_per_1000py <- ifelse(out$person_years > 0, 1000 * out$events / out$person_years, NA_real_)   # from the released (rounded) counts
  out
}

# Tabulate events, person-years and crude rates per 1000 person-years by age in whole years, for the sample
# size calculation (sample_size_incidence.R). Each day at risk counts at the adult's age on that day, from day
# 31 after the index date to the end of follow-up (follow_up() in derive.R), and an event at the age on its
# date. adults_followed counts the adults with any follow-up at that age, so an adult followed across a
# birthday counts at both ages, and the ages' counts add up to more than the last row's, which covers all
# ages together and counts each adult once. Events and person-years add up across the ages.
incidence_by_age <- function(co, fu, nation, study, rule) {
  k <- fu$at_risk
  dob <- co$dob[k]
  first <- co$index_date[k] + study$induction_days + 1L        # the first day at risk
  last <- first + fu$time_days[k] - 1L                           # the last, which is the event's date for an event
  event <- fu$event[k] == 1L
  ages <- if (any(k)) seq(min(age_on(dob, first)), max(age_on(dob, last))) else integer(0)
  rows <- lapply(ages, function(a) {
    days <- pmax(0, as.numeric(pmin(last, add_years(dob, a + 1L) - 1L) - pmax(first, add_years(dob, a)) + 1))
    data.frame(age = as.character(a), adults_followed = sum(days > 0), events = sum(event & age_on(dob, last) == a),
               person_years = sum(days) / 365.25)
  })
  rows[[length(rows) + 1]] <- data.frame(age = if (length(ages)) paste(range(ages), collapse = " to ") else "All",
                                         adults_followed = sum(k), events = sum(event), person_years = sum(fu$time_days[k]) / 365.25)
  out <- cbind(nation = nation, do.call(rbind, rows))
  out$person_years <- sdc_stat(round(out$person_years), out$adults_followed, rule)
  out$adults_followed <- sdc_count(out$adults_followed, rule)
  out$events <- sdc_count(out$events, rule)
  out$rate_per_1000py <- ifelse(out$person_years > 0, 1000 * out$events / out$person_years, NA_real_)   # from the released (rounded) counts
  out
}

# Tabulate RZV uptake by birth month, the proportion whose jump at the cutoff is delta (SAP 6.5.1). Receipt
# is at least one recorded RZV dose by the end of the adult's follow-up for the outcome (follow_up() in
# derive.R), so it can differ between dementia and HZ. run_nation.R passes the dementia follow-up, which
# gives the uptake behind the dementia delta.
uptake_by_bin <- function(co, fu, nation, rule) {
  k <- fu$at_risk
  if (!any(k)) return(data.frame(nation = character(0), bin_start = as.Date(character(0)),
                                 bin_middle = as.Date(character(0)), eligible = integer(0),
                                 n = numeric(0), vaccinated = numeric(0), pct_vaccinated = numeric(0)))
  agg <- rowsum(cbind(n = 1, vacc = fu$vacc_by_end[k]), as.character(co$month_start[k]))
  bins <- as.Date(rownames(agg))
  out <- data.frame(nation = nation, bin_start = bins, bin_middle = month_middle(bins),
                    eligible = as.integer(bins >= STUDY$cutoff),
                    n = sdc_count(agg[, "n"], rule), vaccinated = sdc_count(agg[, "vacc"], rule))
  out$pct_vaccinated <- 100 * out$vaccinated / out$n
  rownames(out) <- NULL
  out
}
