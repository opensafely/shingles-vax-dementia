# The SAP 5.2.1 cohort criteria, applied in order, and the selection tables (SAP 6.3). One cohort, selected
# with the dementia follow-up (fu_primary), serves both outcomes. Dementia up to day 30 excludes adults
# from both analyses. Outcome dates start on day 31 (README, convention 2), so HZ before then excludes no one.
#
# Author/contact: Samantha Ip (hyi20@cam.ac.uk)

cohort_criteria <- function(d, nation, study, fu_primary) {
  has_dob <- !is.na(d$dob)
  # Count dementia medication on days 0 to 30 as dementia in the induction period, when the extract holds
  # out_date_dementia_med (SAP 5.2.1.2).
  med <- if ("out_date_dementia_med" %in% names(d)) d$out_date_dementia_med else rep(as.Date(NA), nrow(d))
  med_early <- (med >= d$index_date & med <= d$index_date + study$induction_days) %in% TRUE
  # Missing inclusion flags fail. Missing exclusion flags do not exclude.
  list(
    # 5.2.1, the criteria for every adult
    "Birth month or week recorded"                  = has_dob,
    "Sex recorded as female or male"                = d$cov_cat_sex %in% c("Female", "Male"),
    "Alive on the index date"                       = has_dob & (is.na(d$var_date_death) | d$var_date_death > d$index_date),
    "Living in the nation on the index date"        = d$var_bin_residency %in% TRUE,
    "Area recorded (urban or rural, and region)"    = d$cov_cat_urban!='Unknown' & (nation != "England" | d$cov_cat_region!='Unknown'),
    "Registered with a GP on the index date"        = d$var_bin_registration %in% TRUE,
    # 5.2.1.1, inclusion
    "Aged 63 to 66 on 31 August 2023"               = in_bandwidth(d, study, study$bandwidths[[study$primary_bandwidth]]),
    # 5.2.1.2, exclusions. Its first bullet, dementia before the index date or in the induction period, is one
    # step. pre_bin_dementia covers dementia diagnosis or medication up to day 30 after the index date. A
    # dementia medication date on days 0 to 30 also counts, as does a dementia outcome date on or before day 30
    # (an event at time 0 or less), which convention 2 rules out and data_checks.csv counts.
    "No dementia up to day 30 after the index date" =
      !(d$pre_bin_dementia %in% TRUE | (fu_primary$event == 1L & fu_primary$time_days <= 0L) %in% TRUE | med_early),
    "Not immunocompromised"                         = !(d$pre_bin_ic %in% TRUE),
    # the flag counts from the earliest possible 65th birthday (README, "What the extract must contain")
    "No RZV or ZVL before the earliest possible index date" = !(d$pre_bin_vax_shingles %in% TRUE),
    # follow-up starts on day 31 (Table 3), so an adult who dies or leaves by day 30 has none
    "In follow-up on day 31 after the index date"   = fu_primary$at_risk
  )
}

# Return the position of each adult's first failed criterion, or 0 if all criteria are met.
first_failed <- function(criteria) {
  fails <- !do.call(cbind, criteria)
  if (!nrow(fails)) return(integer(0))
  ff <- max.col(fails * 1L, ties.method = "first")
  ff[rowSums(fails) == 0] <- 0L
  ff
}

# Return sequential selection counts: remaining after each criterion and excluded at that criterion. Start
# with the extract count, or with the count remaining after criterion `from`.
flow_counts <- function(ff, labels, from = 0L) {
  k <- from:length(labels)
  data.frame(step = k,
             criterion = c("Extracted", labels)[k + 1],
             remaining = vapply(k, function(i) sum(ff == 0L | ff > i), numeric(1)),
             excluded = vapply(k, function(i) if (i == from) NA_real_ else sum(ff == i), numeric(1)))
}

# Tabulate selection counts overall and by eligibility group.
selection_table <- function(d, ff, labels, nation, rule) {
  group <- ifelse(is.na(d$z), "Birth date missing", ifelse(d$z == 1L, "Eligible", "Ineligible"))
  parts <- lapply(c("All", "Eligible", "Ineligible", "Birth date missing"), function(g) {
    keep <- if (g == "All") rep(TRUE, length(ff)) else group == g
    cbind(group = g, flow_counts(ff[keep], labels))
  })
  out <- do.call(rbind, parts)
  out$remaining <- sdc_count(out$remaining, rule)
  out$excluded <- sdc_count(out$excluded, rule)
  cbind(nation = nation, out)
}

# Tabulate birth-month selection counts within the primary bandwidth.
selection_by_bin <- function(d, ff, labels, nation, study, rule) {
  band <- which(labels == "Aged 63 to 66 on 31 August 2023")
  inside <- in_bandwidth(d, study, study$bandwidths[[study$primary_bandwidth]])
  bins <- sort(unique(d$month_start[inside]))
  parts <- lapply(seq_along(bins), function(i) {
    b <- bins[i]
    keep <- inside & d$month_start == b
    cbind(bin_start = b, eligible = d$z[keep][1], flow_counts(ff[keep], labels, from = band))
  })
  out <- if (length(parts)) do.call(rbind, parts) else
    data.frame(bin_start = as.Date(character(0)), eligible = integer(0), step = integer(0),
               criterion = character(0), remaining = numeric(0), excluded = numeric(0))
  out$remaining <- sdc_count(out$remaining, rule)
  out$excluded <- sdc_count(out$excluded, rule)
  cbind(nation = rep(nation, nrow(out)), out)
}
